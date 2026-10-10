# =============================================================================
# Classifier Contract & Accounting Adapter (Minimap2 & Kraken2/Bracken)
# =============================================================================

suppressMessages({
  library(digest)
})

`%||%` <- function(x, fallback) {
  if (is.null(x) || length(x) == 0L) fallback else x
}

RANKS_7_KRAKEN <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
RANKS_8_ABUNDANCE <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")

is_canonical_taxid <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^(0|[1-9][0-9]{0,15})$", x) &&
    (nchar(x) < 16L || x <= "9007199254740991")
}

is_tsv_safe_identifier <- function(x) {
  is.character(x) & !is.na(x) & nzchar(x) & x == trimws(x) &
    !grepl("[[:cntrl:]]", x, perl = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Parse Raw Kraken 2 Report (6 columns)
# -----------------------------------------------------------------------------
# Columns:
# 1. Percentage of fragments covered by clade rooted at this taxon
# 2. Number of fragments covered by clade rooted at this taxon
# 3. Number of fragments assigned directly to this taxon
# 4. Rank code (U, -/R, D, K, P, C, O, F, G, S, etc.)
# 5. NCBI taxonomy ID
# 6. Indented scientific name
# -----------------------------------------------------------------------------

parse_kraken_report <- function(path, sample_id = NULL) {
  sample_label <- if (!is.null(sample_id)) sprintf(" for sample '%s'", sample_id) else ""
  if (is.null(path) || !is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop(sprintf("Kraken2 report%s not found: '%s'", sample_label, path %||% "NULL"), call. = FALSE)
  }

  lines <- readLines(path, warn = FALSE)
  nonempty <- nzchar(lines)
  if (!any(nonempty)) {
    stop(sprintf("Kraken2 report%s is empty: '%s'", sample_label, path), call. = FALSE)
  }

  lines <- lines[nonempty]
  # Appending sentinel tab preserves trailing empty fields
  fields <- strsplit(paste0(lines, "\t"), "\t", fixed = TRUE)
  field_counts <- lengths(fields)
  if (any(field_counts != 6L)) {
    bad_line <- which(field_counts != 6L)[1]
    stop(sprintf(
      "Kraken2 report%s has %d fields at line %d; expected exactly 6.",
      sample_label, field_counts[bad_line], bad_line
    ), call. = FALSE)
  }

  pct_str <- vapply(fields, `[[`, character(1), 1L)
  clade_str <- vapply(fields, `[[`, character(1), 2L)
  direct_str <- vapply(fields, `[[`, character(1), 3L)
  rank_code <- vapply(fields, `[[`, character(1), 4L)
  taxid_str <- vapply(fields, `[[`, character(1), 5L)
  name_raw <- vapply(fields, `[[`, character(1), 6L)

  # Validate TaxIDs
  invalid_taxid <- !vapply(taxid_str, is_canonical_taxid, logical(1))
  if (any(invalid_taxid)) {
    bad_line <- which(invalid_taxid)[1]
    stop(sprintf("Kraken2 report%s has non-canonical TaxID '%s' at line %d.",
                 sample_label, taxid_str[bad_line], bad_line), call. = FALSE)
  }
  if (anyDuplicated(taxid_str)) {
    dup_id <- taxid_str[duplicated(taxid_str)][1]
    stop(sprintf("Kraken2 report%s contains duplicate TaxID: '%s'.",
                 sample_label, dup_id), call. = FALSE)
  }

  # Validate counts
  clade_num <- suppressWarnings(as.numeric(clade_str))
  direct_num <- suppressWarnings(as.numeric(direct_str))
  pct_num <- suppressWarnings(as.numeric(pct_str))

  if (anyNA(clade_num) || any(!is.finite(clade_num)) || any(clade_num < 0) ||
      any(clade_num != floor(clade_num)) || any(clade_num > 9007199254740991)) {
    bad_line <- which(is.na(clade_num) | !is.finite(clade_num) | clade_num < 0 |
                        clade_num != floor(clade_num) | clade_num > 9007199254740991)[1]
    stop(sprintf("Kraken2 report%s has invalid clade count '%s' at line %d.",
                 sample_label, clade_str[bad_line], bad_line), call. = FALSE)
  }

  if (anyNA(direct_num) || any(!is.finite(direct_num)) || any(direct_num < 0) ||
      any(direct_num != floor(direct_num)) || any(direct_num > 9007199254740991)) {
    bad_line <- which(is.na(direct_num) | !is.finite(direct_num) | direct_num < 0 |
                        direct_num != floor(direct_num) | direct_num > 9007199254740991)[1]
    stop(sprintf("Kraken2 report%s has invalid direct count '%s' at line %d.",
                 sample_label, direct_str[bad_line], bad_line), call. = FALSE)
  }

  if (any(direct_num > clade_num)) {
    bad_line <- which(direct_num > clade_num)[1]
    stop(sprintf("Kraken2 report%s has direct count (%g) exceeding clade count (%g) at line %d.",
                 sample_label, direct_num[bad_line], clade_num[bad_line], bad_line), call. = FALSE)
  }

  # Locate unclassified and root rows
  u_idx <- which(rank_code == "U" & taxid_str == "0")
  if (length(u_idx) == 0L) {
    # Check if there is an unclassified row by rank U or TaxID 0
    u_idx <- which(taxid_str == "0")
  }
  if (length(u_idx) > 1L) {
    stop(sprintf("Kraken2 report%s contains multiple unclassified rows.", sample_label), call. = FALSE)
  }
  raw_u <- if (length(u_idx) == 1L) clade_num[u_idx] else 0

  # Root clade: TaxID 1 or rank R
  root_idx <- which(taxid_str == "1" | rank_code == "R")
  if (length(root_idx) == 0L) {
    # If no root row exists, root clade count C is sum of all non-zero direct counts or top-level clades
    # In standard Kraken2 reports, root TaxID 1 exists whenever reads are classified
    raw_c <- 0
  } else {
    raw_c <- clade_num[root_idx[1]]
  }

  raw_total <- raw_c + raw_u

  # Direct counts dictionary
  direct_map <- stats::setNames(direct_num, taxid_str)
  clade_map <- stats::setNames(clade_num, taxid_str)

  # Check that sum of all direct counts equals raw_total
  sum_direct <- sum(direct_num)
  if (sum_direct != raw_total) {
    stop(sprintf(
      "Kraken2 report%s direct count sum (%g) does not equal root + unclassified total (%g).",
      sample_label, sum_direct, raw_total
    ), call. = FALSE)
  }

  report_df <- data.frame(
    pct = pct_num,
    clade_reads = clade_num,
    reads_clade = clade_num,
    direct_reads = direct_num,
    reads_direct = direct_num,
    rank_code = rank_code,
    rank = rank_code,
    taxid = taxid_str,
    name = trimws(name_raw),
    raw_name = name_raw,
    stringsAsFactors = FALSE
  )

  list(
    total_reads = raw_total,
    classified_reads = raw_c,
    unclassified_reads = raw_u,
    raw_total = raw_total,
    raw_classified = raw_c,
    raw_unclassified = raw_u,
    direct_counts = direct_map,
    clade_counts = clade_map,
    report_df = report_df,
    df = report_df,
    path = path,
    sample_id = sample_id
  )
}

# -----------------------------------------------------------------------------
# 2. Parse Headered Bracken Abundance Report (7 columns)
# -----------------------------------------------------------------------------
# Columns:
# name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads,
# new_est_reads, fraction_total_reads
# -----------------------------------------------------------------------------

parse_bracken_abundance <- function(path, sample_id = NULL) {
  sample_label <- if (!is.null(sample_id)) sprintf(" for sample '%s'", sample_id) else ""
  if (is.null(path) || !is.character(path) || length(path) != 1L || !file.exists(path)) {
    stop(sprintf("Bracken report%s not found: '%s'", sample_label, path %||% "NULL"), call. = FALSE)
  }

  header_lines <- readLines(path, n = 2L, warn = FALSE)
  if (length(header_lines) == 0L) {
    stop(sprintf("Bracken report%s is empty: '%s'", sample_label, path), call. = FALSE)
  }

  header <- strsplit(header_lines[1], "\t", fixed = TRUE)[[1]]
  expected_header <- c("name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads",
                       "added_reads", "new_est_reads", "fraction_total_reads")
  if (!identical(header, expected_header)) {
    stop(sprintf(
      "Bracken report%s header violation. Expected: %s. Found: %s.",
      sample_label, paste(expected_header, collapse = "\t"), paste(header, collapse = "\t")
    ), call. = FALSE)
  }

  df <- read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(df) == 0L) {
    return(list(
      species_estimate_sum = 0,
      rows = df,
      taxid_to_est = stats::setNames(numeric(0), character(0)),
      taxid_to_name = stats::setNames(character(0), character(0)),
      path = path,
      sample_id = sample_id
    ))
  }

  taxids <- as.character(df$taxonomy_id)
  invalid_taxids <- !vapply(taxids, is_canonical_taxid, logical(1)) | taxids == "0"
  if (any(invalid_taxids)) {
    bad_line <- which(invalid_taxids)[1]
    stop(sprintf("Bracken report%s contains non-positive or invalid TaxID '%s' at data row %d.",
                 sample_label, taxids[bad_line], bad_line), call. = FALSE)
  }
  if (anyDuplicated(taxids)) {
    dup_id <- taxids[duplicated(taxids)][1]
    stop(sprintf("Bracken report%s contains duplicate TaxID: '%s'.", sample_label, dup_id), call. = FALSE)
  }

  if (any(df$taxonomy_lvl != "S")) {
    bad_line <- which(df$taxonomy_lvl != "S")[1]
    stop(sprintf("Bracken report%s has taxonomy_lvl '%s' at data row %d; expected species rank 'S'.",
                 sample_label, df$taxonomy_lvl[bad_line], bad_line), call. = FALSE)
  }

  count_cols <- c("kraken_assigned_reads", "added_reads", "new_est_reads")
  for (col in count_cols) {
    vals <- df[[col]]
    if (!is.numeric(vals)) vals <- suppressWarnings(as.numeric(vals))
    if (anyNA(vals) || any(!is.finite(vals)) || any(vals < 0) ||
        any(vals != floor(vals)) || any(vals > 9007199254740991)) {
      stop(sprintf("Bracken report%s column '%s' contains non-finite, negative, or non-integer counts.",
                   sample_label, col), call. = FALSE)
    }
    df[[col]] <- vals
  }

  est_sum <- sum(df$new_est_reads)
  if (est_sum > 9007199254740991) {
    stop(sprintf("Bracken report%s new_est_reads sum exceeds max safe integer.", sample_label), call. = FALSE)
  }

  taxid_to_est <- stats::setNames(df$new_est_reads, taxids)
  taxid_to_name <- stats::setNames(as.character(df$name), taxids)

  list(
    species_estimate_sum = est_sum,
    rows = df,
    df = df,
    taxid_to_est = taxid_to_est,
    taxid_to_name = taxid_to_name,
    path = path,
    sample_id = sample_id
  )
}

# -----------------------------------------------------------------------------
# 3. TaxID-to-Lineage Mapping and Narrow Placeholder Normalization
# -----------------------------------------------------------------------------

normalize_placeholder_rank <- function(rank_str) {
  # Map "unclassified <Group> <rank_type>" to "<Group>_Incertae_sedis"
  # Exactly as established from producer evidence.
  sub("^unclassified (.+) (family|genus|order|class|phylum)$", "\\1_Incertae_sedis", rank_str)
}

normalize_placeholder_species <- function(rank_str) {
  normalize_placeholder_rank(rank_str)
}

build_taxid_lineage_map_from_assignments <- function(assignments_df) {
  c_reads <- assignments_df[assignments_df$status == "C" & assignments_df$taxid != "0", , drop = FALSE]
  if (nrow(c_reads) == 0L) return(list(taxid_to_lineage = character(0), conflicts = list()))

  # Aggregate by taxid
  split_by_taxid <- split(c_reads$lineage, c_reads$taxid)
  taxid_map <- character(length(split_by_taxid))
  names(taxid_map) <- names(split_by_taxid)

  conflicts <- list()
  for (tid in names(split_by_taxid)) {
    lins <- split_by_taxid[[tid]]
    tab <- table(lins)
    max_count <- max(tab)
    winners <- names(tab)[tab == max_count]
    # Deterministic winner: alphabetically first if tie
    winner <- sort(winners)[1]
    taxid_map[tid] <- winner
    if (length(tab) > 1L) {
      conflicts[[tid]] <- list(taxid = tid, winner = winner, counts = as.list(tab))
    }
  }

  list(taxid_to_lineage = taxid_map, conflicts = conflicts)
}

# Ranks matcher: checks if an 8-rank abundance row matches a 7-rank assignment lineage
# by omitting kingdom (index 2) and allowing narrow placeholder normalization.
match_abundance_to_assignment_lineage <- function(abundance_lineage, assignment_lineage) {
  ab_parts <- strsplit(abundance_lineage, ";", fixed = TRUE)[[1]]
  as_parts <- strsplit(assignment_lineage, "|", fixed = TRUE)[[1]]

  if (length(ab_parts) != 8L || length(as_parts) != 7L) return(FALSE)

  # 8 ranks: superkingdom, kingdom, phylum, class, order, family, genus, species
  # 7 ranks: superkingdom, phylum, class, order, family, genus, species
  ab_7 <- ab_parts[c(1L, 3L:8L)]

  for (i in 1:7) {
    ab_r <- ab_7[i]
    as_r <- as_parts[i]
    if (identical(ab_r, as_r)) next
    if (identical(ab_r, normalize_placeholder_rank(as_r))) next
    return(FALSE)
  }
  TRUE
}

# -----------------------------------------------------------------------------
# 4. Match Bracken Rows Against Abundance Table
# -----------------------------------------------------------------------------

match_bracken_to_abundance <- function(bracken_data, count_matrix, taxonomy_df,
                                       assignments_df = NULL, taxid_map = NULL,
                                       sample_id = NULL) {
  sample_label <- if (!is.null(sample_id)) sprintf(" for sample '%s'", sample_id) else ""

  if (is.null(taxid_map) && !is.null(assignments_df)) {
    taxid_map <- build_taxid_lineage_map_from_assignments(assignments_df)$taxid_to_lineage
  }
  if (is.null(taxid_map)) {
    taxid_map <- character(0)
  }

  # All positive taxa in abundance matrix for this sample
  unclass_mask <- taxonomy_df$superkingdom == "Unclassified"
  positive_paths <- rownames(count_matrix)[!unclass_mask]
  sample_counts <- count_matrix[positive_paths, sample_id]

  # Only taxa with positive counts in this sample
  active_ab_paths <- positive_paths[sample_counts > 0]
  active_ab_counts <- sample_counts[sample_counts > 0]

  bracken_df <- bracken_data$rows
  if (nrow(bracken_df) == 0L) {
    if (length(active_ab_paths) > 0L) {
      stop(sprintf("Bracken report%s has 0 rows but abundance table has %d positive species.",
                   sample_label, length(active_ab_paths)), call. = FALSE)
    }
    return(TRUE)
  }

  # Verify every Bracken row matches an abundance row
  matched_ab_paths <- character(0)
  for (i in seq_len(nrow(bracken_df))) {
    row <- bracken_df[i, ]
    tid <- as.character(row$taxonomy_id)
    b_name <- as.character(row$name)
    b_est <- as.numeric(row$new_est_reads)

    # Resolve assignment lineage for this TaxID
    as_lineage <- if (!is.null(taxid_map) && tid %in% names(taxid_map)) taxid_map[[tid]] else NULL
    if (is.null(as_lineage) || !nzchar(as_lineage)) {
      # When assignments are absent or TaxID not in assignments,
      # attempt exact species match in taxonomy_df if unique
      matching_idx <- which(taxonomy_df$species == b_name & !unclass_mask)
      if (length(matching_idx) == 1L) {
        matched_path <- taxonomy_df$TaxonPath[matching_idx]
      } else {
        stop(sprintf(
          "Bracken row TaxID '%s' ('%s')%s cannot be resolved to an abundance row without assignment or taxonomy evidence.",
          tid, b_name, sample_label
        ), call. = FALSE)
      }
    } else {
      # Match using 7-rank comparison with placeholder normalization
      matching_paths <- active_ab_paths[vapply(active_ab_paths, function(p) {
        match_abundance_to_assignment_lineage(p, as_lineage)
      }, logical(1))]

      if (length(matching_paths) == 0L) {
        # Try matching against all positive paths even if count is 0
        matching_paths <- positive_paths[vapply(positive_paths, function(p) {
          match_abundance_to_assignment_lineage(p, as_lineage)
        }, logical(1))]
      }

      if (length(matching_paths) != 1L) {
        stop(sprintf(
          "Bracken TaxID '%s' ('%s')%s matches %d abundance rows; expected exactly 1.",
          tid, b_name, sample_label, length(matching_paths)
        ), call. = FALSE)
      }
      matched_path <- matching_paths[1]
    }

    if (matched_path %in% matched_ab_paths) {
      stop(sprintf("Multiple Bracken rows matched abundance row '%s'%s.",
                   matched_path, sample_label), call. = FALSE)
    }
    matched_ab_paths <- c(matched_ab_paths, matched_path)

    # Verify count equality
    ab_count <- as.numeric(count_matrix[matched_path, sample_id])
    if (ab_count != b_est) {
      stop(sprintf(
        "Count mismatch%s for TaxID '%s' ('%s'): Bracken new_est_reads (%.0f) != abundance count (%.0f).",
        sample_label, tid, b_name, b_est, ab_count
      ), call. = FALSE)
    }
  }

  # Verify no extra positive species in abundance table
  extra_ab <- setdiff(active_ab_paths, matched_ab_paths)
  if (length(extra_ab) > 0L) {
    stop(sprintf(
      "Abundance table%s contains %d species not accounted for in Bracken report (e.g. '%s').",
      sample_label, length(extra_ab), extra_ab[1]
    ), call. = FALSE)
  }

  TRUE
}

# -----------------------------------------------------------------------------
# 5. Population Accounting Model
# -----------------------------------------------------------------------------

compute_kraken_sample_accounting <- function(kraken_report, bracken_report,
                                            abundance_unclassified, abundance_species_sum,
                                            assignment_data = NULL, sample_id = NULL) {
  sample_label <- if (!is.null(sample_id)) sprintf(" for sample '%s'", sample_id) else ""

  N <- kraken_report$raw_total
  C <- kraken_report$raw_classified
  U <- kraken_report$raw_unclassified
  E <- bracken_report$species_estimate_sum
  T_ab <- E + abundance_unclassified
  D <- C - E

  # Equalities to enforce:
  # 1. 0 <= E <= C
  if (E < 0 || E > C) {
    stop(sprintf(
      "Accounting failure%s: Bracken species estimates sum E (%.0f) must satisfy 0 <= E <= C (%.0f).",
      sample_label, E, C
    ), call. = FALSE)
  }

  # 2. Abundance unclassified must equal U
  if (abundance_unclassified != U) {
    stop(sprintf(
      "Accounting failure%s: abundance unclassified (%.0f) != Kraken report unclassified U (%.0f).",
      sample_label, as.numeric(abundance_unclassified), as.numeric(U)
    ), call. = FALSE)
  }

  # 3. Abundance classified species sum must equal E
  if (abundance_species_sum != E) {
    stop(sprintf(
      "Accounting failure%s: abundance classified species sum (%.0f) != Bracken estimate sum E (%.0f).",
      sample_label, as.numeric(abundance_species_sum), as.numeric(E)
    ), call. = FALSE)
  }

  assignments_checked <- FALSE
  if (!is.null(assignment_data)) {
    n_assign <- nrow(assignment_data)
    n_c_assign <- sum(assignment_data$status == "C")
    n_u_assign <- sum(assignment_data$status == "U")
    c0 <- sum(assignment_data$status == "C" & assignment_data$taxid == "0")

    if (n_assign != N) {
      stop(sprintf(
        "Accounting failure%s: assignment rows (%.0f) != Kraken report total N (%.0f).",
        sample_label, as.numeric(n_assign), as.numeric(N)
      ), call. = FALSE)
    }
    if (n_c_assign != C) {
      stop(sprintf(
        "Accounting failure%s: assignment status C rows (%.0f) != Kraken report classified C (%.0f).",
        sample_label, as.numeric(n_c_assign), as.numeric(C)
      ), call. = FALSE)
    }
    if (n_u_assign != U) {
      stop(sprintf(
        "Accounting failure%s: assignment status U rows (%.0f) != Kraken report unclassified U (%.0f).",
        sample_label, as.numeric(n_u_assign), as.numeric(U)
      ), call. = FALSE)
    }
    if (c0 != 0L) {
      stop(sprintf(
        "Accounting failure%s: assignment file contains %d status-C reads with TaxID 0.",
        sample_label, c0
      ), call. = FALSE)
    }

    # Verify per-TaxID direct counts between assignments and raw report
    assign_counts <- table(assignment_data$taxid)
    report_directs <- kraken_report$direct_counts[kraken_report$direct_counts > 0]
    assign_pos <- assign_counts[names(assign_counts) != "0"]
    report_pos <- report_directs[names(report_directs) != "0"]

    if (!setequal(names(assign_pos), names(report_pos))) {
      diff_taxids <- setdiff(union(names(assign_pos), names(report_pos)),
                             intersect(names(assign_pos), names(report_pos)))
      stop(sprintf(
        "Direct count taxid set mismatch%s between assignments and Kraken report: %s",
        sample_label, paste(head(diff_taxids, 5), collapse = ", ")
      ), call. = FALSE)
    }
    for (tid in names(assign_pos)) {
      if (as.numeric(assign_pos[[tid]]) != as.numeric(report_pos[[tid]])) {
        stop(sprintf(
          "Direct count mismatch%s for TaxID '%s': assignment (%.0f) != Kraken report (%.0f).",
          sample_label, tid, as.numeric(assign_pos[[tid]]), as.numeric(report_pos[[tid]])
        ), call. = FALSE)
      }
    }
    assignments_checked <- TRUE
  }

  list(
    N = N,
    C = C,
    U = U,
    E = E,
    T = T_ab,
    D = D,
    assignments_checked = assignments_checked,
    reconciliation_pass = TRUE
  )
}

# -----------------------------------------------------------------------------
# 6. Compare Upstream Params Between Cohort Samples
# -----------------------------------------------------------------------------

compare_upstream_params <- function(rep_params, sample_params, sample_id) {
  check_fields <- c("classifier", "database_set", "taxonomic_rank", "min_len", "max_len",
                    "min_read_qual", "abundance_threshold")
  classifier <- rep_params$classifier
  if (identical(classifier, "minimap2")) {
    check_fields <- c(check_fields, "min_percent_identity", "min_ref_coverage")
  } else if (identical(classifier, "kraken2")) {
    check_fields <- c(check_fields, "kraken2_confidence", "bracken_threshold", "bracken_length")
  }

  for (field in check_fields) {
    rep_val <- rep_params[[field]]
    sample_val <- sample_params[[field]]
    if (!identical(rep_val, sample_val)) {
      stop(sprintf(
        "Upstream params contract mismatch for sample '%s' on field '%s': representative '%s' != sample '%s'.",
        sample_id, field, as.character(rep_val %||% "NULL"), as.character(sample_val %||% "NULL")
      ), call. = FALSE)
    }
  }

  # Compare database resources
  db_set <- rep_params$database_set
  rep_db <- rep_params$database_sets[[db_set]]
  sample_db <- sample_params$database_sets[[db_set]]
  for (res in c("reference", "database", "ref2taxid", "taxonomy")) {
    if (!identical(as.character(rep_db[[res]]), as.character(sample_db[[res]]))) {
      stop(sprintf(
        "Upstream database resource mismatch for sample '%s' on '%s.%s'.",
        sample_id, db_set, res
      ), call. = FALSE)
    }
  }

  invisible(TRUE)
}
