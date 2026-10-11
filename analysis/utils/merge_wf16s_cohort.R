#!/usr/bin/env Rscript
# =============================================================================
# wf-16s Cohort Merger Utility
# =============================================================================
# Combines multiple single-sample wf-16s run outputs into a unified cohort bundle
# with validated compatibility, full taxonomic outer-join, neutral placeholder
# metadata, and hashed per-source provenance.
# =============================================================================

get_merge_script_dir <- function() {
  if (!is.null(getOption("wf16s.pipeline_root"))) {
    cand <- file.path(getOption("wf16s.pipeline_root"), "analysis", "utils")
    if (dir.exists(cand)) return(normalizePath(cand, winslash = "/"))
  }
  for (i in seq_len(sys.nframe())) {
    fr <- sys.frame(i)
    if (!is.null(fr$ofile)) {
      cand <- dirname(normalizePath(fr$ofile, winslash = "/"))
      if (dir.exists(cand) && file.exists(file.path(cand, "classifier_contract.R"))) {
        return(cand)
      }
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    cand <- dirname(sub("^--file=", "", file_arg[1]))
    if (dir.exists(cand) && file.exists(file.path(cand, "classifier_contract.R"))) {
      return(normalizePath(cand, winslash = "/"))
    }
  }
  cur <- normalizePath(getwd(), winslash = "/")
  while (nchar(cur) > 3L) {
    cand <- file.path(cur, "analysis", "utils")
    if (dir.exists(cand) && file.exists(file.path(cand, "classifier_contract.R"))) {
      return(normalizePath(cand, winslash = "/"))
    }
    cur <- dirname(cur)
  }
  normalizePath(file.path(getwd(), "analysis", "utils"), winslash = "/")
}

locate_producer_output <- function(input_path) {
  if (!dir.exists(input_path)) {
    stop(sprintf("E_COHORT_INPUT_NOT_FOUND: input directory '%s' does not exist.", input_path), call. = FALSE)
  }
  norm_path <- normalizePath(input_path, winslash = "/", mustWork = TRUE)

  # 1. Direct output root
  if (file.exists(file.path(norm_path, "abundance_table_species.tsv"))) {
    return(norm_path)
  }
  # 2. Check norm_path/output
  if (file.exists(file.path(norm_path, "output", "abundance_table_species.tsv"))) {
    return(normalizePath(file.path(norm_path, "output"), winslash = "/"))
  }
  # 3. Check single nested subdirectory containing output/abundance_table_species.tsv
  subdirs <- list.dirs(norm_path, recursive = FALSE, full.names = TRUE)
  candidates <- c(subdirs, file.path(subdirs, "output"))
  candidates <- unique(candidates[file.exists(file.path(candidates, "abundance_table_species.tsv"))])
  if (length(candidates) == 1L) return(normalizePath(candidates, winslash = "/"))
  if (length(candidates) > 1L) stop("E_COHORT_INPUT_AMBIGUOUS: multiple producer outputs found.", call. = FALSE)
  stop(sprintf("E_COHORT_INPUT_INVALID: could not locate abundance_table_species.tsv in '%s'.", input_path), call. = FALSE)
}

locate_params_json <- function(output_dir) {
  if (file.exists(file.path(output_dir, "params.json"))) {
    return(normalizePath(file.path(output_dir, "params.json"), winslash = "/"))
  }
  parent <- dirname(output_dir)
  if (file.exists(file.path(parent, "params.json"))) {
    return(normalizePath(file.path(parent, "params.json"), winslash = "/"))
  }
  stop(sprintf("E_COHORT_INPUT_INVALID: could not locate params.json for output directory '%s'.", output_dir), call. = FALSE)
}

locate_sample_file <- function(search_dirs, filenames) {
  candidates <- unique(unlist(lapply(search_dirs, function(d) file.path(d, filenames))))
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) > 1L) stop("E_COHORT_INPUT_AMBIGUOUS: multiple sample files found.", call. = FALSE)
  if (length(candidates)) normalizePath(candidates, winslash = "/") else NULL
}

locate_assignments_file <- function(output_dir, sample_id) {
  search_dirs <- c(
    file.path(output_dir, "reads_assignments"),
    file.path(output_dir, "assignments"),
    output_dir
  )
  names <- paste0(sample_id, c(".assignments.tsv", "_lineages.kraken2.assignments.tsv",
                              "_lineages.minimap2.assignments.tsv"))
  locate_sample_file(search_dirs, c(names, paste0(names, ".gz")))
}

locate_kraken_report_file <- function(output_dir, sample_id) {
  search_dirs <- c(
    file.path(output_dir, "kraken2"),
    output_dir
  )
  locate_sample_file(search_dirs, paste0(sample_id, c(".kraken2.report.txt", ".kreport")))
}

locate_bracken_report_file <- function(output_dir, sample_id) {
  search_dirs <- c(
    file.path(output_dir, "bracken"),
    output_dir
  )
  locate_sample_file(search_dirs, paste0(sample_id, c(".kraken2_bracken.report", ".bracken.report")))
}

parse_merge_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  inputs <- character(0)
  output_dir <- NULL
  project_name <- NULL
  validate_only <- FALSE
  overwrite <- FALSE
  taxonomy_cache <- NULL

  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("-i", "--input", "--inputs")) {
      if (i + 1L > length(args)) stop(sprintf("Option '%s' requires an argument.", arg), call. = FALSE)
      val <- args[i + 1L]
      inputs <- c(inputs, if (arg == "--inputs") strsplit(val, "[,;]")[[1]] else val)
      i <- i + 2L
    } else if (startsWith(arg, "--input=")) {
      val <- sub("^--input=", "", arg)
      inputs <- c(inputs, val)
      i <- i + 1L
    } else if (startsWith(arg, "--inputs=")) {
      val <- sub("^--inputs=", "", arg)
      inputs <- c(inputs, strsplit(val, "[,;]")[[1]])
      i <- i + 1L
    } else if (arg %in% c("-o", "--output-dir")) {
      if (i + 1L > length(args)) stop(sprintf("Option '%s' requires an argument.", arg), call. = FALSE)
      output_dir <- args[i + 1L]
      i <- i + 2L
    } else if (startsWith(arg, "--output-dir=")) {
      output_dir <- sub("^--output-dir=", "", arg)
      i <- i + 1L
    } else if (arg %in% c("-p", "--project-name")) {
      if (i + 1L > length(args)) stop(sprintf("Option '%s' requires an argument.", arg), call. = FALSE)
      project_name <- args[i + 1L]
      i <- i + 2L
    } else if (startsWith(arg, "--project-name=")) {
      project_name <- sub("^--project-name=", "", arg)
      i <- i + 1L
    } else if (arg == "--validate-only") {
      validate_only <- TRUE
      i <- i + 1L
    } else if (arg == "--overwrite") {
      overwrite <- TRUE
      i <- i + 1L
    } else if (arg == "--taxonomy-cache") {
      if (i + 1L > length(args)) stop(sprintf("Option '%s' requires an argument.", arg), call. = FALSE)
      taxonomy_cache <- args[i + 1L]
      i <- i + 2L
    } else if (startsWith(arg, "--taxonomy-cache=")) {
      taxonomy_cache <- sub("^--taxonomy-cache=", "", arg)
      i <- i + 1L
    } else if (arg %in% c("-h", "--help")) {
      cat("Usage: Rscript merge_wf16s_cohort.R --input <dir1> --input <dir2> --output-dir <out> [options]\n")
      cat("Options:\n")
      cat("  -i, --input <dir>        Input wf-16s run output directory (repeatable, literal path)\n")
      cat("  --inputs <dirs>          Comma/semicolon-separated inputs; use --input for paths containing these characters\n")
      cat("  -o, --output-dir <dir>   Destination directory for merged cohort bundle\n")
      cat("  -p, --project-name <str> Combined cohort project name\n")
      cat("  --validate-only          Perform compatibility and join validation without writing outputs\n")
      cat("  --overwrite              Overwrite existing output directory\n")
      cat("  -h, --help               Show this help message\n")
      quit(status = 0)
    } else {
      stop(sprintf("Unrecognized option: '%s'", arg), call. = FALSE)
    }
  }

  inputs <- trimws(inputs)
  inputs <- inputs[nzchar(inputs)]

  list(
    inputs = inputs,
    output_dir = output_dir,
    project_name = project_name,
    taxonomy_cache = taxonomy_cache,
    validate_only = validate_only,
    overwrite = overwrite
  )
}

merge_wf16s_cohort <- function(input_dirs, output_dir, project_name = NULL,
                               taxonomy_cache = NULL,
                               validate_only = FALSE, overwrite = FALSE) {
  script_dir <- get_merge_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")

  # Source dependencies if not loaded
  if (!exists("compare_upstream_params", mode = "function")) {
    source(file.path(script_dir, "classifier_contract.R"))
  }
  if (!exists("read_upstream_params", mode = "function")) {
    source(file.path(script_dir, "io.R"))
  }
  if (!exists("compute_file_hash", mode = "function")) {
    source(file.path(script_dir, "atomic_io.R"))
  }
  if (!exists("read_pipeline_version", mode = "function")) {
    source(file.path(script_dir, "version.R"))
  }

  pipeline_version <- tryCatch(
    read_pipeline_version(file.path(repo_root, "VERSION")),
    error = function(e) "0.5.0"
  )

  if (length(input_dirs) < 2L) {
    stop("E_COHORT_INSUFFICIENT_INPUTS: cohort merge requires at least 2 input directories.", call. = FALSE)
  }
  if (is.null(output_dir) || !nzchar(trimws(output_dir))) {
    stop("E_COHORT_MISSING_OUTPUT: --output-dir must be specified.", call. = FALSE)
  }

  project_name <- if (!is.null(project_name) && nzchar(trimws(project_name))) {
    trimws(project_name)
  } else {
    basename(normalizePath(output_dir, winslash = "/", mustWork = FALSE))
  }

  # 1. Resolve and validate input paths
  resolved_roots <- character(length(input_dirs))
  for (idx in seq_along(input_dirs)) {
    resolved_roots[idx] <- locate_producer_output(input_dirs[idx])
  }

  if (anyDuplicated(tolower(resolved_roots)) > 0L) {
    stop("E_COHORT_DUPLICATE_SOURCE: duplicate input path provided.", call. = FALSE)
  }

  # 2. Guard against output directory overlap
  norm_dest <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  for (r_root in unique(c(resolved_roots, normalizePath(input_dirs, winslash = "/")))) {
    if (identical(tolower(norm_dest), tolower(r_root)) ||
        startsWith(tolower(norm_dest), paste0(tolower(r_root), "/")) ||
        startsWith(tolower(r_root), paste0(tolower(norm_dest), "/"))) {
      stop(sprintf("E_COHORT_OUTPUT_OVERLAP: output directory '%s' overlaps with input directory '%s'.",
                   output_dir, r_root), call. = FALSE)
    }
  }

  if (file.exists(norm_dest) && !dir.exists(norm_dest)) stop("E_COHORT_DESTINATION_EXISTS: destination is a file.", call. = FALSE)
  if (dir.exists(norm_dest)) {
    if (!isTRUE(overwrite)) stop("E_COHORT_DESTINATION_EXISTS: destination already exists.", call. = FALSE)
    ownership_path <- file.path(norm_dest, "cohort_merge_manifest.json")
    ownership <- tryCatch(jsonlite::fromJSON(ownership_path, simplifyVector = FALSE), error = function(e) NULL)
    owned <- names(ownership$outputs)
    existing <- list.files(norm_dest, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    if (is.null(ownership) || !length(owned) ||
        any(grepl("[/\\\\]|^[.]|[[:cntrl:]]", owned)) ||
        !setequal(existing, c(owned, "cohort_merge_manifest.json")) ||
        any(vapply(owned, function(p) !identical(compute_file_hash(file.path(norm_dest, p)),
                                               ownership$outputs[[p]]$sha256), logical(1)))) {
      stop("E_COHORT_DESTINATION_UNOWNED: overwrite requires an unchanged merger bundle with no unrelated files.", call. = FALSE)
    }
  }
  if (!is.null(taxonomy_cache)) {
    if (!file.exists(taxonomy_cache)) stop("E_COHORT_CACHE_MISSING: specified taxonomy cache does not exist.", call. = FALSE)
    if (startsWith(tolower(normalizePath(taxonomy_cache, winslash = "/")), paste0(tolower(norm_dest), "/"))) {
      stop("E_COHORT_OUTPUT_OVERLAP: taxonomy cache is inside destination.", call. = FALSE)
    }
  }
  identity_cache <- if (is.null(taxonomy_cache)) list() else jsonlite::fromJSON(taxonomy_cache, simplifyVector = FALSE)
  source_paths <- unique(c(unlist(lapply(resolved_roots, function(p) {
    list.files(p, pattern = "(params[.]json|abundance_table_species[.]tsv|assignments[.]tsv([.]gz)?|[.]report([.]txt)?|[.]kreport)$",
               full.names = TRUE, recursive = TRUE)
  })), taxonomy_cache, vapply(resolved_roots, locate_params_json, character(1))))
  source_hashes <- vapply(source_paths, compute_file_hash, character(1))
  assert_merge_sources <- function() {
    if (!identical(vapply(source_paths, compute_file_hash, character(1)), source_hashes)) {
      stop("E_COHORT_INPUT_CHANGED: source changed during merge.", call. = FALSE)
    }
  }

  # 3. Read params.json and extract upstream contract for each input
  params_list <- list()
  contracts <- list()
  sample_ids <- character(length(resolved_roots))
  params_paths <- character(length(resolved_roots))

  for (idx in seq_along(resolved_roots)) {
    r_root <- resolved_roots[idx]
    p_path <- locate_params_json(r_root)
    params_paths[idx] <- p_path
    p_data <- read_upstream_params(p_path)
    s_id <- p_data$sample
    if (is.null(s_id) || !nzchar(trimws(s_id))) {
      stop(sprintf("E_COHORT_INPUT_INVALID: missing sample ID in '%s'.", p_path), call. = FALSE)
    }
    sample_ids[idx] <- s_id
    params_list[[s_id]] <- p_data
    contracts[[s_id]] <- extract_upstream_contract(p_data)
  }

  if (anyDuplicated(sample_ids) > 0L) {
    dup <- sample_ids[duplicated(sample_ids)][1]
    stop(sprintf("E_COHORT_DUPLICATE_SAMPLE: duplicate sample ID detected: '%s'.", dup), call. = FALSE)
  }
  validate_sample_ids(sample_ids)

  rep_id <- sample_ids[1]
  rep_params <- params_list[[rep_id]]
  rep_contract <- contracts[[rep_id]]
  is_kraken <- identical(rep_contract$classifier, "kraken2")

  # 4. Validate contract compatibility across all samples
  for (i in 2:length(sample_ids)) {
    s_id <- sample_ids[i]
    s_params <- params_list[[s_id]]
    s_contract <- contracts[[s_id]]

    if (!identical(rep_contract$classifier, s_contract$classifier)) {
      stop(sprintf("E_COHORT_INCOMPATIBLE: classifier mismatch between '%s' (%s) and '%s' (%s).",
                   rep_id, rep_contract$classifier, s_id, s_contract$classifier), call. = FALSE)
    }

    tryCatch(
      compare_upstream_params(rep_params, s_params, s_id),
      error = function(e) {
        stop(sprintf("E_COHORT_INCOMPATIBLE: %s", e$message), call. = FALSE)
      }
    )
  }

  # 5. Read abundance tables and validate
  abundance_dfs <- list()
  abundance_paths <- character(length(resolved_roots))

  for (idx in seq_along(resolved_roots)) {
    r_root <- resolved_roots[idx]
    s_id <- sample_ids[idx]
    ab_path <- file.path(r_root, "abundance_table_species.tsv")
    abundance_paths[idx] <- normalizePath(ab_path, winslash = "/", mustWork = TRUE)

    validated <- read_abundance_table(ab_path)
    if (!identical(validated$samples, s_id)) {
      stop(sprintf("E_COHORT_MALFORMED_TABLE: expected only sample '%s' in '%s'.", s_id, ab_path), call. = FALSE)
    }
    df <- utils::read.delim(ab_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    if (!"tax" %in% names(df) || !"total" %in% names(df)) {
      stop(sprintf("E_COHORT_MALFORMED_TABLE: abundance table '%s' missing 'tax' or 'total' column.", ab_path), call. = FALSE)
    }
    non_meta <- setdiff(names(df), c("tax", "total"))
    if (!s_id %in% non_meta) {
      stop(sprintf("E_COHORT_MALFORMED_TABLE: sample '%s' column not found in '%s'.", s_id, ab_path), call. = FALSE)
    }
    if (nrow(df) == 0L || !startsWith(df$tax[1], "Unclassified")) {
      stop(sprintf("E_COHORT_MALFORMED_TABLE: row 1 in '%s' must be Unclassified.", ab_path), call. = FALSE)
    }
    if (anyDuplicated(df$tax) > 0L) {
      stop(sprintf("E_COHORT_DUPLICATE_TAXA: duplicate taxa in '%s'.", ab_path), call. = FALSE)
    }
    abundance_dfs[[s_id]] <- df
  }

  # 6. Locate assignment files and report files
  assignments_map <- list()
  kraken_reports_map <- list()
  bracken_reports_map <- list()

  for (idx in seq_along(resolved_roots)) {
    r_root <- resolved_roots[idx]
    s_id <- sample_ids[idx]

    as_file <- locate_assignments_file(r_root, s_id)
    if (!is.null(as_file)) {
      assignments_map[[s_id]] <- as_file
    }

    if (is_kraken) {
      kr_file <- locate_kraken_report_file(r_root, s_id)
      if (is.null(kr_file)) {
        stop(sprintf("E_COHORT_MISSING_REPORT: could not locate Kraken report for sample '%s' in '%s'.", s_id, r_root), call. = FALSE)
      }
      kraken_reports_map[[s_id]] <- kr_file

      br_file <- locate_bracken_report_file(r_root, s_id)
      if (is.null(br_file)) {
        stop(sprintf("E_COHORT_MISSING_REPORT: could not locate Bracken report for sample '%s' in '%s'.", s_id, r_root), call. = FALSE)
      }
      bracken_reports_map[[s_id]] <- br_file
    }
  }

  # Use the same ingestion and accounting contracts as the pipeline before
  # claiming validate-only success or publishing a bundle.
  sample_accounting <- list()
  for (idx in seq_along(sample_ids)) {
    s_id <- sample_ids[idx]
    ab <- read_abundance_table(abundance_paths[idx])
    u <- ab$count_matrix[ab$unclass_index, s_id]
    e <- sum(ab$count_matrix[-ab$unclass_index, s_id])
    assignments <- if (is.null(assignments_map[[s_id]])) NULL else read_assignments_file(
      assignments_map[[s_id]], s_id, classifier = rep_contract$classifier,
      expected_total = if (is_kraken) NULL else u + e,
      expected_classified = if (is_kraken) NULL else e, expected_unclassified = u)
    if (is_kraken) {
      kr <- parse_kraken_report(kraken_reports_map[[s_id]], s_id)
      br <- parse_bracken_abundance(bracken_reports_map[[s_id]], s_id)
      match_bracken_to_abundance(br, ab$count_matrix, ab$taxonomy, assignments_df = assignments,
                                 sample_id = s_id, taxonomy_cache = identity_cache)
      sample_accounting[[s_id]] <- compute_kraken_sample_accounting(kr, br, u, e, assignments, s_id)
    }
  }

  # 7. Perform outer join across full tax paths
  unclass_tax <- abundance_dfs[[1]]$tax[1]
  classified_taxa <- unique(unlist(lapply(abundance_dfs, function(df) df$tax[-1])))
  classified_taxa <- sort(classified_taxa, method = "radix")
  all_taxa <- c(unclass_tax, classified_taxa)

  merged_counts <- matrix(0.0, nrow = length(all_taxa), ncol = length(sample_ids),
                          dimnames = list(all_taxa, sample_ids))

  for (s_id in sample_ids) {
    df <- abundance_dfs[[s_id]]
    merged_counts[df$tax[1], s_id] <- as.numeric(df[1, s_id])
    if (nrow(df) > 1L) {
      tax_subset <- df$tax[-1]
      counts_subset <- as.numeric(df[-1, s_id])
      merged_counts[tax_subset, s_id] <- counts_subset
    }
  }

  row_totals <- rowSums(merged_counts)
  if (any(!is.finite(row_totals) | row_totals > 9007199254740991)) {
    stop("E_COHORT_COUNT_OVERFLOW: merged aggregate exceeds max safe integer.", call. = FALSE)
  }

  # Format abundance table data frame
  merged_df <- data.frame(
    tax = all_taxa,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (s_id in sample_ids) {
    merged_df[[s_id]] <- sprintf("%.1f", merged_counts[, s_id])
  }
  merged_df[["total"]] <- sprintf("%.1f", row_totals)

  # 8. Build metadata table (neutral placeholder)
  metadata_df <- data.frame(
    SampleID = sample_ids,
    Group = rep("Unspecified", length(sample_ids)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # 9. Build params_by_sample map
  params_by_sample_map <- list()
  for (idx in seq_along(sample_ids)) {
    params_by_sample_map[[sample_ids[idx]]] <- params_paths[idx]
  }

  # Determine ITS bounds and cache naming
  db_name <- tolower(rep_contract$database_set %||% "")
  is_its <- grepl("its", db_name)
  qc_min_len <- as.integer(rep_contract$min_len)
  qc_max_len <- as.integer(rep_contract$max_len)
  cache_basename <- if (!is.null(taxonomy_cache) && nzchar(trimws(taxonomy_cache))) {
    basename(taxonomy_cache)
  } else {
    if (is_its) "taxonomy_cache_its.json" else "taxonomy_cache.json"
  }
  if (cache_basename %in% c("config.yml", "metadata.tsv", "abundance_table_species.tsv", "cohort_merge_manifest.json")) {
    stop("E_COHORT_CACHE_INVALID: cache filename collides with bundle artifacts.", call. = FALSE)
  }

  # 10. Build cohort config
  cohort_cfg <- list(
    schema_version = 1L,
    project_name = project_name,
    mode = "cohort",
    seed = 42L,
    input = list(
      abundance_table = "abundance_table_species.tsv",
      cohort_merge_manifest = "cohort_merge_manifest.json",
      metadata = "metadata.tsv",
      params_json = params_paths[1],
      params_by_sample = params_by_sample_map,
      tax_column = "tax",
      aggregate_columns = list("total"),
      include_samples = NULL,
      assignments = if (length(assignments_map)) assignments_map else NULL,
      kraken2_reports = if (is_kraken) kraken_reports_map else NULL,
      bracken_reports = if (is_kraken) bracken_reports_map else NULL,
      phylogenetic_tree = NULL,
      phylogenetic_tip_map = NULL
    ),
    output = list(
      base_dir = "output"
    ),
    qc = list(
      display_min_length = qc_min_len,
      display_max_length = qc_max_len,
      target_min_length = NULL,
      target_max_length = NULL
    ),
    alpha = list(
      rarefaction_points = 25L,
      resample_depth = 50000L,
      resample_fraction_cap = 0.90,
      resample_iterations = 100L,
      hill_orders = list(0L, 1L, 2L),
      renyi_orders = list(1L)
    ),
    composition = list(
      top_n_taxa = 15L,
      stacked_bar_ranks = list("phylum", "family", "genus"),
      stacked_bar_max_taxa = 15L,
      stacked_bar_min_taxa = 5L,
      stacked_bar_min_mean_relative = 0.005,
      heatmap_ranks = list("phylum", "family", "genus"),
      heatmap_top_n_taxa = 10L,
      heatmap_include_other = TRUE,
      heatmap_transform = "log10_relative"
    ),
    faprotax = list(
      top_n_functions = 20L
    ),
    krona = list(
      enabled = FALSE,
      render_html = TRUE,
      html_renderer = "builtin",
      executable = "ktImportText"
    ),
    pavian = list(
      enabled = FALSE,
      render_html = TRUE,
      sankey = list(
        enabled = TRUE,
        render_html = TRUE,
        ranks = list("D", "K", "P", "C", "O", "F", "G", "S"),
        max_taxa_per_rank = 10L
      )
    ),
    beta = list(
      distances = list("bray", "jaccard"),
      primary_distance = "bray",
      permutations = 999L,
      minimum_count = 1L,
      strata_column = NULL,
      resampling = list(
        enabled = FALSE,
        iterations = 100L,
        depth_fraction_of_minimum = 0.75,
        minimum_success_fraction = 0.50
      )
    ),
    shared_taxa = list(
      rank = "species",
      minimum_count = 1L,
      group_prevalence = 0.5
    ),
    taxonomy = list(
      cache = cache_basename,
      network_mode = "cache_only",
      unresolved_policy = "warn",
      email_env = "NCBI_EMAIL",
      api_key_env = "NCBI_API_KEY"
    )
  )

  # 11. If validate_only, return validation summary without touching disk
  assert_merge_sources()
  if (isTRUE(validate_only)) {
    cat(sprintf("[merge_wf16s_cohort] Validation successful: %d samples, %d unique taxa joined.\n",
                length(sample_ids), length(all_taxa)))
    return(invisible(list(
      status = "validated",
      samples = sample_ids,
      taxa_count = length(all_taxa),
      classifier = rep_contract$classifier,
      database_set = rep_contract$database_set
    )))
  }

  # 12. Atomic staging and publication
  norm_dest <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (dir.exists(norm_dest) && !isTRUE(overwrite)) {
    stop(sprintf("E_COHORT_DESTINATION_EXISTS: destination directory '%s' already exists; pass --overwrite to replace.",
                 output_dir), call. = FALSE)
  }

  parent_dir <- dirname(norm_dest)
  if (!dir.exists(parent_dir)) {
    dir.create(parent_dir, recursive = TRUE)
  }

  stage_id <- sprintf(".__stage_%s_%s_%05d", basename(norm_dest), format(Sys.time(), "%Y%m%d_%H%M%S"), sample(10000:99999, 1))
  stage_dir <- file.path(parent_dir, stage_id)
  dir.create(stage_dir, recursive = TRUE)

  cleanup_needed <- TRUE
  on.exit({
    if (cleanup_needed && dir.exists(stage_dir)) {
      unlink(stage_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  # Write abundance_table_species.tsv
  out_ab_path <- file.path(stage_dir, "abundance_table_species.tsv")
  utils::write.table(merged_df, out_ab_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  # Write metadata.tsv
  out_meta_path <- file.path(stage_dir, "metadata.tsv")
  utils::write.table(metadata_df, out_meta_path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  # Write config.yml
  out_cfg_path <- file.path(stage_dir, "config.yml")
  cfg_header <- c(
    "# =============================================================================",
    "# Merged wf-16s Cohort Configuration (Generated by merge_wf16s_cohort.R)",
    "# NOTE: Cohort descriptive and experimental metadata unconfirmed.",
    "# ============================================================================="
  )
  cfg_yaml <- yaml::as.yaml(cohort_cfg)
  writeLines(c(cfg_header, cfg_yaml), out_cfg_path)

  # Write or copy taxonomy cache
  out_cache_path <- file.path(stage_dir, cache_basename)
  if (!is.null(taxonomy_cache) && file.exists(taxonomy_cache)) {
    file.copy(taxonomy_cache, out_cache_path, overwrite = TRUE)
  } else {
    writeLines("{}\n", out_cache_path)
  }

  # Build source records for manifest
  sources_records <- list()
  for (idx in seq_along(sample_ids)) {
    s_id <- sample_ids[idx]
    s_rec <- list(
      sample_id = s_id,
      output_root = resolved_roots[idx],
      abundance_table = list(path = abundance_paths[idx], sha256 = compute_file_hash(abundance_paths[idx])),
      params_json = list(path = params_paths[idx], sha256 = compute_file_hash(params_paths[idx]))
    )
    if (!is.null(assignments_map[[s_id]])) {
      s_rec$assignments <- list(path = assignments_map[[s_id]], sha256 = compute_file_hash(assignments_map[[s_id]]))
    }
    if (is_kraken) {
      s_rec$kraken2_report <- list(path = kraken_reports_map[[s_id]], sha256 = compute_file_hash(kraken_reports_map[[s_id]]))
      s_rec$bracken_report <- list(path = bracken_reports_map[[s_id]], sha256 = compute_file_hash(bracken_reports_map[[s_id]]))
    }
    sources_records[[s_id]] <- s_rec
  }

  output_records <- list(
    "abundance_table_species.tsv" = list(
      size_bytes = as.numeric(file.info(out_ab_path)$size),
      sha256 = compute_file_hash(out_ab_path)
    ),
    "metadata.tsv" = list(
      size_bytes = as.numeric(file.info(out_meta_path)$size),
      sha256 = compute_file_hash(out_meta_path)
    ),
    "config.yml" = list(
      size_bytes = as.numeric(file.info(out_cfg_path)$size),
      sha256 = compute_file_hash(out_cfg_path)
    ),
    taxonomy_cache_record = list(
      size_bytes = as.numeric(file.info(out_cache_path)$size),
      sha256 = compute_file_hash(out_cache_path)
    )
  )
  names(output_records)[length(output_records)] <- cache_basename

  # Write cohort_merge_manifest.json
  manifest_data <- list(
    manifest_version = "1.0",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    pipeline_version = pipeline_version,
    project_name = project_name,
    compatibility = list(
      classifier = rep_contract$classifier,
      database_set = rep_contract$database_set,
      taxonomic_rank = rep_contract$taxonomic_rank,
      upstream_contract = rep_contract,
      status = "compatible"
    ),
    samples = sample_ids,
    sample_accounting = sample_accounting,
    total_taxa = length(all_taxa),
    sources = unname(sources_records),
    outputs = output_records
  )
  out_manifest_path <- file.path(stage_dir, "cohort_merge_manifest.json")
  jsonlite::write_json(manifest_data, out_manifest_path, auto_unbox = TRUE, pretty = TRUE)

  # Publish atomically
  assert_merge_sources()
  backup_dir <- NULL
  if (dir.exists(norm_dest)) {
    backup_id <- sprintf(".__backup_%s_%s_%05d", basename(norm_dest), format(Sys.time(), "%Y%m%d_%H%M%S"), sample(10000:99999, 1))
    backup_dir <- file.path(parent_dir, backup_id)
    if (!file.rename(norm_dest, backup_dir)) {
      stop(sprintf("E_COHORT_PUBLISH_FAILED: could not stage backup for existing '%s'.", norm_dest), call. = FALSE)
    }
  }

  if (!file.rename(stage_dir, norm_dest)) {
    if (!is.null(backup_dir) && dir.exists(backup_dir)) {
      file.rename(backup_dir, norm_dest)
    }
    stop(sprintf("E_COHORT_PUBLISH_FAILED: could not publish merged cohort to '%s'.", norm_dest), call. = FALSE)
  }

  if (!is.null(backup_dir) && dir.exists(backup_dir)) {
    unlink(backup_dir, recursive = TRUE, force = TRUE)
  }

  cleanup_needed <- FALSE
  cat(sprintf("[merge_wf16s_cohort] Successfully merged %d samples (%d taxa) into '%s'.\n",
              length(sample_ids), length(all_taxa), norm_dest))

  invisible(list(
    status = "merged",
    output_dir = norm_dest,
    samples = sample_ids,
    taxa_count = length(all_taxa)
  ))
}

is_direct_merge_cli_call <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0L) return(FALSE)
  script_name <- basename(sub("^--file=", "", file_arg[1]))
  identical(script_name, "merge_wf16s_cohort.R")
}

# CLI Entrypoint
if (is_direct_merge_cli_call()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) {
    args <- c("--help")
  }
  opts <- parse_merge_cli_args(args)
  res <- merge_wf16s_cohort(
    input_dirs = opts$inputs,
    output_dir = opts$output_dir,
    project_name = opts$project_name,
    taxonomy_cache = opts$taxonomy_cache,
    validate_only = opts$validate_only,
    overwrite = opts$overwrite
  )
}
