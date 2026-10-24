# Read-only independent comparison of the legacy microeco FAPROTAX and
# official FAPROTAX-mapper 2 taxonomy products from one completed run.
# Per-sample, feature, and function evidence belongs only under output/.

compare_faprotax_backends <- function(run_root, private_evidence_dir) {
  fail <- function(message) stop(message, call. = FALSE)
  same <- function(actual, expected, message) {
    if (!isTRUE(identical(actual, expected))) fail(message)
  }
  close_enough <- function(actual, expected, message, tolerance = 0) {
    if (length(actual) != length(expected) ||
        !isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tolerance,
                          check.attributes = FALSE))) fail(message)
  }
  strings <- function(x, label, allow_empty = FALSE) {
    if (is.null(x)) return(character())
    if (is.character(x)) out <- unname(x) else if (is.list(x) && is.null(names(x))) {
      out <- vapply(x, function(item) {
        if (!is.character(item) || length(item) != 1L || is.na(item)) {
          fail(paste("Invalid string array item:", label))
        }
        item
      }, character(1))
    } else fail(paste("Expected JSON string array:", label))
    if (!allow_empty && (!length(out) || anyNA(out) || any(!nzchar(out)))) {
      fail(paste("Empty or missing string array:", label))
    }
    unname(out)
  }
  sha256 <- function(path) digest::digest(file = path, algo = "sha256")
  norm <- function(path) normalizePath(path, winslash = "/", mustWork = TRUE)
  read_tsv <- function(path, label) {
    if (!file.exists(path) || dir.exists(path)) fail(paste("Missing file:", label))
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = "", fileEncoding = "UTF-8",
               na.strings = "NA")
  }
  to_matrix <- function(tab, id_col, label) {
    if (ncol(tab) < 2L || !identical(names(tab)[1L], id_col) ||
        anyNA(tab[[1L]]) || any(!nzchar(as.character(tab[[1L]]))) ||
        anyDuplicated(tab[[1L]]) || anyDuplicated(names(tab))) {
      fail(paste("Invalid matrix table:", label))
    }
    values <- tab[, -1L, drop = FALSE]
    if (!all(vapply(values, function(x) is.numeric(x) ||
                    (is.logical(x) && all(is.na(x))), logical(1)))) {
      fail(paste("Nonnumeric matrix column:", label))
    }
    out <- as.matrix(values)
    storage.mode(out) <- "numeric"
    rownames(out) <- as.character(tab[[1L]])
    out
  }
  read_manifest <- function(path) {
    if (!file.exists(path)) fail("Completed run_manifest.json is missing.")
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
  output_path <- function(module, filename, root, manifest) {
    paths <- strings(manifest$modules[[module]]$outputs,
                     paste0("modules.", module, ".outputs"), allow_empty = TRUE)
    paths <- paths[basename(gsub("\\\\", "/", paths)) == filename]
    if (length(paths) != 1L) fail(paste("Expected one", filename, "output from", module))
    path <- paths[[1L]]
    if (!grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]", path)) path <- file.path(root, path)
    path <- norm(path)
    if (!startsWith(tolower(path), paste0(tolower(root), "/"))) {
      fail("Module output escapes the supplied run root.")
    }
    path
  }
  verify_input_records <- function(inputs) {
    checked <- 0L
    walk <- function(x) {
      if (!is.list(x)) return(invisible(NULL))
      if (!is.null(x$path) && !is.null(x$sha256) && !is.null(x$size_bytes)) {
        path <- as.character(x$path)
        if (!file.exists(path) || dir.exists(path) ||
            !identical(sha256(path), as.character(x$sha256)) ||
            as.numeric(file.info(path)$size) != as.numeric(x$size_bytes)) {
          fail("An original input/resource differs from the recorded run manifest.")
        }
        checked <<- checked + 1L
        return(invisible(NULL))
      }
      lapply(x, walk)
      invisible(NULL)
    }
    walk(inputs)
    checked
  }

  run_root <- norm(run_root)
  # The CLI sets the working directory to the repository root. When sourced by
  # another verifier, callers likewise run from the repository root.
  repo_root <- norm(getwd())
  output_root <- norm(file.path(repo_root, "output"))
  private_candidate <- if (grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]", private_evidence_dir)) {
    normalizePath(private_evidence_dir, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(repo_root, private_evidence_dir), winslash = "/", mustWork = FALSE)
  }
  if (!startsWith(tolower(private_candidate), paste0(tolower(output_root), "/"))) {
    fail("Private comparison evidence must be written under the ignored repository output/ directory.")
  }
  if (identical(tolower(private_candidate), tolower(run_root)) ||
      startsWith(tolower(private_candidate), paste0(tolower(run_root), "/"))) {
    fail("Private comparison evidence cannot be written into or beneath the completed run root.")
  }
  dir.create(private_candidate, recursive = TRUE, showWarnings = FALSE)
  private_evidence_dir <- norm(private_candidate)
  if (!startsWith(tolower(private_evidence_dir), paste0(tolower(output_root), "/"))) {
    fail("Private comparison evidence must be written under the ignored repository output/ directory.")
  }
  if (identical(tolower(private_evidence_dir), tolower(run_root)) ||
      startsWith(tolower(private_evidence_dir), paste0(tolower(run_root), "/"))) {
    fail("Private comparison evidence cannot be written into or beneath the completed run root.")
  }

  manifest <- read_manifest(file.path(run_root, "run_manifest.json"))
  if (!identical(as.character(manifest$run_status), "completed")) {
    fail("Backend comparison requires a completed pipeline run.")
  }
  if (!identical(as.character(manifest$modules$faprotax$status), "completed") ||
      !identical(as.character(manifest$modules$faprotax2$status), "completed")) {
    fail("Both FAPROTAX modules must have completed; structured skips are not comparable.")
  }
  samples <- strings(manifest$samples, "manifest.samples")
  if (anyDuplicated(samples)) fail("Run manifest sample IDs are duplicated.")
  if (!requireNamespace("digest", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("yaml", quietly = TRUE)) {
    fail("Comparison requires the already provisioned digest, jsonlite, and yaml packages.")
  }

  # Independent raw FAPROTAX2 validation, then direct reads for cross-backend checks.
  helper_path <- file.path(repo_root, "tests", "verify_faprotax2_run.R")
  if (!file.exists(helper_path)) fail("Independent FAPROTAX2 oracle is missing.")
  source(helper_path, local = TRUE, encoding = "UTF-8")
  verify_faprotax2_run(run_root, manifest)
  input_records_checked <- verify_input_records(manifest$inputs)
  expected_basis <- as.character(manifest$count_basis)
  f1_prov <- jsonlite::fromJSON(output_path("faprotax", "faprotax_provenance.json", run_root, manifest),
                                simplifyVector = FALSE)
  f2_prov <- jsonlite::fromJSON(output_path("faprotax2", "faprotax2_provenance.json", run_root, manifest),
                                simplifyVector = FALSE)
  same(as.character(f1_prov$count_basis), expected_basis, "Legacy count basis differs from the run manifest.")
  same(as.character(f2_prov$count_basis), expected_basis, "FAPROTAX2 count basis differs from the run manifest.")
  same(as.character(f2_prov$mode), "taxonomy", "This comparison is only defined for FAPROTAX2 taxonomy mode.")
  same(as.character(f2_prov$normalize), "none", "FAPROTAX2 run did not use normalize=none.")
  same(isTRUE(f2_prov$probability_weighting), FALSE, "FAPROTAX2 probability weighting was enabled.")

  config <- yaml::read_yaml(file.path(run_root, "resolved_config.yml"))
  tax_column <- as.character(config$input$tax_column)
  abundance_path <- as.character(manifest$inputs$abundance_table$path)
  source_table <- read_tsv(abundance_path, "original abundance table")
  if (!tax_column %in% names(source_table) || !all(samples %in% names(source_table))) {
    fail("Original abundance table does not contain the configured taxonomy/sample columns.")
  }
  tax_paths <- as.character(source_table[[tax_column]])
  if (anyNA(tax_paths) || any(!nzchar(tax_paths)) || anyDuplicated(tax_paths)) {
    fail("Original abundance table has invalid or duplicate TaxonPath values.")
  }
  source_counts <- as.matrix(source_table[, samples, drop = FALSE])
  storage.mode(source_counts) <- "numeric"
  if (anyNA(source_counts) || any(!is.finite(source_counts)) || any(source_counts < 0) ||
      any(abs(source_counts - round(source_counts)) > 1e-8)) {
    fail("Original abundance table does not contain valid integer counts for this count basis.")
  }
  rank_fields <- strsplit(paste0(tax_paths, ";"), ";", fixed = TRUE)
  if (any(lengths(rank_fields) != 8L)) fail("Original TaxonPath values do not contain eight ranks.")
  superkingdom <- vapply(rank_fields, `[[`, character(1), 1L)
  eligible <- superkingdom %in% c("Bacteria", "Archaea") & rowSums(source_counts) > 0
  eligible_indices <- which(eligible)
  eligible_ids <- sprintf("Taxon_%06d", eligible_indices)
  if (!length(eligible_ids) || anyDuplicated(eligible_ids)) fail("No unique eligible source taxon universe was recovered.")
  eligible_counts <- source_counts[eligible, , drop = FALSE]
  rownames(eligible_counts) <- eligible_ids

  # FAPROTAX2 mapper input must be the same ordered eligible rows and counts.
  f2_input <- read_tsv(output_path("faprotax2", "mapper_input.tsv", run_root, manifest), "FAPROTAX2 mapper input")
  if (!identical(names(f2_input), c("OTU", samples, "taxonomy")) ||
      !identical(as.character(f2_input$OTU), eligible_ids)) {
    fail("FAPROTAX2 mapper input IDs/order differ from the original eligible count universe.")
  }
  f2_input_counts <- as.matrix(f2_input[, samples, drop = FALSE])
  storage.mode(f2_input_counts) <- "numeric"
  if (!isTRUE(all.equal(unname(f2_input_counts), unname(eligible_counts), tolerance = 0))) {
    fail("FAPROTAX2 mapper input counts differ from the original abundance table.")
  }
  identity <- read_tsv(output_path("faprotax2", "taxon_identity.tsv", run_root, manifest), "FAPROTAX2 taxon identity")
  if (!identical(names(identity), c("FeatureID", "TaxonPath", "MapperTaxonomy")) ||
      !identical(as.character(identity$FeatureID), eligible_ids) ||
      !identical(as.character(identity$TaxonPath), tax_paths[eligible])) {
    fail("FAPROTAX2 identity crosswalk differs from original TaxonPath rows.")
  }
  expected_mapper_taxonomy <- vapply(rank_fields[eligible], function(fields) {
    mapped <- fields[c(1L, 3L:8L)]
    mapped[tolower(trimws(mapped)) %in% c("unknown", "unclassified")] <- ""
    paste(mapped, collapse = ";")
  }, character(1))
  same(as.character(identity$MapperTaxonomy), unname(expected_mapper_taxonomy),
       "FAPROTAX2 seven-rank taxonomy translation differs from the original input.")

  # Independent FAPROTAX1 reconstruction from its feature/function assignments.
  assignment <- read_tsv(output_path("faprotax", "faprotax_taxon_function_assignments.tsv", run_root, manifest),
                         "legacy taxon-function assignments")
  if (!identical(names(assignment), c("FeatureID", "Function", "TaxonPath")) ||
      anyNA(assignment) || any(!nzchar(as.matrix(assignment)))) {
    fail("Legacy taxon-function assignment schema or values are invalid.")
  }
  if (any(!assignment$FeatureID %in% eligible_ids) ||
      anyDuplicated(paste(assignment$FeatureID, assignment$Function, sep = "\r"))) {
    fail("Legacy assignments contain an unknown eligible feature or duplicate feature/function pair.")
  }
  assignment_index <- match(assignment$FeatureID, eligible_ids)
  if (nrow(assignment) && any(assignment$TaxonPath != tax_paths[eligible][assignment_index])) {
    fail("Legacy assignment TaxonPath does not identify its original abundance row.")
  }
  f1_labels <- sort(unique(as.character(assignment$Function)), method = "radix")
  f1_binary <- matrix(0, nrow = length(eligible_ids), ncol = length(f1_labels),
                      dimnames = list(eligible_ids, f1_labels))
  if (nrow(assignment)) {
    f1_binary[cbind(match(assignment$FeatureID, eligible_ids),
                    match(assignment$Function, f1_labels))] <- 1
  }
  f1_function_counts <- if (length(f1_labels)) t(f1_binary) %*% eligible_counts else
    matrix(numeric(), nrow = 0L, ncol = length(samples), dimnames = list(character(), samples))
  f1_function_table <- read_tsv(output_path("faprotax", "faprotax_function_abundance.tsv", run_root, manifest),
                                "legacy function abundance")
  required_f1_function <- c("SampleID", "Function", "FunctionReadCount",
                            "PctEligibleProkaryoticReads", "PctClassifiedReads", "PctTotalReads")
  if (!identical(names(f1_function_table), required_f1_function) ||
      anyDuplicated(paste(f1_function_table$SampleID, f1_function_table$Function, sep = "\r"))) {
    fail("Legacy function abundance does not match its output contract.")
  }
  f1_active <- sort(unique(as.character(f1_function_table$Function)), method = "radix")
  same(f1_active, f1_labels, "Legacy assignment labels and function abundance labels differ.")
  if (!setequal(f1_function_table$SampleID, samples) ||
      !setequal(unique(f1_function_table$Function), f1_labels) ||
      nrow(f1_function_table) != length(samples) * length(f1_labels)) {
    fail("Legacy function abundance omits or duplicates a sample/function pair.")
  }
  eligible_total <- colSums(eligible_counts)
  all_total <- colSums(source_counts)
  unclassified_index <- which(superkingdom == "Unclassified")
  if (length(unclassified_index) != 1L) fail("Original abundance table must have one canonical unclassified row.")
  classified_total <- all_total - source_counts[unclassified_index, ]
  f1_coverage <- read_tsv(output_path("faprotax", "faprotax_mapping_coverage.tsv", run_root, manifest),
                          "legacy mapping coverage")
  same(as.character(f1_coverage$SampleID), samples, "Legacy coverage sample order differs from the manifest.")
  for (field in c("TotalReads", "UnclassifiedReads", "ClassifiedReads", "EligibleProkaryoticReads")) {
    expected <- switch(field,
      TotalReads = all_total,
      UnclassifiedReads = source_counts[unclassified_index, ],
      ClassifiedReads = classified_total,
      EligibleProkaryoticReads = eligible_total)
    close_enough(f1_coverage[[field]], expected, paste("Legacy coverage disagrees with source counts:", field))
  }
  close_enough(f1_coverage$EligibleProkaryoticTaxa,
               colSums(eligible_counts > 0), "Legacy eligible taxon counts differ from the original table.")
  f1_mapped_rows <- rowSums(f1_binary) > 0
  f1_mapped_reads <- colSums(eligible_counts[f1_mapped_rows, , drop = FALSE])
  f1_mapped_taxa <- colSums(eligible_counts[f1_mapped_rows, , drop = FALSE] > 0)
  close_enough(f1_coverage$FunctionMappedReads, f1_mapped_reads,
               "Legacy mapped read arithmetic differs from reconstructed binary assignments.")
  close_enough(f1_coverage$FunctionUnmappedEligibleReads, eligible_total - f1_mapped_reads,
               "Legacy unmapped eligible read arithmetic is inconsistent.")
  close_enough(f1_coverage$FunctionMappedTaxa, f1_mapped_taxa,
               "Legacy mapped taxon arithmetic differs from reconstructed binary assignments.")
  for (i in seq_len(nrow(f1_function_table))) {
    sample <- as.character(f1_function_table$SampleID[i])
    fun <- as.character(f1_function_table$Function[i])
    expected <- f1_function_counts[fun, sample]
    close_enough(f1_function_table$FunctionReadCount[i], expected,
                 "Legacy function abundance differs from assignments multiplied by source counts.")
    pct <- function(n, d) if (d > 0) 100 * n / d else 0
    for (field in c("PctEligibleProkaryoticReads", "PctClassifiedReads", "PctTotalReads")) {
      denominator <- switch(field, PctEligibleProkaryoticReads = eligible_total[[sample]],
                            PctClassifiedReads = classified_total[[sample]], PctTotalReads = all_total[[sample]])
      close_enough(f1_function_table[[field]][i], pct(expected, denominator),
                   paste("Legacy percentage arithmetic differs for", field), tolerance = 1e-12)
    }
  }

  # Independently recompute FAPROTAX2 state products and compare the long summary.
  states_table <- read_tsv(output_path("faprotax2", "mapper_feature_states.tsv", run_root, manifest),
                           "FAPROTAX2 feature states")
  states <- to_matrix(states_table, "OTU", "FAPROTAX2 feature states")
  if (!identical(rownames(states), eligible_ids)) {
    fail("FAPROTAX2 state matrix feature order differs from its input.")
  }
  if (any(!is.na(states) & !states %in% c(0, 1))) fail("FAPROTAX2 states are not 0, 1, or NA.")
  states_zeroed <- states
  states_zeroed[is.na(states_zeroed)] <- 0
  f2_function_counts <- t(states_zeroed) %*% eligible_counts
  f2_function_table <- read_tsv(output_path("faprotax2", "faprotax2_function_abundance.tsv", run_root, manifest),
                                "FAPROTAX2 function abundance")
  if (!identical(names(f2_function_table), c("Function", "SampleID", "EstimatedFunctionAbundance", "CountBasis", "Unit")) ||
      anyDuplicated(paste(f2_function_table$Function, f2_function_table$SampleID, sep = "\r"))) {
    fail("FAPROTAX2 function abundance does not match its output contract.")
  }
  same(sort(unique(as.character(f2_function_table$SampleID)), method = "radix"),
       sort(samples, method = "radix"), "FAPROTAX2 summary samples differ from the manifest.")
  same(as.character(f2_function_table$CountBasis), rep(expected_basis, nrow(f2_function_table)),
       "FAPROTAX2 summary count basis differs from the manifest.")
  f2_labels <- colnames(states)
  f2_active <- sort(rownames(f2_function_counts)[rowSums(f2_function_counts) > 0], method = "radix")
  if (!setequal(f2_function_table$Function, f2_labels) ||
      nrow(f2_function_table) != length(f2_labels) * length(samples)) {
    fail("FAPROTAX2 summary does not include every catalog function/sample pair.")
  }
  for (i in seq_len(nrow(f2_function_table))) {
    sample <- as.character(f2_function_table$SampleID[i])
    fun <- as.character(f2_function_table$Function[i])
    expected <- f2_function_counts[fun, sample]
    close_enough(f2_function_table$EstimatedFunctionAbundance[i], expected,
                 "FAPROTAX2 function abundance differs from the independent unweighted state/count product.")
  }
  f2_predicted <- rowSums(states == 1, na.rm = TRUE) > 0
  f2_unknown_without_presence <- !f2_predicted & rowSums(is.na(states)) > 0
  f2_known_absence_only <- !f2_predicted & !f2_unknown_without_presence
  f2_mapped_reads <- colSums(eligible_counts[f2_predicted, , drop = FALSE])
  f2_mapped_taxa <- colSums(eligible_counts[f2_predicted, , drop = FALSE] > 0)
  f2_coverage <- read_tsv(output_path("faprotax2", "faprotax2_coverage.tsv", run_root, manifest),
                          "FAPROTAX2 coverage")
  same(as.character(f2_coverage$SampleID), samples, "FAPROTAX2 coverage sample order differs from the manifest.")
  close_enough(f2_coverage$EligibleTaxonAbundance, eligible_total,
               "FAPROTAX2 eligible abundance differs from the original count table.")
  close_enough(f2_coverage$TaxaWithPredictedPresence, f2_mapped_taxa,
               "FAPROTAX2 predicted taxon counts differ from feature states.")
  close_enough(f2_coverage$PredictedPresenceAbundance, f2_mapped_reads,
               "FAPROTAX2 predicted presence abundance differs from feature states.")
  close_enough(f2_coverage$UnknownWithoutPresenceAbundance,
               colSums(eligible_counts[f2_unknown_without_presence, , drop = FALSE]),
               "FAPROTAX2 unknown-without-presence abundance differs from feature states.")
  close_enough(f2_coverage$KnownAbsenceOnlyAbundance,
               colSums(eligible_counts[f2_known_absence_only, , drop = FALSE]),
               "FAPROTAX2 known-absence abundance differs from feature states.")

  # The F1 per-run table exposes only active labels; F2's raw state matrix has its
  # full 101-column configured catalog. Call the set difference descriptive.
  label_union <- sort(union(f1_labels, f2_labels), method = "radix")
  catalog <- data.frame(
    Function = label_union,
    LegacyObservedAssignment = label_union %in% f1_labels,
    FAPROTAX2Catalog = label_union %in% f2_labels,
    LegacyPositiveAbundance = label_union %in% f1_active,
    FAPROTAX2PositiveAbundance = label_union %in% f2_active,
    stringsAsFactors = FALSE
  )
  catalog$ObservedComparison <- ifelse(catalog$LegacyPositiveAbundance & catalog$FAPROTAX2PositiveAbundance,
    "positive in both", ifelse(catalog$LegacyPositiveAbundance, "legacy only in this run",
    ifelse(catalog$FAPROTAX2PositiveAbundance, "FAPROTAX2 only in this run", "no positive abundance")))
  utils::write.table(catalog, file.path(private_evidence_dir, "catalog_comparison.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

  function_rows <- expand.grid(Function = label_union, SampleID = samples,
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  function_rows$LegacyFunctionReadCount <- NA_real_
  function_rows$FAPROTAX2FunctionAbundance <- NA_real_
  for (i in seq_len(nrow(function_rows))) {
    fun <- function_rows$Function[i]; sample <- function_rows$SampleID[i]
    if (fun %in% f1_labels) function_rows$LegacyFunctionReadCount[i] <- f1_function_counts[fun, sample]
    if (fun %in% f2_labels) function_rows$FAPROTAX2FunctionAbundance[i] <- f2_function_counts[fun, sample]
  }
  utils::write.table(function_rows, file.path(private_evidence_dir, "function_sample_counts.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

  f1_unknown <- !f1_mapped_rows
  sample_summary <- data.frame(
    SampleID = samples,
    EligibleTaxonAbundance = as.numeric(eligible_total[samples]),
    EligibleTaxa = as.numeric(colSums(eligible_counts > 0)[samples]),
    LegacyPositiveTaxa = as.numeric(f1_mapped_taxa[samples]),
    LegacyNoPositiveFunctionTaxa = as.numeric(colSums(eligible_counts[f1_unknown, , drop = FALSE] > 0)[samples]),
    LegacyPositiveFunctionAbundance = as.numeric(f1_mapped_reads[samples]),
    FAPROTAX2PositiveTaxa = as.numeric(f2_mapped_taxa[samples]),
    FAPROTAX2NoPresenceTaxa = as.numeric(colSums(eligible_counts[!f2_predicted, , drop = FALSE] > 0)[samples]),
    FAPROTAX2UnknownWithoutPresenceTaxa = as.numeric(colSums(eligible_counts[f2_unknown_without_presence, , drop = FALSE] > 0)[samples]),
    FAPROTAX2KnownAbsenceOnlyTaxa = as.numeric(colSums(eligible_counts[f2_known_absence_only, , drop = FALSE] > 0)[samples]),
    FAPROTAX2PositiveFunctionAbundance = as.numeric(f2_mapped_reads[samples]),
    LegacySummedFunctionAbundance = as.numeric(colSums(f1_function_counts)[samples]),
    FAPROTAX2SummedFunctionAbundance = as.numeric(colSums(f2_function_counts)[samples]),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  sample_summary$LegacyFunctionSumExceedsEligible <-
    sample_summary$LegacySummedFunctionAbundance > sample_summary$EligibleTaxonAbundance
  sample_summary$FAPROTAX2FunctionSumExceedsEligible <-
    sample_summary$FAPROTAX2SummedFunctionAbundance > sample_summary$EligibleTaxonAbundance
  utils::write.table(sample_summary, file.path(private_evidence_dir, "sample_metrics.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

  feature_evidence <- data.frame(
    FeatureID = eligible_ids, TaxonPath = tax_paths[eligible],
    LegacyPositiveFunction = rowSums(f1_binary) > 0,
    FAPROTAX2PositiveFunction = f2_predicted,
    FAPROTAX2UnknownWithoutPresence = f2_unknown_without_presence,
    FAPROTAX2KnownAbsenceOnly = f2_known_absence_only,
    stringsAsFactors = FALSE
  )
  feature_evidence <- cbind(feature_evidence, as.data.frame(eligible_counts, check.names = FALSE))
  utils::write.table(feature_evidence, file.path(private_evidence_dir, "feature_taxon_evidence.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")

  result <- list(
    verified = TRUE, mode = "taxonomy", count_basis = expected_basis,
    samples = length(samples), eligible_taxa = length(eligible_ids),
    eligible_abundance = as.numeric(eligible_total),
    legacy_observed_labels = f1_labels, faprotax2_catalog_labels = f2_labels,
    legacy_active_labels = f1_active, faprotax2_active_labels = f2_active,
    shared_active_labels = intersect(f1_active, f2_active),
    f2_catalog_not_observed_in_legacy = setdiff(f2_labels, f1_labels),
    legacy_observed_not_in_f2_catalog = setdiff(f1_labels, f2_labels),
    input_records_verified = input_records_checked,
    probabilities_weighted = FALSE, normalized = FALSE,
    private_evidence_dir = private_evidence_dir
  )
  class(result) <- c("faprotax_backend_comparison", "list")
  invisible(result)
}

# Standalone CLI: Rscript tests/compare_faprotax_backends.R RUN_ROOT PRIVATE_OUTPUT_DIR
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
  if (identical(basename(script_path), "compare_faprotax_backends.R")) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 2L) {
      stop("Usage: Rscript tests/compare_faprotax_backends.R RUN_ROOT PRIVATE_OUTPUT_DIR", call. = FALSE)
    }
    repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
    setwd(repo_root)
    result <- compare_faprotax_backends(args[[1L]], args[[2L]])
    cat(sprintf(paste0("Verified taxonomy comparison: %d samples, %d eligible taxa, ",
                       "%d shared active labels; %d input/resource fingerprints unchanged; details under output/.\n"),
                result$samples, result$eligible_taxa, length(result$shared_active_labels),
                result$input_records_verified))
  }
}
