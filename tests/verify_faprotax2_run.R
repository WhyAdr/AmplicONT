# Independent raw-output oracle for completed FAPROTAX2 taxonomy runs.
# This file intentionally does not source analysis/09_faprotax2.R or its validators.

verify_faprotax2_run <- function(root, manifest = NULL) {
  fail <- function(message) stop(message, call. = FALSE)
  assert <- function(condition, message) {
    if (!isTRUE(condition)) fail(message)
  }
  as_strings <- function(value, label, allow_null = FALSE) {
    if (is.null(value)) {
      if (allow_null) return(character())
      fail(paste("Missing string array:", label))
    }
    if (is.character(value)) return(unname(value))
    if (!is.list(value) || !is.null(names(value))) fail(paste("Expected JSON string array:", label))
    result <- vapply(value, function(item) {
      if (!is.character(item) || length(item) != 1L || is.na(item)) {
        fail(paste("Invalid string array item:", label))
      }
      item
    }, character(1))
    unname(result)
  }
  sha256 <- function(path) digest::digest(file = path, algo = "sha256")
  normalize_path <- function(path, must_work = TRUE) {
    normalizePath(path, winslash = "/", mustWork = must_work)
  }
  read_tsv <- function(path, label) {
    if (!file.exists(path) || dir.exists(path)) fail(paste("Missing regular file:", label, path))
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
               quote = "", comment.char = "", fileEncoding = "UTF-8",
               na.strings = "NA")
  }
  read_matrix <- function(path, key, label) {
    table <- read_tsv(path, label)
    if (ncol(table) < 2L || !identical(names(table)[1L], key) ||
        anyDuplicated(names(table)) || anyNA(table[[1L]]) ||
        any(!nzchar(table[[1L]])) || anyDuplicated(table[[1L]])) {
      fail(paste("Invalid matrix schema or identifiers:", label))
    }
    values <- table[, -1L, drop = FALSE]
    valid <- vapply(values, function(column) {
      is.numeric(column) || (is.logical(column) && all(is.na(column)))
    }, logical(1))
    if (!all(valid)) fail(paste("Nonnumeric matrix column:", label))
    matrix <- as.matrix(values)
    storage.mode(matrix) <- "numeric"
    rownames(matrix) <- as.character(table[[1L]])
    matrix
  }
  json_record_path <- function(path, output_root) {
    normalized_root <- normalize_path(output_root)
    if (grepl("^[A-Za-z]:[/\\\\]|^[/\\\\]", path)) {
      normalized <- normalize_path(path)
      prefix <- paste0(tolower(normalized_root), "/")
      if (!startsWith(tolower(normalized), prefix)) fail(paste("Manifest output is outside its root:", path))
      substring(normalized, nchar(normalized_root) + 2L)
    } else {
      gsub("\\\\", "/", path)
    }
  }
  verify_ownership <- function(module_record, run_manifest, output_root, expected_status) {
    module_outputs <- as_strings(module_record$outputs, "modules.faprotax2.outputs", allow_null = TRUE)
    owned <- tolower(gsub("\\\\", "/", as_strings(run_manifest$owned_outputs, "owned_outputs")))
    artifacts <- run_manifest$artifacts
    if (!is.list(artifacts)) fail("Manifest artifacts must be an array.")
    artifact_path <- vapply(artifacts, function(item) as.character(item$relative_path), character(1))
    producer <- vapply(artifacts, function(item) {
      if (is.null(item$producer_module)) "" else as.character(item$producer_module)
    }, character(1))
    module_relative <- vapply(module_outputs, json_record_path, character(1),
                              output_root = output_root)
    if (any(!tolower(gsub("\\\\", "/", module_relative)) %in% owned)) {
      fail("An FAPROTAX2 module output is absent from owned_outputs.")
    }
    artifact_index <- match(tolower(gsub("\\\\", "/", module_relative)),
                            tolower(gsub("\\\\", "/", artifact_path)))
    if (anyNA(artifact_index) || any(producer[artifact_index] != "faprotax2")) {
      fail("Every FAPROTAX2 module output must have a faprotax2 producer artifact.")
    }
    prefixed <- startsWith(tolower(gsub("\\\\", "/", artifact_path)), "09_faprotax2/")
    if (any(producer[prefixed] != "faprotax2")) {
      fail("A 09_FAPROTAX2 artifact is not owned by faprotax2.")
    }
    invisible(list(outputs = module_outputs, relative = module_relative))
  }
  verify_resources <- function(run_manifest, provenance = NULL) {
    records <- run_manifest$inputs$faprotax2_resources
    expected <- c(
      mapper_script = "8bd6ff97aa53f17fb0122b98c436275c230e34629f7a83c446b2b3854e53a64f",
      "default_config.tsv" = "d9c597304f966a2b2ab96a5db8de7a14a894240d074c3ab5ddcad93c27584c37",
      "reference_alignments_sanitized_subset.fasta.gz" = "3a21836125cad9cc9bc0ee1cf1e62bda2a5cf98c54ea43723eeaff2f74dcbbe2",
      "reference_alignments_sanitized.fasta.gz" = "0f30fd59650c738494a901c2651a022b51d058c56488990c510f584400cca5b9",
      "reference2species_and_traits.tsv" = "a8b3dc971470c224a7f33874021622bf24270e695076740cb2f151dd7a151f96",
      "tree_model.txt" = "a3c3382d0c2d33c1a4662f3af93c28ba0fe4f91c5d49cd500930faac33b7a7b9",
      "tree.tre" = "1b40b2e3f4165b6eea49290064aadaff3c74b4744c16e91c35506afcfdfa977d"
    )
    if (!is.list(records) || !setequal(names(records), names(expected))) {
      fail("Manifest does not contain the exact audited FAPROTAX2 resource fingerprints.")
    }
    for (name in names(expected)) {
      record <- records[[name]]
      path <- as.character(record$path)
      if (!file.exists(path) || dir.exists(path) ||
          !identical(as.character(record$sha256), unname(expected[name])) ||
          !identical(sha256(path), unname(expected[name])) ||
          !is.numeric(record$size_bytes) ||
          as.numeric(file.info(path)$size) != as.numeric(record$size_bytes)) {
        fail(paste("FAPROTAX2 resource fingerprint changed:", name))
      }
      if (!is.null(provenance)) {
        provenance_record <- provenance$resources[[name]]
        if (is.null(provenance_record) ||
            !identical(as.character(provenance_record$path), path) ||
            !identical(as.character(provenance_record$sha256), unname(expected[name])) ||
            as.numeric(provenance_record$size_bytes) != as.numeric(record$size_bytes)) {
          fail(paste("FAPROTAX2 provenance resource differs from manifest:", name))
        }
      }
    }
    invisible(records)
  }

  root <- normalize_path(root)
  if (is.null(manifest)) {
    manifest_path <- file.path(root, "run_manifest.json")
    if (!file.exists(manifest_path)) fail("run_manifest.json is missing.")
    manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  }
  module <- manifest$modules$faprotax2
  if (is.null(module)) {
    return(invisible(list(status = "not_requested", verified = TRUE)))
  }
  status <- as.character(module$status)
  module_outputs <- as_strings(module$outputs, "modules.faprotax2.outputs", allow_null = TRUE)
  if (identical(status, "not_run")) {
    if (length(module_outputs) != 0L) fail("Unrequested FAPROTAX2 must not declare outputs.")
    return(invisible(list(status = "not_requested", verified = TRUE)))
  }
  if (identical(status, "skipped")) {
    if (!identical(as.character(module$reason), "No eligible prokaryotic taxa")) {
      fail("FAPROTAX2 skipped with an unexpected reason.")
    }
    owned <- verify_ownership(module, manifest, root, status)
    if (length(owned$relative) != 1L ||
        !identical(unname(tolower(gsub("\\\\", "/", owned$relative))),
                   "09_faprotax2/faprotax2_skipped.tsv")) {
      fail(paste("FAPROTAX2 no-prokaryote skip must own only faprotax2_skipped.tsv; found:",
                 paste(owned$relative, collapse = ", ")))
    }
    skip <- read_tsv(file.path(root, owned$relative[1L]), "FAPROTAX2 skip record")
    if (!identical(names(skip), c("Status", "Reason")) || nrow(skip) != 1L ||
        !identical(as.character(skip$Status), "skipped") ||
        !identical(as.character(skip$Reason), "No positive-count Bacteria or Archaea rows")) {
      fail("FAPROTAX2 skip output does not match the no-prokaryote contract.")
    }
    verify_resources(manifest)
    input_record <- manifest$inputs$abundance_table
    if (is.null(input_record$path) || is.null(input_record$sha256) ||
        !file.exists(as.character(input_record$path)) ||
        !identical(sha256(as.character(input_record$path)), as.character(input_record$sha256)) ||
        as.numeric(file.info(as.character(input_record$path))$size) != as.numeric(input_record$size_bytes)) {
      fail("Skipped FAPROTAX2 original abundance fingerprint changed.")
    }
    skipped_config <- yaml::read_yaml(file.path(root, "resolved_config.yml"))
    skipped_source <- read_tsv(as.character(input_record$path), "original abundance table")
    skipped_tax_column <- as.character(skipped_config$input$tax_column)
    skipped_samples <- as_strings(manifest$samples, "samples")
    if (!skipped_tax_column %in% names(skipped_source) || !all(skipped_samples %in% names(skipped_source))) {
      fail("Skipped FAPROTAX2 source table lacks taxonomy or sample fields.")
    }
    skipped_paths <- strsplit(as.character(skipped_source[[skipped_tax_column]]), ";", fixed = TRUE)
    skipped_counts <- as.matrix(skipped_source[, skipped_samples, drop = FALSE])
    storage.mode(skipped_counts) <- "numeric"
    if (any(lengths(skipped_paths) != 8L) || anyNA(skipped_counts) ||
        any(!is.finite(skipped_counts)) || any(skipped_counts < 0) ||
        any(vapply(skipped_paths, function(x) x[1L] %in% c("Bacteria", "Archaea"), logical(1)) &
              rowSums(skipped_counts) > 0)) {
      fail("Manifest claims a no-prokaryote skip, but eligible source rows exist or inputs are invalid.")
    }
    return(invisible(list(status = "skipped", verified = TRUE, reason = module$reason)))
  }
  if (!identical(status, "completed")) fail(paste("FAPROTAX2 module is not completed:", status))

  ownership <- verify_ownership(module, manifest, root, status)
  paths <- stats::setNames(file.path(root, ownership$relative), basename(ownership$relative))
  required_files <- c("mapper_input.tsv", "taxon_identity.tsv", "mapper_function_abundance.tsv",
    "mapper_feature_states.tsv", "mapper_annotation_probabilities.tsv", "mapper_report.txt",
    "mapper_stdout.txt", "mapper_stderr.txt", "faprotax2_function_abundance.tsv",
    "faprotax2_coverage.tsv", "faprotax2_feature_function_assignments.tsv",
    "faprotax2_top_functions.png", "faprotax2_provenance.json")
  if (!all(required_files %in% names(paths))) fail("Completed FAPROTAX2 outputs are incomplete.")
  if (!all(file.exists(paths[required_files]))) fail("A declared FAPROTAX2 output is missing.")

  provenance <- jsonlite::fromJSON(paths[["faprotax2_provenance.json"]], simplifyVector = FALSE)
  if (!identical(as.character(provenance$mapper), "FAPROTAX-mapper") ||
      !identical(as.character(provenance$mapper_version), "2.0 beta") ||
      !identical(as.character(provenance$mode), "taxonomy") ||
      !identical(provenance$mapper_modified, FALSE) ||
      !identical(as.character(provenance$archive_sha256),
                 "859b99133b55f1e4f637e6ac14a804f3321ed27f76ee03d95e3f215382206bac") ||
      !identical(as.character(provenance$normalize), "none") ||
      !identical(provenance$probability_weighting, FALSE) ||
      !identical(as.character(provenance$hsp_algorithm), "config") ||
      !identical(as.character(provenance$unknown_ref_states), "config")) {
    fail("FAPROTAX2 provenance does not describe the audited unweighted taxonomy run.")
  }
  invocation <- as_strings(provenance$arguments, "faprotax2_provenance.arguments")
  option_value <- function(flag) {
    index <- match(flag, invocation)
    if (is.na(index) || index == length(invocation)) return(NA_character_)
    invocation[[index + 1L]]
  }
  if (any(invocation == "--weigh_by_probability") ||
      !identical(option_value("--normalize"), "none") ||
      !identical(option_value("--hsp_algorithm"), "config") ||
      !identical(option_value("--unknown_ref_states"), "config") ||
      !identical(option_value("--Nthreads"), as.character(provenance$nthreads))) {
    fail("Recorded FAPROTAX2 invocation does not match unweighted normalize-none settings.")
  }
  if (!identical(as.character(provenance$count_basis), as.character(manifest$count_basis))) {
    fail("FAPROTAX2 count basis differs from the original run manifest.")
  }
  verify_resources(manifest, provenance)

  abundance_record <- manifest$inputs$abundance_table
  if (is.null(abundance_record$path) || is.null(abundance_record$sha256) ||
      !file.exists(as.character(abundance_record$path)) ||
      !identical(sha256(as.character(abundance_record$path)), as.character(abundance_record$sha256)) ||
      as.numeric(file.info(as.character(abundance_record$path))$size) != as.numeric(abundance_record$size_bytes)) {
    fail("Original abundance-table fingerprint does not match the run manifest.")
  }
  if (!identical(sha256(paths[["mapper_input.tsv"]]), as.character(provenance$input_sha256))) {
    fail("Mapper-input hash does not match FAPROTAX2 provenance.")
  }

  config <- yaml::read_yaml(file.path(root, "resolved_config.yml"))
  tax_column <- as.character(config$input$tax_column)
  sample_ids <- as_strings(manifest$samples, "samples")
  source <- read_tsv(as.character(abundance_record$path), "original abundance table")
  if (!tax_column %in% names(source) || !all(sample_ids %in% names(source))) {
    fail("Original abundance table lacks the configured taxonomy or sample columns.")
  }
  source_paths <- as.character(source[[tax_column]])
  if (anyNA(source_paths) || any(!nzchar(source_paths)) || anyDuplicated(source_paths)) {
    fail("Original abundance taxonomy paths must be present and unique.")
  }
  source_counts <- as.matrix(source[, sample_ids, drop = FALSE])
  storage.mode(source_counts) <- "numeric"
  if (anyNA(source_counts) || any(!is.finite(source_counts)) || any(source_counts < 0)) {
    fail("Original abundance counts are not finite nonnegative numbers.")
  }
  source_parts <- strsplit(source_paths, ";", fixed = TRUE)
  if (any(lengths(source_parts) != 8L)) fail("Original taxonomy rows do not contain exactly eight ranks.")
  prokaryotic <- vapply(source_parts, function(x) x[1L] %in% c("Bacteria", "Archaea"), logical(1))
  eligible <- which(prokaryotic & rowSums(source_counts) > 0)
  expected_ids <- sprintf("Taxon_%06d", eligible)
  expected_paths <- source_paths[eligible]
  expected_counts <- source_counts[eligible, , drop = FALSE]
  if (length(expected_ids) == 0L) fail("Completed FAPROTAX2 run has no eligible source features.")

  identity <- read_tsv(paths[["taxon_identity.tsv"]], "FAPROTAX2 source taxonomy mapping")
  if (!identical(names(identity), c("FeatureID", "TaxonPath", "MapperTaxonomy")) ||
      !identical(as.character(identity$FeatureID), expected_ids) ||
      !identical(as.character(identity$TaxonPath), expected_paths)) {
    fail("Original source TaxonPath values or feature ordering were not retained exactly.")
  }
  expected_mapper_taxonomy <- vapply(source_parts[eligible], function(parts) {
    selected <- parts[c(1L, 3L, 4L, 5L, 6L, 7L, 8L)]
    missing <- is.na(selected) | tolower(trimws(selected)) %in% c("unknown", "unclassified")
    selected[missing] <- ""
    paste(selected, collapse = ";")
  }, character(1))
  if (!identical(as.character(identity$MapperTaxonomy), unname(expected_mapper_taxonomy))) {
    fail("The eight-rank TaxonPath to seven-rank mapper taxonomy translation differs.")
  }

  mapper_input <- read_tsv(paths[["mapper_input.tsv"]], "raw mapper input")
  expected_input_names <- c("OTU", sample_ids, "taxonomy")
  if (!identical(names(mapper_input), expected_input_names) ||
      !identical(as.character(mapper_input$OTU), expected_ids) ||
      !identical(as.character(mapper_input$taxonomy), unname(expected_mapper_taxonomy))) {
    fail("Raw mapper input IDs, samples, or translated taxonomy differ from source input.")
  }
  mapper_counts <- as.matrix(mapper_input[, sample_ids, drop = FALSE])
  storage.mode(mapper_counts) <- "numeric"
  rownames(mapper_counts) <- as.character(mapper_input$OTU)
  if (!isTRUE(all.equal(unname(mapper_counts), unname(expected_counts), tolerance = 0))) {
    fail("Raw mapper counts differ from the eligible original abundance rows.")
  }

  states <- read_matrix(paths[["mapper_feature_states.tsv"]], "OTU", "raw feature-function states")
  probabilities <- read_matrix(paths[["mapper_annotation_probabilities.tsv"]], "OTU",
                               "raw annotation probabilities")
  abundance <- read_matrix(paths[["mapper_function_abundance.tsv"]], "function",
                           "raw function abundance")
  default_config <- read_tsv(as.character(manifest$inputs$faprotax2_resources[["default_config.tsv"]]$path),
                            "official default function config")
  if (!"function" %in% names(default_config) || anyDuplicated(default_config[["function"]])) {
    fail("Official default function catalog is invalid.")
  }
  functions <- as.character(default_config[["function"]])
  if (!identical(rownames(states), expected_ids) ||
      !identical(rownames(probabilities), expected_ids) ||
      !identical(colnames(probabilities), colnames(states)) ||
      !identical(rownames(abundance), functions) ||
      !identical(colnames(abundance), sample_ids) ||
      !identical(colnames(states), functions)) {
    fail("Raw mapper feature, sample, and function catalogs or orders disagree.")
  }
  if (any(!is.na(states) & !states %in% c(0, 1))) fail("Raw mapper states contain values outside 0/1/NA.")
  if (any(!is.na(probabilities) &
          (!is.finite(probabilities) | probabilities < 0 | probabilities > 1))) {
    fail("Raw annotation probabilities contain values outside [0,1].")
  }
  if (anyNA(abundance) || any(!is.finite(abundance)) || any(abundance < 0)) {
    fail("Raw function abundance contains missing, nonfinite, or negative values.")
  }
  if (any((states == 1 & probabilities < 0.50000001) |
          (states == 0 & probabilities > 0.49999999), na.rm = TRUE)) {
    fail("Observed 0/1 states disagree with the official configured probability thresholds.")
  }
  binary <- states
  binary[is.na(binary)] <- 0
  expected_abundance <- t(binary) %*% mapper_counts
  if (!isTRUE(all.equal(unname(abundance), unname(expected_abundance), tolerance = 0))) {
    fail("Raw function table is not the exact unnormalized state × abundance product.")
  }

  summary <- read_tsv(paths[["faprotax2_function_abundance.tsv"]], "long function abundance")
  summary_columns <- c("Function", "SampleID", "EstimatedFunctionAbundance", "CountBasis", "Unit")
  if (!identical(names(summary), summary_columns) ||
      anyDuplicated(paste(summary$Function, summary$SampleID, sep = "\r")) ||
      nrow(summary) != length(functions) * length(sample_ids) ||
      !setequal(as.character(summary$Function), functions) ||
      !setequal(as.character(summary$SampleID), sample_ids) ||
      any(as.character(summary$CountBasis) != as.character(manifest$count_basis)) ||
      any(as.character(summary$Unit) != "sum_of_input_taxon_abundances_with_predicted_presence")) {
    fail("Long function summary schema, catalog, or count-basis metadata is invalid.")
  }
  function_index <- match(as.character(summary$Function), functions)
  sample_index <- match(as.character(summary$SampleID), sample_ids)
  expected_long <- abundance[cbind(function_index, sample_index)]
  if (!isTRUE(all.equal(as.numeric(summary$EstimatedFunctionAbundance),
                        as.numeric(expected_long), tolerance = 0))) {
    fail("Long summary abundances differ from the independent raw-matrix result.")
  }

  assignments <- read_tsv(paths[["faprotax2_feature_function_assignments.tsv"]],
                          "feature-function assignment summary")
  assignment_columns <- c("FeatureID", "Function", "PredictedState", "AnnotationProbability")
  if (!identical(names(assignments), assignment_columns)) fail("Assignment summary has an unexpected schema.")
  positive_pairs <- which(states == 1, arr.ind = TRUE)
  expected_pair_key <- paste(rownames(states)[positive_pairs[, 1L]],
                             colnames(states)[positive_pairs[, 2L]], sep = "\r")
  actual_pair_key <- paste(as.character(assignments$FeatureID),
                           as.character(assignments$Function), sep = "\r")
  if (anyDuplicated(actual_pair_key) || !setequal(actual_pair_key, expected_pair_key) ||
      nrow(assignments) != nrow(positive_pairs) ||
      any(as.numeric(assignments$PredictedState) != 1)) {
    fail("Assignment rows do not equal the unique positive feature-function states.")
  }
  pair_row <- match(actual_pair_key, expected_pair_key)
  expected_probability <- probabilities[positive_pairs]
  if (!isTRUE(all.equal(as.numeric(assignments$AnnotationProbability),
                        as.numeric(expected_probability[pair_row]), tolerance = 0))) {
    fail("Assignment probabilities differ from the raw mapper probability matrix.")
  }

  coverage <- read_tsv(paths[["faprotax2_coverage.tsv"]], "taxonomy coverage summary")
  coverage_columns <- c("SampleID", "InputTaxonAbundance", "ExcludedInputAbundance",
    "EligibleTaxonAbundance", "PredictedPresenceAbundance", "NoPredictedPresenceAbundance",
    "UnknownWithoutPresenceAbundance", "KnownAbsenceOnlyAbundance", "EligibleTaxa",
    "TaxaWithPredictedPresence", "CountBasis")
  if (!identical(names(coverage), coverage_columns) ||
      !identical(as.character(coverage$SampleID), sample_ids) ||
      any(as.character(coverage$CountBasis) != as.character(manifest$count_basis))) {
    fail("Coverage summary schema, sample order, or count basis is invalid.")
  }
  predicted <- rowSums(states == 1, na.rm = TRUE) > 0
  unknown_without_presence <- !predicted & rowSums(is.na(states)) > 0
  known_absence_only <- !predicted & !unknown_without_presence
  expected_coverage <- data.frame(
    InputTaxonAbundance = colSums(source_counts),
    ExcludedInputAbundance = colSums(source_counts) - colSums(mapper_counts),
    EligibleTaxonAbundance = colSums(mapper_counts),
    PredictedPresenceAbundance = colSums(mapper_counts[predicted, , drop = FALSE]),
    NoPredictedPresenceAbundance = colSums(mapper_counts[!predicted, , drop = FALSE]),
    UnknownWithoutPresenceAbundance = colSums(mapper_counts[unknown_without_presence, , drop = FALSE]),
    KnownAbsenceOnlyAbundance = colSums(mapper_counts[known_absence_only, , drop = FALSE]),
    EligibleTaxa = colSums(mapper_counts > 0),
    TaxaWithPredictedPresence = colSums(mapper_counts[predicted, , drop = FALSE] > 0)
  )
  for (column in names(expected_coverage)) {
    if (!isTRUE(all.equal(as.numeric(coverage[[column]]),
                          as.numeric(expected_coverage[[column]]), tolerance = 0))) {
      fail(paste("Coverage arithmetic differs for", column))
    }
  }
  partition_total <- coverage$PredictedPresenceAbundance +
    coverage$UnknownWithoutPresenceAbundance + coverage$KnownAbsenceOnlyAbundance
  if (!isTRUE(all.equal(as.numeric(partition_total),
                        as.numeric(coverage$EligibleTaxonAbundance), tolerance = 0))) {
    fail("Presence/unknown/known-absence abundance categories do not sum to eligible abundance.")
  }

  if (!identical(as.character(provenance$taxonomy_translation),
      "NCBI superkingdom,phylum,class,order,family,genus,species; intervening kingdom omitted; Unknown/Unclassified/NA become empty fields; original TaxonPath retained") ||
      !grepl("no normalization applied", as.character(provenance$denominator_definition), fixed = TRUE)) {
    fail("Provenance taxonomy translation/denominator description is inconsistent.")
  }
  if (!grepl("Overlapping", as.character(provenance$interpretation), fixed = TRUE) ||
      !grepl("without probability weighting", as.character(provenance$interpretation), ignore.case = TRUE)) {
    fail("Provenance omits overlap or probability-weighting semantics.")
  }

  result <- list(
    status = "completed", verified = TRUE, features = nrow(states),
    functions = ncol(states), samples = sample_ids,
    state_counts = c(present = sum(states == 1, na.rm = TRUE),
                     absent = sum(states == 0, na.rm = TRUE), unknown = sum(is.na(states))),
    input_abundance = colSums(source_counts),
    eligible_abundance = colSums(mapper_counts),
    summed_function_abundance = colSums(abundance),
    overlapping_function_total_exceeds_eligible =
      any(colSums(abundance) > colSums(mapper_counts))
  )
  invisible(result)
}

# Direct Rscript entry point. Sourcing this helper from another verifier has no side effects.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE)
  if (identical(basename(script_path), "verify_faprotax2_run.R")) {
    script_args <- commandArgs(trailingOnly = TRUE)
    if (length(script_args) != 1L) {
      stop("Usage: Rscript tests/verify_faprotax2_run.R OUTPUT_DIR", call. = FALSE)
    }
    result <- verify_faprotax2_run(script_args[[1L]])
    if (is.null(result)) cat("FAPROTAX2: no-op\n") else {
      cat(sprintf("FAPROTAX2 %s verified", result$status))
      if (!is.null(result$features)) cat(sprintf(": %d features x %d functions x %d samples",
                                                result$features, result$functions,
                                                length(result$samples)))
      cat("\n")
    }
  }
}
