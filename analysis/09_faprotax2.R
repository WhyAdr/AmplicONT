# Official FAPROTAX-mapper 2.0 adapter. No upstream software/data are bundled.
FAPROTAX2_ARCHIVE_SHA256 <- "859b99133b55f1e4f637e6ac14a804f3321ed27f76ee03d95e3f215382206bac"
FAPROTAX2_RELEASE_HASHES <- c(
  mapper_script = "8bd6ff97aa53f17fb0122b98c436275c230e34629f7a83c446b2b3854e53a64f",
  default_config.tsv = "d9c597304f966a2b2ab96a5db8de7a14a894240d074c3ab5ddcad93c27584c37",
  reference_alignments_sanitized_subset.fasta.gz = "3a21836125cad9cc9bc0ee1cf1e62bda2a5cf98c54ea43723eeaff2f74dcbbe2",
  reference_alignments_sanitized.fasta.gz = "0f30fd59650c738494a901c2651a022b51d058c56488990c510f584400cca5b9",
  reference2species_and_traits.tsv = "a8b3dc971470c224a7f33874021622bf24270e695076740cb2f151dd7a151f96",
  tree_model.txt = "a3c3382d0c2d33c1a4662f3af93c28ba0fe4f91c5d49cd500930faac33b7a7b9",
  tree.tre = "1b40b2e3f4165b6eea49290064aadaff3c74b4744c16e91c35506afcfdfa977d"
)

faprotax2_asset_paths <- function(cfg) {
  f2 <- cfg$faprotax2
  if (!identical(f2$mode, "taxonomy")) {
    stop("FAPROTAX2 sequence mode has not passed the native workflow gate; see patch/faprotax2-sequence-feasibility.md.", call. = FALSE)
  }
  if (is.null(f2$mapper_path) || is.null(f2$reference_dir)) {
    stop("FAPROTAX2 requires explicit mapper_path and reference_dir; no automatic downloads are performed.", call. = FALSE)
  }
  c(mapper_script = f2$mapper_path,
    stats::setNames(file.path(f2$reference_dir, names(FAPROTAX2_RELEASE_HASHES)[-1]),
                    names(FAPROTAX2_RELEASE_HASHES)[-1]))
}

faprotax2_resource_inventory <- function(cfg) {
  paths <- faprotax2_asset_paths(cfg)
  records <- lapply(names(paths), function(role) {
    record <- fingerprint_file(paths[[role]], paste0("faprotax2.", role))
    if (!identical(record$sha256, FAPROTAX2_RELEASE_HASHES[[role]])) {
      stop(sprintf("FAPROTAX2 release hash mismatch for '%s'; expected the audited mapper 2.0 archive.", role), call. = FALSE)
    }
    record
  })
  stats::setNames(records, names(paths))
}

faprotax2_child_environment <- function(cfg) {
  env <- Sys.getenv()
  libraries <- unique(c(cfg$faprotax2$mapper_library, .libPaths()))
  libraries <- normalizePath(libraries, winslash = "/", mustWork = TRUE)
  env["R_LIBS"] <- paste(libraries, collapse = .Platform$path.sep)
  env["RENV_CONFIG_AUTOLOADER_ENABLED"] <- "FALSE"
  env
}

validate_faprotax2_runtime <- function(cfg) {
  resources <- faprotax2_resource_inventory(cfg)
  if (!is.null(cfg$faprotax2$mapper_library) && !dir.exists(cfg$faprotax2$mapper_library)) {
    stop("FAPROTAX2 mapper_library does not exist.", call. = FALSE)
  }
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  env <- faprotax2_child_environment(cfg)
  code <- paste(
    'pkgs <- c("castor","naturalsort","optparse","stringr","Rcpp","Matrix","RSpectra","jsonlite");',
    'for(p in pkgs) if(!requireNamespace(p,quietly=TRUE)) stop(paste("Missing mapper dependency",p));',
    'if(as.character(packageVersion("castor")) != "1.8.7" || as.character(packageVersion("naturalsort")) != "0.1.3") stop("Expected audited castor 1.8.7 and naturalsort 0.1.3");',
    'records <- lapply(pkgs,function(p) list(package=p,version=as.character(packageVersion(p)),path=find.package(p),description_sha256=digest::digest(file=file.path(find.package(p),"DESCRIPTION"),algo="sha256")));',
    'cat(jsonlite::toJSON(list(r=R.version.string,platform=R.version$platform,packages=records),auto_unbox=TRUE))'
  )
  result <- processx::run(rscript, c("--vanilla", "-e", code), env = env,
                          timeout = 30000, error_on_status = FALSE, cleanup_tree = TRUE)
  if (result$status != 0L) {
    stop(sprintf("FAPROTAX2 dependency preflight failed: %s", trimws(result$stderr)), call. = FALSE)
  }
  list(rscript = rscript, env = env, resources = resources,
       dependency_runtime = jsonlite::fromJSON(result$stdout, simplifyVector = FALSE))
}

prepare_faprotax2_taxonomy <- function(context) {
  counts <- context$count_matrix[, context$samples, drop = FALSE]
  indices <- setdiff(seq_len(nrow(counts)), context$unclass_index)
  source <- context$taxonomy[indices, , drop = FALSE]
  keep <- !is.na(source$superkingdom) & source$superkingdom %in% c("Bacteria", "Archaea") &
    rowSums(counts[indices, , drop = FALSE]) > 0
  if (!any(keep)) return(NULL)
  indices <- indices[keep]
  source <- source[keep, , drop = FALSE]
  # Mapper references use domain, phylum, class, order, family, genus, species.
  ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  taxonomy <- as.matrix(source[, ranks, drop = FALSE])
  taxonomy[is.na(taxonomy) | tolower(trimws(taxonomy)) %in% c("unknown", "unclassified")] <- ""
  ids <- sprintf("Taxon_%06d", indices)
  counts <- counts[indices, , drop = FALSE]
  rownames(counts) <- ids
  list(counts = counts, taxonomy = stats::setNames(apply(taxonomy, 1, paste, collapse = ";"), ids),
       feature_ids = ids, source_taxonomy = source)
}

read_faprotax2_matrix <- function(path, id_column) {
  table <- utils::read.delim(path, check.names = FALSE, quote = "", comment.char = "",
                            fileEncoding = "UTF-8", na.strings = "NA", stringsAsFactors = FALSE)
  if (ncol(table) < 2L || !identical(names(table)[1], id_column) ||
      anyDuplicated(names(table)) || anyNA(table[[1]]) ||
      any(!nzchar(table[[1]])) || anyDuplicated(table[[1]])) {
    stop(sprintf("FAPROTAX2 invalid matrix schema: %s", basename(path)), call. = FALSE)
  }
  values <- table[, -1, drop = FALSE]
  if (!all(vapply(values, is.numeric, logical(1)))) {
    # A column consisting entirely of NA is read as logical by read.delim.
    valid <- vapply(values, function(x) is.numeric(x) || (is.logical(x) && all(is.na(x))), logical(1))
    if (!all(valid)) stop("FAPROTAX2 output contains a nonnumeric column.", call. = FALSE)
  }
  result <- as.matrix(values)
  storage.mode(result) <- "numeric"
  rownames(result) <- table[[1]]
  result
}

validate_faprotax2_matrices <- function(states, probabilities, abundance, counts) {
  valid_names <- function(x) !is.null(x) && length(x) > 0L && !anyNA(x) &&
    all(nzchar(x)) && !anyDuplicated(x)
  matrices <- list(states, probabilities, abundance, counts)
  if (!all(vapply(matrices, function(x) is.matrix(x) && is.numeric(x) &&
                  valid_names(rownames(x)) && valid_names(colnames(x)), logical(1)))) {
    stop("FAPROTAX2 matrices require unique sample, feature and function identifiers with nonempty names.", call. = FALSE)
  }
  if (!identical(rownames(states), rownames(counts)) ||
      !identical(dimnames(probabilities), dimnames(states)) ||
      !identical(rownames(abundance), colnames(states)) ||
      !identical(colnames(abundance), colnames(counts)) ||
      anyDuplicated(rownames(states)) || anyDuplicated(colnames(states)) || !ncol(states)) {
    stop("FAPROTAX2 matrix feature, function or sample identifiers/order disagree.", call. = FALSE)
  }
  if (any(!is.na(states) & !states %in% c(0, 1)) ||
      any(!is.na(probabilities) & (!is.finite(probabilities) | probabilities < 0 | probabilities > 1)) ||
      anyNA(abundance) || any(!is.finite(abundance) | abundance < 0) ||
      anyNA(counts) || any(!is.finite(counts) | counts < 0)) {
    stop("FAPROTAX2 matrices contain invalid states, probabilities or abundances.", call. = FALSE)
  }
  binary <- states
  binary[is.na(binary)] <- 0
  expected <- t(binary) %*% counts
  if (!isTRUE(all.equal(unname(abundance), unname(expected), tolerance = 1e-12))) {
    stop("FAPROTAX2 abundance differs from unweighted states multiplied by input counts.", call. = FALSE)
  }
  invisible(TRUE)
}

summarize_faprotax2_coverage <- function(states, counts, context) {
  present <- rowSums(states == 1, na.rm = TRUE) > 0
  unknown_without_presence <- !present & rowSums(is.na(states)) > 0
  data.frame(SampleID = colnames(counts),
    EligibleTaxonAbundance = colSums(counts),
    PredictedPresenceAbundance = colSums(counts[present, , drop = FALSE]),
    NoPredictedPresenceAbundance = colSums(counts[!present, , drop = FALSE]),
    UnknownWithoutPresenceAbundance = colSums(counts[unknown_without_presence, , drop = FALSE]),
    KnownAbsenceOnlyAbundance = colSums(counts[!present & !unknown_without_presence, , drop = FALSE]),
    EligibleTaxa = colSums(counts > 0),
    TaxaWithPredictedPresence = colSums(counts[present, , drop = FALSE] > 0),
    CountBasis = context$count_basis, stringsAsFactors = FALSE, check.names = FALSE)
}

run_faprotax2 <- function(context, runner = processx::run, runtime = NULL) {
  cfg <- context$config
  runtime <- runtime %||% validate_faprotax2_runtime(cfg)
  out <- cfg$output$dirs$faprotax2
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  prepared <- prepare_faprotax2_taxonomy(context)
  if (is.null(prepared)) {
    path <- file.path(out, "faprotax2_skipped.tsv")
    write.table(data.frame(Status = "skipped", Reason = "No positive-count Bacteria or Archaea rows"),
                path, sep = "\t", row.names = FALSE, quote = FALSE)
    return(list(status = "skipped", reason = "No eligible prokaryotic taxa", outputs = path))
  }
  input <- data.frame(OTU = prepared$feature_ids, prepared$counts, taxonomy = unname(prepared$taxonomy),
                      check.names = FALSE, stringsAsFactors = FALSE)
  input_path <- file.path(out, "mapper_input.tsv")
  write.table(input, input_path, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  source_path <- file.path(out, "taxon_identity.tsv")
  write.table(data.frame(FeatureID = prepared$feature_ids,
                         TaxonPath = prepared$source_taxonomy$TaxonPath,
                         MapperTaxonomy = unname(prepared$taxonomy)),
                source_path, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  raw <- file.path(out, c("mapper_function_abundance.tsv", "mapper_feature_states.tsv",
                         "mapper_annotation_probabilities.tsv", "mapper_report.txt"))
  log_paths <- file.path(out, c("mapper_stdout.txt", "mapper_stderr.txt"))
  args <- c("--vanilla", cfg$faprotax2$mapper_path,
    "-d", cfg$faprotax2$reference_dir, "-i", input_path, "-t", "column:taxonomy", "-k", "OTU",
    "-o", raw[1], "-e", raw[2], "-q", raw[3], "-r", raw[4],
    "--normalize", "none", "--hsp_algorithm", "config", "--unknown_ref_states", "config",
    "--placement_failure", "error", "--Nthreads", as.character(cfg$faprotax2$nthreads))
  process <- runner(runtime$rscript, args, env = runtime$env,
                    timeout = cfg$faprotax2$timeout_seconds * 1000,
                    error_on_status = FALSE, cleanup_tree = TRUE)
  writeLines(process$stdout, log_paths[1], useBytes = TRUE)
  writeLines(process$stderr, log_paths[2], useBytes = TRUE)
  if (!identical(process$status, 0L)) {
    stop(sprintf("FAPROTAX2 mapper exited %s: %s", process$status,
                 paste(tail(strsplit(process$stderr, "\n", fixed = TRUE)[[1]], 8), collapse = "\n")), call. = FALSE)
  }
  if (!all(file.exists(raw))) stop("FAPROTAX2 mapper omitted required raw outputs.", call. = FALSE)
  states <- read_faprotax2_matrix(raw[2], "OTU")
  probabilities <- read_faprotax2_matrix(raw[3], "OTU")
  abundance <- read_faprotax2_matrix(raw[1], "function")
  validate_faprotax2_matrices(states, probabilities, abundance, prepared$counts)
  coverage <- summarize_faprotax2_coverage(states, prepared$counts, context)
  long <- expand.grid(Function = rownames(abundance), SampleID = colnames(abundance),
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  long$EstimatedFunctionAbundance <- as.numeric(abundance)
  long$CountBasis <- context$count_basis
  long$Unit <- "sum_of_input_taxon_abundances_with_predicted_presence"
  summary_paths <- file.path(out, c("faprotax2_function_abundance.tsv", "faprotax2_coverage.tsv"))
  write.table(long, summary_paths[1], sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  write.table(coverage, summary_paths[2], sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  provenance_path <- file.path(out, "faprotax2_provenance.json")
  jsonlite::write_json(list(schema_version = 1L, mapper = "FAPROTAX-mapper", mapper_version = "2.0 beta",
    mode = "taxonomy", mapper_modified = FALSE, archive_sha256 = FAPROTAX2_ARCHIVE_SHA256,
    resources = runtime$resources, runtime = runtime$dependency_runtime,
    executable = runtime$rscript, arguments = json_array(args),
    count_basis = context$count_basis, normalize = "none", probability_weighting = FALSE,
    hsp_algorithm = "config", unknown_ref_states = "config", nthreads = cfg$faprotax2$nthreads,
    input_sha256 = compute_file_hash(input_path),
    interpretation = paste("Overlapping taxon-based phenotype predictions; no evidence of genes, expression or measured activity.",
      "Unrecognized taxonomies may receive broad fallback predictions. Coverage means positive function states, not successful taxonomic placement.",
      "Annotation probabilities are unitless; estimated abundances sum the input taxon-count basis without probability weighting.")),
    provenance_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  list(status = "completed", outputs = c(input_path, source_path, raw, log_paths, summary_paths, provenance_path),
       warnings = "FAPROTAX2 predicts overlapping phenotypes from taxonomic placement; unknown lineages may receive broad fallback predictions.")
}
