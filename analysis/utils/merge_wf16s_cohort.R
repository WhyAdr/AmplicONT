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

canonicalize_merge_destination <- function(path) {
  # normalizePath() can leave a nonexistent Windows path unchanged, including
  # backslashes and short-name ancestors. Resolve the existing ancestor first.
  path <- gsub("\\\\", "/", path)
  if (file.exists(path) || dir.exists(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  parent <- dirname(path)
  if (identical(parent, path)) stop("E_COHORT_DESTINATION_INVALID: cannot resolve destination.", call. = FALSE)
  file.path(canonicalize_merge_destination(parent), basename(path))
}

cohort_path_identity <- function(path) {
  normalized <- gsub("\\\\", "/", path)
  if (identical(.Platform$OS.type, "windows")) tolower(normalized) else normalized
}

cohort_transaction_paths <- function(destination, transaction_id = NULL) {
  canonical <- canonicalize_merge_destination(destination)
  identity <- cohort_path_identity(canonical)
  # Keep transaction paths usable with Windows path-length limits. A 96-bit
  # destination key and 128-bit per-attempt suffix retain ample collision
  # resistance while avoiding paths that R cannot open after creating them.
  key <- substr(digest::digest(identity, algo = "sha256", serialize = FALSE), 1L, 24L)
  journal <- file.path(dirname(canonical), paste0(".wf16s-cohort-", key, ".json"))
  if (is.null(transaction_id)) return(list(journal = journal, key = key))
  list(
    journal = journal, key = key,
    stage = file.path(dirname(canonical), paste0(".wf16s-cohort-stage-", key, "-", transaction_id)),
    backup = file.path(dirname(canonical), paste0(".wf16s-cohort-backup-", key, "-", transaction_id))
  )
}

new_cohort_transaction_id <- function() {
  entropy <- paste(Sys.time(), Sys.getpid(), tempfile("cohort_tx_"), sep = "|")
  paste0("ct-", substr(digest::digest(entropy, algo = "sha256", serialize = FALSE), 1L, 32L))
}

cohort_directory_snapshot <- function(path) {
  if (!dir.exists(path)) stop(sprintf("Cohort directory is missing: '%s'.", path), call. = FALSE)
  children <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE,
                        recursive = FALSE, include.dirs = TRUE)
  if (!length(children)) {
    entries <- list()
  } else {
    is_dir <- dir.exists(children)
    links <- vapply(children, function(candidate) {
      link <- tryCatch(Sys.readlink(candidate), error = function(e) "")
      length(link) == 1L && nzchar(link)
    }, logical(1))
    if (any(is_dir) || any(links)) {
      stop("Cohort bundle contains a directory or symlink entry.", call. = FALSE)
    }
    children <- children[order(basename(children), method = "radix")]
    entries <- lapply(children, function(candidate) {
      info <- file.info(candidate)
      if (is.na(info$size) || is.na(info$isdir) || info$isdir) {
        stop(sprintf("Cohort bundle entry is not a regular file: '%s'.", candidate), call. = FALSE)
      }
      list(name = basename(candidate), size_bytes = as.numeric(info$size),
           sha256 = compute_file_hash(candidate))
    })
  }
  lines <- if (length(entries)) vapply(entries, function(entry) {
    paste(entry$name, format(entry$size_bytes, scientific = FALSE, trim = TRUE),
          entry$sha256, sep = "\t")
  }, character(1)) else character(0)
  list(entries = entries,
       sha256 = digest::digest(paste(lines, collapse = "\n"), algo = "sha256", serialize = FALSE))
}

cohort_snapshots_equal <- function(left, right) {
  if (is.null(left) || is.null(right) || !is.list(left$entries) || !is.list(right$entries) ||
      length(left$entries) != length(right$entries)) return(FALSE)
  canonical <- function(snapshot) {
    if (!length(snapshot$entries)) return(character(0))
    vapply(snapshot$entries, function(entry) {
      paste(as.character(entry$name),
            format(as.numeric(entry$size_bytes), scientific = FALSE, trim = TRUE),
            as.character(entry$sha256), sep = "\t")
    }, character(1))
  }
  identical(canonical(left), canonical(right)) && identical(as.character(left$sha256), as.character(right$sha256))
}

cohort_basename_is_safe <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value) &&
    !value %in% c(".", "..") && !startsWith(value, ".") &&
    !grepl("[/\\\\[:cntrl:]]", value, perl = TRUE) && identical(basename(value), value)
}

cohort_basenames_are_portably_unique <- function(values) {
  if (!length(values) || !all(vapply(values, cohort_basename_is_safe, logical(1)))) return(FALSE)
  !anyDuplicated(tolower(enc2utf8(values)))
}

validate_cohort_bundle <- function(path) {
  snapshot <- cohort_directory_snapshot(path)
  manifest_path <- file.path(path, "cohort_merge_manifest.json")
  manifest <- tryCatch(jsonlite::fromJSON(manifest_path, simplifyVector = FALSE),
                       error = function(e) NULL)
  if (is.null(manifest) || !identical(as.character(manifest$manifest_version), "1.0") ||
      !is.list(manifest$outputs) || is.null(names(manifest$outputs)) ||
      !length(manifest$outputs)) {
    stop("E_COHORT_BUNDLE_INVALID: cohort manifest is missing or malformed.", call. = FALSE)
  }
  output_names <- names(manifest$outputs)
  all_names <- c(output_names, "cohort_merge_manifest.json")
  if (!cohort_basenames_are_portably_unique(all_names)) {
    stop("E_COHORT_BUNDLE_INVALID: output names are unsafe or case-colliding.", call. = FALSE)
  }
  physical_names <- vapply(snapshot$entries, `[[`, character(1), "name")
  if (!setequal(tolower(physical_names), tolower(all_names))) {
    stop("E_COHORT_BUNDLE_INVALID: physical file census differs from the manifest.", call. = FALSE)
  }
  for (name in output_names) {
    record <- manifest$outputs[[name]]
    actual <- snapshot$entries[[match(name, physical_names)]]
    if (is.null(record$sha256) || is.null(record$size_bytes) ||
        !identical(tolower(as.character(record$sha256)), tolower(actual$sha256)) ||
        !identical(as.numeric(record$size_bytes), as.numeric(actual$size_bytes))) {
      stop(sprintf("E_COHORT_BUNDLE_INVALID: output fingerprint mismatch for '%s'.", name), call. = FALSE)
    }
  }
  config <- tryCatch(yaml::read_yaml(file.path(path, "config.yml")), error = function(e) NULL)
  if (is.null(config) || !identical(config$mode, "cohort") || is.null(config$input)) {
    stop("E_COHORT_BUNDLE_INVALID: generated config.yml is not a parseable cohort config.", call. = FALSE)
  }
  required <- c(abundance_table = "abundance_table_species.tsv", metadata = "metadata.tsv",
                cohort_merge_manifest = "cohort_merge_manifest.json")
  for (key in names(required)) {
    if (!identical(as.character(config$input[[key]]), unname(required[[key]])) ||
        !required[[key]] %in% output_names && key != "cohort_merge_manifest") {
      stop(sprintf("E_COHORT_BUNDLE_INVALID: config input.%s does not resolve within the bundle.", key),
           call. = FALSE)
    }
  }
  cache_name <- as.character(config$taxonomy$cache %||% "")
  if (!cohort_basename_is_safe(cache_name) || !cache_name %in% output_names) {
    stop("E_COHORT_BUNDLE_INVALID: config taxonomy.cache does not resolve to a bundle output.",
         call. = FALSE)
  }
  list(snapshot = snapshot, manifest = manifest, config = config)
}

cohort_path_matches_snapshot <- function(path, expected, validate_bundle = FALSE) {
  if (!dir.exists(path) || is.null(expected)) return(FALSE)
  actual <- tryCatch({
    if (validate_bundle) validate_cohort_bundle(path)$snapshot else cohort_directory_snapshot(path)
  }, error = function(e) NULL)
  !is.null(actual) && cohort_snapshots_equal(actual, expected)
}

assert_cohort_transaction_path <- function(candidate, expected, label) {
  candidate_norm <- canonicalize_merge_destination(candidate)
  expected_norm <- canonicalize_merge_destination(expected)
  if (!identical(cohort_path_identity(candidate_norm), cohort_path_identity(expected_norm)) ||
      !identical(cohort_path_identity(dirname(candidate_norm)),
                 cohort_path_identity(dirname(canonicalize_merge_destination(expected_norm))))) {
    stop(sprintf("E_COHORT_RECOVERY_REQUIRED: journal %s escapes its transaction paths.", label),
         call. = FALSE)
  }
  invisible(candidate_norm)
}

write_cohort_publication_journal <- function(path, journal) {
  atomic_write_json(journal, path, auto_unbox = TRUE)
  if (!file.exists(path)) stop("E_COHORT_JOURNAL_WRITE: publication journal was not persisted.", call. = FALSE)
  invisible(path)
}

remove_cohort_journal_checked <- function(path) {
  if (file.exists(path) && (unlink(path, force = TRUE) != 0L || file.exists(path))) {
    stop(sprintf("E_COHORT_CLEANUP_FAILED: publication journal retained at '%s'.", path), call. = FALSE)
  }
  invisible(TRUE)
}

remove_cohort_directory_checked <- function(path, expected = NULL, allowed_names = NULL) {
  if (!dir.exists(path)) return(invisible(TRUE))
  snapshot <- cohort_directory_snapshot(path)
  if (!is.null(expected) && !cohort_snapshots_equal(snapshot, expected)) {
    stop(sprintf("E_COHORT_RECOVERY_REQUIRED: refusing to remove changed transaction directory '%s'.", path),
         call. = FALSE)
  }
  if (is.null(expected)) {
    names <- if (length(snapshot$entries)) vapply(snapshot$entries, `[[`, character(1), "name") else character(0)
    if (is.null(allowed_names) || any(!names %in% allowed_names)) {
      stop(sprintf("E_COHORT_RECOVERY_REQUIRED: unexpected files prevent cleanup of '%s'.", path),
           call. = FALSE)
    }
  }
  if (unlink(path, recursive = TRUE, force = TRUE) != 0L || dir.exists(path)) {
    stop(sprintf("E_COHORT_CLEANUP_FAILED: could not remove transaction directory '%s'.", path),
         call. = FALSE)
  }
  invisible(TRUE)
}

recover_cohort_publication <- function(destination) {
  canonical <- canonicalize_merge_destination(destination)
  paths <- cohort_transaction_paths(canonical)
  journal_path <- paths$journal
  if (!file.exists(journal_path)) return(invisible(FALSE))
  journal <- tryCatch(jsonlite::fromJSON(journal_path, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(journal) || !identical(journal$journal_schema_version, 1L) ||
      !is.character(journal$transaction_id) || length(journal$transaction_id) != 1L ||
      !grepl("^ct-[0-9a-f]{32}$", journal$transaction_id) ||
      !identical(cohort_path_identity(as.character(journal$destination)), cohort_path_identity(canonical)) ||
      !is.logical(journal$had_prior) || length(journal$had_prior) != 1L || is.na(journal$had_prior) ||
      !journal$phase %in% c("prepared", "ready_to_publish", "prior_moved", "stage_published") ||
      !is.list(journal$allowed_stage_basenames)) {
    stop(sprintf("E_COHORT_RECOVERY_REQUIRED: malformed publication journal at '%s'.", journal_path),
         call. = FALSE)
  }
  expected_paths <- cohort_transaction_paths(canonical, journal$transaction_id)
  assert_cohort_transaction_path(journal$stage, expected_paths$stage, "stage")
  assert_cohort_transaction_path(journal$backup, expected_paths$backup, "backup")
  if (!identical(cohort_path_identity(journal_path), cohort_path_identity(expected_paths$journal))) {
    stop("E_COHORT_RECOVERY_REQUIRED: journal identity mismatch.", call. = FALSE)
  }
  old_snapshot <- journal$prior_snapshot
  new_snapshot <- journal$expected_new_snapshot
  final_old <- isTRUE(journal$had_prior) && cohort_path_matches_snapshot(canonical, old_snapshot, TRUE)
  final_new <- !is.null(new_snapshot) && cohort_path_matches_snapshot(canonical, new_snapshot, TRUE)
  backup_exists <- dir.exists(journal$backup)
  backup_old <- isTRUE(journal$had_prior) && backup_exists &&
    cohort_path_matches_snapshot(journal$backup, old_snapshot, TRUE)
  stage_exists <- dir.exists(journal$stage)
  stage_new <- !is.null(new_snapshot) && stage_exists &&
    cohort_path_matches_snapshot(journal$stage, new_snapshot, TRUE)
  allowed <- vapply(journal$allowed_stage_basenames, as.character, character(1))

  if (final_new) {
    if (backup_exists) {
      if (!backup_old) stop(sprintf(
        "E_COHORT_RECOVERY_REQUIRED: published bundle is valid but changed backup was retained at '%s'.",
        journal$backup), call. = FALSE)
      remove_cohort_directory_checked(journal$backup, expected = old_snapshot)
    }
    if (stage_exists) remove_cohort_directory_checked(journal$stage, expected = new_snapshot)
    remove_cohort_journal_checked(journal_path)
    return(invisible(TRUE))
  }
  if (final_old && !backup_exists) {
    if (stage_exists) {
      if (stage_new) remove_cohort_directory_checked(journal$stage, expected = new_snapshot) else
        remove_cohort_directory_checked(journal$stage, allowed_names = allowed)
    }
    remove_cohort_journal_checked(journal_path)
    return(invisible(TRUE))
  }
  if (!dir.exists(canonical) && backup_exists && isTRUE(journal$had_prior)) {
    if (!backup_old) {
      if (!file.exists(canonical) && file.rename(journal$backup, canonical)) {
        stop(sprintf("E_COHORT_RECOVERY_REQUIRED: unexpected backup was restored to '%s'; inspect it before retrying.",
                     canonical), call. = FALSE)
      }
      stop(sprintf("E_COHORT_RECOVERY_REQUIRED: backup differs from its snapshot and remains at '%s'.",
                   journal$backup), call. = FALSE)
    }
    if (!file.rename(journal$backup, canonical) ||
        !cohort_path_matches_snapshot(canonical, old_snapshot, TRUE)) {
      stop(sprintf("E_COHORT_RECOVERY_REQUIRED: failed to restore verified backup '%s' to '%s'.",
                   journal$backup, canonical), call. = FALSE)
    }
    if (stage_exists) {
      if (stage_new) remove_cohort_directory_checked(journal$stage, expected = new_snapshot) else
        remove_cohort_directory_checked(journal$stage, allowed_names = allowed)
    }
    remove_cohort_journal_checked(journal_path)
    return(invisible(TRUE))
  }
  if (!isTRUE(journal$had_prior) && !dir.exists(canonical) && !file.exists(canonical) && stage_new) {
    if (!file.rename(journal$stage, canonical) ||
        !cohort_path_matches_snapshot(canonical, new_snapshot, TRUE)) {
      stop(sprintf("E_COHORT_RECOVERY_REQUIRED: failed to publish verified stage '%s'.", journal$stage),
           call. = FALSE)
    }
    remove_cohort_journal_checked(journal_path)
    return(invisible(TRUE))
  }
  if (!isTRUE(journal$had_prior) && !dir.exists(canonical) && !file.exists(canonical) &&
      !stage_new && stage_exists && identical(journal$phase, "prepared")) {
    remove_cohort_directory_checked(journal$stage, allowed_names = allowed)
    remove_cohort_journal_checked(journal_path)
    return(invisible(TRUE))
  }
  stop(sprintf("E_COHORT_RECOVERY_REQUIRED: ambiguous publication state for '%s'; stage and backup retained.",
               canonical), call. = FALSE)
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
                               validate_only = FALSE, overwrite = FALSE,
                               publication_hook = NULL, lock_timeout_ms = 10000) {
  script_dir <- get_merge_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")

  # Source dependencies if not loaded
  if (!exists("compare_upstream_params", mode = "function")) {
    source(file.path(script_dir, "classifier_contract.R"))
  }
  if (!exists("read_upstream_params", mode = "function")) {
    source(file.path(script_dir, "io.R"))
  }
  if (!exists("compute_file_hash", mode = "function") ||
      !exists("acquire_output_lock", mode = "function") ||
      !exists("atomic_write_json", mode = "function")) {
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
  norm_dest <- canonicalize_merge_destination(output_dir)
  for (r_root in unique(c(resolved_roots, normalizePath(input_dirs, winslash = "/")))) {
    if (identical(tolower(norm_dest), tolower(r_root)) ||
        startsWith(tolower(norm_dest), paste0(tolower(r_root), "/")) ||
        startsWith(tolower(r_root), paste0(tolower(norm_dest), "/"))) {
      stop(sprintf("E_COHORT_OUTPUT_OVERLAP: output directory '%s' overlaps with input directory '%s'.",
                   output_dir, r_root), call. = FALSE)
    }
  }

  destination_lock <- NULL
  publication_started <- FALSE
  prior_destination_snapshot <- NULL
  transaction_journal_path <- cohort_transaction_paths(norm_dest)$journal
  if (!isTRUE(validate_only)) {
    destination_lock <- acquire_output_lock(norm_dest, timeout_ms = lock_timeout_ms)
    on.exit({
      if (!publication_started) {
        tryCatch(recover_cohort_publication(norm_dest), error = function(e) {
          warning(conditionMessage(e), call. = FALSE)
        })
      }
      release_output_lock(destination_lock)
    }, add = TRUE)
    recover_cohort_publication(norm_dest)
  } else if (file.exists(transaction_journal_path)) {
    stop(sprintf("E_COHORT_RECOVERY_REQUIRED: pending publication journal at '%s'; validate-only cannot mutate it.",
                 transaction_journal_path), call. = FALSE)
  }

  if (file.exists(norm_dest) && !dir.exists(norm_dest)) stop("E_COHORT_DESTINATION_EXISTS: destination is a file.", call. = FALSE)
  if (dir.exists(norm_dest)) {
    if (!isTRUE(overwrite)) stop("E_COHORT_DESTINATION_EXISTS: destination already exists.", call. = FALSE)
    prior_bundle <- tryCatch(validate_cohort_bundle(norm_dest), error = function(e) NULL)
    if (is.null(prior_bundle)) {
      stop("E_COHORT_DESTINATION_UNOWNED: overwrite requires an unchanged merger bundle with no unrelated files.", call. = FALSE)
    }
    prior_destination_snapshot <- prior_bundle$snapshot
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
  reserved_names <- c("config.yml", "metadata.tsv", "abundance_table_species.tsv",
                      "cohort_merge_manifest.json")
  if (!cohort_basename_is_safe(cache_basename) ||
      !cohort_basenames_are_portably_unique(c(reserved_names, cache_basename))) {
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

  # 12. Journaled staging and publication under the destination lock.
  norm_dest <- canonicalize_merge_destination(output_dir)
  parent_dir <- dirname(norm_dest)
  if (!dir.exists(parent_dir) && !dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE) &&
      !dir.exists(parent_dir)) {
    stop(sprintf("E_COHORT_STAGE_FAILED: could not create destination parent '%s'.", parent_dir),
         call. = FALSE)
  }
  transaction_id <- new_cohort_transaction_id()
  tx_paths <- cohort_transaction_paths(norm_dest, transaction_id)
  stage_dir <- tx_paths$stage
  backup_dir <- tx_paths$backup
  stage_allowed_names <- c("abundance_table_species.tsv", "metadata.tsv", "config.yml",
                           cache_basename, "cohort_merge_manifest.json")
  had_prior <- !is.null(prior_destination_snapshot)
  if (had_prior && !cohort_path_matches_snapshot(norm_dest, prior_destination_snapshot, TRUE)) {
    stop("E_COHORT_DESTINATION_CHANGED: destination changed after ownership validation; no files were moved.",
         call. = FALSE)
  }
  if (!had_prior && (dir.exists(norm_dest) || file.exists(norm_dest))) {
    stop("E_COHORT_DESTINATION_CHANGED: a destination appeared after validation; it was left untouched.",
         call. = FALSE)
  }
  journal <- list(
    journal_schema_version = 1L, transaction_id = transaction_id,
    destination = norm_dest, stage = stage_dir,
    backup = backup_dir, had_prior = had_prior,
    prior_snapshot = prior_destination_snapshot, expected_new_snapshot = NULL,
    allowed_stage_basenames = as.list(stage_allowed_names), phase = "prepared"
  )
  write_cohort_publication_journal(tx_paths$journal, journal)
  if (!dir.create(stage_dir, recursive = FALSE, showWarnings = FALSE)) {
    stop(sprintf("E_COHORT_STAGE_FAILED: could not create transaction stage '%s'.", stage_dir),
         call. = FALSE)
  }
  run_publication_hook <- function(phase) {
    if (is.function(publication_hook)) publication_hook(phase, norm_dest, stage_dir, backup_dir)
  }

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
    copied <- file.copy(taxonomy_cache, out_cache_path, overwrite = FALSE)
    if (!isTRUE(copied) || !identical(compute_file_hash(out_cache_path), compute_file_hash(taxonomy_cache))) {
      stop("E_COHORT_CACHE_COPY_FAILED: taxonomy cache was not copied byte-identically.", call. = FALSE)
    }
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

  staged_bundle <- tryCatch(validate_cohort_bundle(stage_dir), error = function(e) {
    stop(sprintf("E_COHORT_STAGE_INVALID: %s", conditionMessage(e)), call. = FALSE)
  })
  journal$expected_new_snapshot <- staged_bundle$snapshot
  journal$phase <- "ready_to_publish"
  write_cohort_publication_journal(tx_paths$journal, journal)
  publication_started <- TRUE
  assert_merge_sources()
  run_publication_hook("before_prior_rename")
  if (had_prior) {
    if (!cohort_path_matches_snapshot(norm_dest, prior_destination_snapshot, TRUE)) {
      stop("E_COHORT_DESTINATION_CHANGED: destination changed before publication; it was left untouched.",
           call. = FALSE)
    }
    if (!file.rename(norm_dest, backup_dir)) {
      stop(sprintf("E_COHORT_PUBLISH_FAILED: could not preserve prior destination at '%s'.", backup_dir),
           call. = FALSE)
    }
    run_publication_hook("after_prior_rename_before_journal")
    journal$phase <- "prior_moved"
    write_cohort_publication_journal(tx_paths$journal, journal)
    run_publication_hook("after_prior_rename")
    if (!cohort_path_matches_snapshot(backup_dir, prior_destination_snapshot, TRUE)) {
      if (!dir.exists(norm_dest) && !file.exists(norm_dest) && file.rename(backup_dir, norm_dest)) {
        stop("E_COHORT_DESTINATION_CHANGED: destination changed at the backup boundary; the changed files were restored and retained.",
             call. = FALSE)
      }
      stop(sprintf("E_COHORT_RECOVERY_REQUIRED: changed prior destination retained at '%s'.", backup_dir),
           call. = FALSE)
    }
  }

  run_publication_hook("before_stage_rename")
  if (dir.exists(norm_dest) || file.exists(norm_dest)) {
    stop("E_COHORT_DESTINATION_CHANGED: a new destination appeared before publish; it was left untouched.",
         call. = FALSE)
  }
  assert_merge_sources()
  if (!file.rename(stage_dir, norm_dest)) {
    if (had_prior && !dir.exists(norm_dest) && !file.exists(norm_dest) &&
        cohort_path_matches_snapshot(backup_dir, prior_destination_snapshot, TRUE)) {
      if (!file.rename(backup_dir, norm_dest) ||
          !cohort_path_matches_snapshot(norm_dest, prior_destination_snapshot, TRUE)) {
        stop(sprintf("E_COHORT_RECOVERY_REQUIRED: publication failed and prior bundle could not be verified after restore; backup or destination retained at '%s'.",
                     if (dir.exists(backup_dir)) backup_dir else norm_dest), call. = FALSE)
      }
      recover_cohort_publication(norm_dest)
    }
    stop(sprintf("E_COHORT_PUBLISH_FAILED: could not publish staged cohort to '%s'.", norm_dest),
         call. = FALSE)
  }
  run_publication_hook("after_stage_rename_before_journal")
  journal$phase <- "stage_published"
  write_cohort_publication_journal(tx_paths$journal, journal)
  run_publication_hook("after_stage_rename")
  if (!cohort_path_matches_snapshot(norm_dest, staged_bundle$snapshot, TRUE)) {
    stop(sprintf("E_COHORT_RECOVERY_REQUIRED: published bundle changed; recovery material retained at '%s'.",
                 if (dir.exists(backup_dir)) backup_dir else norm_dest), call. = FALSE)
  }

  if (had_prior && dir.exists(backup_dir)) {
    run_publication_hook("before_backup_cleanup")
    if (!cohort_path_matches_snapshot(backup_dir, prior_destination_snapshot, TRUE)) {
      stop(sprintf("E_COHORT_RECOVERY_REQUIRED: changed backup was retained at '%s'.", backup_dir),
           call. = FALSE)
    }
    remove_cohort_directory_checked(backup_dir, expected = prior_destination_snapshot)
  }
  remove_cohort_journal_checked(tx_paths$journal)

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
