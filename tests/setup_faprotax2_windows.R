fail <- function(message) stop(message, call. = FALSE)
default_if_null <- function(value, fallback) if (is.null(value)) fallback else value

usage <- function() {
  cat(paste0(
    "Prepare the pinned FAPROTAX-mapper 2.0 native-Windows runtime.\n\n",
    "Usage:\n",
    "  Rscript --vanilla tests/setup_faprotax2_windows.R [--locked-library PATH] [--mapper-library PATH]\n",
    "  Rscript --vanilla tests/setup_faprotax2_windows.R --help\n\n",
    "Paths are absolute or relative to the repository root. Defaults:\n",
    "  --locked-library  renv/library/windows/R-4.5/x86_64-w64-mingw32\n",
    "  --mapper-library  output/faprotax2-runtime/library\n\n",
    "Only native Windows R 4.5.3 is supported. Pinned mapper assets and binary packages\n",
    "are SHA-256 checked before use. Packages install only into --mapper-library.\n"
  ))
}

parse_arguments <- function(args) {
  if (identical(args, "--help")) {
    usage()
    quit(save = "no", status = 0L, runLast = FALSE)
  }
  parsed <- list(locked_library = NULL, mapper_library = NULL)
  names_by_flag <- c("--locked-library" = "locked_library",
                     "--mapper-library" = "mapper_library")
  i <- 1L
  while (i <= length(args)) {
    flag <- args[[i]]
    if (!flag %in% names(names_by_flag)) {
      fail(sprintf("Unsupported argument: %s\nRun with --help for usage.", flag))
    }
    if (i == length(args) || startsWith(args[[i + 1L]], "--") ||
        !nzchar(args[[i + 1L]])) {
      fail(sprintf("Argument %s requires a non-empty path.", flag))
    }
    key <- unname(names_by_flag[[flag]])
    if (!is.null(parsed[[key]])) fail(sprintf("Argument %s may be supplied only once.", flag))
    parsed[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  parsed
}

absolute_from_repo <- function(path, repo_root) {
  normalized <- gsub("\\\\", "/", path)
  is_absolute <- grepl("^[A-Za-z]:/", normalized) || startsWith(normalized, "//")
  if (!is_absolute) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

resolve_existing_directory <- function(path, repo_root, label) {
  candidate <- absolute_from_repo(path, repo_root)
  if (!dir.exists(candidate)) fail(sprintf("%s does not exist as a directory.", label))
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

canonical_candidate_directory <- function(path, repo_root, label) {
  candidate <- absolute_from_repo(path, repo_root)
  if (file.exists(candidate) && !dir.exists(candidate)) {
    fail(sprintf("%s exists but is not a directory.", label))
  }
  if (dir.exists(candidate)) return(normalizePath(candidate, winslash = "/", mustWork = TRUE))

  cursor <- candidate
  suffix <- character()
  while (!dir.exists(cursor)) {
    if (file.exists(cursor)) fail(sprintf("An existing path component for %s is not a directory.", label))
    parent <- dirname(cursor)
    if (identical(parent, cursor)) fail(sprintf("Could not resolve the parent directory for %s.", label))
    suffix <- c(basename(cursor), suffix)
    cursor <- parent
  }
  resolved <- normalizePath(cursor, winslash = "/", mustWork = TRUE)
  if (length(suffix)) resolved <- paste(c(resolved, suffix), collapse = "/")
  sub("/+$", "", resolved)
}

path_key <- function(path) tolower(sub("/+$", "", gsub("\\\\", "/", path)))

same_or_descendant <- function(path, root) {
  p <- path_key(path)
  r <- path_key(root)
  identical(p, r) || startsWith(p, paste0(r, "/"))
}

paths_overlap <- function(left, right) {
  same_or_descendant(left, right) || same_or_descendant(right, left)
}

sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

ensure_cached_archive <- function(path, url, expected_sha256, label) {
  if (dir.exists(path)) fail(sprintf("%s cache target is a directory; refusing to overwrite it.", label))
  if (file.exists(path)) {
    observed <- sha256_file(path)
    if (!identical(tolower(observed), tolower(expected_sha256))) {
      fail(sprintf("%s cache checksum mismatch; refusing to replace the existing file.", label))
    }
    return(invisible(list(path = path, sha256 = observed, downloaded = FALSE)))
  }

  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
    fail(sprintf("Could not create the cache directory for %s.", label))
  }
  temporary <- tempfile(pattern = ".faprotax2-download-", tmpdir = parent)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  status <- tryCatch(
    utils::download.file(url, destfile = temporary, mode = "wb", method = "libcurl", quiet = TRUE),
    error = function(e) fail(sprintf("Could not download %s from its pinned URL: %s", label, e$message))
  )
  if (!identical(as.integer(status), 0L) || !file.exists(temporary) || dir.exists(temporary)) {
    fail(sprintf("Download failed for %s.", label))
  }
  observed <- sha256_file(temporary)
  if (!identical(tolower(observed), tolower(expected_sha256))) {
    fail(sprintf("Downloaded %s checksum mismatch; the file was not cached.", label))
  }
  if (file.exists(path) || dir.exists(path) || !file.rename(temporary, path)) {
    fail(sprintf("Could not place the verified %s into its cache.", label))
  }
  invisible(list(path = path, sha256 = observed, downloaded = TRUE))
}

archive_entries <- function(tar_exe, archive, process_runner = processx::run) {
  result <- process_runner(tar_exe, c("-tzf", archive), timeout = 60,
                           error_on_status = FALSE, cleanup_tree = TRUE)
  if (isTRUE(result$timeout)) fail("Timed out listing the pinned mapper archive with native tar.exe.")
  if (!identical(result$status, 0L)) fail("Native tar.exe could not list the pinned mapper archive.")
  entries <- strsplit(result$stdout, "\\r?\\n", perl = TRUE)[[1L]]
  entries <- entries[nzchar(entries)]
  entries <- sub("^\\./", "", gsub("\\\\", "/", entries))
  if (!length(entries)) fail("The pinned mapper archive has no entries.")
  unsafe <- grepl("^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)", entries, perl = TRUE)
  if (any(unsafe)) fail("The pinned mapper archive contains an unsafe path; extraction was refused.")
  top <- vapply(strsplit(entries, "/", fixed = TRUE), `[[`, character(1), 1L)
  # The official macOS-created tar has one root-level AppleDouble sidecar.
  allowed <- top == "FAPROTAX-mapper_v2.0" | entries == "._FAPROTAX-mapper_v2.0"
  if (any(!allowed)) {
    fail("The pinned mapper archive has an unexpected top-level path; extraction was refused.")
  }
  entries
}

validate_mapper_resource_dir <- function(mapper_dir) {
  cfg <- list(faprotax2 = list(
    mode = "taxonomy",
    mapper_path = file.path(mapper_dir, "faprotax.R"),
    reference_dir = file.path(mapper_dir, "FAPROTAX_ref_data")
  ))
  resources <- tryCatch(faprotax2_resource_inventory(cfg), error = function(e) {
    fail(sprintf("Production mapper resource validation failed: %s", e$message))
  })
  if (length(resources) != 7L || !identical(names(resources), names(FAPROTAX2_RELEASE_HASHES))) {
    fail("Production resource validation did not return the mapper and six pinned reference assets.")
  }
  invisible(resources)
}

extract_mapper_archive <- function(archive, unpacked_root, mapper_dir,
                                   process_runner = processx::run) {
  if (file.exists(mapper_dir) && !dir.exists(mapper_dir)) {
    fail("The expected unpacked mapper path exists but is not a directory.")
  }
  if (dir.exists(mapper_dir)) {
    validate_mapper_resource_dir(mapper_dir)
    return(invisible(FALSE))
  }

  tar_exe <- unname(Sys.which("tar"))
  if (!nzchar(tar_exe)) fail("Native tar.exe is required to unpack the mapper archive.")
  tar_exe <- normalizePath(tar_exe, winslash = "/", mustWork = TRUE)
  if (!tolower(basename(tar_exe)) %in% c("tar.exe", "tar")) {
    fail("The resolved tar command is not the native tar.exe program.")
  }
  archive_entries(tar_exe, archive, process_runner)

  if (!dir.exists(unpacked_root) && !dir.create(unpacked_root, recursive = TRUE, showWarnings = FALSE)) {
    fail("Could not create the mapper extraction directory.")
  }
  canonical_root <- normalizePath(unpacked_root, winslash = "/", mustWork = TRUE)
  expected_final <- file.path(canonical_root, "FAPROTAX-mapper_v2.0")
  if (!identical(path_key(normalizePath(mapper_dir, winslash = "/", mustWork = FALSE)),
                 path_key(expected_final))) {
    fail("The mapper extraction destination must be the pinned archive's direct child directory.")
  }
  if (file.exists(expected_final) || dir.exists(expected_final)) {
    fail("The mapper destination appeared during extraction setup; it was left untouched.")
  }

  stage_root <- tempfile(pattern = ".faprotax2-extract-", tmpdir = canonical_root)
  if (!dir.create(stage_root, showWarnings = FALSE)) {
    fail("Could not create a unique temporary mapper extraction directory.")
  }
  stage_owned <- TRUE
  on.exit({
    if (stage_owned && (file.exists(stage_root) || dir.exists(stage_root))) {
      canonical_stage <- normalizePath(stage_root, winslash = "/", mustWork = FALSE)
      is_owned_child <- startsWith(path_key(canonical_stage),
                                   paste0(path_key(canonical_root), "/")) &&
        startsWith(basename(canonical_stage), ".faprotax2-extract-")
      if (is_owned_child) {
        removed <- unlink(stage_root, recursive = TRUE, force = TRUE)
        if (!identical(as.integer(removed), 0L)) {
          warning("Could not remove the owned temporary mapper extraction directory.", call. = FALSE)
        }
      } else {
        warning("Temporary mapper extraction cleanup refused an unexpected path.", call. = FALSE)
      }
    }
  }, add = TRUE)

  result <- process_runner(tar_exe, c("-xzf", archive, "-C", stage_root),
                           timeout = 120, error_on_status = FALSE, cleanup_tree = TRUE)
  if (isTRUE(result$timeout)) fail("Timed out extracting the pinned mapper archive with native tar.exe.")
  if (!identical(result$status, 0L)) fail("Native tar.exe failed to extract the pinned mapper archive.")
  stage_mapper <- file.path(stage_root, "FAPROTAX-mapper_v2.0")
  if (!dir.exists(stage_mapper)) fail("The expected mapper directory was not created by extraction.")
  validate_mapper_resource_dir(stage_mapper)
  if (file.exists(expected_final) || dir.exists(expected_final)) {
    fail("The mapper destination appeared before promotion; it was left untouched.")
  }
  if (!file.rename(stage_mapper, expected_final)) {
    fail("Could not promote the validated mapper directory into its final location.")
  }
  validate_mapper_resource_dir(expected_final)
  invisible(TRUE)
}

lock_versions_from_installed <- function(packages, libraries) {
  versions <- setNames(rep(NA_character_, length(packages)), names(packages))
  sources <- setNames(rep(NA_character_, length(packages)), names(packages))
  for (library in libraries) {
    if (!dir.exists(library)) next
    installed <- utils::installed.packages(lib.loc = library)
    if (!nrow(installed)) next
    names_here <- installed[, "Package"]
    wanted <- intersect(names(packages), names_here)
    unresolved <- wanted[is.na(versions[wanted])]
    if (!length(unresolved)) next
    rows <- match(unresolved, names_here)
    versions[unresolved] <- installed[rows, "Version"]
    sources[unresolved] <- library
  }
  list(versions = versions, sources = sources)
}

verify_locked_library <- function(repo_root, locked_library) {
  lock_path <- file.path(repo_root, "renv.lock")
  if (!file.exists(lock_path) || dir.exists(lock_path)) fail("The project renv.lock is missing.")
  lock <- tryCatch(jsonlite::fromJSON(lock_path, simplifyVector = FALSE),
                   error = function(e) fail(sprintf("Could not read renv.lock: %s", e$message)))
  packages <- lock[["Packages"]]
  if (!is.list(packages) || length(packages) != 69L || is.null(names(packages))) {
    fail("Expected exactly 69 package entries in the committed renv.lock.")
  }
  lock_r <- as.character(lock[["R"]][["Version"]])
  if (!identical(lock_r, "4.5.3")) fail("renv.lock does not target R 4.5.3.")
  expected <- vapply(packages, function(record) as.character(record[["Version"]]), character(1))
  if (anyNA(expected) || any(!nzchar(expected))) fail("renv.lock contains a package without a version.")

  baseline_libraries <- unique(c(locked_library, .Library.site, .Library))
  baseline_libraries <- baseline_libraries[!is.na(baseline_libraries) & nzchar(baseline_libraries) &
                                             dir.exists(baseline_libraries)]
  actual <- lock_versions_from_installed(packages, baseline_libraries)
  missing <- names(expected)[is.na(actual$versions)]
  mismatched <- names(expected)[!is.na(actual$versions) & actual$versions != expected]
  if (length(missing) || length(mismatched)) {
    issues <- c(if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
                if (length(mismatched)) paste0("version mismatch: ", paste(mismatched, collapse = ", ")))
    fail(paste("Locked project/base libraries do not match renv.lock;", paste(issues, collapse = "; ")))
  }
  list(packages = packages, expected = expected, actual = actual,
       project_package_count = nrow(utils::installed.packages(lib.loc = locked_library)),
       baseline_libraries = baseline_libraries)
}

expected_mapper_binaries <- list(
  list(package = "naturalsort", archive_version = "0.1.3", description_version = "0.1.3",
       package_version = "0.1.3", filename = "naturalsort_0.1.3.zip",
       sha256 = "aa835431b76b3e979ee2e99ae36e538dce53a3144cac0ecfc0360d1c6fa10259"),
  list(package = "castor", archive_version = "1.8.7", description_version = "1.8.7",
       package_version = "1.8.7", filename = "castor_1.8.7.zip",
       sha256 = "dbcc5d83d9bcf63cd0b4de9b94d6079553435a42616e01fffdb4d4085c5df847"),
  list(package = "RSpectra", archive_version = "0.16-2", description_version = "0.16-2",
       package_version = "0.16.2", filename = "RSpectra_0.16-2.zip",
       sha256 = "d3e22c9a67d04d834a772856fa5372f2b1c04a1f8ff6cea32f6e1669bfb211f1")
)

validate_optional_package <- function(package_record, library) {
  package_dir <- file.path(library, package_record$package)
  if (!file.exists(package_dir) && !dir.exists(package_dir)) return(FALSE)
  if (!dir.exists(package_dir)) {
    fail(sprintf("Existing %s entry in mapper_library is not a package directory; refusing to overwrite it.",
                 package_record$package))
  }
  canonical_package_dir <- normalizePath(package_dir, winslash = "/", mustWork = TRUE)
  if (!same_or_descendant(canonical_package_dir, library)) {
    fail(sprintf("Existing %s package folder resolves outside mapper_library.", package_record$package))
  }
  description <- file.path(package_dir, "DESCRIPTION")
  if (!file.exists(description) || dir.exists(description)) {
    fail(sprintf("Existing %s package folder has no valid DESCRIPTION; refusing to overwrite it.",
                 package_record$package))
  }
  fields <- tryCatch(read.dcf(description, fields = c("Package", "Version")),
                     error = function(e) fail(sprintf("Existing %s DESCRIPTION is invalid; refusing to overwrite it.",
                                                      package_record$package)))
  if (nrow(fields) != 1L || !identical(unname(fields[1L, "Package"]), package_record$package) ||
      !identical(unname(fields[1L, "Version"]), package_record$description_version)) {
    fail(sprintf("Existing %s package version differs from the pin; refusing to overwrite it.",
                 package_record$package))
  }
  installed <- utils::installed.packages(lib.loc = library)
  rows <- which(installed[, "Package"] == package_record$package)
  if (length(rows) != 1L || !identical(installed[rows, "Version"], package_record$description_version)) {
    fail(sprintf("Existing %s folder is not a valid installed package; refusing to overwrite it.",
                 package_record$package))
  }
  version <- tryCatch(as.character(utils::packageVersion(package_record$package, lib.loc = library)),
                      error = function(e) NA_character_)
  if (!identical(version, package_record$package_version)) {
    fail(sprintf("Installed %s package version differs from the pin; refusing to overwrite it.",
                 package_record$package))
  }
  TRUE
}

install_optional_package <- function(package_record, archive, library) {
  package_dir <- file.path(library, package_record$package)
  lock_dirs <- c(file.path(library, "00LOCK"), file.path(library, paste0("00LOCK-", package_record$package)))
  if (any(file.exists(lock_dirs) | dir.exists(lock_dirs))) {
    fail(sprintf("An installation lock exists for %s; remove or resolve it manually before setup.",
                 package_record$package))
  }
  present <- validate_optional_package(package_record, library)
  if (!present) {
    if (file.exists(package_dir) || dir.exists(package_dir)) {
      fail(sprintf("An invalid %s package folder already exists; refusing to overwrite it.",
                   package_record$package))
    }
    utils::install.packages(archive, repos = NULL, type = "win.binary", lib = library,
                            dependencies = FALSE, quiet = FALSE)
  }
  if (!validate_optional_package(package_record, library)) {
    fail(sprintf("Pinned %s binary did not install into mapper_library.", package_record$package))
  }
  invisible(!present)
}

main <- function() {
  original_lib_paths <- .libPaths()
  on.exit(.libPaths(original_lib_paths), add = TRUE)

  parsed <- parse_arguments(commandArgs(trailingOnly = TRUE))
  if (!identical(.Platform$OS.type, "windows") || !identical(as.character(getRversion()), "4.5.3") ||
      !identical(R.version[["platform"]], "x86_64-w64-mingw32")) {
    fail(sprintf("Requires native Windows R 4.5.3 x86_64-w64-mingw32; found %s (%s).",
                 R.version[["version.string"]], R.version[["platform"]]))
  }

  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) != 1L) fail("Could not determine the helper script path.")
  script_path <- normalizePath(sub("^--file=", "", file_argument), winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
  setwd(repo_root)
  options(wf16s.pipeline_root = repo_root, timeout = max(120, getOption("timeout")))

  default_locked <- file.path(repo_root, "renv", "library", "windows", "R-4.5",
                              "x86_64-w64-mingw32")
  locked_library <- resolve_existing_directory(default_if_null(parsed$locked_library, default_locked),
                                               repo_root, "Locked project library")
  default_mapper <- file.path(repo_root, "output", "faprotax2-runtime", "library")
  mapper_library <- canonical_candidate_directory(default_if_null(parsed$mapper_library, default_mapper),
                                                  repo_root, "Mapper library")

  # Keep package installation out of every active/project/base library, including descendants.
  protected <- unique(c(locked_library, .libPaths(), .Library.site, .Library))
  protected <- protected[!is.na(protected) & nzchar(protected) & dir.exists(protected)]
  protected <- vapply(protected, normalizePath, character(1), winslash = "/", mustWork = TRUE)
  if (any(vapply(protected, function(path) paths_overlap(mapper_library, path), logical(1)))) {
    fail("--mapper-library must be separate from the project, base, and active .libPaths() libraries.")
  }

  output_root <- canonical_candidate_directory(file.path(repo_root, "output"), repo_root,
                                               "Output directory")
  if (identical(path_key(mapper_library), path_key(output_root))) {
    fail("--mapper-library must name a dedicated library directory beneath output/, not output/ itself.")
  }
  if (same_or_descendant(mapper_library, repo_root) && !same_or_descendant(mapper_library, output_root)) {
    fail("A mapper library inside the checkout must be under output/; tracked source paths are not valid targets.")
  }
  if (same_or_descendant(repo_root, mapper_library)) {
    fail("--mapper-library cannot be an ancestor of the repository checkout.")
  }

  baseline_paths <- unique(c(locked_library, .Library.site, .Library))
  baseline_paths <- baseline_paths[!is.na(baseline_paths) & nzchar(baseline_paths) & dir.exists(baseline_paths)]
  .libPaths(baseline_paths)
  lock_state <- verify_locked_library(repo_root, locked_library)

  source(file.path(repo_root, "analysis", "utils", "io.R"), local = .GlobalEnv)
  source(file.path(repo_root, "analysis", "utils", "preflight.R"), local = .GlobalEnv)
  source(file.path(repo_root, "analysis", "09_faprotax2.R"), local = .GlobalEnv)

  upstream_root <- file.path(repo_root, "output", "faprotax2-audit", "upstream")
  archive_path <- file.path(upstream_root, "FAPROTAX-mapper_v2.0.tar.gz")
  archive_url <- paste0(
    "https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX2/",
    "SECTION_Download/MODULE_Downloads/CLASS_Latest%20release/",
    "UNIT_A_FAPROTAX_v2.0_mapper/FAPROTAX-mapper_v2.0.tar.gz"
  )
  mapper_archive <- ensure_cached_archive(archive_path, archive_url,
                                          FAPROTAX2_ARCHIVE_SHA256,
                                          "Official FAPROTAX-mapper 2.0 archive")
  unpacked_root <- file.path(upstream_root, "unpacked")
  mapper_dir <- file.path(unpacked_root, "FAPROTAX-mapper_v2.0")
  extract_mapper_archive(archive_path, unpacked_root, mapper_dir)

  cfg <- list(faprotax2 = list(
    mode = "taxonomy",
    mapper_path = file.path(mapper_dir, "faprotax.R"),
    reference_dir = file.path(mapper_dir, "FAPROTAX_ref_data"),
    mapper_library = mapper_library
  ))
  resources <- faprotax2_resource_inventory(cfg)
  if (length(resources) != 7L || !identical(names(resources), names(FAPROTAX2_RELEASE_HASHES))) {
    fail("Production resource validation did not return the mapper and six pinned reference assets.")
  }

  binary_cache <- file.path(repo_root, "output", "faprotax2-runtime", "downloads")
  binary_archives <- setNames(vector("list", length(expected_mapper_binaries)),
                              vapply(expected_mapper_binaries, `[[`, character(1), "package"))
  for (record in expected_mapper_binaries) {
    cache_path <- file.path(binary_cache, record$filename)
    url <- paste0("https://cloud.r-project.org/bin/windows/contrib/4.5/", record$filename)
    check <- ensure_cached_archive(cache_path, url, record$sha256,
                                   paste(record$package, record$archive_version, "Windows binary"))
    binary_archives[[record$package]] <- check$path
  }

  if (file.exists(mapper_library) && !dir.exists(mapper_library)) fail("Mapper library target is not a directory.")
  if (!dir.exists(mapper_library) && !dir.create(mapper_library, recursive = TRUE, showWarnings = FALSE)) {
    fail("Could not create the isolated mapper library.")
  }
  mapper_library <- normalizePath(mapper_library, winslash = "/", mustWork = TRUE)
  .libPaths(unique(c(mapper_library, locked_library, .Library.site, .Library)))

  installed_now <- character()
  for (record in expected_mapper_binaries) {
    did_install <- install_optional_package(record, binary_archives[[record$package]], mapper_library)
    if (isTRUE(did_install)) installed_now <- c(installed_now, record$package)
  }

  baseline_after <- lock_versions_from_installed(lock_state$packages, lock_state$baseline_libraries)
  baseline_mismatch <- names(lock_state$expected)[is.na(baseline_after$versions) |
                                                   baseline_after$versions != lock_state$expected]
  if (length(baseline_mismatch)) {
    fail("The locked project/base package versions changed during setup; inspect the local installation log.")
  }

  runtime <- validate_faprotax2_runtime(cfg)
  dependency_records <- runtime$dependency_runtime$packages
  dependency_names <- vapply(dependency_records, function(record) as.character(record$package), character(1))
  required <- c("castor", "naturalsort", "optparse", "stringr", "Rcpp", "Matrix",
                "RSpectra", "jsonlite", "digest")
  required_rows <- match(required, dependency_names)
  if (anyNA(required_rows)) fail("Production runtime validation omitted a required dependency version.")
  dependency_versions <- vapply(dependency_records[required_rows],
                                 function(record) as.character(record$version), character(1))
  names(dependency_versions) <- required
  if (!identical(unname(dependency_versions[["castor"]]), "1.8.7") ||
      !identical(unname(dependency_versions[["naturalsort"]]), "0.1.3") ||
      !identical(unname(dependency_versions[["RSpectra"]]), "0.16.2")) {
    fail("The production preflight resolved a mapper binary at an unexpected version.")
  }

  cat("FAPROTAX-mapper 2.0 native-Windows setup completed.\n")
  cat(sprintf("R: %s (%s); renv.lock packages checked: %d exact versions.\n",
              runtime$dependency_runtime$r, runtime$dependency_runtime$platform,
              length(lock_state$expected)))
  cat(sprintf("Locked project library entries: %d; binary overlay packages installed: %s.\n",
              lock_state$project_package_count,
              if (length(installed_now)) paste(installed_now, collapse = ", ") else "none (pins already present)"))
  cat("Production resource validator: PASS (mapper script plus six reference assets).\n")
  cat(sprintf("Official archive SHA-256: %s\n", mapper_archive$sha256))
  cat("Resolved mapper dependency versions:\n")
  for (package in required) cat(sprintf("  %s %s\n", package, dependency_versions[[package]]))
  cat("Resource SHA-256 values:\n")
  for (role in names(resources)) {
    record <- resources[[role]]
    cat(sprintf("  %s %s (%d bytes)\n", role, record$sha256, as.integer(record$size_bytes)))
  }
  cat("Review the upstream FAPROTAX download and license terms:\n")
  cat("  https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX2/lib/php/index.php?section=Download\n")
  cat("  https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX2/lib/php/index.php?section=License\n")
  invisible(NULL)
}

main()
