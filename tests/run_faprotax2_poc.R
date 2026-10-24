fail <- function(x) stop(x, call. = FALSE)
assert <- function(x, message) if (!isTRUE(x)) fail(message)
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) != 10L) {
  fail("Usage: Rscript tests/run_faprotax2_poc.R --mapper PATH --reference-dir PATH --mapper-library PATH --locked-library PATH --output-dir PATH")
}
args <- as.list(argv[seq.int(2L, 10L, by = 2L)])
names(args) <- substring(argv[seq.int(1L, 9L, by = 2L)], 3L)
required <- c("mapper", "reference-dir", "mapper-library", "locked-library", "output-dir")
if (!setequal(names(args), required)) fail("Provide all five named options exactly once")

repo <- normalizePath(getwd(), winslash = "\\", mustWork = TRUE)
mapper <- normalizePath(args[["mapper"]], winslash = "\\", mustWork = TRUE)
ref <- normalizePath(args[["reference-dir"]], winslash = "\\", mustWork = TRUE)
mapper_lib <- normalizePath(args[["mapper-library"]], winslash = "\\", mustWork = TRUE)
locked_lib <- normalizePath(args[["locked-library"]], winslash = "\\", mustWork = TRUE)
out <- normalizePath(args[["output-dir"]], winslash = "\\", mustWork = FALSE)
assert(.Platform$OS.type == "windows" && as.character(getRversion()) == "4.5.3",
       paste("Requires native Windows R 4.5.3; found", R.version.string))
assert(file.exists(mapper) && file.exists(file.path(ref, "default_config.tsv")),
       "Mapper script or official reference directory is missing")
assert(dir.exists(mapper_lib) && dir.exists(locked_lib), "Both R libraries must exist")
assert(grepl(" ", out, fixed = TRUE) && grepl("Rényi", out, fixed = TRUE),
       "Output path must contain spaces and Rényi")
assert(!file.exists(out), "Use a new output directory for each POC run")
dir.create(out, recursive = TRUE)
.libPaths(unique(c(mapper_lib, locked_lib, .Library.site, .Library)))
for (pkg in c("processx", "digest", "ps")) {
  assert(requireNamespace(pkg, quietly = TRUE), paste("Missing harness dependency:", pkg))
}

fixture_dir <- file.path(repo, "tests", "fixtures", "faprotax2")
input_path <- file.path(fixture_dir, "taxonomy_known_species.tsv")
unknown_path <- file.path(fixture_dir, "taxonomy_unknown_root.tsv")
invalid_path <- file.path(fixture_dir, "taxonomy_invalid_count.tsv")
input <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE,
                    quote = "", comment.char = "", fileEncoding = "UTF-8")
samples <- c("Rényi_α", "Control_β")
assert(identical(names(input), c("OTU", samples, "taxonomy")), "Fixture labels/schema changed")
assert(!anyDuplicated(input$OTU), "Fixture feature IDs must be unique")
counts <- as.matrix(input[, samples, drop = FALSE])
storage.mode(counts) <- "numeric"
rownames(counts) <- input$OTU
assert(all(is.finite(counts) & counts >= 0 & counts == floor(counts)),
       "Fixture counts must be nonnegative integers")

profile <- file.path(out, "mapper isolated library.Rprofile")
profile_lines <- sprintf('.libPaths(c("%s", "%s", .Library.site, .Library))',
                         normalizePath(mapper_lib, winslash = "/"),
                         normalizePath(locked_lib, winslash = "/"))
writeLines(profile_lines, profile, useBytes = TRUE)
child_env <- Sys.getenv()
child_env["RENV_CONFIG_AUTOLOADER_ENABLED"] <- "FALSE"
child_env["R_LIBS"] <- paste(mapper_lib, locked_lib, sep = .Platform$path.sep)
child_env["R_LIBS_USER"] <- mapper_lib
child_env["R_PROFILE_USER"] <- profile
child_env["LC_ALL"] <- "English_United States.utf8"
rscript <- file.path(R.home("bin"), "Rscript.exe")
preflight <- processx::run(
  rscript, c("--no-save", "--no-restore", "-e",
             "cat(R.version.string,'\\n'); for(p in c('castor','naturalsort','optparse','stringr')) cat(p,as.character(packageVersion(p)),'\\n',sep='\\t'); cat(paste(.libPaths(),collapse='\\n'))"),
  env = child_env, wd = repo, timeout = 30000, error_on_status = FALSE
)
writeLines(preflight$stdout, file.path(out, "runtime preflight.txt"), useBytes = TRUE)
writeLines(preflight$stderr, file.path(out, "runtime preflight stderr.txt"), useBytes = TRUE)
assert(identical(preflight$status, 0L), "Child R runtime/dependency preflight failed")
assert(grepl("castor\\t1\\.8\\.7", preflight$stdout), "Expected isolated castor 1.8.7")
assert(grepl("naturalsort\\t0\\.1\\.3", preflight$stdout), "Expected isolated naturalsort 0.1.3")

run_mapper <- function(label, table = input_path, failure = "error", only_function = NULL) {
  d <- file.path(out, label)
  dir.create(d)
  files <- c(
    function_table = file.path(d, "function abundance.tsv"),
    states = file.path(d, "feature function states.tsv"),
    probs = file.path(d, "annotation probabilities.tsv"),
    overlaps = file.path(d, "function overlaps.tsv"),
    report = file.path(d, "report.txt"),
    log = file.path(d, "mapper log.txt"),
    stdout = file.path(d, "stdout.txt"),
    stderr = file.path(d, "stderr.txt")
  )
  a <- c("--no-save", "--no-restore", mapper,
         "-d", ref, "-i", table, "-t", "column:taxonomy", "-k", "OTU",
         "-o", files["function_table"], "-e", files["states"], "-q", files["probs"],
         "-g", files["overlaps"], "-r", files["report"], "-l", files["log"],
         "--normalize", "none", "--hsp_algorithm", "config", "--unknown_ref_states", "config",
         "--placement_failure", failure, "--Nthreads", "1")
  if (!is.null(only_function)) a <- c(a, "--only_functions", only_function)
  writeLines(paste(c(rscript, a), collapse = "\n"), file.path(d, "argument vector.txt"), useBytes = TRUE)
  p <- processx::process$new(rscript, a, env = child_env, wd = repo,
                             stdout = files["stdout"], stderr = files["stderr"],
                             cleanup = TRUE, cleanup_tree = TRUE)
  start <- Sys.time()
  peak <- NA_real_
  repeat {
    if (!p$is_alive()) break
    memory <- tryCatch(ps::ps_memory_info(ps::ps_handle(p$get_pid())), error = function(e) NULL)
    if (!is.null(memory) && "rss" %in% names(memory))
      peak <- max(peak, unname(memory[["rss"]]), na.rm = TRUE)
    if (as.numeric(difftime(Sys.time(), start, units = "secs")) > 600) {
      p$kill()
      p$wait()
      fail(paste("Mapper timed out after 600 seconds:", label))
    }
    Sys.sleep(0.25)
  }
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  status <- p$get_exit_status()
  writeLines(c(paste("exit_status", status, sep = "\t"),
               paste("elapsed_seconds", sprintf("%.3f", elapsed), sep = "\t"),
               paste("peak_rss_bytes", ifelse(is.finite(peak), format(peak, scientific = FALSE), "unavailable"), sep = "\t")),
             file.path(d, "process metrics.txt"), useBytes = TRUE)
  list(dir = d, files = files, args = a, status = status, elapsed = elapsed, peak = peak)
}

read_table <- function(path, id_name) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "", fileEncoding = "UTF-8", na.strings = "NA")
  assert(nrow(d) > 0L && names(d)[1L] == id_name, paste("Unexpected schema:", path))
  ids <- d[[1L]]
  assert(!anyDuplicated(ids) && !anyNA(ids) && all(nzchar(ids)), paste("Duplicate/empty IDs:", path))
  m <- as.matrix(d[, -1L, drop = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- ids
  m
}
verify <- function(r) {
  assert(identical(r$status, 0L), paste("Mapper returned", r$status, "for", r$dir))
  for (f in r$files[1:6]) assert(file.exists(f) && file.info(f)$size > 0, paste("Missing output:", f))
  states <- read_table(r$files["states"], "OTU")
  probs <- read_table(r$files["probs"], "OTU")
  functions <- read_table(r$files["function_table"], "function")
  assert(identical(rownames(states), input$OTU), "Feature IDs/order changed")
  assert(identical(rownames(probs), rownames(states)) && identical(colnames(probs), colnames(states)),
         "Probability table is not synchronized to feature-state table")
  assert(identical(rownames(functions), colnames(states)), "Function catalog differs across mapper outputs")
  assert(identical(colnames(functions), samples), "Unicode sample labels changed")
  assert(all(is.na(states) | states %in% c(0, 1)), "States contain values outside 0/1/NA")
  assert(any(states == 0, na.rm = TRUE) && any(states == 1, na.rm = TRUE), "Did not observe both 0 and 1 states")
  assert(all(is.na(probs) | (is.finite(probs) & probs >= 0 & probs <= 1)),
         "Probabilities outside [0,1]")
  assert(!any((states == 1 & probs < 0.50000001) | (states == 0 & probs > 0.49999999), na.rm = TRUE),
         "State disagrees with configured probability thresholds")
  expected <- t(ifelse(is.na(states), 0, states)) %*% counts
  assert(isTRUE(all.equal(functions, expected, check.attributes = FALSE, tolerance = 0)),
         "Function abundance is not exactly state(NA->0) multiplied by input counts")
  report <- paste(readLines(r$files["report"], warn = FALSE), collapse = "\n")
  assert(grepl("Considered 101 FAPROTAX functions", report, fixed = TRUE), "Default function config not used")
  list(states = states, probs = probs, functions = functions)
}

one <- run_mapper("run one")
v1 <- verify(one)
two <- run_mapper("Rényi repeat two")
v2 <- verify(two)
assert(identical(v1$states, v2$states), "Repeated state output differs")
assert(identical(v1$probs, v2$probs), "Repeated probability output differs")
assert(identical(v1$functions, v2$functions), "Repeated function output differs")
hash <- function(r) vapply(r$files[c("function_table", "states", "probs", "overlaps")],
                           function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE), character(1))
h1 <- hash(one)
h2 <- hash(two)

unknown <- run_mapper("unknown taxonomy boundary", table = unknown_path,
                      only_function = "aerobic_chemoheterotrophy")
assert(identical(unknown$status, 0L), "Unknown taxonomy probe did not complete")
unknown_report <- paste(readLines(unknown$files["report"], warn = FALSE), collapse = "\n")
assert(grepl("Number of processed OTUs, i.e., after all OTU filtering: 1", unknown_report, fixed = TRUE),
       "Unknown taxonomy fixture was not retained as one processed feature")
unknown_states <- read_table(unknown$files["states"], "OTU")
assert(identical(rownames(unknown_states), "feature_unknown_taxonomy"),
       "Unknown taxonomy feature ID was not preserved")

negative <- run_mapper("invalid numeric input", table = invalid_path, failure = "error")
assert(!identical(negative$status, 0L), "Invalid numeric input unexpectedly succeeded")
assert(!file.exists(negative$files["function_table"]), "Failed mapper unexpectedly wrote a function table")
failure_text <- paste(c(readLines(negative$files["stdout"], warn = FALSE),
                        readLines(negative$files["stderr"], warn = FALSE)), collapse = "\n")
assert(nzchar(failure_text), "Mapper returned nonzero without preserving diagnostic output")

timeout_dir <- file.path(out, "timeout cleanup")
dir.create(timeout_dir)
marker <- file.path(timeout_dir, "completion marker.txt")
timeout_env <- child_env
timeout_env["POC_MARKER"] <- marker
slow <- processx::process$new(
  rscript, c("--no-save", "--no-restore", "-e",
             "Sys.sleep(30); writeLines('finished', Sys.getenv('POC_MARKER'))"),
  env = timeout_env, wd = repo, stdout = file.path(timeout_dir, "stdout.txt"),
  stderr = file.path(timeout_dir, "stderr.txt"), cleanup = TRUE, cleanup_tree = TRUE
)
Sys.sleep(0.5)
assert(slow$is_alive(), "Timeout probe exited before interruption")
slow$kill()
slow$wait()
assert(!slow$is_alive() && !file.exists(marker), "Killed child remained active or wrote completion marker")
writeLines("Killed and reaped before the 30-second completion marker.", file.path(timeout_dir, "result.txt"))

out_lines <- c(
  paste("R:", R.version.string),
  paste("mapper SHA-256:", digest::digest(file = mapper, algo = "sha256", serialize = FALSE)),
  paste("input SHA-256:", digest::digest(file = input_path, algo = "sha256", serialize = FALSE)),
  paste("features x samples:", nrow(counts), "x", ncol(counts)),
  paste("functions:", nrow(v1$functions)),
  paste("observed state cells present/absent/NA:",
        sum(v1$states == 1, na.rm = TRUE), sum(v1$states == 0, na.rm = TRUE), sum(is.na(v1$states))),
  paste("raw output table hashes identical across changed paths:", identical(h1, h2)),
  paste("unknown taxonomy feature exit / retained:", unknown$status,
        grepl("Number of processed OTUs, i.e., after all OTU filtering: 1", unknown_report, fixed = TRUE)),
  paste("invalid numeric input exit:", negative$status),
  paste("run one seconds / peak RSS bytes:", sprintf("%.3f", one$elapsed), ifelse(is.finite(one$peak), format(one$peak, scientific = FALSE), "unavailable")),
  paste("run two seconds / peak RSS bytes:", sprintf("%.3f", two$elapsed), ifelse(is.finite(two$peak), format(two$peak, scientific = FALSE), "unavailable")),
  "Repeat reused local mapper/reference; no installs or network inputs. Network access was not disabled."
)
writeLines(out_lines, file.path(out, "POC summary.txt"), useBytes = TRUE)
writeLines(capture.output(dput(list(run_one = h1, run_two = h2))),
           file.path(out, "table SHA-256.txt"), useBytes = TRUE)
cat(paste(out_lines, collapse = "\n"), "\n")
