# Independent contract tests for the official FAPROTAX-mapper 2.0 adapter.
source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "preflight.R"))
source(file.path("..", "..", "analysis", "utils", "manifest.R"))
source(file.path("..", "..", "analysis", "09_faprotax2.R"))

make_faprotax2_test_context <- function(output_dir, include_prokaryotes = TRUE) {
  samples <- c("Rényi_α", "Control_β")
  lineages <- rbind(
    c("Unclassified", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"),
    c("Bacteria", "Bacillati", "Pseudomonadota", "Gammaproteobacteria", "Nitrosomonadales",
      "Nitrosomonadaceae", "Nitrosomonas", "Nitrosomonas europaea"),
    c("Archaea", "Archaea", "Thermoproteota", "Nitrososphaeria", "Nitrososphaerales",
      "Nitrososphaeraceae", "Nitrososphaera", "Nitrososphaera viennensis"),
    c("Bacteria", "Bacillati", "Bacteroidota", "Bacteroidia", "Bacteroidales",
      "Bacteroidaceae", "Bacteroides", "Bacteroides testii"),
    c("Bacteria", "Bacillati", "Firmicutes", "Bacilli", "Bacillales",
      "Bacillaceae", "Bacillus", "Bacillus testii"),
    c("Bacteria", "Bacillati", "Firmicutes", "Bacilli", "Bacillales",
      "Bacillaceae", "Bacillus", "Bacillus zero"),
    c("Eukaryota", "Opisthokonta", "Ascomycota", "Saccharomycetes", "Saccharomycetales",
      "Saccharomycetaceae", "Saccharomyces", "Saccharomyces cerevisiae")
  )
  colnames(lineages) <- RANKS_8
  if (!include_prokaryotes) {
    lineages <- lineages[c(1L, 6L, 7L), , drop = FALSE]
  }
  paths <- apply(lineages, 1L, paste, collapse = ";")
  rownames(lineages) <- paths
  counts <- if (include_prokaryotes) {
    rbind(c(2, 1), c(4, 0), c(0, 3), c(2, 1), c(1, 2), c(0, 0), c(9, 1))
  } else {
    rbind(c(2, 1), c(0, 0), c(9, 1))
  }
  rownames(counts) <- paths
  colnames(counts) <- samples
  taxonomy <- as.data.frame(lineages, stringsAsFactors = FALSE)
  taxonomy$TaxonPath <- paths

  list(
    config = list(
      output = list(dirs = list(faprotax2 = output_dir)),
      faprotax2 = list(
        mode = "taxonomy", mapper_path = "mapper path with spaces/faprotax.R",
        reference_dir = "reference path with spaces", mapper_library = NULL,
        nthreads = 1L, timeout_seconds = 10L
      )
    ),
    samples = samples,
    count_matrix = counts,
    taxonomy = taxonomy,
    unclass_index = 1L,
    count_basis = "minimap2_alignments"
  )
}

make_faprotax2_matrix_fixture <- function() {
  states <- rbind(
    feature_a = c(1, 0),
    feature_b = c(NA_real_, 1),
    feature_c = c(NA_real_, 0),
    feature_d = c(0, 0)
  )
  colnames(states) <- c("nitrification", "methanogenesis")
  probabilities <- rbind(
    feature_a = c(0.9, 0.1),
    feature_b = c(NA_real_, 0.8),
    feature_c = c(NA_real_, 0.2),
    feature_d = c(0.1, 0.3)
  )
  dimnames(probabilities) <- dimnames(states)
  counts <- rbind(feature_a = c(4, 0), feature_b = c(0, 3),
                  feature_c = c(2, 1), feature_d = c(1, 2))
  colnames(counts) <- c("Rényi_α", "Control_β")
  binary <- states
  binary[is.na(binary)] <- 0
  abundance <- t(binary) %*% counts
  list(states = states, probabilities = probabilities,
       abundance = abundance, counts = counts)
}

test_that("taxonomy preparation uses seven mapper ranks and only positive-count prokaryotes", {
  prepared <- prepare_faprotax2_taxonomy(make_faprotax2_test_context(tempdir()))

  expect_identical(prepared$feature_ids, c("Taxon_000002", "Taxon_000003", "Taxon_000004", "Taxon_000005"))
  expect_identical(rownames(prepared$counts), prepared$feature_ids)
  expect_identical(unname(prepared$taxonomy[[1L]]),
    "Bacteria;Pseudomonadota;Gammaproteobacteria;Nitrosomonadales;Nitrosomonadaceae;Nitrosomonas;Nitrosomonas europaea")
  expect_identical(unname(prepared$taxonomy[[2L]]),
    "Archaea;Thermoproteota;Nitrososphaeria;Nitrososphaerales;Nitrososphaeraceae;Nitrososphaera;Nitrososphaera viennensis")
  expect_false(any(grepl("Bacillati", prepared$taxonomy, fixed = TRUE)))
  expect_false(any(grepl("Eukaryota", prepared$taxonomy, fixed = TRUE)))
  expect_false("Taxon_000006" %in% prepared$feature_ids)
  expect_identical(colnames(prepared$counts), c("Rényi_α", "Control_β"))
})

test_that("taxonomy preparation turns unknown ranks into empty mapper fields", {
  context <- make_faprotax2_test_context(tempdir())
  context$taxonomy$genus[2L] <- "Unknown"
  context$taxonomy$species[2L] <- NA_character_

  prepared <- prepare_faprotax2_taxonomy(context)
  expect_identical(unname(prepared$taxonomy[[1L]]),
    "Bacteria;Pseudomonadota;Gammaproteobacteria;Nitrosomonadales;Nitrosomonadaceae;;")
})

test_that("taxonomy preparation returns no rows when positive-count prokaryotes are absent", {
  prepared <- prepare_faprotax2_taxonomy(
    make_faprotax2_test_context(tempdir(), include_prokaryotes = FALSE))
  expect_null(prepared)
})

test_that("mapper matrices preserve NA states and verify unweighted abundance by direct multiplication", {
  fixture <- make_faprotax2_matrix_fixture()
  expect_invisible(validate_faprotax2_matrices(
    fixture$states, fixture$probabilities, fixture$abundance, fixture$counts))

  fractional_counts <- fixture$counts
  fractional_counts["feature_a", "Rényi_α"] <- 2.5
  binary <- fixture$states
  binary[is.na(binary)] <- 0
  fractional_abundance <- t(binary) %*% fractional_counts
  expect_invisible(validate_faprotax2_matrices(
    fixture$states, fixture$probabilities, fractional_abundance, fractional_counts))
  expect_false(all(fractional_abundance == floor(fractional_abundance)))
})

test_that("mapper matrix validation rejects invalid states, probabilities and abundance arithmetic", {
  fixture <- make_faprotax2_matrix_fixture()

  altered <- fixture
  altered$states[1L, 1L] <- 0.5
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "invalid states")

  altered <- fixture
  altered$probabilities[1L, 1L] <- 1.01
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "invalid states, probabilities")

  altered <- fixture
  altered$probabilities[1L, 1L] <- NA_real_
  expect_invisible(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts))

  altered <- fixture
  altered$abundance[1L, 1L] <- altered$abundance[1L, 1L] + 0.25
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "differs from unweighted states")
})

test_that("mapper matrix validation rejects feature, function and sample identifier mismatches", {
  fixture <- make_faprotax2_matrix_fixture()

  altered <- fixture
  rownames(altered$counts) <- rev(rownames(altered$counts))
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "identifiers/order disagree")

  altered <- fixture
  rownames(altered$probabilities) <- rev(rownames(altered$probabilities))
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "identifiers/order disagree")

  altered <- fixture
  rownames(altered$abundance) <- rev(rownames(altered$abundance))
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "identifiers/order disagree")

  altered <- fixture
  colnames(altered$abundance) <- rev(colnames(altered$abundance))
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "identifiers/order disagree")

  altered <- fixture
  colnames(altered$counts) <- c("Rényi_α", "Rényi_α")
  colnames(altered$abundance) <- c("Rényi_α", "Rényi_α")
  expect_error(validate_faprotax2_matrices(altered$states, altered$probabilities,
    altered$abundance, altered$counts), "unique sample")
})

test_that("coverage reports input count basis and separates predicted presence from unknowns", {
  fixture <- make_faprotax2_matrix_fixture()
  context <- list(count_basis = "bracken_species_estimates")
  coverage <- summarize_faprotax2_coverage(fixture$states, fixture$counts, context)

  expect_identical(names(coverage), c(
    "SampleID", "EligibleTaxonAbundance", "PredictedPresenceAbundance",
    "NoPredictedPresenceAbundance", "UnknownWithoutPresenceAbundance",
    "KnownAbsenceOnlyAbundance", "EligibleTaxa", "TaxaWithPredictedPresence", "CountBasis"
  ))
  expect_false(any(grepl("placement|mapped_reads|functionmappedreads", names(coverage), ignore.case = TRUE)))
  expect_identical(coverage$SampleID, c("Rényi_α", "Control_β"))
  expect_equal(coverage$EligibleTaxonAbundance, c(7, 6))
  expect_equal(coverage$PredictedPresenceAbundance, c(4, 3))
  expect_equal(coverage$NoPredictedPresenceAbundance, c(3, 3))
  expect_equal(coverage$UnknownWithoutPresenceAbundance, c(2, 1))
  expect_equal(coverage$KnownAbsenceOnlyAbundance, c(1, 2))
  expect_equal(coverage$EligibleTaxa, c(3, 3))
  expect_equal(coverage$TaxaWithPredictedPresence, c(1, 1))
  expect_identical(coverage$CountBasis, rep("bracken_species_estimates", 2L))
})

test_that("matrix reader preserves Unicode names and rejects malformed identifiers", {
  path <- tempfile(fileext = ".tsv")
  writeLines(c("OTU\tRényi_α\tControl_β", "feature_a\t1\tNA", "feature_b\t0\t2"),
             path, useBytes = TRUE)
  matrix <- read_faprotax2_matrix(path, "OTU")
  expect_identical(rownames(matrix), c("feature_a", "feature_b"))
  expect_identical(colnames(matrix), c("Rényi_α", "Control_β"))
  expect_equal(matrix["feature_b", "Rényi_α"], 0)
  expect_true(is.na(matrix["feature_a", "Control_β"]))

  writeLines(c("wrong_id\tS1", "feature_a\t1"), path)
  expect_error(read_faprotax2_matrix(path, "OTU"), "invalid matrix schema")
  writeLines(c("OTU\tS1", "feature_a\t1", "feature_a\t2"), path)
  expect_error(read_faprotax2_matrix(path, "OTU"), "invalid matrix schema")
  writeLines(c("OTU\tS1", "feature_a\tnot-a-number"), path)
  expect_error(read_faprotax2_matrix(path, "OTU"), "nonnumeric column")
})

test_that("no eligible prokaryotes produce a structured skip without launching the mapper", {
  context <- make_faprotax2_test_context(
    file.path(tempdir(), "faprotax2 no prokaryotes"), include_prokaryotes = FALSE)
  called <- FALSE
  runner <- function(command, args, ...) {
    called <<- TRUE
    stop("mapper must not run for an inapplicable input")
  }
  runtime <- list(rscript = "unused", env = character(), resources = list(),
                  dependency_runtime = list())

  result <- run_faprotax2(context, runner = runner, runtime = runtime)
  expect_identical(result$status, "skipped")
  expect_identical(result$reason, "No eligible prokaryotic taxa")
  expect_true(file.exists(result$outputs))
  expect_false(called)
  skipped <- read.delim(result$outputs, check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(skipped$Status, "skipped")
  expect_identical(skipped$Reason, "No positive-count Bacteria or Archaea rows")
})

test_that("taxonomy adapter passes argument vectors and validates independent mapper outputs", {
  output_dir <- file.path(tempdir(), "faprotax2 output Rényi with spaces")
  context <- make_faprotax2_test_context(output_dir)
  expected_ids <- c("Taxon_000002", "Taxon_000003", "Taxon_000004", "Taxon_000005")
  states <- rbind(c(1, 1), c(NA_real_, 1), c(NA_real_, 0), c(0, 0))
  rownames(states) <- expected_ids
  colnames(states) <- c("nitrification", "methanogenesis")
  probabilities <- rbind(
    c(0.9, 0.9), c(NA_real_, 0.8), c(NA_real_, 0.2), c(0.1, 0.3)
  )
  rownames(probabilities) <- expected_ids
  colnames(probabilities) <- colnames(states)
  input_counts <- prepare_faprotax2_taxonomy(context)$counts
  binary <- states
  binary[is.na(binary)] <- 0
  expected_abundance <- t(binary) %*% input_counts
  observed <- list()
  runtime <- list(
    rscript = "Rscript path with spaces/Rscript.exe", env = c(R_LIBS = "isolated library"),
    resources = list(mapper_script = list(path = "official mapper", sha256 = "fixture-hash")),
    dependency_runtime = list(r = "R 4.5.3", packages = list())
  )
  runner <- function(command, args, env, timeout, error_on_status, cleanup_tree) {
    observed$command <<- command
    observed$args <<- args
    observed$env <<- env
    observed$timeout <<- timeout
    observed$error_on_status <<- error_on_status
    observed$cleanup_tree <<- cleanup_tree
    expect_identical(command, runtime$rscript)
    expect_true(all(c("--vanilla", "--normalize", "none", "--Nthreads", "1") %in% args))
    expect_identical(env, runtime$env)

    value_after <- function(flag) args[match(flag, args) + 1L]
    input <- read.delim(value_after("-i"), check.names = FALSE, stringsAsFactors = FALSE)
    expect_identical(input$OTU, expected_ids)
    expect_identical(names(input), c("OTU", "Rényi_α", "Control_β", "taxonomy"))
    expect_false(any(grepl("Bacillati", input$taxonomy, fixed = TRUE)))

    write.table(data.frame(OTU = expected_ids, states, check.names = FALSE), value_after("-e"),
      sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data.frame(OTU = expected_ids, probabilities, check.names = FALSE), value_after("-q"),
      sep = "\t", row.names = FALSE, quote = FALSE)
    function_table <- data.frame(rownames(expected_abundance), expected_abundance,
                                 check.names = FALSE)
    names(function_table)[1L] <- "function"
    write.table(function_table, value_after("-o"), sep = "\t", row.names = FALSE, quote = FALSE)
    writeLines("fake mapper report", value_after("-r"))
    list(status = 0L, stdout = "mapper stdout", stderr = "")
  }

  result <- run_faprotax2(context, runner = runner, runtime = runtime)
  expect_identical(result$status, "completed")
  expect_true(all(file.exists(result$outputs)))
  expect_true(observed$timeout == 10000)
  expect_identical(observed$error_on_status, FALSE)
  expect_identical(observed$cleanup_tree, TRUE)

  function_abundance <- read.delim(file.path(output_dir, "faprotax2_function_abundance.tsv"),
                                   check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(function_abundance$Function,
                   rep(c("nitrification", "methanogenesis"), times = 2L))
  expect_identical(function_abundance$SampleID,
                   rep(c("Rényi_α", "Control_β"), each = 2L))
  expect_equal(function_abundance$EstimatedFunctionAbundance, c(4, 4, 0, 3))
  expect_gt(sum(function_abundance$EstimatedFunctionAbundance[
    function_abundance$SampleID == "Rényi_α"]), 7)
  expect_equal(function_abundance$EstimatedFunctionAbundance[
    function_abundance$SampleID == "Rényi_α"], c(4, 4))
  expect_identical(function_abundance$CountBasis,
                   rep("minimap2_alignments", 4L))
  expect_true(all(grepl("sum_of_input_taxon_abundances", function_abundance$Unit, fixed = TRUE)))

  coverage <- read.delim(file.path(output_dir, "faprotax2_coverage.tsv"),
                         check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(coverage$EligibleTaxonAbundance, c(7, 6))
  expect_equal(coverage$PredictedPresenceAbundance, c(4, 3))
  expect_equal(coverage$UnknownWithoutPresenceAbundance, c(2, 1))

  provenance <- jsonlite::fromJSON(file.path(output_dir, "faprotax2_provenance.json"),
                                   simplifyVector = FALSE)
  expect_identical(provenance$mapper_modified, FALSE)
  expect_identical(provenance$probability_weighting, FALSE)
  expect_identical(provenance$count_basis, "minimap2_alignments")
  expect_identical(provenance$archive_sha256, FAPROTAX2_ARCHIVE_SHA256)
})

test_that("a nonzero mapper exit is surfaced with its stderr", {
  context <- make_faprotax2_test_context(
    file.path(tempdir(), "faprotax2 failed process"))
  runtime <- list(rscript = "Rscript", env = character(), resources = list(),
                  dependency_runtime = list())
  runner <- function(command, args, ...) {
    list(status = 1L, stdout = "partial stdout", stderr = "intentional fixture failure")
  }

  expect_error(run_faprotax2(context, runner = runner, runtime = runtime),
               "exited 1: intentional fixture failure")
  expect_true(file.exists(file.path(context$config$output$dirs$faprotax2, "mapper_stderr.txt")))
  expect_false(file.exists(file.path(context$config$output$dirs$faprotax2,
                                     "faprotax2_function_abundance.tsv")))
})
