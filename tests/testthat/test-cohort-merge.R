# =============================================================================
# Unit Tests: wf-16s Cohort Merger Utility
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "classifier_contract.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "atomic_io.R"))
source(file.path("..", "..", "analysis", "utils", "merge_wf16s_cohort.R"))

make_synthetic_sample_bundle <- function(dir, sample_id, taxa_counts, unclass_count = 10,
                                         classifier = "kraken2", database_set = "ncbi_16s_18s_28s_ITS",
                                         min_len = 250, max_len = 1000) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- file.path(dir, "output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # 1. abundance_table_species.tsv
  unclass_row <- "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown"
  taxa <- names(taxa_counts)
  all_tax <- c(unclass_row, taxa)
  all_counts <- c(unclass_count, as.numeric(taxa_counts))

  ab_df <- data.frame(
    tax = all_tax,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  ab_df[[sample_id]] <- sprintf("%.1f", all_counts)
  ab_df[["total"]] <- sprintf("%.1f", all_counts)

  utils::write.table(ab_df, file.path(out_dir, "abundance_table_species.tsv"),
                     sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  # 2. params.json
  params_data <- list(
    sample = sample_id,
    classifier = classifier,
    database_set = database_set,
    taxonomic_rank = "S",
    min_len = min_len,
    max_len = max_len,
    min_read_qual = 15,
    abundance_threshold = 0,
    database_sets = load_upstream_database_contracts()
  )
  if (classifier == "kraken2") {
    params_data$kraken2_confidence <- 0
    params_data$bracken_threshold <- 10
    params_data$bracken_length <- 1000L
  } else {
    params_data$min_percent_identity <- 95
    params_data$min_ref_coverage <- 90
  }
  jsonlite::write_json(params_data, file.path(out_dir, "params.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null")

  # Real producer-shaped fixtures: assignments, raw report and estimates agree.
  assign_dir <- file.path(out_dir, "reads_assignments")
  dir.create(assign_dir, recursive = TRUE, showWarnings = FALSE)
  tids <- as.character(100L + match(sub(".*;", "", taxa), c("TaxonA", "TaxonB", "TaxonC")))
  lineage <- vapply(strsplit(taxa, ";", fixed = TRUE), function(x) paste(x[c(1L, 3L:8L)], collapse = "|"), character(1))
  assignment_rows <- character(0)
  for (i in seq_along(taxa)) {
    for (j in seq_len(taxa_counts[i])) {
      assignment_rows <- c(assignment_rows, if (classifier == "kraken2")
        paste("C", paste0("r", length(assignment_rows)+1), tids[i], "500", paste0(tids[i], ":1"), lineage[i], sep="\t") else
        paste("C", paste0("r", length(assignment_rows)+1), tids[i], "500", lineage[i], sep="\t"))
    }
  }
  for (j in seq_len(unclass_count)) assignment_rows <- c(assignment_rows, if (classifier == "kraken2")
    paste("U", paste0("u", j), "0", "500", "0:1", "", sep="\t") else
    paste("U", paste0("u", j), "0", "500", "", sep="\t"))
  writeLines(assignment_rows, file.path(assign_dir, sprintf("%s.assignments.tsv", sample_id)))
  if (classifier == "kraken2") {
    total <- sum(all_counts)
    classified <- sum(taxa_counts)
    report_row <- function(n, direct, rank, tid, name) paste(sprintf("%.2f", 100*n/total), n, direct, rank, tid, name, sep="\t")
    rows <- c(report_row(unclass_count, unclass_count, "U", "0", "unclassified"),
              report_row(classified, 0, "R", "1", "root"))
    for (i in seq_along(taxa)) rows <- c(rows, report_row(taxa_counts[i], taxa_counts[i], "S", tids[i], paste0("  ", sub(".*;", "", taxa[i]))))
    k_dir <- file.path(out_dir, "kraken2")
    dir.create(k_dir, recursive = TRUE)
    writeLines(rows, file.path(k_dir, sprintf("%s.kraken2.report.txt", sample_id)))
    b_dir <- file.path(out_dir, "bracken")
    dir.create(b_dir, recursive = TRUE)
    header <- "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
    rows <- vapply(seq_along(taxa), function(i) paste(sub(".*;", "", taxa[i]), tids[i], "S", taxa_counts[i], 0, taxa_counts[i], taxa_counts[i]/classified, sep="\t"), character(1))
    writeLines(c(header, rows), file.path(b_dir, sprintf("%s.bracken.report", sample_id)))
  }

  out_dir
}

run_cohort_merge_child <- function(repo_root, bundles, output_dir, crash_phase = "",
                                   lock_timeout_ms = 10000L) {
  launcher <- tempfile("cohort_merge_child_", fileext = ".R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "source(file.path(args[[1]], 'analysis', 'utils', 'merge_wf16s_cohort.R'))",
    paste0("merge_wf16s_cohort(args[2:3], args[[4]], overwrite = TRUE, ",
           "lock_timeout_ms = as.integer(args[[6]]), "),
    "  publication_hook = function(phase, destination, stage, backup) {",
    "    if (nzchar(args[[5]]) && identical(phase, args[[5]])) quit(save = 'no', status = 86L)",
    "  })"
  ), launcher)
  child_env <- Sys.getenv()
  child_env[["R_LIBS_USER"]] <- .libPaths()[1]
  child_env[["RENV_CONFIG_AUTOLOADER_ENABLED"]] <- "FALSE"
  processx::run(Sys.which("Rscript"),
    c(launcher, repo_root, bundles[[1]], bundles[[2]], output_dir,
      crash_phase, as.character(lock_timeout_ms)),
    wd = repo_root, env = child_env, timeout = 120000,
    error_on_status = FALSE, echo = FALSE
  )
}

test_that("merge_wf16s_cohort rejects fewer than 2 inputs", {
  expect_error(
    merge_wf16s_cohort(input_dirs = c("dir1"), output_dir = "out"),
    "E_COHORT_INSUFFICIENT_INPUTS"
  )
})

test_that("merge_wf16s_cohort rejects non-existent input directory", {
  expect_error(
    merge_wf16s_cohort(input_dirs = c("nonexistent_dir1", "nonexistent_dir2"), output_dir = "out"),
    "E_COHORT_INPUT_NOT_FOUND"
  )
})

test_that("merge_wf16s_cohort rejects duplicate source directories", {
  root <- tempfile("cohort_test_dup_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle1), output_dir = file.path(root, "merged")),
    "E_COHORT_DUPLICATE_SOURCE"
  )
})

test_that("merge_wf16s_cohort rejects output directory overlapping with input", {
  root <- tempfile("cohort_test_overlap_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(bundle1, "nested_out")),
    "E_COHORT_OUTPUT_OVERLAP"
  )
})

test_that("merge destinations resolve existing ancestors before overlap checks", {
  root <- tempfile("cohort_destination_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  destination <- file.path(bundle1, "missing", "nested_out")
  if (.Platform$OS.type == "windows") {
    destination <- gsub("/", "\\\\", file.path(shortPathName(bundle1), "missing", "nested_out"))
  }
  expect_identical(canonicalize_merge_destination(destination),
                   file.path(normalizePath(bundle1, winslash = "/"), "missing", "nested_out"))
  expect_error(merge_wf16s_cohort(c(bundle1, bundle2), destination, validate_only = TRUE),
               "E_COHORT_OUTPUT_OVERLAP")
  expect_false(dir.exists(file.path(bundle1, "missing")))
})

test_that("merge_wf16s_cohort rejects colliding sample IDs", {
  root <- tempfile("cohort_test_same_sample_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "SameSampleID", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "SameSampleID", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(root, "merged")),
    "E_COHORT_DUPLICATE_SAMPLE"
  )
})

test_that("merge_wf16s_cohort rejects mixed classifiers across inputs", {
  root <- tempfile("cohort_test_mixed_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10), classifier = "kraken2")
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20), classifier = "minimap2")

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(root, "merged")),
    "E_COHORT_INCOMPATIBLE"
  )
})

test_that("merge_wf16s_cohort rejects threshold mismatch across inputs", {
  root <- tempfile("cohort_test_thresh_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10), min_len = 250)
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20), min_len = 300)

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(root, "merged")),
    "E_COHORT_INCOMPATIBLE"
  )
})

test_that("merge_wf16s_cohort --validate-only succeeds without mutating filesystem", {
  root <- tempfile("cohort_test_val_only_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  out_dir <- file.path(root, "merged_val_only")

  res <- merge_wf16s_cohort(
    input_dirs = c(bundle1, bundle2),
    output_dir = out_dir,
    validate_only = TRUE
  )
  expect_identical(res$status, "validated")
  expect_false(dir.exists(out_dir))
})

test_that("merge_wf16s_cohort performs full outer join, neutral metadata, config, and hashed manifest", {
  root <- tempfile("cohort_test_merge_")
  dir.create(root)
  # Sample 1: TaxonA (10), TaxonB (20), Unclassified (5)
  # Sample 2: TaxonB (30), TaxonC (40), Unclassified (8)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1",
                                         c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10, "Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20), unclass_count = 5)
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2",
                                         c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 30, "Eukaryota;Fungi;P;C;O;F;G;TaxonC" = 40), unclass_count = 8)
  out_dir <- file.path(root, "merged_cohort")

  res <- merge_wf16s_cohort(
    input_dirs = c(bundle1, bundle2),
    output_dir = out_dir,
    project_name = "Synthetic_ITS_Cohort"
  )
  expect_identical(res$status, "merged")
  expect_true(dir.exists(out_dir))

  # 1. Verify abundance_table_species.tsv
  ab_path <- file.path(out_dir, "abundance_table_species.tsv")
  expect_true(file.exists(ab_path))
  ab_df <- utils::read.delim(ab_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(names(ab_df), c("tax", "Sample_1", "Sample_2", "total"))
  expect_equal(nrow(ab_df), 4) # Unclassified + TaxonA + TaxonB + TaxonC

  # Row 1 must be Unclassified
  expect_true(startsWith(ab_df$tax[1], "Unclassified"))
  expect_equal(as.numeric(ab_df$Sample_1[1]), 5.0)
  expect_equal(as.numeric(ab_df$Sample_2[1]), 8.0)
  expect_equal(as.numeric(ab_df$total[1]), 13.0)

  # TaxonA: present in S1 (10), absent in S2 (0)
  row_a <- ab_df[ab_df$tax == "Eukaryota;Fungi;P;C;O;F;G;TaxonA", ]
  expect_equal(nrow(row_a), 1)
  expect_equal(as.numeric(row_a$Sample_1), 10.0)
  expect_equal(as.numeric(row_a$Sample_2), 0.0)
  expect_equal(as.numeric(row_a$total), 10.0)

  # TaxonB: present in both S1 (20) and S2 (30)
  row_b <- ab_df[ab_df$tax == "Eukaryota;Fungi;P;C;O;F;G;TaxonB", ]
  expect_equal(nrow(row_b), 1)
  expect_equal(as.numeric(row_b$Sample_1), 20.0)
  expect_equal(as.numeric(row_b$Sample_2), 30.0)
  expect_equal(as.numeric(row_b$total), 50.0)

  # TaxonC: absent in S1 (0), present in S2 (40)
  row_c <- ab_df[ab_df$tax == "Eukaryota;Fungi;P;C;O;F;G;TaxonC", ]
  expect_equal(nrow(row_c), 1)
  expect_equal(as.numeric(row_c$Sample_1), 0.0)
  expect_equal(as.numeric(row_c$Sample_2), 40.0)
  expect_equal(as.numeric(row_c$total), 40.0)

  # 2. Verify metadata.tsv
  meta_path <- file.path(out_dir, "metadata.tsv")
  expect_true(file.exists(meta_path))
  meta_df <- utils::read.delim(meta_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  expect_identical(names(meta_df), c("SampleID", "Group"))
  expect_identical(meta_df$SampleID, c("Sample_1", "Sample_2"))
  expect_identical(meta_df$Group, c("Unspecified", "Unspecified"))

  # 3. Verify config.yml
  cfg_path <- file.path(out_dir, "config.yml")
  expect_true(file.exists(cfg_path))
  cfg <- yaml::yaml.load_file(cfg_path)
  expect_identical(cfg$project_name, "Synthetic_ITS_Cohort")
  expect_identical(cfg$mode, "cohort")
  expect_identical(cfg$input$abundance_table, "abundance_table_species.tsv")
  expect_identical(cfg$input$metadata, "metadata.tsv")
  expect_named(cfg$input$params_by_sample, c("Sample_1", "Sample_2"))
  expect_named(cfg$input$kraken2_reports, c("Sample_1", "Sample_2"))
  expect_named(cfg$input$bracken_reports, c("Sample_1", "Sample_2"))
  expect_equal(cfg$qc$display_min_length, 250)
  expect_equal(cfg$qc$display_max_length, 1000)

  # 4. Verify cohort_merge_manifest.json
  manifest_path <- file.path(out_dir, "cohort_merge_manifest.json")
  expect_true(file.exists(manifest_path))
  mf <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  expect_identical(mf$manifest_version, "1.0")
  expect_identical(mf$project_name, "Synthetic_ITS_Cohort")
  expect_identical(mf$compatibility$status, "compatible")
  expect_identical(mf$compatibility$classifier, "kraken2")
  expect_length(mf$sources, 2)
  expect_true("abundance_table_species.tsv" %in% names(mf$outputs))
  expect_true("metadata.tsv" %in% names(mf$outputs))
  expect_true("config.yml" %in% names(mf$outputs))
})

test_that("merge_wf16s_cohort refuses existing destination without overwrite and succeeds with overwrite", {
  root <- tempfile("cohort_test_overwrite_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  out_dir <- file.path(root, "merged_dest")

  merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = out_dir)
  expect_true(dir.exists(out_dir))

  # Refuses without overwrite
  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = out_dir, overwrite = FALSE),
    "E_COHORT_DESTINATION_EXISTS"
  )

  # Overwrites cleanly with overwrite = TRUE
  res <- merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = out_dir, overwrite = TRUE)
  expect_identical(res$status, "merged")
})

test_that("parse_merge_cli_args parses repeated and comma-separated options", {
  args1 <- c("-i", "dir1", "-i", "dir2", "-o", "my_out", "-p", "MyProject", "--overwrite")
  parsed1 <- parse_merge_cli_args(args1)
  expect_identical(parsed1$inputs, c("dir1", "dir2"))
  expect_identical(parsed1$output_dir, "my_out")
  expect_identical(parsed1$project_name, "MyProject")
  expect_true(parsed1$overwrite)
  expect_false(parsed1$validate_only)

  args2 <- c("--inputs=dirA,dirB", "--output-dir=out2", "--validate-only")
  parsed2 <- parse_merge_cli_args(args2)
  expect_identical(parsed2$inputs, c("dirA", "dirB"))
  expect_identical(parsed2$output_dir, "out2")
  expect_true(parsed2$validate_only)
})

test_that("merger validates report accounting before validate-only success", {
  root <- tempfile(); dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root,"a"),"A",c("Eukaryota;Fungi;P;C;O;F;G;TaxonA"=10))
  b <- make_synthetic_sample_bundle(file.path(root,"b"),"B",c("Eukaryota;Fungi;P;C;O;F;G;TaxonB"=20))
  report <- file.path(a,"kraken2","A.kraken2.report.txt")
  writeLines("100.00\t10\t10\tU\t0\tunclassified",report)
  out <- file.path(root,"merged")
  expect_error(merge_wf16s_cohort(c(a,b),out,validate_only=TRUE),"Accounting failure")
  expect_false(dir.exists(out))
})

test_that("merger protects unrelated files and honors actual configured length bounds", {
  root <- tempfile(); dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root,"a"),"A",c("Eukaryota;Fungi;P;C;O;F;G;TaxonA"=10),min_len=300,max_len=900)
  b <- make_synthetic_sample_bundle(file.path(root,"b"),"B",c("Eukaryota;Fungi;P;C;O;F;G;TaxonB"=20),min_len=300,max_len=900)
  out <- file.path(root,"merged")
  merge_wf16s_cohort(c(a,b),out)
  cfg <- yaml::read_yaml(file.path(out,"config.yml"))
  expect_equal(cfg$qc$display_min_length,300)
  expect_equal(cfg$qc$display_max_length,900)
  cfg <- load_config(file.path(out,"config.yml"))
  expect_length(validate_cohort_merge_input(cfg),2)
  source_table <- file.path(a,"abundance_table_species.tsv")
  writeLines(c(readLines(source_table),""),source_table)
  expect_error(validate_cohort_merge_input(cfg),"E_COHORT_PROVENANCE")
  unrelated <- file.path(out,"notes.txt"); writeLines("keep",unrelated)
  expect_error(merge_wf16s_cohort(c(a,b),out,overwrite=TRUE),"E_COHORT_DESTINATION_UNOWNED")
  expect_equal(readLines(unrelated),"keep")
  expect_null(locate_kraken_report_file(a,"unrelated_sample"))
})

test_that("portable cache names cannot collide with generated bundle artifacts", {
  root <- tempfile("cohort_cache_name_identity_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  reserved_variants <- c("CONFIG.YML", "Metadata.tsv", "ABUNDANCE_TABLE_SPECIES.tsv",
                         "COHORT_MERGE_MANIFEST.JSON")
  for (name in reserved_variants) {
    cache <- file.path(root, name)
    writeLines("{}", cache)
    expect_error(
      merge_wf16s_cohort(c(a, b), file.path(root, paste0("rejected_", name)),
                         taxonomy_cache = cache),
      "E_COHORT_CACHE_INVALID"
    )
  }
  cache <- file.path(root, "alternate_taxonomy_cache.json")
  writeLines('{"fixture":true}', cache)
  output <- file.path(root, "accepted")
  expect_identical(merge_wf16s_cohort(c(a, b), output, taxonomy_cache = cache)$status, "merged")
  expect_identical(compute_file_hash(file.path(output, basename(cache))), compute_file_hash(cache))
  staged_config <- yaml::read_yaml(file.path(output, "config.yml"))
  expect_identical(staged_config$taxonomy$cache, basename(cache))
})

test_that("cohort overwrite detects destination changes at every publication boundary", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  scenarios <- c("before_prior_rename", "after_prior_rename", "before_stage_rename",
                 "after_stage_rename", "before_backup_cleanup")
  for (phase_to_change in scenarios) {
    for (mutation in c("add", "remove", "modify")) {
      root <- tempfile(paste0("cohort_boundary_", phase_to_change, "_", mutation, "_"))
      dir.create(root)
      a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
      b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
      output <- file.path(root, "merged")
      merge_wf16s_cohort(c(a, b), output)
      changed_path <- NULL
      inject <- function(phase, destination, stage, backup) {
        if (!identical(phase, phase_to_change)) return(invisible(NULL))
        target <- if (phase == "after_prior_rename" || phase == "before_backup_cleanup") backup else destination
        if (phase == "before_stage_rename" && !dir.exists(target)) dir.create(target)
        if (identical(mutation, "add") || phase == "before_stage_rename") {
          changed_path <<- file.path(target, "late-user-note.txt")
          writeLines(paste("preserve", phase), changed_path)
          if (identical(mutation, "remove")) unlink(changed_path)
          if (identical(mutation, "modify")) writeLines("modified late note", changed_path)
        } else {
          changed_path <<- file.path(target, "abundance_table_species.tsv")
          if (identical(mutation, "remove")) {
            unlink(changed_path)
          } else {
            writeLines("modified abundance artifact", changed_path)
          }
        }
      }
      expect_error(
        merge_wf16s_cohort(c(a, b), output, overwrite = TRUE, publication_hook = inject),
        "E_COHORT_DESTINATION_CHANGED|E_COHORT_RECOVERY_REQUIRED",
        info = paste(phase_to_change, mutation)
      )
      preserved_path <- if (phase_to_change == "after_prior_rename") {
        file.path(output, basename(changed_path))
      } else changed_path
      expect_identical(file.exists(preserved_path), mutation != "remove",
                       info = paste(phase_to_change, mutation))
      if (identical(mutation, "modify")) {
        expected_content <- if (phase_to_change == "before_stage_rename") {
          "modified late note"
        } else "modified abundance artifact"
        expect_identical(readLines(preserved_path, warn = FALSE), expected_content,
                         info = paste(phase_to_change, mutation))
      }
      if (phase_to_change %in% c("before_prior_rename", "after_prior_rename")) {
        expect_true(dir.exists(output), info = paste(phase_to_change, mutation))
      } else if (phase_to_change %in% c("before_stage_rename", "after_stage_rename",
                                        "before_backup_cleanup")) {
        sibling_artifacts <- list.files(root, pattern = "^\\.wf16s-cohort-backup-",
                                        full.names = TRUE, all.files = TRUE)
        expect_true(length(sibling_artifacts) >= 1L, info = paste(phase_to_change, mutation))
      }
    }
  }
})

test_that("cohort rollback retains backup when a concurrent destination blocks restoration", {
  root <- tempfile("cohort_restore_conflict_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  output <- file.path(root, "merged")
  merge_wf16s_cohort(c(a, b), output)
  backup_path <- NULL
  occupant <- file.path(output, "concurrent-destination.txt")
  inject <- function(phase, destination, stage, backup) {
    if (!identical(phase, "after_prior_rename")) return(invisible(NULL))
    backup_path <<- backup
    writeLines("late addition", file.path(backup, "late-user-note.txt"))
    dir.create(destination)
    writeLines("concurrent occupant", occupant)
  }

  expect_error(
    merge_wf16s_cohort(c(a, b), output, overwrite = TRUE, publication_hook = inject),
    "E_COHORT_RECOVERY_REQUIRED"
  )
  expect_true(file.exists(occupant))
  expect_true(file.exists(file.path(backup_path, "abundance_table_species.tsv")))
  expect_true(file.exists(file.path(backup_path, "late-user-note.txt")))
  expect_true(file.exists(cohort_transaction_paths(output)$journal))
  expect_error(recover_cohort_publication(output), "E_COHORT_RECOVERY_REQUIRED")
  expect_true(file.exists(occupant))
  expect_true(file.exists(file.path(backup_path, "late-user-note.txt")))
})

test_that("a destination lock held by another process prevents cohort publication", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  root <- tempfile("cohort_cross_process_lock_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  output <- file.path(root, "merged")
  held <- acquire_output_lock(output, timeout_ms = 1000)
  on.exit(release_output_lock(held), add = TRUE)
  child <- run_cohort_merge_child(repo_root, c(a, b), output, lock_timeout_ms = 50L)
  expect_gt(child$status, 0L)
  expect_match(paste(child$stderr, child$stdout), "E_OUTPUT_BUSY")
  expect_false(dir.exists(output))
  expect_false(file.exists(cohort_transaction_paths(output)$journal))
})

test_that("cohort publication recovers durable states after child-process termination", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  root <- tempfile("cohort_crash_recovery_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  output <- file.path(root, "merged")
  merge_wf16s_cohort(c(a, b), output)

  before_prior <- run_cohort_merge_child(repo_root, c(a, b), output, "before_prior_rename")
  expect_equal(before_prior$status, 86L)
  expect_true(file.exists(cohort_transaction_paths(output)$journal))
  expect_true(cohort_path_matches_snapshot(output,
    validate_cohort_bundle(output)$snapshot, validate_bundle = TRUE))
  expect_identical(merge_wf16s_cohort(c(a, b), output, overwrite = TRUE)$status, "merged")
  expect_false(file.exists(cohort_transaction_paths(output)$journal))

  after_prior_before_record <- run_cohort_merge_child(
    repo_root, c(a, b), output, "after_prior_rename_before_journal"
  )
  expect_equal(after_prior_before_record$status, 86L)
  prior_gap_journal <- jsonlite::fromJSON(
    cohort_transaction_paths(output)$journal, simplifyVector = FALSE
  )
  expect_identical(prior_gap_journal$phase, "ready_to_publish")
  expect_true(cohort_path_matches_snapshot(prior_gap_journal$backup,
    prior_gap_journal$prior_snapshot, validate_bundle = TRUE))
  expect_identical(merge_wf16s_cohort(c(a, b), output, overwrite = TRUE)$status, "merged")
  expect_false(file.exists(cohort_transaction_paths(output)$journal))

  before_stage <- run_cohort_merge_child(repo_root, c(a, b), output, "before_stage_rename")
  expect_equal(before_stage$status, 86L)
  expect_false(dir.exists(output))
  expect_true(file.exists(cohort_transaction_paths(output)$journal))
  expect_identical(merge_wf16s_cohort(c(a, b), output, overwrite = TRUE)$status, "merged")
  expect_false(file.exists(cohort_transaction_paths(output)$journal))

  after_prior <- run_cohort_merge_child(repo_root, c(a, b), output, "after_prior_rename")
  expect_equal(after_prior$status, 86L)
  expect_false(dir.exists(output))
  expect_true(file.exists(cohort_transaction_paths(output)$journal))
  expect_identical(merge_wf16s_cohort(c(a, b), output, overwrite = TRUE)$status, "merged")
  expect_true(isTRUE(validate_cohort_bundle(output)$manifest$compatibility$status == "compatible"))
  expect_false(file.exists(cohort_transaction_paths(output)$journal))

  after_publish_before_record <- run_cohort_merge_child(
    repo_root, c(a, b), output, "after_stage_rename_before_journal"
  )
  expect_equal(after_publish_before_record$status, 86L)
  published_gap_journal <- jsonlite::fromJSON(
    cohort_transaction_paths(output)$journal, simplifyVector = FALSE
  )
  expect_identical(published_gap_journal$phase, "prior_moved")
  expect_true(cohort_path_matches_snapshot(output,
    published_gap_journal$expected_new_snapshot, validate_bundle = TRUE))
  expect_true(cohort_path_matches_snapshot(published_gap_journal$backup,
    published_gap_journal$prior_snapshot, validate_bundle = TRUE))
  expect_identical(merge_wf16s_cohort(c(a, b), output, overwrite = TRUE)$status, "merged")
  expect_false(file.exists(cohort_transaction_paths(output)$journal))

  after_publish <- run_cohort_merge_child(repo_root, c(a, b), output, "after_stage_rename")
  expect_equal(after_publish$status, 86L)
  expect_true(dir.exists(output))
  expect_true(file.exists(cohort_transaction_paths(output)$journal))
  expect_identical(merge_wf16s_cohort(c(a, b), output, overwrite = TRUE)$status, "merged")
  expect_true(isTRUE(validate_cohort_bundle(output)$manifest$compatibility$status == "compatible"))
  expect_false(file.exists(cohort_transaction_paths(output)$journal))
})

test_that("malformed cohort journals cannot redirect recovery outside transaction paths", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  root <- tempfile("cohort_bad_journal_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  output <- file.path(root, "merged")
  merge_wf16s_cohort(c(a, b), output)
  crashed <- run_cohort_merge_child(repo_root, c(a, b), output, "after_prior_rename")
  expect_equal(crashed$status, 86L)
  journal_path <- cohort_transaction_paths(output)$journal
  journal <- jsonlite::fromJSON(journal_path, simplifyVector = FALSE)
  outside <- file.path(root, "outside")
  dir.create(outside)
  sentinel <- file.path(outside, "must-survive.txt")
  writeLines("preserve", sentinel)
  journal$stage <- outside
  jsonlite::write_json(journal, journal_path, auto_unbox = TRUE, pretty = TRUE)
  expect_error(recover_cohort_publication(output), "E_COHORT_RECOVERY_REQUIRED")
  expect_true(file.exists(sentinel))
  expect_false(dir.exists(output))
})

test_that("cohort transaction paths fit conservative Windows path limits", {
  root <- tempfile("cohort_transaction_path_length_")
  dir.create(root)
  destination <- file.path(root, "merged_cohort_destination")
  paths <- cohort_transaction_paths(destination, new_cohort_transaction_id())
  artifact <- file.path(paths$stage, "abundance_table_species.tsv")
  expect_lt(nchar(artifact), 240L)
})

test_that("a fresh prepared journal recovers before mkdir and after interrupted cleanup", {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  root <- tempfile("cohort_empty_prepared_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  output <- file.path(root, "merged")
  child <- run_cohort_merge_child(repo_root, c(a, b), output, "after_prepare_journal")
  expect_equal(child$status, 86L)
  journal_path <- cohort_transaction_paths(output)$journal
  journal <- jsonlite::fromJSON(journal_path, simplifyVector = FALSE)
  expect_false(dir.exists(journal$stage))
  expect_true(recover_cohort_publication(output))
  expect_false(file.exists(journal_path))

  # Reproduce a failure between successful stage removal and journal removal.
  write_cohort_publication_journal(journal_path, journal)
  dir.create(journal$stage)
  writeLines("partial", file.path(journal$stage, "metadata.tsv"))
  env <- new.env(parent = environment(recover_cohort_publication))
  env$recover <- recover_cohort_publication
  environment(env$recover) <- env
  env$remove_cohort_journal_checked <- function(path) stop("injected journal removal interruption")
  expect_error(env$recover(output), "injected journal removal interruption")
  expect_false(dir.exists(journal$stage))
  expect_true(file.exists(journal_path))
  expect_true(recover_cohort_publication(output))
  expect_false(file.exists(journal_path))
  expect_identical(merge_wf16s_cohort(c(a, b), output)$status, "merged")
})

test_that("failed stage mkdir is retryable without manual journal removal", {
  root <- tempfile("cohort_mkdir_failure_")
  dir.create(root)
  a <- make_synthetic_sample_bundle(file.path(root, "a"), "A", c("Eukaryota;Fungi;P;C;O;F;G;TaxonA" = 10))
  b <- make_synthetic_sample_bundle(file.path(root, "b"), "B", c("Eukaryota;Fungi;P;C;O;F;G;TaxonB" = 20))
  output <- file.path(root, "merged")
  env <- new.env(parent = environment(merge_wf16s_cohort))
  env$merge <- merge_wf16s_cohort
  environment(env$merge) <- env
  env$dir.create <- function(path, ...) {
    if (startsWith(basename(path), ".wf16s-cohort-stage-")) return(FALSE)
    base::dir.create(path, ...)
  }
  expect_error(env$merge(c(a, b), output), "E_COHORT_STAGE_FAILED")
  expect_false(file.exists(cohort_transaction_paths(output)$journal))
  expect_identical(merge_wf16s_cohort(c(a, b), output)$status, "merged")
})

test_that("empty prepared recovery retains unexpected occupants and later phases", {
  root <- tempfile("cohort_prepared_occupants_")
  dir.create(root)
  output <- file.path(root, "merged")
  transaction_id <- new_cohort_transaction_id()
  paths <- cohort_transaction_paths(output, transaction_id)
  journal <- list(journal_schema_version = 1L, transaction_id = transaction_id,
    destination = canonicalize_merge_destination(output), stage = paths$stage,
    backup = paths$backup, had_prior = FALSE, prior_snapshot = NULL,
    expected_new_snapshot = NULL, allowed_stage_basenames = list("metadata.tsv"), phase = "prepared")
  for (target in c(output, paths$stage, paths$backup)) {
    for (kind in c("file", "directory")) {
      write_cohort_publication_journal(paths$journal, journal)
      if (kind == "file") writeLines("preserve", target) else {
        dir.create(target)
        writeLines("preserve", file.path(target, "user.txt"))
      }
      expect_error(recover_cohort_publication(output), "E_COHORT_RECOVERY_REQUIRED")
      expect_true(file.exists(paths$journal))
      expect_identical(readLines(if (kind == "file") target else file.path(target, "user.txt")), "preserve")
      unlink(target, recursive = TRUE)
    }
  }
  for (phase in c("ready_to_publish", "prior_moved", "stage_published")) {
    journal$phase <- phase
    write_cohort_publication_journal(paths$journal, journal)
    expect_error(recover_cohort_publication(output), "E_COHORT_RECOVERY_REQUIRED")
    expect_true(file.exists(paths$journal))
  }
})
