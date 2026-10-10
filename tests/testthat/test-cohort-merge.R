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

  # 3. assignments
  assign_dir <- file.path(out_dir, "reads_assignments")
  dir.create(assign_dir, recursive = TRUE, showWarnings = FALSE)
  assign_file <- file.path(assign_dir, sprintf("%s.assignments.tsv", sample_id))
  writeLines(c("read_id\tstatus\ttaxid", sprintf("r1\tC\t123")), assign_file)

  # 4. reports if kraken2
  if (classifier == "kraken2") {
    k_dir <- file.path(out_dir, "kraken2")
    dir.create(k_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(c(" 100.00\t10\t10\tU\t0\tunclassified"),
               file.path(k_dir, sprintf("%s.kraken2.report.txt", sample_id)))

    b_dir <- file.path(out_dir, "bracken")
    dir.create(b_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(c("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"),
               file.path(b_dir, sprintf("%s.bracken.report", sample_id)))
  }

  out_dir
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
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Fungi;TaxonA" = 10))

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle1), output_dir = file.path(root, "merged")),
    "E_COHORT_DUPLICATE_SOURCE"
  )
})

test_that("merge_wf16s_cohort rejects output directory overlapping with input", {
  root <- tempfile("cohort_test_overlap_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Fungi;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Fungi;TaxonB" = 20))

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(bundle1, "nested_out")),
    "E_COHORT_OUTPUT_OVERLAP"
  )
})

test_that("merge_wf16s_cohort rejects colliding sample IDs", {
  root <- tempfile("cohort_test_same_sample_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "SameSampleID", c("Fungi;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "SameSampleID", c("Fungi;TaxonB" = 20))

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(root, "merged")),
    "E_COHORT_DUPLICATE_SAMPLE"
  )
})

test_that("merge_wf16s_cohort rejects mixed classifiers across inputs", {
  root <- tempfile("cohort_test_mixed_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Fungi;TaxonA" = 10), classifier = "kraken2")
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Fungi;TaxonB" = 20), classifier = "minimap2")

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(root, "merged")),
    "E_COHORT_INCOMPATIBLE"
  )
})

test_that("merge_wf16s_cohort rejects threshold mismatch across inputs", {
  root <- tempfile("cohort_test_thresh_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Fungi;TaxonA" = 10), min_len = 250)
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Fungi;TaxonB" = 20), min_len = 300)

  expect_error(
    merge_wf16s_cohort(input_dirs = c(bundle1, bundle2), output_dir = file.path(root, "merged")),
    "E_COHORT_INCOMPATIBLE"
  )
})

test_that("merge_wf16s_cohort --validate-only succeeds without mutating filesystem", {
  root <- tempfile("cohort_test_val_only_")
  dir.create(root)
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Fungi;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Fungi;TaxonB" = 20))
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
                                         c("Fungi;TaxonA" = 10, "Fungi;TaxonB" = 20), unclass_count = 5)
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2",
                                         c("Fungi;TaxonB" = 30, "Fungi;TaxonC" = 40), unclass_count = 8)
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
  row_a <- ab_df[ab_df$tax == "Fungi;TaxonA", ]
  expect_equal(nrow(row_a), 1)
  expect_equal(as.numeric(row_a$Sample_1), 10.0)
  expect_equal(as.numeric(row_a$Sample_2), 0.0)
  expect_equal(as.numeric(row_a$total), 10.0)

  # TaxonB: present in both S1 (20) and S2 (30)
  row_b <- ab_df[ab_df$tax == "Fungi;TaxonB", ]
  expect_equal(nrow(row_b), 1)
  expect_equal(as.numeric(row_b$Sample_1), 20.0)
  expect_equal(as.numeric(row_b$Sample_2), 30.0)
  expect_equal(as.numeric(row_b$total), 50.0)

  # TaxonC: absent in S1 (0), present in S2 (40)
  row_c <- ab_df[ab_df$tax == "Fungi;TaxonC", ]
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
  bundle1 <- make_synthetic_sample_bundle(file.path(root, "s1"), "Sample_1", c("Fungi;TaxonA" = 10))
  bundle2 <- make_synthetic_sample_bundle(file.path(root, "s2"), "Sample_2", c("Fungi;TaxonB" = 20))
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
