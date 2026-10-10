# =============================================================================
# Unit Tests: Kraken2 / Bracken Classifier Contract & Adapter
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "classifier_contract.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))

test_that("parse_kraken_report parses standard 6-column reports and computes totals", {
  lines <- c(
    " 40.00\t20\t20\tU\t0\tunclassified",
    " 60.00\t30\t0\tR\t1\troot",
    " 60.00\t30\t0\tD\t2759\tEukaryota",
    " 60.00\t30\t0\tK\t4751\tFungi",
    " 60.00\t30\t0\tP\t4890\tAscomycota",
    " 60.00\t30\t0\tC\t147537\tSaccharomycetes",
    " 60.00\t30\t0\tO\t4891\tSaccharomycetales",
    " 60.00\t30\t0\tF\t4892\tSaccharomycetaceae",
    " 60.00\t30\t0\tG\t4930\tSaccharomyces",
    " 56.00\t28\t28\tS\t4932\t  Saccharomyces cerevisiae",
    "  4.00\t2\t2\tS\t4934\t  Saccharomyces pastorianus"
  )
  path <- tempfile("kraken_rep_", fileext = ".txt")
  writeLines(lines, path)

  rep <- parse_kraken_report(path)
  expect_equal(rep$total_reads, 50)
  expect_equal(rep$classified_reads, 30)
  expect_equal(rep$unclassified_reads, 20)
  expect_equal(nrow(rep$df), 11)
  expect_true(all(c("pct", "clade_reads", "direct_reads", "rank_code", "taxid", "name") %in% names(rep$df)))
})

test_that("parse_bracken_abundance parses 7-column report and checks column format", {
  lines <- c(
    "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads",
    "Saccharomyces cerevisiae\t4932\tS\t28\t0\t28\t1.00000"
  )
  path <- tempfile("bracken_rep_", fileext = ".report")
  writeLines(lines, path)

  br <- parse_bracken_abundance(path)
  expect_equal(br$species_estimate_sum, 28)
  expect_equal(nrow(br$df), 1)
  expect_equal(br$df$taxonomy_id[1], 4932)
  expect_equal(br$df$new_est_reads[1], 28)
  expect_equal(br$df$taxonomy_lvl[1], "S")

  # Malformed header rejection
  bad_lines <- c("bad\theader", "Saccharomyces cerevisiae\t4932")
  bad_path <- tempfile("bracken_bad_", fileext = ".report")
  writeLines(bad_lines, bad_path)
  expect_error(parse_bracken_abundance(bad_path), "Bracken report.*header violation")
})

test_that("placeholder normalization standardizes unclassified names", {
  expect_equal(normalize_placeholder_rank("unclassified Boletales family"), "Boletales_Incertae_sedis")
  expect_equal(normalize_placeholder_rank("unclassified Dipodascales genus"), "Dipodascales_Incertae_sedis")
  expect_equal(normalize_placeholder_rank("Boletus edulis"), "Boletus edulis")
  expect_equal(normalize_placeholder_species("unclassified Boletales family"), "Boletales_Incertae_sedis")
})

test_that("build_taxid_lineage_map_from_assignments aggregates lineages deterministically", {
  assignments <- data.frame(
    status = c("C", "C", "C", "U"),
    taxid = c("4932", "4932", "4934", "0"),
    lineage = c(
      "Eukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Saccharomyces|Saccharomyces cerevisiae",
      "Eukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Saccharomyces|Saccharomyces cerevisiae",
      "Eukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Saccharomyces|Saccharomyces pastorianus",
      ""
    ),
    stringsAsFactors = FALSE
  )
  res <- build_taxid_lineage_map_from_assignments(assignments)
  expect_equal(length(res$taxid_to_lineage), 2)
  expect_true("4932" %in% names(res$taxid_to_lineage))
  expect_true(grepl("Saccharomyces cerevisiae", res$taxid_to_lineage[["4932"]]))
  expect_equal(length(res$conflicts), 0)
})

test_that("match_abundance_to_assignment_lineage bridges 8-rank table and 7-rank lineages", {
  ab_lineage <- "Eukaryota;Fungi;Ascomycota;Saccharomycetes;Saccharomycetales;Saccharomycetaceae;Saccharomyces;Saccharomyces cerevisiae"
  as_lineage <- "Eukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Saccharomyces|Saccharomyces cerevisiae"
  expect_true(match_abundance_to_assignment_lineage(ab_lineage, as_lineage))

  # Test placeholder rank normalization
  ab_placeholder <- "Eukaryota;Fungi;Chytridiomycota;Chytridiomycetes;Polychytriales;Polychytriales_Incertae_sedis;Polychytrium;Polychytrium aggregatum"
  as_placeholder <- "Eukaryota|Chytridiomycota|Chytridiomycetes|Polychytriales|unclassified Polychytriales family|Polychytrium|Polychytrium aggregatum"
  expect_true(match_abundance_to_assignment_lineage(ab_placeholder, as_placeholder))

  # Non-matching lineages
  as_mismatch <- "Eukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Candida|Candida albicans"
  expect_false(match_abundance_to_assignment_lineage(ab_lineage, as_mismatch))
})

test_that("compute_kraken_sample_accounting validates N == C + U and D == C - E", {
  kraken_rep <- list(
    raw_total = 50,
    raw_classified = 30,
    raw_unclassified = 20,
    direct_counts = c("0" = 20, "1" = 2, "4932" = 28),
    report_df = data.frame()
  )
  bracken_rep <- list(
    species_estimate_sum = 28,
    rows = data.frame()
  )
  assignment_data <- data.frame(
    status = c(rep("C", 30), rep("U", 20)),
    taxid = c(rep("4932", 28), rep("1", 2), rep("0", 20)),
    stringsAsFactors = FALSE
  )

  acc <- compute_kraken_sample_accounting(
    kraken_report = kraken_rep,
    bracken_report = bracken_rep,
    abundance_unclassified = 20,
    abundance_species_sum = 28,
    assignment_data = assignment_data,
    sample_id = "sample1"
  )

  expect_equal(acc$N, 50)
  expect_equal(acc$C, 30)
  expect_equal(acc$U, 20)
  expect_equal(acc$E, 28)
  expect_equal(acc$T, 48)
  expect_equal(acc$D, 2)
  expect_equal(acc$N, acc$C + acc$U)
  expect_equal(acc$D, acc$C - acc$E)
  expect_true(acc$assignments_checked)
  expect_true(acc$reconciliation_pass)
})

test_that("scan_assignments_file validates 6-field kraken2 assignments format", {
  good_lines <- c(
    "U\tread1\t0\t500\t0:50\t",
    "C\tread2\t4932\t600\t4932:60\tEukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Saccharomyces|Saccharomyces cerevisiae"
  )
  good_path <- tempfile("k2_assign_good_", fileext = ".tsv")
  writeLines(good_lines, good_path)

  scanned <- scan_assignments_file(good_path, classifier = "kraken2")
  expect_equal(scanned, 2)

  reads <- read_assignments_file(good_path, "sample1", expected_unclassified = 1, classifier = "kraken2")
  expect_equal(nrow(reads), 2)
  expect_equal(sum(reads$status == "C"), 1)
  expect_equal(sum(reads$status == "U"), 1)

  # Check invalid hit token
  bad_hit_lines <- c(
    "C\tread2\t4932\t600\tinvalid_hit\tEukaryota|Ascomycota|Saccharomycetes|Saccharomycetales|Saccharomycetaceae|Saccharomyces|Saccharomyces cerevisiae"
  )
  bad_hit_path <- tempfile("k2_bad_hit_", fileext = ".tsv")
  writeLines(bad_hit_lines, bad_hit_path)
  expect_error(scan_assignments_file(bad_hit_path, classifier = "kraken2"), "malformed hit list")

  # Check invalid status U with non-zero taxid
  bad_u_lines <- c("U\tread1\t4932\t500\t4932:50\t")
  bad_u_path <- tempfile("k2_bad_u_", fileext = ".tsv")
  writeLines(bad_u_lines, bad_u_path)
  expect_error(scan_assignments_file(bad_u_path, classifier = "kraken2"), "status U with positive TaxID")
})

test_that("compare_upstream_params checks cohort consistency across kraken parameters", {
  rep_params <- list(
    classifier = "kraken2",
    database_set = "ncbi_16s_18s_28s_ITS",
    taxonomic_rank = "S",
    min_len = 250,
    max_len = 1000,
    min_read_qual = 15,
    abundance_threshold = 1,
    kraken2_confidence = 0,
    bracken_threshold = 10,
    bracken_length = NULL
  )
  mismatched <- rep_params
  mismatched$bracken_threshold <- 5

  expect_silent(compare_upstream_params(rep_params, rep_params, "sample1"))
  expect_error(compare_upstream_params(rep_params, mismatched, "sample2"), "contract mismatch.*bracken_threshold")
})

test_that("build_context loads synthetic kraken2 fixture cleanly", {
  fixture_dir <- file.path("..", "..", "tests", "fixtures", "synthetic_kraken2")
  cfg_path <- file.path(fixture_dir, "config.yml")
  cfg <- load_config(cfg_path)
  expect_equal(cfg$project_name, "synthetic_kraken2_test")
  expect_true(!is.null(cfg$input$kraken2_reports))
  expect_true(!is.null(cfg$input$bracken_reports))

  ctx <- build_context(cfg)
  expect_equal(ctx$upstream_contract$classifier, "kraken2")
  expect_equal(ctx$upstream_contract$abundance_source, "bracken_species_estimates")
  expect_false(ctx$upstream_contract$identity_coverage_active)
  expect_true("sample1" %in% names(ctx$sample_accounting))
  acc <- ctx$sample_accounting[["sample1"]]
  expect_equal(acc$N, 50)
  expect_equal(acc$C, 30)
  expect_equal(acc$U, 20)
  expect_equal(acc$E, 28)
  expect_equal(acc$D, 2)
})
