# =============================================================================
# Unit Tests: Cohort Modules & Synthetic Multi-Sample Gating
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "metrics.R"))
source(file.path("..", "..", "analysis", "utils", "alpha_phylogeny.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "02_alpha_diversity.R"))
source(file.path("..", "..", "analysis", "03_beta_diversity.R"))
source(file.path("..", "..", "analysis", "04_taxa_composition.R"))
source(file.path("..", "..", "analysis", "05_ordination.R"))
source(file.path("..", "..", "analysis", "06_shared_taxa.R"))

test_that("Single-sample context skips beta, ordination, and shared_taxa gracefully", {
  tmp <- tempdir()
  ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = c("S1"))

  cfg <- get_default_config()
  cfg$input$abundance_table <- ab_file
  cfg$input$params_json <- create_temp_params(tmp)
  cfg$output$base_dir <- file.path(tmp, "out_single")
  cfg$output$dirs <- list(
    beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"),
    alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"),
    ordination = file.path(cfg$output$base_dir, "05_Ordination"),
    shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa")
  )

  context <- build_context(cfg)
  expect_equal(context$mode, "single")

  # Run beta
  res_beta <- run_beta(context)
  expect_equal(res_beta$status, "skipped")
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "beta_diversity_skipped.tsv")))

  # Run ordination
  res_ord <- run_ordination(context)
  expect_equal(res_ord$status, "skipped")
  expect_true(file.exists(file.path(cfg$output$dirs$ordination, "ordination_skipped.tsv")))

  # Run shared taxa
  res_shared <- run_shared_taxa(context)
  expect_equal(res_shared$status, "skipped")
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "shared_taxa_skipped.tsv")))
})

test_that("Synthetic cohort (2 groups x 3 replicates) passes cohort gates and pairs PERMANOVA with betadisper", {
  tmp <- tempdir()
  sample_names <- c("Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3")
  ab_file <- create_temp_abundance(tmp, n_species = 20, sample_names = sample_names)
  meta_file <- create_temp_metadata(tmp, sample_names = sample_names,
                                    groups = c("Control", "Control", "Control", "Treated", "Treated", "Treated"))

  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- ab_file
  cfg$input$metadata <- meta_file
  cfg$input$params_json <- create_temp_params(tmp)
  cfg$output$base_dir <- file.path(tmp, "out_cohort")
  cfg$output$dirs <- list(
    alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"),
    beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"),
    ordination = file.path(cfg$output$base_dir, "05_Ordination"),
    shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa")
  )

  context <- build_context(cfg)
  expect_equal(context$mode, "cohort")
  expect_equal(length(context$samples), 6)

  # Run cohort alpha, including the metadata-aware group test path.
  res_alpha <- run_alpha(context)
  expect_equal(res_alpha$status, "completed")
  alpha_out <- read.delim(file.path(cfg$output$dirs$alpha, "alpha_diversity.tsv"), check.names = FALSE)
  expect_equal(alpha_out$SampleID, sample_names)
  expect_true("Group" %in% names(alpha_out))
  expect_false(any(c("Group.x", "Group.y") %in% names(alpha_out)))

  # Run beta
  res_beta <- run_beta(context)
  expect_equal(res_beta$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "permanova.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "betadisper.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "pcoa_scores_bray.tsv")))

  # Run ordination
  res_ord <- run_ordination(context)
  expect_equal(res_ord$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$ordination, "pca_scores.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$ordination, "pca_variance.tsv")))
  nmds_diag <- file.path(cfg$output$dirs$ordination, "nmds_diagnostics.tsv")
  expect_true(file.exists(nmds_diag))
  nmds <- read.delim(nmds_diag, check.names = FALSE)
  expect_equal(nmds$Converged[1], nmds$BestSolutionRepetitions[1] > 0)

  # Run shared taxa
  res_shared <- run_shared_taxa(context)
  expect_equal(res_shared$status, "completed")
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "core_taxa.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "unique_taxa.tsv")))
  expect_true(file.exists(file.path(cfg$output$dirs$shared_taxa, "group_prevalence.tsv")))
})

test_that("Forced cohort mode rejects a one-sample table", {
  tmp <- tempfile("forced_cohort_")
  dir.create(tmp)
  ab_file <- create_temp_abundance(tmp, n_species = 5, sample_names = "Only1")
  meta_file <- create_temp_metadata(tmp, sample_names = "Only1", groups = "Control")
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- ab_file
  cfg$input$metadata <- meta_file
  cfg$input$params_json <- create_temp_params(tmp)
  expect_error(build_context(cfg), "at least 2")
})

test_that("Under-replicated cohort skips PERMANOVA with explicit reason", {
  tmp <- tempdir()
  # 2 groups with only 1 sample each
  sample_names <- c("Ctrl1", "Trt1")
  ab_file <- create_temp_abundance(tmp, n_species = 10, sample_names = sample_names)
  meta_file <- create_temp_metadata(tmp, sample_names = sample_names, groups = c("Control", "Treated"))

  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- ab_file
  cfg$input$metadata <- meta_file
  cfg$input$params_json <- create_temp_params(tmp)
  cfg$output$base_dir <- file.path(tmp, "out_underrep")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))

  context <- build_context(cfg)
  res_beta <- run_beta(context)

  perm_tsv <- file.path(cfg$output$dirs$beta, "permanova.tsv")
  expect_true(file.exists(perm_tsv))
  perm_df <- read.delim(perm_tsv)
  expect_equal(perm_df$Status[1], "Skipped")
  expect_match(perm_df$Reason[1], "at least 2 groups with >= 2 samples")
})

test_that("Cohort distance results are invariant to sample and metadata order", {
  root <- tempfile("cohort_order_")
  dir.create(root)
  samples <- c("Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3")
  abundance_a <- create_temp_abundance(root, n_species = 12, sample_names = samples)
  abundance_df <- read.delim(abundance_a, check.names = FALSE)
  permuted <- rev(samples)
  abundance_b <- file.path(root, "synthetic_abundance_permuted.tsv")
  write.table(
    abundance_df[, c("tax", permuted, "total")], abundance_b,
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  groups <- c(Ctrl1 = "Control", Ctrl2 = "Control", Ctrl3 = "Control",
              Trt1 = "Treated", Trt2 = "Treated", Trt3 = "Treated")
  metadata_a <- create_temp_metadata(root, samples, unname(groups[samples]))
  metadata_b <- file.path(root, "metadata_permuted.tsv")
  write.table(
    data.frame(SampleID = permuted, Group = unname(groups[permuted])), metadata_b,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  run_beta_fixture <- function(abundance, metadata, output_name) {
    cfg <- get_default_config()
    cfg$mode <- "cohort"
    cfg$input$abundance_table <- abundance
    cfg$input$metadata <- metadata
    cfg$input$params_json <- create_temp_params(root)
    cfg$beta$permutations <- 19L
    cfg$output$base_dir <- file.path(root, output_name)
    cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
    context <- build_context(cfg)
    expect_equal(run_beta(context)$status, "completed")
    result <- read.delim(
      file.path(cfg$output$dirs$beta, "distance_bray.tsv"),
      check.names = FALSE
    )
    result <- result[order(result$SampleID), c("SampleID", sort(samples))]
    rownames(result) <- NULL
    result
  }

  expect_equal(
    run_beta_fixture(abundance_a, metadata_a, "first"),
    run_beta_fixture(abundance_b, metadata_b, "second"),
    tolerance = 1e-12
  )
})

test_that("PERMANOVA retains the complete exchangeability design and matches vegan oracles", {
  groups <- factor(c("A", "A", "B", "B"))
  design <- build_group_permutation_design(groups, n_permutations = 23L, seed = 42L)
  expect_true(design$testable)
  expect_equal(nrow(design$matrix), 23L)
  label_changing <- apply(design$matrix, 1L, function(index) {
    !identical(unname(as.character(groups[index])), unname(as.character(groups)))
  })
  expect_true(any(!label_changing))

  oracle_distance <- stats::dist(matrix(c(0, 1, 10, 11), ncol = 1L))
  oracle <- vegan::adonis2(oracle_distance ~ groups, permutations = design$matrix)
  expect_equal(oracle$`Pr(>F)`[1], 1 / 3, tolerance = 1e-12)
  dispersion <- vegan::permutest(
    vegan::betadisper(oracle_distance, groups), permutations = design$matrix
  )
  expect_true(is.finite(dispersion$tab$`Pr(>F)`[1]))

  blocked_groups <- factor(c("A", "A", "A", "B", "B", "B"))
  blocked_strata <- factor(c("A-only", "A-only", "mixed", "mixed", "B-only", "B-only"))
  blocked <- build_group_permutation_design(
    blocked_groups, n_permutations = 100L, strata = blocked_strata, seed = 7L
  )
  expect_true(blocked$testable)
  expect_equal(nrow(blocked$matrix), 7L)
  blocked_changes <- apply(blocked$matrix, 1L, function(index) {
    !identical(unname(as.character(blocked_groups[index])), unname(as.character(blocked_groups)))
  })
  expect_true(any(!blocked_changes))
  confounded <- assess_group_exchangeability(
    factor(c("A", "A", "B", "B")), factor(c("a", "a", "b", "b"))
  )
  expect_false(confounded$testable)
})

test_that("Pipeline PERMANOVA and betadisper use every row of the reported design", {
  root <- tempfile("beta_permutation_oracle_")
  dir.create(root)
  samples <- c("A1", "A2", "B1", "B2")
  abundance <- data.frame(
    tax = c("Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown",
            "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;A",
            "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;B"),
    A1 = c(1, 10, 0), A2 = c(1, 10, 0),
    B1 = c(1, 0, 10), B2 = c(1, 0, 10), check.names = FALSE
  )
  abundance$total <- rowSums(abundance[, samples])
  abundance_path <- file.path(root, "abundance.tsv")
  write.table(abundance, abundance_path, sep = "\t", row.names = FALSE, quote = FALSE)
  metadata_path <- create_temp_metadata(root, samples, c("A", "A", "B", "B"))
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- abundance_path
  cfg$input$metadata <- metadata_path
  cfg$input$params_json <- create_temp_params(root)
  cfg$beta$distances <- "bray"
  cfg$beta$primary_distance <- "bray"
  cfg$beta$permutations <- 23L
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
  context <- build_context(cfg)
  expect_equal(run_beta(context)$status, "completed")

  beta_dir <- cfg$output$dirs$beta
  distance_frame <- read.delim(file.path(beta_dir, "distance_bray.tsv"), check.names = FALSE)
  distance <- stats::as.dist(as.matrix(distance_frame[, samples]))
  design <- build_group_permutation_design(context$metadata$Group, 23L, seed = cfg$seed)
  direct_permanova <- vegan::adonis2(distance ~ Group, data = context$metadata,
                                    permutations = design$matrix)
  reported_permanova <- read.delim(file.path(beta_dir, "permanova.tsv"), check.names = FALSE)
  group_row <- reported_permanova[reported_permanova$Term == "Model", , drop = FALSE]
  expect_equal(nrow(group_row), 1L)
  expect_equal(group_row$EffectivePermutations, nrow(design$matrix))
  expect_equal(group_row[["Pr(>F)"]], direct_permanova$`Pr(>F)`[1], tolerance = 1e-12)

  direct_dispersion <- vegan::permutest(
    vegan::betadisper(distance, context$metadata$Group), permutations = design$matrix
  )
  # The zero-dispersion fixture can produce an undefined p-value on Linux.
  # Preserve the numeric column type even when its only value is missing.
  reported_dispersion <- read.delim(file.path(beta_dir, "betadisper.tsv"),
                                    check.names = FALSE, colClasses = c(P_Value = "numeric"))
  expect_equal(reported_dispersion$EffectivePermutations, nrow(design$matrix))
  expect_equal(reported_dispersion$P_Value, direct_dispersion$tab$`Pr(>F)`[1], tolerance = 1e-12)
})

test_that("Threshold-empty Jaccard profiles are skipped without invalidating Bray distances", {
  root <- tempfile("jaccard_empty_profiles_")
  dir.create(root)
  samples <- c("C1", "C2", "T1", "T2")
  lineages <- c(
    "Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;A",
    "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;B"
  )
  abundance <- data.frame(tax = lineages, check.names = FALSE)
  abundance$C1 <- c(1, 101, 0)
  abundance$C2 <- c(1, 101, 0)
  abundance$T1 <- c(1, 50, 2)
  abundance$T2 <- c(1, 50, 2)
  abundance$total <- rowSums(abundance[, samples])
  abundance_path <- file.path(root, "abundance.tsv")
  write.table(abundance, abundance_path, sep = "\t", row.names = FALSE, quote = FALSE)
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- abundance_path
  cfg$input$metadata <- create_temp_metadata(root, samples, c("C", "C", "T", "T"))
  cfg$input$params_json <- create_temp_params(root)
  cfg$beta$minimum_count <- 100L
  cfg$beta$primary_distance <- "jaccard"
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
  result <- run_beta(build_context(cfg))
  expect_true(any(grepl("E_EMPTY_THRESHOLD_PROFILE", result$warnings, fixed = TRUE)))
  expect_true(file.exists(file.path(cfg$output$dirs$beta, "distance_bray.tsv")))
  expect_false(file.exists(file.path(cfg$output$dirs$beta, "distance_jaccard.tsv")))
  skipped <- read.delim(file.path(cfg$output$dirs$beta, "distance_jaccard_skipped.tsv"),
                        check.names = FALSE)
  expect_identical(skipped$ReasonCode, "E_EMPTY_THRESHOLD_PROFILE")
  expect_match(skipped$AffectedSamples, "T1")
  permanova <- read.delim(file.path(cfg$output$dirs$beta, "permanova.tsv"), check.names = FALSE)
  expect_identical(permanova$DesignStatus, "primary_distance_undefined")
})

test_that("Jaccard handles all-empty, single-empty, and finite-identical profiles", {
  run_profile_case <- function(label, classified, minimum_count) {
    root <- tempfile(paste0("jaccard_profile_", label, "_"))
    dir.create(root)
    samples <- colnames(classified)
    abundance <- data.frame(
      tax = c("Unclassified;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown",
              "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;A",
              "Bacteria;Bacillati;Bacillota;Bacilli;Bacillales;Bacillaceae;Bacillus;B"),
      check.names = FALSE
    )
    for (index in seq_along(samples)) {
      abundance[[samples[[index]]]] <- c(1, classified[1, index], classified[2, index])
    }
    abundance$total <- rowSums(abundance[, samples, drop = FALSE])
    abundance_path <- file.path(root, "abundance.tsv")
    write.table(abundance, abundance_path, sep = "\t", row.names = FALSE, quote = FALSE)
    cfg <- get_default_config()
    cfg$mode <- "cohort"
    cfg$input$abundance_table <- abundance_path
    cfg$input$metadata <- create_temp_metadata(root, samples, c("A", "A", "B", "B"))
    cfg$input$params_json <- create_temp_params(root)
    cfg$beta$minimum_count <- minimum_count
    cfg$beta$primary_distance <- "jaccard"
    cfg$output$base_dir <- file.path(root, "output")
    cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
    result <- run_beta(build_context(cfg))
    list(config = cfg, result = result)
  }

  samples <- c("A1", "A2", "B1", "B2")
  all_empty <- run_profile_case(
    "all_empty", matrix(c(10, 20, 20, 10, 20, 10, 10, 20), nrow = 2L,
                         dimnames = list(NULL, samples)), 100L
  )
  expect_true(any(grepl("E_EMPTY_THRESHOLD_PROFILE", all_empty$result$warnings, fixed = TRUE)))
  all_empty_diag <- read.delim(file.path(all_empty$config$output$dirs$beta,
                                         "distance_jaccard_skipped.tsv"), check.names = FALSE)
  expect_match(all_empty_diag$AffectedSamples, "A1.*A2.*B1.*B2")

  one_empty <- run_profile_case(
    "one_empty", matrix(c(101, 0, 101, 0, 50, 2, 101, 0), nrow = 2L,
                         dimnames = list(NULL, samples)), 100L
  )
  one_empty_diag <- read.delim(file.path(one_empty$config$output$dirs$beta,
                                         "distance_jaccard_skipped.tsv"), check.names = FALSE)
  expect_identical(one_empty_diag$AffectedSamples, "B1")
  expect_false(file.exists(file.path(one_empty$config$output$dirs$beta, "distance_jaccard.tsv")))

  identical_profiles <- run_profile_case(
    "identical", matrix(c(101, 0, 101, 0, 101, 0, 101, 0), nrow = 2L,
                         dimnames = list(NULL, samples)), 100L
  )
  distance_path <- file.path(identical_profiles$config$output$dirs$beta, "distance_jaccard.tsv")
  expect_true(file.exists(distance_path))
  distance <- read.delim(distance_path, check.names = FALSE)[, samples]
  expect_true(all(is.finite(as.matrix(distance))))
  expect_equal(unname(as.matrix(distance)), matrix(0, 4L, 4L), tolerance = 0)
  expect_false(file.exists(file.path(identical_profiles$config$output$dirs$beta,
                                     "distance_jaccard_skipped.tsv")))
})

test_that("Shared taxa preserves numeric-like sample IDs and arbitrary group labels", {
  root <- tempfile("shared_identity_")
  dir.create(root)
  samples <- c("01", "02")
  abundance <- create_temp_abundance(root, n_species = 4, sample_names = samples)
  metadata <- create_temp_metadata(root, samples, c("0", "Group-one"))
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- abundance
  cfg$input$metadata <- metadata
  cfg$input$params_json <- create_temp_params(root)
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(shared_taxa = file.path(cfg$output$base_dir, "06_Shared_Taxa"))
  context <- build_context(cfg)
  expect_identical(context$metadata$SampleID, samples)
  expect_identical(context$metadata$Group, c("0", "Group-one"))
  expect_equal(run_shared_taxa(context)$status, "completed")
  prevalence <- read.delim(file.path(cfg$output$dirs$shared_taxa, "group_prevalence.tsv"),
                           check.names = FALSE)
  expect_true(all(c("0", "Group-one") %in% names(prevalence)))
})

test_that("One-row composition heatmaps do not attempt row clustering", {
  root <- tempfile("one_row_heatmap_")
  dir.create(root)
  samples <- c("S1", "S2")
  abundance <- create_temp_abundance(root, n_species = 3, sample_names = samples)
  metadata <- create_temp_metadata(root, samples, c("Control", "Treated"))
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- abundance
  cfg$input$metadata <- metadata
  cfg$input$params_json <- create_temp_params(root)
  cfg$composition$top_n_taxa <- 1L
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(composition = file.path(cfg$output$base_dir, "04_Taxa_Composition"))
  context <- build_context(cfg)
  expect_equal(run_taxa_composition(context)$status, "completed")
  heatmap <- file.path(cfg$output$dirs$composition, "04_heatmap_genus.png")
  expect_true(file.exists(heatmap))
  expect_gt(file.info(heatmap)$size, 0)
})

test_that("NMDS repetition helper follows vegan metaMDS semantics", {
  expect_false(nmds_solution_repeated(0L))
  expect_true(nmds_solution_repeated(1L))
  expect_true(nmds_solution_repeated(2L))
  expect_false(nmds_solution_repeated(NULL))
  expect_false(nmds_solution_repeated(NA))
  expect_false(nmds_solution_repeated("not-a-number"))
})

test_that("Beta rarefaction stability records isolated Procrustes failures", {
  root <- tempfile("beta_resampling_failure_")
  dir.create(root)
  samples <- c("Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3")
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- create_temp_abundance(root, n_species = 12, sample_names = samples)
  cfg$input$metadata <- create_temp_metadata(
    root, samples, c("Control", "Control", "Control", "Treated", "Treated", "Treated")
  )
  cfg$input$params_json <- create_temp_params(root)
  cfg$beta$permutations <- 19L
  cfg$beta$resampling$enabled <- TRUE
  cfg$beta$resampling$iterations <- 3L
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "03_Beta_Diversity"))
  context <- build_context(cfg)

  calls <- 0L
  fail_once <- function(X, Y) {
    calls <<- calls + 1L
    if (calls == 1L) stop("injected Procrustes failure")
    vegan::procrustes(X, Y)
  }
  expect_equal(run_beta(context, procrustes_fn = fail_once)$status, "completed")
  diagnostics <- read.delim(
    file.path(cfg$output$dirs$beta, "pcoa_rarefaction_diagnostics.tsv"), check.names = FALSE
  )
  stability <- read.delim(
    file.path(cfg$output$dirs$beta, "pcoa_rarefaction_stability.tsv"), check.names = FALSE
  )
  expect_equal(diagnostics$SuccessfulIterations, 2L)
  expect_identical(as.character(diagnostics$FailedIterations), "1")
  expect_setequal(unique(stability$Iteration), c(2L, 3L))

  always_failed <- run_beta(context, procrustes_fn = function(X, Y) stop("always fail"))
  expect_equal(always_failed$status, "completed")
  skipped_diagnostics <- read.delim(
    file.path(cfg$output$dirs$beta, "pcoa_rarefaction_diagnostics.tsv"), check.names = FALSE
  )
  expect_equal(skipped_diagnostics$Status, "Skipped")
  expect_identical(skipped_diagnostics$ReasonCode, "E_MINIMUM_SUCCESS_NOT_MET")
  expect_equal(skipped_diagnostics$MinimumSuccessfulIterations, 2L)
  expect_equal(skipped_diagnostics$SuccessfulIterations, 0L)
  expect_false(file.exists(file.path(cfg$output$dirs$beta, "pcoa_rarefaction_stability.tsv")))
})
