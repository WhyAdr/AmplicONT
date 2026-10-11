# =============================================================================
# Workload preflight regression tests
# =============================================================================

source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "metrics.R"))
source(file.path("..", "..", "analysis", "utils", "alpha_phylogeny.R"))
source(file.path("..", "..", "analysis", "03_beta_diversity.R"))
source(file.path("..", "..", "analysis", "05_ordination.R"))
source(file.path("..", "..", "analysis", "utils", "preflight.R"))

test_that("workload preflight covers beta, ordination, permutation, and resampling", {
  root <- tempfile("workload_preflight_")
  dir.create(root)
  samples <- c("C1", "C2", "C3", "T1", "T2", "T3")
  cfg <- get_default_config()
  cfg$mode <- "cohort"
  cfg$input$abundance_table <- create_temp_abundance(root, 20L, samples)
  cfg$input$metadata <- create_temp_metadata(
    root, samples, c("C", "C", "C", "T", "T", "T")
  )
  cfg$input$params_json <- create_temp_params(root)
  cfg$output$base_dir <- file.path(root, "output")
  cfg$output$dirs <- list(beta = file.path(cfg$output$base_dir, "beta"),
                          ordination = file.path(cfg$output$base_dir, "ordination"))
  context <- build_context(cfg)

  beta_context <- context
  beta_context$config$resource_budgets$max_distance_cells <- 1L
  expect_error(assert_module_workload(beta_context, "beta"),
               "E_RESOURCE_BUDGET.*distance_matrix_cells")
  expect_error(run_module_preflight(beta_context, "beta"),
               "E_RESOURCE_BUDGET.*distance_matrix_cells")
  expect_error(run_beta(beta_context), "E_RESOURCE_BUDGET.*distance_matrix_cells")
  expect_false(dir.exists(cfg$output$base_dir))

  ordination_context <- context
  ordination_context$config$resource_budgets$max_distance_work <- 1L
  expect_error(assert_module_workload(ordination_context, "ordination"),
               "E_RESOURCE_BUDGET.*distance_work")
  expect_error(run_module_preflight(ordination_context, "ordination"),
               "E_RESOURCE_BUDGET.*distance_work")
  expect_error(run_ordination(ordination_context), "E_RESOURCE_BUDGET.*distance_work")
  expect_false(dir.exists(cfg$output$base_dir))

  cfg$beta$resampling$enabled <- TRUE
  cfg$beta$resampling$iterations <- 10L
  workload_20 <- estimate_module_workload(6L, 20L, "beta", cfg, permutation_rows = 719L)
  workload_40 <- estimate_module_workload(6L, 40L, "beta", cfg, permutation_rows = 719L)
  cells_20 <- workload_20$Value[workload_20$Metric == "rarefaction_cells"]
  cells_40 <- workload_40$Value[workload_40$Metric == "rarefaction_cells"]
  expect_equal(cells_20, 6 * 20 * 10)
  expect_equal(cells_40, 6 * 40 * 10)
  expect_gt(cells_40, cells_20)
  expect_equal(permutation_rows_within_limit(6L, 999L), 719)
  exhaustive <- build_group_permutation_design(
    factor(c("A", "A", "A", "B", "B", "B")), n_permutations = 999L
  )
  expect_equal(nrow(exhaustive$matrix), 719L)
})
