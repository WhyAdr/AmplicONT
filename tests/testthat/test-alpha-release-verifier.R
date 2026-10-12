# Exercise the same independent checker invoked by verify_release_run.R.
source(file.path("..", "verify_alpha_outputs.R"))
source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "metrics.R"))
source(file.path("..", "..", "analysis", "utils", "alpha_phylogeny.R"))
source(file.path("..", "..", "analysis", "utils", "plotting.R"))
source(file.path("..", "..", "analysis", "02_alpha_diversity.R"))

for (phylo in c(FALSE, TRUE)) {
  for (domain in c("ordinary", "unsupported", "mixed")) {
    test_that(sprintf("release alpha gate validates %s counts with phylogeny %s", domain, phylo), {
      root <- tempfile("alpha_release_gate_")
      dir.create(root)
      cfg <- get_default_config()
      cfg$input$abundance_table <- create_temp_abundance(root, n_species = 5L, sample_names = c("A", "B"))
      cfg$input$params_json <- create_temp_params(root)
      cfg$input$metadata <- create_temp_metadata(root, sample_names = c("A", "B"))
      input <- read.delim(cfg$input$abundance_table, check.names = FALSE)
      classified <- !grepl("^Unclassified", input[[1]])
      # Test both a cell exceeding the integer limit and a total exceeding it
      # despite every individual count fitting in an integer.
      input$A[classified] <- if (domain == "unsupported") c(3e9, 10, 20, 30, 40) else c(100, 10, 20, 30, 40)
      input$B[classified] <- if (domain != "ordinary") rep(5e8, 5L) else c(50, 20, 30, 40, 60)
      input$total <- input$A + input$B
      write.table(input, cfg$input$abundance_table, sep = "\t", quote = FALSE, row.names = FALSE)
      if (phylo) {
        cfg$input$phylogenetic_tree <- file.path(root, "tree.nwk")
        cfg$input$phylogenetic_tip_map <- file.path(root, "map.tsv")
        writeLines("((t1:1,t2:1):1,(t3:1,(t4:1,t5:1):1):1);", cfg$input$phylogenetic_tree)
        write.table(data.frame(TaxonPath = input[[1]][classified], TipLabel = paste0("t", 1:5)),
                    cfg$input$phylogenetic_tip_map, sep = "\t", quote = FALSE, row.names = FALSE)
      }
      cfg$alpha$resample_depth <- 50L
      cfg$alpha$resample_iterations <- 2L
      cfg$alpha$rarefaction_points <- 3L
      cfg$output$base_dir <- file.path(root, "out")
      cfg$output$dirs <- list(alpha = file.path(cfg$output$base_dir, "02_Alpha_Diversity"))
      context <- build_context(cfg)
      suppressWarnings(run_alpha(context))
      alpha_path <- function(name) file.path(cfg$output$dirs$alpha, name)
      provenance <- read.delim(alpha_path("alpha_phylogeny_status.tsv"), check.names = FALSE)
      manifest <- list(alpha = list(phylogeny = as.list(provenance[1, ])), inputs = list())
      for (key in c("phylogenetic_tree", "phylogenetic_tip_map")) {
        if (!is.null(cfg$input[[key]])) manifest$inputs[[key]] <- list(
          sha256 = digest::digest(file = cfg$input[[key]], algo = "sha256"))
      }
      counts <- context$count_matrix[-context$unclass_index, , drop = FALSE]
      rownames(counts) <- context$taxonomy$TaxonPath[-context$unclass_index]
      verify <- function() verify_alpha_outputs(cfg$output$base_dir, manifest, cfg, c("A", "B"), counts)
      expect_error(verify(), NA)
      mutate_table <- function(name, transform) {
        path <- alpha_path(name)
        original <- readBin(path, "raw", n = file.info(path)$size)
        on.exit(writeBin(original, path))
        table <- read.delim(path, check.names = FALSE)
        write.table(transform(table), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
        expect_error(verify())
      }
      remove_output <- function(name) {
        path <- alpha_path(name)
        original <- readBin(path, "raw", n = file.info(path)$size)
        on.exit(writeBin(original, path))
        unlink(path)
        suppressWarnings(expect_error(verify()))
      }
      if (domain != "ordinary") {
        remove_output("rarefaction_availability.tsv")
        mutate_table("rarefaction_availability.tsv", function(x) x[-1, ])
        mutate_table("rarefaction_availability.tsv", function(x) rbind(x, x[1, ]))
        mutate_table("rarefaction_availability.tsv", function(x) { x$ReasonCode[1] <- "bogus"; x })
        mutate_table("rarefaction_availability.tsv", function(x) { x$Reason[1] <- ""; x })
        mutate_table("rarefaction_availability.tsv", function(x) { x$SampleID[1] <- "absent"; x })
      } else {
        write.table(data.frame(Operation = "resampling", SampleID = "A", Status = "Skipped",
          ReasonCode = "E_RAREFACTION_NUMERIC_DOMAIN", Reason = "false skip"),
          alpha_path("rarefaction_availability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        expect_error(verify())
        unlink(alpha_path("rarefaction_availability.tsv"))
      }
      if (domain != "unsupported") {
        remove_output("rarefaction_resamples_long.tsv")
        remove_output("alpha_plot_registry.tsv")
        mutate_table("rarefaction_resamples_long.tsv", function(x) x[-1, ])
        mutate_table("rarefaction_resamples_long.tsv", function(x) { x$SampleID[1] <- "absent"; x })
        remove_output("02b_resample_boxplots.png")
      } else {
        writeLines("stale output", alpha_path("rarefaction_resamples_long.tsv"))
        expect_error(verify())
        unlink(alpha_path("rarefaction_resamples_long.tsv"))
      }
      mutate_table("alpha_metric_status.tsv", function(x) { x$Value[1] <- -999; x })
      mutate_table("alpha_diversity.tsv", function(x) {
        x[["Observed species richness (S)"]] <- NULL
        x
      })
      expect_error(verify(), NA)
    })
  }
}
