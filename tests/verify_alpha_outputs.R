# Independent alpha artifact checks shared by the release gate and regression tests.
verify_alpha_outputs <- function(root, manifest, resolved_config, release_samples, classified_counts) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  sha256_file <- function(path) digest::digest(file = path, algo = "sha256")
  assert_unique_key <- function(data, columns, label) {
    if (anyDuplicated(data[columns])) stop(sprintf("Duplicate key in %s.", label))
  }
  q_token <- function(q) {
    value <- sprintf("%.17g", as.numeric(q))
    value <- sub("[.]0+([eE]|$)", "\\1", value)
    value <- sub("([.][0-9]*[1-9])0+([eE]|$)", "\\1\\2", value)
    value <- gsub("[.]", "p", value)
    value <- gsub("[eE][+]", "e", value)
    value <- gsub("[eE]", "e", value)
    gsub("-", "m", value)
  }
  q_label <- function(q) sprintf("%.17g", as.numeric(q))
  hill_orders <- as.numeric(unlist(resolved_config$alpha$hill_orders))
  renyi_orders <- as.numeric(unlist(resolved_config$alpha$renyi_orders))
  tree_path <- resolved_config$input$phylogenetic_tree %||% NULL
  tip_map_path <- resolved_config$input$phylogenetic_tip_map %||% NULL
  phylo_enabled <- !is.null(tree_path)
  expected_keys <- c(
    "richness", "shannon", "ens", "simpson", "invsimpson", "fisher_alpha",
    "chao1", "ace", "pielou", "berger_parker", "heip", "evar", "mcintosh",
    paste0("renyi_q", vapply(renyi_orders, q_token, character(1))),
    paste0("hill_q", vapply(hill_orders, q_token, character(1))),
    if (phylo_enabled) c("faith_pd", "psr", "pse") else character(0)
  )
  expected_labels <- c(
    "Richness (S)", "Shannon (H)", "ENS (e^H)", "Simpson (D)",
    "Inv. Simpson", "Fisher's alpha", "Chao1", "ACE",
    "Pielou's evenness (J)", "Berger-Parker index", "Heip's evenness",
    "Smith-Wilson Evar", "McIntosh diversity index",
    sprintf("R\u00e9nyi entropy (q=%s)", vapply(renyi_orders, q_label, character(1))),
    sprintf("Hill number (q=%s)", vapply(hill_orders, q_label, character(1))),
    if (phylo_enabled) c("Faith PD", "PSR", "PSE") else character(0)
  )
  names(expected_labels) <- expected_keys

  phylo_contract <- NULL
  if (phylo_enabled) {
    stopifnot(file.exists(tree_path), !dir.exists(tree_path),
              !is.null(manifest$inputs$phylogenetic_tree),
              identical(sha256_file(tree_path), manifest$inputs$phylogenetic_tree$sha256))
    tree <- ape::read.tree(tree_path)
    stopifnot(inherits(tree, "phylo"), ape::is.rooted(tree),
              length(tree$tip.label) > 0L, !anyNA(tree$tip.label),
              all(nzchar(tree$tip.label)), !anyDuplicated(tree$tip.label),
              !is.null(tree$edge.length), length(tree$edge.length) == nrow(tree$edge),
              all(is.finite(tree$edge.length)), all(tree$edge.length >= 0))
    positive_paths <- rownames(classified_counts)[rowSums(classified_counts) > 0]
    if (is.null(tip_map_path)) {
      map <- data.frame(TaxonPath = positive_paths, TipLabel = positive_paths,
                        stringsAsFactors = FALSE)
      stopifnot(is.null(manifest$inputs$phylogenetic_tip_map))
    } else {
      stopifnot(file.exists(tip_map_path), !dir.exists(tip_map_path),
                !is.null(manifest$inputs$phylogenetic_tip_map),
                identical(sha256_file(tip_map_path),
                          manifest$inputs$phylogenetic_tip_map$sha256))
      header <- strsplit(readLines(tip_map_path, n = 1L, warn = FALSE), "\t",
                         fixed = TRUE)[[1]]
      stopifnot(identical(header, c("TaxonPath", "TipLabel")))
      map <- read.delim(tip_map_path, sep = "\t", quote = "", comment.char = "",
                        check.names = FALSE, stringsAsFactors = FALSE,
                        colClasses = "character")
      stopifnot(nrow(map) > 0L, !anyNA(map), all(nzchar(map$TaxonPath)),
                all(nzchar(map$TipLabel)), !anyDuplicated(map$TaxonPath),
                !anyDuplicated(map$TipLabel),
                !any(grepl("^Unclassified(?:;|$)", map$TaxonPath, perl = TRUE)))
    }
    stopifnot(all(positive_paths %in% map$TaxonPath))
    mapping <- stats::setNames(map$TipLabel, map$TaxonPath)
    retained_tips <- unique(unname(mapping[positive_paths]))
    stopifnot(!anyNA(retained_tips), all(retained_tips %in% tree$tip.label))
    original_tip_count <- length(tree$tip.label)
    extra_tips <- setdiff(tree$tip.label, retained_tips)
    phylo_contract <- list(
      tree = tree, mapping = mapping, original_tip_count = original_tip_count,
      retained_tip_count = length(retained_tips), pruned_tip_count = length(extra_tips)
    )
  }

  independent_phylo <- function(counts) {
    counts <- counts[counts > 0]
    tips <- unname(phylo_contract$mapping[names(counts)])
    stopifnot(length(tips) == length(counts), !anyNA(tips),
              all(tips %in% phylo_contract$tree$tip.label))
    tree <- phylo_contract$tree
    active <- match(tips, tree$tip.label)
    included <- rep(FALSE, nrow(tree$edge))
    repeat {
      hit <- which(tree$edge[, 2] %in% active & !included)
      if (!length(hit)) break
      included[hit] <- TRUE
      active <- unique(c(active, tree$edge[hit, 1]))
    }
    faith <- sum(tree$edge.length[included])
    S <- length(counts)
    psr <- pse <- NA_real_
    if (S > 1L) {
      covariance <- ape::vcv.phylo(tree, corr = TRUE)[tips, tips, drop = FALSE]
      psv <- (S * sum(diag(covariance)) - sum(covariance)) / (S * (S - 1))
      psr <- psv * S
      abundance <- as.numeric(counts)
      N <- sum(abundance)
      numerator <- N * sum(diag(covariance) * abundance) -
        as.numeric(t(abundance) %*% covariance %*% abundance)
      denominator <- N^2 - N * mean(abundance)
      if (is.finite(denominator) && denominator > 0) pse <- numerator / denominator
    }
    c(faith_pd = faith, psr = psr, pse = pse)
  }

  independent_alpha <- function(counts) {
    counts <- round(counts[counts > 0])
    N <- sum(counts)
    S <- length(counts)
    if (!N) return(stats::setNames(rep(NA_real_, length(expected_keys)), expected_keys))
    p <- counts / N
    H <- -sum(p * log(p))
    estimate <- tryCatch(suppressWarnings(vegan::estimateR(counts)),
                         error = function(e) rep(NA_real_, 5L))
    fisher <- tryCatch(suppressWarnings(as.numeric(vegan::fisher.alpha(counts))),
                       error = function(e) NA_real_)
    value <- c(
      richness = S, shannon = H, ens = exp(H), simpson = 1 - sum(p^2),
      invsimpson = 1 / sum(p^2), fisher_alpha = fisher,
      chao1 = unname(estimate["S.chao1"]), ace = unname(estimate["S.ACE"]),
      pielou = if (S > 1L) H / log(S) else NA_real_,
      berger_parker = max(counts) / N,
      heip = if (S > 1L) (exp(H) - 1) / (S - 1) else NA_real_,
      evar = if (S > 1L) 1 - (2 / pi) * atan(mean((log(counts) - mean(log(counts)))^2)) else NA_real_,
      mcintosh = if (N > 1L) (N - sqrt(sum(counts^2))) / (N - sqrt(N)) else NA_real_
    )
    renyi <- function(q) {
      if (q == 1) return(H)
      epsilon <- q - 1
      log_p <- log(p)
      if (abs(epsilon) < 1e-4) {
        mean_log <- sum(p * log_p)
        centered <- log_p - mean_log
        return(-mean_log - epsilon * sum(p * centered^2) / 2 -
                 epsilon^2 * sum(p * centered^3) / 6)
      }
      max_log <- max(log_p)
      (q / (1 - q)) * max_log + log(sum(exp(q * (log_p - max_log)))) / (1 - q)
    }
    for (q in renyi_orders) value[[paste0("renyi_q", q_token(q))]] <- renyi(q)
    for (q in hill_orders) value[[paste0("hill_q", q_token(q))]] <- exp(renyi(q))
    if (phylo_enabled) value <- c(value, independent_phylo(counts))
    value[expected_keys]
  }

  alpha_table <- read.delim(file.path(root, "02_Alpha_Diversity", "alpha_diversity.tsv"),
                            check.names = FALSE, stringsAsFactors = FALSE)
  status_table <- read.delim(file.path(root, "02_Alpha_Diversity", "alpha_metric_status.tsv"),
                             check.names = FALSE, stringsAsFactors = FALSE)
  definitions <- read.delim(file.path(root, "02_Alpha_Diversity", "alpha_metric_definitions.tsv"),
                            check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(definitions$MetricKey),
            identical(definitions$MetricKey[seq_along(expected_keys)], expected_keys),
            identical(as.character(definitions$MetricLabel[seq_along(expected_keys)]),
                      unname(expected_labels)))
  assert_unique_key(status_table, c("SampleID", "MetricKey"), "full-depth alpha status")
  stopifnot(setequal(status_table$SampleID, release_samples),
            setequal(expected_keys, status_table$MetricKey),
            nrow(status_table) == length(release_samples) * length(expected_keys),
            !any(is.nan(status_table$Value)),
            !any(is.infinite(status_table$Value), na.rm = TRUE))
  full_column_labels <- c(
    richness = "Observed species richness (S)",
    chao1 = "Chao1 (estimated richness)", shannon = "Shannon (H)",
    ens = "Effective number of species (e^H)",
    simpson = "Simpson's D (1-sum p^2)", invsimpson = "Inverse Simpson",
    pielou = "Pielou's evenness (J)", fisher_alpha = "Fisher's alpha",
    berger_parker = "Berger-Parker dominance", ace = "ACE",
    heip = "Heip's evenness", evar = "Smith-Wilson Evar",
    mcintosh = "McIntosh diversity index"
  )
  required_labels <- expected_labels
  required_labels[names(full_column_labels)] <- full_column_labels
  stopifnot(all(unname(required_labels) %in% names(alpha_table)),
            setequal(alpha_table$SampleID, release_samples), !anyDuplicated(alpha_table$SampleID))
  for (sample_id in release_samples) {
    expected <- independent_alpha(classified_counts[, sample_id])
    actual_status <- status_table[status_table$SampleID == sample_id, ]
    actual <- stats::setNames(actual_status$Value, actual_status$MetricKey)
    for (key in expected_keys) {
      if (is.finite(expected[[key]])) {
        stopifnot(isTRUE(actual_status$Valid[actual_status$MetricKey == key]),
                  isTRUE(all.equal(actual[[key]], expected[[key]], tolerance = 1e-10)))
      } else {
        row <- actual_status[actual_status$MetricKey == key, ]
        stopifnot(nrow(row) == 1L, !isTRUE(row$Valid), is.na(row$Value), nzchar(row$Reason))
      }
      column <- if (key %in% names(full_column_labels)) {
        full_column_labels[[key]]
      } else {
        expected_labels[[key]]
      }
      if (column %in% names(alpha_table)) {
        observed <- alpha_table[alpha_table$SampleID == sample_id, column]
        if (is.finite(expected[[key]])) stopifnot(isTRUE(all.equal(observed, expected[[key]], tolerance = 1e-10)))
        else stopifnot(is.na(observed))
      }
    }
  }

  # Derive eligibility from the input counts, independently of producer helpers
  # and the skip ledger. A supported cell is not enough: totals must fit too.
  eligible <- release_samples[vapply(release_samples, function(sample_id) {
    counts <- classified_counts[, sample_id]
    all(counts <= .Machine$integer.max) && sum(counts) <= .Machine$integer.max
  }, logical(1))]
  unsupported <- setdiff(release_samples, eligible)
  alpha_path <- function(name) file.path(root, "02_Alpha_Diversity", name)
  ledger_path <- alpha_path("rarefaction_availability.tsv")
  if (length(unsupported)) {
    ledger <- read.delim(ledger_path, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(identical(names(ledger), c("Operation", "SampleID", "Status", "ReasonCode", "Reason")),
              !anyNA(ledger), !anyDuplicated(ledger[c("Operation", "SampleID")]),
              all(ledger$Status == "Skipped"), all(nzchar(trimws(ledger$Reason))),
              all(ledger$ReasonCode == "E_RAREFACTION_NUMERIC_DOMAIN"))
    expected_skips <- expand.grid(Operation = c("analytical_curve", "resampling"),
                                  SampleID = unsupported, stringsAsFactors = FALSE)
    stopifnot(setequal(paste(ledger$Operation, ledger$SampleID, sep = "\r"),
                      paste(expected_skips$Operation, expected_skips$SampleID, sep = "\r")))
  } else stopifnot(!file.exists(ledger_path))
  curve <- read.delim(alpha_path("rarefaction_curve.tsv"), check.names = FALSE)
  stopifnot(identical(names(curve), c("depth", "mean_richness", "sd_richness", "SampleID")),
            setequal(curve$SampleID, eligible))
  has_resamples <- length(eligible) > 0L
  stopifnot(file.exists(alpha_path("02a_rarefaction_curve.png")) == has_resamples,
            file.exists(alpha_path("rarefaction_curve_skipped.tsv")) == !has_resamples)
  if (!has_resamples) {
    skipped <- read.delim(alpha_path("rarefaction_curve_skipped.tsv"))
    stopifnot(nrow(skipped) == 1L, skipped$Status == "Skipped",
              !is.na(skipped$Reason), nzchar(trimws(skipped$Reason)),
              !any(file.exists(alpha_path(c("rarefaction_resamples.tsv",
                "rarefaction_resamples_long.tsv", "alpha_plot_registry.tsv")))),
              !length(list.files(alpha_path(""), pattern = "^02[b-f]_.*[.]png$")))
  }
  plot_registry <- data.frame(FigureKey = character(), Status = character())
  if (has_resamples) {
    wide <- read.delim(alpha_path("rarefaction_resamples.tsv"), check.names = FALSE)
    stopifnot(setequal(wide$SampleID, eligible),
              !anyDuplicated(wide[c("SampleID", "iteration")]),
              setequal(wide$iteration, seq_len(resolved_config$alpha$resample_iterations)),
              nrow(wide) == length(eligible) * resolved_config$alpha$resample_iterations)
    long <- read.delim(file.path(root, "02_Alpha_Diversity", "rarefaction_resamples_long.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
    expected_long_columns <- c("SampleID", "iteration", "subsample_depth", "MetricKey",
                               "MetricLabel", "Parameter", "Value", "Valid", "Reason")
    stopifnot(identical(names(long), expected_long_columns),
              !anyDuplicated(long[c("SampleID", "iteration", "MetricKey")]),
              setequal(long$SampleID, eligible),
              identical(unique(as.integer(long$subsample_depth)),
                        as.integer(min(resolved_config$alpha$resample_depth,
                          max(1L, floor(min(colSums(classified_counts)) *
                                        resolved_config$alpha$resample_fraction_cap))))),
              setequal(unique(long$MetricKey), expected_keys),
              setequal(long$iteration, seq_len(resolved_config$alpha$resample_iterations)),
              nrow(long) == length(eligible) * resolved_config$alpha$resample_iterations *
                length(expected_keys),
              !any(is.nan(long$Value)), !any(is.infinite(long$Value), na.rm = TRUE),
              all(nzchar(long$Reason[!long$Valid])))
    stopifnot(all(as.character(long$MetricLabel) == unname(expected_labels[long$MetricKey])))

    plot_registry <- read.delim(file.path(root, "02_Alpha_Diversity", "alpha_plot_registry.tsv"),
                                check.names = FALSE, stringsAsFactors = FALSE)
    expected_02b <- c("richness", "shannon", "ens", "simpson", "invsimpson", "fisher_alpha")
    expected_02c <- c("chao1", "ace", "pielou", "berger_parker")
    expected_02d <- c("heip", "evar", "mcintosh", "renyi_q1")
    expected_02e <- c(paste0("hill_q", vapply(hill_orders, q_token, character(1))),
                      "faith_pd", "psr", "pse")
    stopifnot(identical(plot_registry$MetricKey[plot_registry$FigureKey == "02b"], expected_02b),
              identical(plot_registry$MetricKey[plot_registry$FigureKey == "02c"], expected_02c),
              identical(plot_registry$MetricKey[plot_registry$FigureKey == "02d"], expected_02d),
              identical(plot_registry$MetricKey[plot_registry$FigureKey == "02e"], expected_02e),
              all(plot_registry$WidthIn == 10), all(plot_registry$HeightIn == 7),
              all(plot_registry$DPI == 150),
              all(plot_registry$Status %in% c("Completed", "Ineligible", "Skipped")))

    png_dimensions <- function(path) {
      connection <- file(path, open = "rb")
      on.exit(close(connection))
      signature <- readBin(connection, what = "raw", n = 8L)
      stopifnot(identical(as.integer(signature),
                          c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
      chunk_length <- readBin(connection, what = "integer", n = 1L, size = 4L,
                              signed = TRUE, endian = "big")
      chunk_type <- rawToChar(readBin(connection, what = "raw", n = 4L))
      stopifnot(chunk_length == 13L, identical(chunk_type, "IHDR"))
      c(width = readBin(connection, what = "integer", n = 1L, size = 4L,
                        signed = TRUE, endian = "big"),
        height = readBin(connection, what = "integer", n = 1L, size = 4L,
                         signed = TRUE, endian = "big"))
    }
    for (figure in unique(plot_registry$FigureFile[plot_registry$Status == "Completed"])) {
      path <- file.path(root, "02_Alpha_Diversity", figure)
      stopifnot(file.exists(path), file.info(path)$size > 0,
                identical(unname(png_dimensions(path)), c(1500L, 1050L)))
    }
  }
  phylo_status <- read.delim(file.path(root, "02_Alpha_Diversity", "alpha_phylogeny_status.tsv"),
                             check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(nrow(phylo_status) == 1L)
  stopifnot(!is.null(manifest$alpha$phylogeny),
            identical(as.character(manifest$alpha$phylogeny$Status),
                      as.character(phylo_status$Status[[1]])))
  if (phylo_enabled) {
    stopifnot(identical(phylo_status$Status[[1]], "Completed"),
              identical(as.character(phylo_status$TreeSHA256[[1]]), sha256_file(tree_path)),
              identical(as.character(manifest$alpha$phylogeny$TreeSHA256),
                        sha256_file(tree_path)),
              as.integer(phylo_status$OriginalTipCount[[1]]) == phylo_contract$original_tip_count,
              as.integer(phylo_status$RetainedTipCount[[1]]) == phylo_contract$retained_tip_count,
              as.integer(phylo_status$PrunedTipCount[[1]]) == phylo_contract$pruned_tip_count,
              identical(as.character(phylo_status$PruningMethod[[1]]),
                        "extra tips excluded from metric traversal; original root retained"),
              !file.exists(file.path(root, "02_Alpha_Diversity", "alpha_phylogeny_skipped.tsv")),
              file.exists(file.path(root, "02_Alpha_Diversity", "02e_resample_boxplots_phylogenetic.png")) == has_resamples)
    if (!is.null(tip_map_path)) {
      stopifnot(identical(as.character(phylo_status$TipMapSHA256[[1]]),
                          sha256_file(tip_map_path)))
    }
  } else {
    stopifnot(identical(phylo_status$Status[[1]], "Skipped"),
              all(plot_registry$Status[plot_registry$FigureKey == "02e"] == "Skipped"))
    stopifnot(file.exists(file.path(root, "02_Alpha_Diversity", "alpha_phylogeny_skipped.tsv")) == has_resamples,
              !file.exists(file.path(root, "02_Alpha_Diversity", "02e_resample_boxplots_phylogenetic.png")))
  }
}
