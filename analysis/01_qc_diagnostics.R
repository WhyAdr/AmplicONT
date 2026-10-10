# =============================================================================
# Module 01: QC Diagnostics (Read Length, Classification Reconciliation)
# =============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

run_qc <- function(context) {
  cfg <- context$config
  qc_dir <- cfg$output$dirs$qc
  assignments_map <- context$assignments
  is_kraken <- identical(context$upstream_contract$classifier, "kraken2")
  classifier <- context$upstream_contract$classifier %||% "minimap2"
  count_basis <- context$count_basis %||% (if (is_kraken) "bracken_species_estimates" else "minimap2_alignments")
  is_16s <- identical(context$params$database_set, "ncbi_16s_18s")

  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  all_outputs <- character(0)
  reconciliation_rows <- list()
  length_summary_rows <- list()
  accounting_rows <- list()
  investigation_rows <- list()
  median_or_na <- function(x) if (length(x)) stats::median(x) else NA_real_

  # Target lengths from params.json or config fallbacks
  min_len_target <- if (!is.null(context$params$min_len)) {
    context$params$min_len
  } else {
    cfg$qc$target_min_length
  }
  max_len_target <- if (!is.null(context$params$max_len)) {
    context$params$max_len
  } else {
    cfg$qc$target_max_length
  }

  display_min <- cfg$qc$display_min_length %||% 1200L
  display_max <- cfg$qc$display_max_length %||% 1800L

  for (sample_id in context$samples) {
    asgn_path <- assignments_map[[sample_id]]
    # Expected counts from abundance context
    stat_row <- context$sample_stats[context$sample_stats$SampleID == sample_id, ]
    exp_total <- stat_row$TotalReads[1]
    exp_class <- stat_row$ClassifiedReads[1]
    exp_unclass <- stat_row$UnclassifiedReads[1]
    bamstats_path <- if (is_kraken) NULL else context$bamstats[[sample_id]]
    bamstats_available <- !is.null(bamstats_path) && !is.na(bamstats_path) && file.exists(bamstats_path)
    assignment_available <- !is.null(asgn_path) && file.exists(asgn_path)

    kraken_acct <- if (is_kraken) context$sample_accounting[[sample_id]] else NULL

    if (is_kraken) {
      accounting_rows[[sample_id]] <- data.frame(
        SampleID = sample_id,
        Classifier = "kraken2",
        CountBasis = count_basis,
        AssignmentAvailable = assignment_available,
        AbundanceTotal = exp_total,
        AbundanceClassified = exp_class,
        AbundanceUnclassified = exp_unclass,
        RawC = if (assignment_available) kraken_acct$C else NA_integer_,
        RawU = if (assignment_available) kraken_acct$U else NA_integer_,
        C_TaxID0 = if (assignment_available) 0L else NA_integer_,
        TaxID_GT0 = if (assignment_available) kraken_acct$C else NA_integer_,
        EffectiveClassifiedPct = 100 * exp_class / exp_total,
        C0ShareOfEffectiveUnclassifiedPct = 0,
        RawTotal_N = kraken_acct$N,
        RawClassified_C = kraken_acct$C,
        RawUnclassified_U = kraken_acct$U,
        BrackenEstimate_E = kraken_acct$E,
        AbundanceSum_T = kraken_acct$T,
        Difference_D = kraken_acct$D,
        ReportAvailable = TRUE,
        BrackenAvailable = TRUE,
        ReportVerified = TRUE,
        BrackenVerified = TRUE,
        AssignmentVerified = if (assignment_available) isTRUE(kraken_acct$assignments_checked) else NA,
        stringsAsFactors = FALSE
      )
    } else {
      accounting_rows[[sample_id]] <- data.frame(
        SampleID = sample_id,
        Classifier = "minimap2",
        CountBasis = count_basis,
        AssignmentAvailable = assignment_available,
        AbundanceTotal = exp_total,
        AbundanceClassified = exp_class,
        AbundanceUnclassified = exp_unclass,
        RawC = NA_integer_,
        RawU = NA_integer_,
        C_TaxID0 = NA_integer_,
        TaxID_GT0 = NA_integer_,
        EffectiveClassifiedPct = 100 * exp_class / exp_total,
        C0ShareOfEffectiveUnclassifiedPct = NA_real_,
        RawTotal_N = NA_integer_,
        RawClassified_C = NA_integer_,
        RawUnclassified_U = NA_integer_,
        BrackenEstimate_E = NA_integer_,
        AbundanceSum_T = exp_total,
        Difference_D = 0L,
        ReportAvailable = FALSE,
        BrackenAvailable = FALSE,
        ReportVerified = NA,
        BrackenVerified = NA,
        AssignmentVerified = NA,
        stringsAsFactors = FALSE
      )
    }

    investigation_rows[[sample_id]] <- data.frame(
      SampleID = sample_id,
      AssignmentAvailable = assignment_available,
      BamstatsAvailable = bamstats_available,
      MedianClassifiedLength = NA_real_,
      MedianC0Length = NA_real_,
      MedianRawULength = NA_real_,
      MinPercentIdentity = if (is_kraken) NA_real_ else context$params$min_percent_identity,
      MinRefCoverage = if (is_kraken) NA_real_ else context$params$min_ref_coverage,
      BamstatsC0Matched = NA_integer_,
      IdentityOnlyFailed = NA_integer_,
      RefCoverageOnlyFailed = NA_integer_,
      BothFailed = NA_integer_,
      stringsAsFactors = FALSE
    )

    if (is_kraken && !assignment_available) {
      reconciliation_rows[[sample_id]] <- data.frame(
        SampleID = sample_id,
        Classifier = "kraken2",
        CountBasis = count_basis,
        TotalReads = exp_total,
        AbundanceClassified = exp_class,
        AbundanceUnclassified = exp_unclass,
        AssignmentAvailable = FALSE,
        AssignmentTotal = NA_integer_,
        AssignmentRawStatusC = NA_integer_,
        AssignmentRawStatusU = NA_integer_,
        AssignmentFailedFilter = NA_integer_,
        EffectiveClassified = NA_integer_,
        EffectiveUnclassified = NA_integer_,
        RawTotal_N = kraken_acct$N,
        RawClassified_C = kraken_acct$C,
        RawUnclassified_U = kraken_acct$U,
        BrackenEstimate_E = kraken_acct$E,
        AbundanceSum_T = kraken_acct$T,
        Difference_D = kraken_acct$D,
        ReportVerified = TRUE,
        BrackenVerified = TRUE,
        AssignmentVerified = NA,
        ReconciliationPass = isTRUE(kraken_acct$reconciliation_pass),
        stringsAsFactors = FALSE
      )
      next
    }
    if (!is_kraken && !assignment_available) next

    sample_out_dir <- file.path(qc_dir, sanitize_filename(sample_id))
    dir.create(sample_out_dir, recursive = TRUE, showWarnings = FALSE)

    reads <- context$assignment_data[[sample_id]]
    if (is.null(reads)) {
      reads <- read_assignments_file(
        path = asgn_path,
        sample_id = sample_id,
        expected_total = if (is_kraken) kraken_acct$N else exp_total,
        expected_classified = if (is_kraken) kraken_acct$C else exp_class,
        expected_unclassified = if (is_kraken) kraken_acct$U else exp_unclass,
        classifier = classifier
      )
    }

    n_status_C <- sum(reads$status == "C")
    n_status_U <- sum(reads$status == "U")
    n_qc_reclass <- sum(reads$status == "C" & reads$taxid == "0")
    n_eff_class <- sum(reads$effective_classified)
    n_eff_unclass <- sum(!reads$effective_classified)
    accounting_rows[[sample_id]][c("RawC", "RawU", "C_TaxID0", "TaxID_GT0")] <-
      list(n_status_C, n_status_U, n_qc_reclass, n_eff_class)

    if (!is_kraken) {
      accounting_rows[[sample_id]]$RawTotal_N <- nrow(reads)
      accounting_rows[[sample_id]]$RawClassified_C <- n_status_C
      accounting_rows[[sample_id]]$RawUnclassified_U <- n_status_U
      accounting_rows[[sample_id]]$AssignmentVerified <- (n_eff_class == exp_class) && (nrow(reads) == exp_total)
      accounting_rows[[sample_id]]$C0ShareOfEffectiveUnclassifiedPct <- if (exp_unclass > 0) {
        100 * n_qc_reclass / exp_unclass
      } else {
        NA_real_
      }
    } else {
      accounting_rows[[sample_id]]$AssignmentVerified <- (n_status_C == kraken_acct$C) &&
        (n_status_U == kraken_acct$U) && (n_qc_reclass == 0L)
    }

    investigation_rows[[sample_id]]$MedianClassifiedLength <-
      median_or_na(reads$read_length[reads$effective_classified])
    investigation_rows[[sample_id]]$MedianC0Length <-
      median_or_na(reads$read_length[reads$status == "C" & reads$taxid == "0"])
    investigation_rows[[sample_id]]$MedianRawULength <-
      median_or_na(reads$read_length[reads$status == "U"])

    if (bamstats_available && !is_kraken) {
      partition <- partition_minimap2_failures(reads, bamstats_path, context$params, sample_id)
      investigation_rows[[sample_id]][c(
        "BamstatsC0Matched", "IdentityOnlyFailed", "RefCoverageOnlyFailed", "BothFailed"
      )] <- list(partition$matched, partition$identity_only, partition$coverage_only, partition$both)
    }

    # 1. Reconciliation Table Row
    if (is_kraken) {
      reconciliation_rows[[sample_id]] <- data.frame(
        SampleID = sample_id,
        Classifier = "kraken2",
        CountBasis = count_basis,
        TotalReads = exp_total,
        AbundanceClassified = exp_class,
        AbundanceUnclassified = exp_unclass,
        AssignmentAvailable = TRUE,
        AssignmentTotal = nrow(reads),
        AssignmentRawStatusC = n_status_C,
        AssignmentRawStatusU = n_status_U,
        AssignmentFailedFilter = n_qc_reclass,
        EffectiveClassified = n_eff_class,
        EffectiveUnclassified = n_eff_unclass,
        RawTotal_N = kraken_acct$N,
        RawClassified_C = kraken_acct$C,
        RawUnclassified_U = kraken_acct$U,
        BrackenEstimate_E = kraken_acct$E,
        AbundanceSum_T = kraken_acct$T,
        Difference_D = kraken_acct$D,
        ReportVerified = TRUE,
        BrackenVerified = TRUE,
        AssignmentVerified = (n_status_C == kraken_acct$C && n_status_U == kraken_acct$U && n_qc_reclass == 0L),
        ReconciliationPass = (n_status_C == kraken_acct$C && n_status_U == kraken_acct$U &&
                                n_qc_reclass == 0L && isTRUE(kraken_acct$reconciliation_pass)),
        stringsAsFactors = FALSE
      )
    } else {
      reconciliation_rows[[sample_id]] <- data.frame(
        SampleID = sample_id,
        Classifier = "minimap2",
        CountBasis = count_basis,
        TotalReads = exp_total,
        AbundanceClassified = exp_class,
        AbundanceUnclassified = exp_unclass,
        AssignmentAvailable = TRUE,
        AssignmentTotal = nrow(reads),
        AssignmentRawStatusC = n_status_C,
        AssignmentRawStatusU = n_status_U,
        AssignmentFailedFilter = n_qc_reclass,
        EffectiveClassified = n_eff_class,
        EffectiveUnclassified = n_eff_unclass,
        RawTotal_N = nrow(reads),
        RawClassified_C = n_status_C,
        RawUnclassified_U = n_status_U,
        BrackenEstimate_E = NA_integer_,
        AbundanceSum_T = exp_total,
        Difference_D = 0L,
        ReportVerified = NA,
        BrackenVerified = NA,
        AssignmentVerified = (n_eff_class == exp_class) && (nrow(reads) == exp_total),
        ReconciliationPass = (n_eff_class == exp_class) && (nrow(reads) == exp_total),
        stringsAsFactors = FALSE
      )
    }

    # 2. Length summary statistics and Figure 0a
    if (is_kraken) {
      reads$status_category <- ifelse(reads$status == "C", "Classified", "Unclassified")
      reads$effective_status <- reads$status_category

      accounting_plot <- data.frame(
        Category = factor(c("Kraken unclassified", "Kraken classified"),
                          levels = c("Kraken unclassified", "Kraken classified")),
        Count = c(n_status_U, n_status_C)
      ) %>%
        mutate(Fraction = Count / sum(Count),
               Label = sprintf("%s\n%s (%.1f%%)", Category, scales::comma(Count), 100 * Fraction))

      p0a <- ggplot(accounting_plot, aes(x = 3, y = Count, fill = Category)) +
        geom_col(width = 1, color = "white", linewidth = 1.2) +
        coord_polar(theta = "y") +
        xlim(c(1, 4)) +
        geom_text(data = accounting_plot[accounting_plot$Count > 0, , drop = FALSE],
                  aes(x = 3.5, label = Label),
                  position = position_stack(vjust = 0.5), size = 3.1) +
        annotate("text", x = 1, y = 0,
                 label = sprintf("%s\nreads", scales::comma(kraken_acct$N)),
                 size = 4.2, fontface = "bold") +
        scale_fill_manual(
          values = c("Kraken unclassified" = "#bdbdbd", "Kraken classified" = "#1b9e77"),
          labels = c(
            "Kraken unclassified" = "Unclassified by Kraken2 against reference database",
            "Kraken classified" = "Classified by Kraken2 (raw read assignment)"
          )
        ) +
        theme_void() +
        labs(title = sprintf("Exact Read Accounting: %s", sample_id), fill = NULL) +
        theme(legend.position = "bottom",
              plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
    } else {
      reads$status_category <- ifelse(
        reads$status == "C" & reads$taxid != "0", "Classified",
        ifelse(reads$status == "C" & reads$taxid == "0", "QC-filtered", "Never aligned")
      )
      reads$effective_status <- ifelse(reads$effective_classified, "Classified", "Unclassified")

      accounting_plot <- data.frame(
        Category = factor(c("Raw U", "C + TaxID 0", "TaxID > 0"),
                          levels = c("Raw U", "C + TaxID 0", "TaxID > 0")),
        Count = c(n_status_U, n_qc_reclass, n_eff_class)
      ) %>%
        mutate(Fraction = Count / sum(Count),
               Label = sprintf("%s\n%s (%.1f%%)", Category, scales::comma(Count), 100 * Fraction))

      p0a <- ggplot(accounting_plot, aes(x = 3, y = Count, fill = Category)) +
        geom_col(width = 1, color = "white", linewidth = 1.2) +
        coord_polar(theta = "y") +
        xlim(c(1, 4)) +
        geom_text(data = accounting_plot[accounting_plot$Count > 0, , drop = FALSE],
                  aes(x = 3.5, label = Label),
                  position = position_stack(vjust = 0.5), size = 3.1) +
        annotate("text", x = 1, y = 0,
                 label = sprintf("%s\nreads", scales::comma(exp_total)),
                 size = 4.2, fontface = "bold") +
        scale_fill_manual(
          values = c("Raw U" = "#bdbdbd", "C + TaxID 0" = "#e6ab02", "TaxID > 0" = "#1b9e77"),
          labels = c(
            "Raw U" = "Never aligned/classified against the reference",
            "C + TaxID 0" = "Aligned, then failed identity and/or reference-coverage QC",
            "TaxID > 0" = "Effectively classified; included in classified denominators"
          )
        ) +
        theme_void() +
        labs(title = sprintf("Exact Read Accounting: %s", sample_id), fill = NULL) +
        theme(legend.position = "bottom",
              plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
    }
    p0a_path <- file.path(sample_out_dir, "00a_read_accounting_donut.png")
    save_plot(p0a_path, p0a, width = 7.5, height = 6)
    all_outputs <- c(all_outputs, p0a_path)

    cats_to_summarize <- if (is_kraken) c("All", "Classified", "Unclassified") else c("All", "Classified", "Unclassified", "QC-filtered")
    for (cat_name in cats_to_summarize) {
      sub_lens <- if (cat_name == "All") {
        reads$read_length
      } else if (cat_name %in% c("Classified", "Unclassified")) {
        reads$read_length[reads$effective_status == cat_name]
      } else {
        reads$read_length[reads$status_category == cat_name]
      }

      if (length(sub_lens) > 0) {
        length_summary_rows[[paste(sample_id, cat_name, sep = "_")]] <- data.frame(
          SampleID = sample_id,
          Category = cat_name,
          Count = length(sub_lens),
          Min = min(sub_lens),
          Q25 = as.numeric(quantile(sub_lens, 0.25)),
          Median = median(sub_lens),
          Mean = mean(sub_lens),
          Q75 = as.numeric(quantile(sub_lens, 0.75)),
          Max = max(sub_lens),
          SD = sd(sub_lens),
          stringsAsFactors = FALSE
        )
      }
    }

    # 3. Figure 1a: Donut Plot
    donut_df <- reads %>%
      count(effective_status) %>%
      mutate(
        frac = n / sum(n),
        ymax = cumsum(frac),
        ymin = c(0, head(ymax, -1)),
        label = sprintf("%s\n%s reads\n(%.1f%%)", effective_status, scales::comma(n), 100 * frac)
      )

    center_donut_count <- if (is_kraken) kraken_acct$N else sum(donut_df$n)
    p1a <- ggplot(donut_df, aes(ymin = ymin, ymax = ymax, xmin = 3, xmax = 4, fill = effective_status)) +
      geom_rect(color = "white", linewidth = 1.2) +
      coord_polar(theta = "y") +
      xlim(c(1, 4)) +
      scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
      annotate("text", x = 1, y = 0, label = sprintf("%s\nreads", scales::comma(center_donut_count)),
               size = 4.2, fontface = "bold") +
      theme_void() +
      labs(title = sprintf("Read Classification Outcome: %s", sample_id), fill = NULL) +
      theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
      geom_text(aes(x = 3.5, y = (ymin + ymax) / 2, label = label), inherit.aes = FALSE,
                data = donut_df, color = "black", size = 3.3, fontface = "bold")

    p1a_path <- file.path(sample_out_dir, "01a_classification_donut.png")
    save_plot(p1a_path, p1a, width = 5.5, height = 5.5)
    all_outputs <- c(all_outputs, p1a_path)

    # 4. Figure 1b: Read Length Histogram
    hist_subtitle <- if (is_16s) {
      if (!is.null(min_len_target) && !is.null(max_len_target)) {
        sprintf("Expected full-length 16S ~1400-1550 bp | Upstream target [%d, %d] bp",
                min_len_target, max_len_target)
      } else {
        "Expected full-length 16S ~1400-1550 bp"
      }
    } else {
      if (!is.null(min_len_target) && !is.null(max_len_target)) {
        sprintf("Target amplicon [%d, %d] bp", min_len_target, max_len_target)
      } else {
        "Amplicon read length distribution"
      }
    }

    p1b <- ggplot(reads, aes(x = read_length, fill = effective_status)) +
      geom_histogram(binwidth = 10, alpha = 0.85, position = "identity") +
      scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
      coord_cartesian(xlim = c(display_min, display_max)) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = sprintf("Read Length Distribution: %s", sample_id),
        subtitle = hist_subtitle,
        x = "Read length (bp)", y = "Read count", fill = "Effective Status"
      ) +
      theme_amplicon() +
      theme(legend.position = "bottom")

    if (!is.null(min_len_target)) {
      p1b <- p1b + geom_vline(xintercept = min_len_target, linetype = "dotted", color = "grey30")
    }
    if (!is.null(max_len_target)) {
      p1b <- p1b + geom_vline(xintercept = max_len_target, linetype = "dotted", color = "grey30")
    }

    p1b_path <- file.path(sample_out_dir, "01b_read_length_distribution.png")
    save_plot(p1b_path, p1b, width = 7.5, height = 5)
    all_outputs <- c(all_outputs, p1b_path)

    # 5. Figure 1c: Diagnostic Bar
    if (is_kraken) {
      qc_diag_df <- data.frame(
        Category = factor(
          c("Kraken classified\n(raw C)",
            "Bracken species\nestimates (E)",
            "Accounting diff\n(D = C - E)",
            "Unclassified\n(raw U)"),
          levels = c("Kraken classified\n(raw C)",
                     "Bracken species\nestimates (E)",
                     "Accounting diff\n(D = C - E)",
                     "Unclassified\n(raw U)")
        ),
        Count = c(kraken_acct$C, kraken_acct$E, kraken_acct$D, kraken_acct$U)
      )

      p1c <- ggplot(qc_diag_df, aes(x = Category, y = Count, fill = Category)) +
        geom_col(width = 0.6) +
        geom_text(aes(label = scales::comma(Count)), vjust = -0.4, size = 3.6, fontface = "bold") +
        scale_fill_manual(values = c(
          "Kraken classified\n(raw C)" = "#1b9e77",
          "Bracken species\nestimates (E)" = "#386cb0",
          "Accounting diff\n(D = C - E)" = "#fdb462",
          "Unclassified\n(raw U)" = "#bdbdbd"
        )) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
        labs(
          title = sprintf("Raw Kraken Classification vs. Bracken Estimates (%s)", sample_id),
          subtitle = sprintf(
            "Raw classified C=%s, Bracken estimates E=%s, Count difference D=%s",
            scales::comma(kraken_acct$C), scales::comma(kraken_acct$E), scales::comma(kraken_acct$D)
          ),
          x = NULL, y = "Count"
        ) +
        theme_amplicon() +
        theme(legend.position = "none")
    } else {
      qc_diag_df <- data.frame(
        Category = factor(
          c("True classified\n(aligned + passed QC)",
            "QC-filtered\n(aligned but failed\ncoverage/identity)",
            "Never aligned\n(status = U)"),
          levels = c("True classified\n(aligned + passed QC)",
                     "QC-filtered\n(aligned but failed\ncoverage/identity)",
                     "Never aligned\n(status = U)")
        ),
        Count = c(n_eff_class, n_qc_reclass, n_status_U)
      )

      p1c <- ggplot(qc_diag_df, aes(x = Category, y = Count, fill = Category)) +
        geom_col(width = 0.6) +
        geom_text(aes(label = scales::comma(Count)), vjust = -0.4, size = 3.6, fontface = "bold") +
        scale_fill_manual(values = c(
          "True classified\n(aligned + passed QC)" = "#1b9e77",
          "QC-filtered\n(aligned but failed\ncoverage/identity)" = "#e6ab02",
          "Never aligned\n(status = U)" = "#bdbdbd"
        )) +
        scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15))) +
        labs(
          title = sprintf("Raw Alignment Status vs. Effective Classification (%s)", sample_id),
          subtitle = sprintf(
            "Alignment QC filter reclassifies %s reads to taxid=0 without altering raw status letter",
            scales::comma(n_qc_reclass)
          ),
          x = NULL, y = "Read count"
        ) +
        theme_amplicon() +
        theme(legend.position = "none")
    }

    p1c_path <- file.path(sample_out_dir, "01c_qc_filter_diagnostic.png")
    save_plot(p1c_path, p1c, width = 7, height = 5.5)
    all_outputs <- c(all_outputs, p1c_path)

    # 6. Figure 1d: Violin & Boxplot of Read Length by Effective Class
    violin_subtitle <- if (is_16s) {
      "Full-length 16S amplicon read length comparison"
    } else {
      "Amplicon read length comparison"
    }

    p1d <- ggplot(reads, aes(x = effective_status, y = read_length, fill = effective_status)) +
      geom_violin(alpha = 0.6, trim = FALSE) +
      geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
      scale_fill_manual(values = c("Classified" = "#1b9e77", "Unclassified" = "#bdbdbd")) +
      coord_cartesian(ylim = c(display_min, display_max)) +
      scale_y_continuous(labels = scales::comma) +
      labs(
        title = sprintf("Read Length Distribution by Classification (%s)", sample_id),
        subtitle = violin_subtitle,
        x = "Effective Status", y = "Read length (bp)"
      ) +
      theme_amplicon() +
      theme(legend.position = "none")

    p1d_path <- file.path(sample_out_dir, "01d_read_length_violin.png")
    save_plot(p1d_path, p1d, width = 6, height = 5)
    all_outputs <- c(all_outputs, p1d_path)
  }

  # Export summary TSVs
  accounting_df <- do.call(rbind, accounting_rows)
  stopifnot(all(accounting_df$AbundanceClassified + accounting_df$AbundanceUnclassified ==
                  accounting_df$AbundanceTotal))
  available <- accounting_df$AssignmentAvailable

  if (!is_kraken) {
    stopifnot(all(accounting_df$RawC[available] + accounting_df$RawU[available] ==
                    accounting_df$AbundanceTotal[available]))
    stopifnot(all(accounting_df$RawU[available] + accounting_df$C_TaxID0[available] ==
                    accounting_df$AbundanceUnclassified[available]))
    stopifnot(all(accounting_df$TaxID_GT0[available] == accounting_df$AbundanceClassified[available]))
  } else {
    stopifnot(all(accounting_df$RawTotal_N == accounting_df$RawClassified_C + accounting_df$RawUnclassified_U))
    stopifnot(all(accounting_df$Difference_D == accounting_df$RawClassified_C - accounting_df$BrackenEstimate_E))
    stopifnot(all(accounting_df$AbundanceTotal == accounting_df$BrackenEstimate_E + accounting_df$RawUnclassified_U))
    if (any(available)) {
      stopifnot(all(accounting_df$RawC[available] == accounting_df$RawClassified_C[available]))
      stopifnot(all(accounting_df$RawU[available] == accounting_df$RawUnclassified_U[available]))
      stopifnot(all(accounting_df$C_TaxID0[available] == 0L))
      stopifnot(all(accounting_df$TaxID_GT0[available] == accounting_df$RawClassified_C[available]))
    }
  }

  accounting_file <- file.path(qc_dir, "00_read_accounting.tsv")
  write.table(accounting_df, accounting_file, sep = "\t", row.names = FALSE, quote = FALSE)
  all_outputs <- c(all_outputs, accounting_file)

  investigation_df <- do.call(rbind, investigation_rows)
  if (!is_kraken) {
    with_partition <- investigation_df$BamstatsAvailable & investigation_df$AssignmentAvailable
    if (any(with_partition)) {
      stopifnot(all(
        investigation_df$BamstatsC0Matched[with_partition] ==
          investigation_df$IdentityOnlyFailed[with_partition] +
          investigation_df$RefCoverageOnlyFailed[with_partition] +
          investigation_df$BothFailed[with_partition]
      ))
    }
  } else {
    stopifnot(all(investigation_df$BamstatsAvailable == FALSE))
    stopifnot(all(is.na(investigation_df$BamstatsC0Matched)))
    stopifnot(all(is.na(investigation_df$IdentityOnlyFailed)))
    stopifnot(all(is.na(investigation_df$RefCoverageOnlyFailed)))
    stopifnot(all(is.na(investigation_df$BothFailed)))
  }

  investigation_file <- file.path(qc_dir, "00_read_investigation.tsv")
  write.table(investigation_df, investigation_file, sep = "\t", row.names = FALSE, quote = FALSE,
              na = "NA")
  all_outputs <- c(all_outputs, investigation_file)

  reconciliation_file <- file.path(qc_dir, "classification_reconciliation.tsv")
  if (length(reconciliation_rows) > 0) {
    reconciliation_df <- do.call(rbind, reconciliation_rows)
    stopifnot(all(reconciliation_df$ReconciliationPass == TRUE))
    write.table(reconciliation_df, reconciliation_file,
                sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, reconciliation_file)
  }

  length_file <- file.path(qc_dir, "read_length_by_status.tsv")
  if (length(length_summary_rows) > 0) {
    write.table(do.call(rbind, length_summary_rows), length_file,
                sep = "\t", row.names = FALSE, quote = FALSE)
    all_outputs <- c(all_outputs, length_file)
  }

  list(
    status = "completed",
    outputs = all_outputs
  )
}
