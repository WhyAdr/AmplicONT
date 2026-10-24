# Offline adapter and real child-process contracts; no mapper assets required.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript tests/run_faprotax2_unit.R")
repo <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/")))
setwd(repo)
result <- testthat::test_dir("tests/testthat", filter = "^faprotax2$", reporter = "summary")
summary <- as.data.frame(result)
cat(sprintf("FAPROTAX2 contracts: %d tests, %d expectations, %d failures, %d errors.\n",
            nrow(summary), sum(summary$nb), sum(summary$failed), sum(summary$error)))
if (any(summary$failed > 0 | summary$error)) quit(status = 1)
