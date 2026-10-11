source(file.path("..", "..", "analysis", "utils", "config.R"))
source(file.path("..", "..", "analysis", "utils", "io.R"))
source(file.path("..", "..", "analysis", "utils", "merge_wf16s_cohort.R"))

test_that("native report percentages, identities and every clade reconcile", {
  lines <- c("40.00\t2\t2\tU\t0\tunclassified", "60.00\t3\t0\tR\t1\troot",
             "60.00\t3\t0\tD\t2759\t  Eukaryota", "60.00\t3\t3\tS\t4932\t    Yeast")
  check <- function(x) { p <- tempfile(); writeLines(x, p); parse_kraken_report(p) }
  expect_equal(check(lines)$raw_total, 5)
  bad <- lines; bad[3] <- sub("3\t0", "2\t0", bad[3], fixed=TRUE)
  expect_error(check(bad), "clade accounting")
  bad <- lines; bad[1] <- sub("40.00", "NaN", bad[1], fixed=TRUE)
  expect_error(check(bad), "percentage")
  bad <- lines; bad[1] <- sub("40.00", "50.00", bad[1], fixed=TRUE)
  expect_error(check(bad), "percentage")
  bad <- lines; bad[1] <- sub("U\t0", "S\t0", bad[1], fixed=TRUE)
  expect_error(check(bad), "identity")
  bad <- lines; bad[4] <- sub("    Yeast", "      Yeast", bad[4], fixed=TRUE)
  expect_error(check(bad), "hierarchy")
})

test_that("Bracken preserves TaxID spelling and validates physical rows and estimates", {
  header <- "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads"
  check <- function(row) { p <- tempfile(); writeLines(c(header,row), p); parse_bracken_abundance(p) }
  expect_error(check("Yeast\t04932\tS\t2\t1\t3\t1"), "invalid TaxID")
  expect_error(check("Yeast\t4932\tS\t2\t1\t4\t1"), "assigned.*added")
  expect_error(check("Yeast\t4932\tS\t2\t1\t3\tNaN"), "fractions")
  expect_error(check("Yeast\t4932\tS\t2\t1\t3\t1\textra"), "exactly 7")
  expect_equal(check(character(0))$species_estimate_sum, 0)
})

test_that("species spelling cannot substitute for TaxID evidence", {
  path <- "Eukaryota;Fungi;P;C;O;F;G;Yeast"
  tax <- as.data.frame(as.list(setNames(strsplit(path,";",fixed=TRUE)[[1]], RANKS_8)))
  tax$TaxonPath <- path
  counts <- matrix(3,1,1,dimnames=list(path,"sample"))
  br <- list(rows=data.frame(taxonomy_id="4932",name="Yeast",new_est_reads=3))
  expect_error(match_bracken_to_abundance(br,counts,tax,sample_id="sample"), "without assignment or taxonomy evidence")
  expect_true(match_bracken_to_abundance(br,counts,tax,sample_id="sample",taxonomy_cache=setNames(list(4932),path)))
  assignments <- data.frame(status=c("C","C"),taxid="4932",lineage=c("a|b|c|d|e|f|g","a|b|c|d|e|f|h"))
  expect_error(match_bracken_to_abundance(br,counts,tax,assignments_df=assignments,sample_id="sample"), "Conflicting")
})

test_that("Kraken assignment lineage cannot contain a trailing extra rank", {
  p <- tempfile(); writeLines("C\tr1\t4932\t500\t4932:1\ta|b|c|d|e|f|g|",p)
  expect_error(read_assignments_file(p,"sample",classifier="kraken2"), "invalid lineage")
  for (hits in c("9007199254740992:1", "4932:9007199254740992", "4932:01")) {
    writeLines(paste("C","r1","4932","500",hits,"a|b|c|d|e|f|g",sep="\t"),p)
    expect_error(read_assignments_file(p,"sample",classifier="kraken2"), "hit list")
  }
})

test_that("repeatable merger input accepts punctuation literally", {
  opts <- parse_merge_cli_args(c("--input","D:/a,b;c","--input=D:/d,e;f","--output-dir","out"))
  expect_identical(opts$inputs,c("D:/a,b;c","D:/d,e;f"))
})
