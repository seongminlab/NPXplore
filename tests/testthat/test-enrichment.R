enrichr_fixture <- function() {
  data.frame(Term = "Example term", P.value = 0.001,
             Adjusted.P.value = 0.02, Combined.Score = 5, Genes = "TP53;EGFR")
}

test_that("GO works without attaching enrichR and restores session options", {
  skip_if_not_installed("enrichR")
  old_options <- options(enrichR.base.address = NULL, enrichR.quiet = NULL,
                enrichR.live = NULL)
  on.exit(options(old_options), add = TRUE)
  local_mocked_bindings(enrichr = function(genes, databases) {
    expect_identical(getOption("enrichR.base.address"), "https://maayanlab.cloud/Enrichr/")
    expect_identical(getOption("enrichR.quiet"), FALSE)
    expect_identical(getOption("enrichR.live"), TRUE)
    expect_identical(genes, c("TP53", "EGFR"))
    setNames(list(enrichr_fixture()), databases)
  }, .package = "enrichR")
  result <- npx_go_enrichment(c(" TP53 ", "EGFR", "", " ", NA, "TP53"), "GO_test")
  expect_named(result, "GO_test")
  expect_null(getOption("enrichR.base.address"))
  expect_null(getOption("enrichR.quiet"))
  expect_null(getOption("enrichR.live"))
})

test_that("GO respects custom Enrichr options and restores them on errors", {
  skip_if_not_installed("enrichR")
  old_options <- options(enrichR.base.address = "https://example.test/Enrichr/",
                enrichR.quiet = TRUE, enrichR.live = TRUE)
  on.exit(options(old_options), add = TRUE)
  local_mocked_bindings(enrichr = function(...) {
    expect_identical(getOption("enrichR.base.address"), "https://example.test/Enrichr/")
    expect_true(getOption("enrichR.quiet"))
    options(enrichR.live = FALSE)
    stop("HTTP 503")
  }, .package = "enrichR")
  expect_error(npx_go_enrichment("TP53", "GO_test"), "HTTP 503")
  expect_true(getOption("enrichR.live"))
})

test_that("GO rejects partial or malformed service responses", {
  skip_if_not_installed("enrichR")
  local_mocked_bindings(enrichr = function(...) list(GO_test = enrichr_fixture()),
                        .package = "enrichR")
  expect_error(npx_go_enrichment("TP53", c("GO_test", "missing")), "missing")
  local_mocked_bindings(enrichr = function(...) list(GO_test = data.frame(error = "bad request")),
                        .package = "enrichR")
  expect_error(npx_go_enrichment("TP53", "GO_test"), "columns")
})

test_that("unnamed GO databases produce results and a workbook", {
  skip_if_not_installed("openxlsx")
  local_mocked_bindings(npx_go_enrichment = function(...) list(GO_test = enrichr_fixture()),
                        npx_go_save_plot = function(...) character())
  result <- NPXplore:::npx_go_analyze_set(data.frame(Assay = c("TP53", "EGFR")),
                                         "example", tempfile(), databases = "GO_test")
  expect_named(result$results, "GO_test")
  expect_equal(result$results$GO_test$GeneRatio, 1)
  expect_true(file.exists(result$files[["excel"]]))
})

test_that("GO plot cutoff selects adjusted p-values", {
  skip_if_not_installed("ggplot2")
  tab <- rbind(enrichr_fixture(), transform(enrichr_fixture(), Term = "Not significant",
                                           Adjusted.P.value = 0.2))
  tab <- NPXplore:::npx_enrichr_table(tab, c("TP53", "EGFR"))
  plot <- NPXplore:::npx_go_dotplot(tab, "GO", p_cutoff = 0.05)
  expect_equal(as.character(plot$data$Term), "Example term")
})

test_that("functional analysis retains failed sets separately from empty sets", {
  local_mocked_bindings(npx_go_enrichment = function(...) stop("service unavailable"))
  analysis <- list(levels = c("Control", "Disease"),
                   deps = data.frame(Assay = "TP53", Log2FC = 1))
  expect_warning(result <- npx_functional_analysis(analysis, output_dir = tempfile(), databases = "GO_test"),
                 "service unavailable")
  expect_identical(result$results$up_regulated$status, "Failed")
  expect_match(result$results$up_regulated$reason, "service unavailable")
  expect_identical(result$results$down_regulated$status, "Skipped")
  expect_identical(result$status, "Failed")
  expect_true(file.exists(tail(result$files, 1)))
})

test_that("pipeline does not impose a one-minute functional stage limit by default", {
  expect_identical(formals(npx_pipeline)$enrichment_timeout, quote(Inf))
  expect_equal(NPXplore:::npx_try_timeout(42, "test", seconds = Inf), 42)
  expect_warning(expect_null(NPXplore:::npx_try_timeout(stop("failure"), "test", Inf)),
                 "failure")
})

test_that("GO set wrappers propagate explicit stage time limits", {
  local_mocked_bindings(npx_go_enrichment = function(...) stop("reached elapsed time limit"))
  expect_error(NPXplore:::npx_go_run_set(data.frame(Assay = "TP53"), "example", tempfile()),
               "reached elapsed time limit")
})

test_that("GO libraries fail independently and retain successful tables", {
  skip_if_not_installed("openxlsx")
  calls <- character()
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) {
    calls <<- c(calls, databases)
    if (databases == "bad") stop("HTTP 503 service unavailable")
    setNames(list(enrichr_fixture()), databases)
  }, npx_go_save_plot = function(...) character())
  expect_warning(result <- NPXplore:::npx_go_analyze_set(data.frame(Assay = "TP53"), "example", tempfile(),
    databases = c(BP = "bad", MF = "good")), "HTTP 503")
  expect_identical(calls, c("bad", "good"))
  expect_identical(result$status, "Partial")
  expect_named(result$results, "MF")
  expect_identical(result$diagnostics$Code, c("QueryError", "Completed"))
  expect_true(file.exists(result$files[["excel"]]))
  expect_true(file.exists(result$files[["diagnostics"]]))
})

test_that("GO distinguishes empty input, empty response, nonsignificance and malformed response", {
  skip_if_not_installed("openxlsx")
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) {
    if (databases == "empty") return(list(empty = enrichr_fixture()[FALSE, ]))
    if (databases == "nonsignificant") return(list(nonsignificant = transform(enrichr_fixture(), Adjusted.P.value = 0.8)))
    list(malformed = data.frame(error = "bad response"))
  }, npx_go_save_plot = function(...) stop("No plot should be attempted"))
  expect_warning(result <- NPXplore:::npx_go_analyze_set(data.frame(Assay = "TP53"), "example", tempfile(),
    databases = c("empty", "nonsignificant", "malformed")), "InvalidResponse")
  expect_identical(result$diagnostics$Code, c("EmptyResult", "NoSignificantTerms", "InvalidResponse"))
  expect_identical(result$diagnostics$Status, c("Success", "Success", "Failed"))
  expect_equal(nrow(result$results$nonsignificant), 1)
  empty <- NPXplore:::npx_go_analyze_set(data.frame(Assay = character()), "empty", tempfile(), databases = "empty")
  expect_identical(empty$status, "Skipped")
  expect_identical(empty$diagnostics$Code, "InputEmpty")
  expect_true(file.exists(empty$files[["diagnostics"]]))
})

test_that("GO preserves tables and processes next library after rendering failure", {
  skip_if_not_installed("openxlsx")
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) setNames(list(enrichr_fixture()), databases),
    npx_go_dotplot = function(result, title, ...) { if (grepl("bad", title)) stop("plot failed"); "plot" },
    npx_go_save_plot = function(...) stop("disk full"))
  result <- NPXplore:::npx_go_analyze_set(data.frame(Assay = "TP53"), "example", tempfile(), databases = c("bad", "good"))
  expect_named(result$results, c("bad", "good"))
  expect_identical(result$diagnostics$Code, c("PlotError", "OutputError"))
  expect_identical(result$status, "Partial")
  expect_true(file.exists(result$files[["excel"]]))
})

test_that("GO returns output errors and in-memory diagnostics when workbook export fails", {
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) setNames(list(enrichr_fixture()), databases),
    npx_go_save_plot = function(...) character(), npx_go_write_workbook = function(...) stop("workbook locked"))
  result <- NPXplore:::npx_go_analyze_set(data.frame(Assay = "TP53"), "example", tempfile(), databases = "GO_test")
  expect_identical(result$status, "Partial")
  expect_equal(nrow(result$results$GO_test), 1)
  expect_true("OutputError" %in% result$diagnostics$Code)
  expect_true(file.exists(result$files[["diagnostics"]]))
})

test_that("GO selects terms by adjusted p-value and Combined Score, then plots by GeneRatio", {
  tab <- data.frame(
    Term = c("Selected lower ratio", "Excluded lower score", "Selected higher ratio", "Boundary selected", "Excluded nonsignificant"),
    P.value = 0.001,
    Adjusted.P.value = c(0.01, 0.01, 0.02, 0.05, 0.06),
    Combined.Score = c(50, 10, 40, 30, 100),
    Genes = c("TP53", "TP53;EGFR;MYC;AKT1", "TP53;EGFR;MYC", "TP53;EGFR", "TP53;EGFR;MYC;AKT1;MTOR")
  )
  plot <- NPXplore:::npx_go_dotplot(NPXplore:::npx_enrichr_table(tab, c("TP53", "EGFR", "MYC", "AKT1")),
                                    "GO", top_n = 3, p_cutoff = 0.05)
  expect_identical(as.character(plot$data$Term), c("Selected higher ratio", "Boundary selected", "Selected lower ratio"))
})

test_that("ANOVA GO reports invalid assignments and uses valid clusters separately", {
  skip_if_not_installed("openxlsx")
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) setNames(list(enrichr_fixture()), databases),
    npx_go_save_plot = function(...) character())
  analysis <- list(method = "anova", levels = c("A", "B", "C"), deps = data.frame(Assay = c("TP53", "EGFR", "MYC")))
  clusters <- list(clusters = data.frame(Assay = c("TP53", "EGFR", "MYC"), cluster = c("1", "2", NA)))
  result <- npx_functional_analysis(analysis, cluster_result = clusters, output_dir = tempfile(), databases = "GO_test")
  expect_named(result$results, c("Unassigned", "Cluster_1", "Cluster_2"))
  expect_identical(result$results$Cluster_1$genes, "TP53")
  expect_identical(result$results$Cluster_2$genes, "EGFR")
  expect_true("UnassignedCluster" %in% result$diagnostics$Code)
  clusters$clusters$cluster <- NA_character_
  invalid <- npx_functional_analysis(analysis, cluster_result = clusters, output_dir = tempfile(), databases = "GO_test")
  expect_identical(invalid$status, "Failed")
  expect_identical(invalid$diagnostics$Code, "ClusteringError")
  expect_true(file.exists(invalid$diagnostics_file))
})

test_that("GO retains successful results when output directory cannot be written", {
  target <- tempfile()
  writeLines("this is a file", target)
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) setNames(list(transform(enrichr_fixture(), Adjusted.P.value = 0.8)), databases))
  result <- suppressWarnings(NPXplore:::npx_go_analyze_set(data.frame(Assay = "TP53"), "example", target, databases = "GO_test"))
  expect_identical(result$status, "Partial")
  expect_equal(nrow(result$results$GO_test), 1)
  expect_true("OutputError" %in% result$diagnostics$Code)
  expect_null(result$diagnostics_file)
})

test_that("GO retains partial plot exports and their diagnostics", {
  skip_if_not_installed("openxlsx")
  local_mocked_bindings(npx_go_enrichment = function(genes, databases) setNames(list(enrichr_fixture()), databases),
    npx_go_save_plot = function(...) structure(c(pdf = "retained.pdf"), errors = "png: disk full"))
  result <- NPXplore:::npx_go_analyze_set(data.frame(Assay = "TP53"), "example", tempfile(), databases = "GO_test")
  expect_identical(result$status, "Partial")
  expect_identical(result$files[["pdf"]], "retained.pdf")
  expect_identical(result$diagnostics$Code, "OutputError")
  expect_match(result$diagnostics$Reason, "disk full")
})


test_that("GO results and diagnostics respect the v0.2.0 output layout", {
  out <- tempfile()
  local_mocked_bindings(
    npx_go_enrichment = function(genes, databases) stats::setNames(list(data.frame(
      Term = "Example", P.value = 0.2, Adjusted.P.value = 0.3, Combined.Score = 5, Genes = "A"
    )), databases)
  )
  result <- npx_functional_analysis(list(method = "ttest_independent", levels = c("Healthy", "Disease"),
    deps = data.frame(Assay = "A", Log2FC = 1)), output_dir = out, databases = "test")
  expect_setequal(list.files(out), c("Results", "log"))
  expect_identical(dirname(result$diagnostics_file), file.path(out, "log"))
  expect_true(file.exists(file.path(out, "Results", "GeneOntology", "GeneOntology_up_regulated.xlsx")))
  expect_true(file.exists(file.path(out, "log", "GeneOntology_up_regulated_diagnostics.csv")))
})
