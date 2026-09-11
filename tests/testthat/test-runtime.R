
test_that("empty results can be exported with the complete schema", {
  expect_equal(nrow(npx_posthoc_table_for_output(data.frame())), 0)
  expect_equal(nrow(npx_statistical_table_for_output(data.frame(), list(method = "anova"))), 0)
})

test_that("real Olink two-group tests and ANOVA run on synthetic data", {
  d <- runtime_fixture()
  two <- d[d$Condition != "Followup", ]
  for (method in c("ttest", "wilcox")) {
    result <- npx_analyze(two, test = method)
    expect_equal(nrow(result$result), 8)
    expect_equal(nrow(result$deps), 8)
  }
  result <- npx_analyze(d)
  expect_identical(result$method, "anova")
  expect_equal(nrow(result$deps), 8)
})

test_that("sparse and constant assays do not abort all t-tests", {
  d <- runtime_fixture()
  d <- d[d$Condition != "Followup", ]
  rows <- which(d$Assay == "A1" & d$Condition == "Healthy")
  d$NPX[rows[-1]] <- NA
  d$NPX[d$Assay == "A2"] <- 1
  result <- npx_ttest(d)
  expect_equal(nrow(result$result), 6)
  expect_setequal(result$excluded_assays, c("OID10001", "OID10002"))
  expect_true(all(nzchar(result$excluded_assay_reasons$Reason)))
})

test_that("paired preflight requires complete variable pairs", {
  d <- runtime_fixture()
  d <- d[d$Condition != "Followup", ]
  d$Subject <- rep(rep(paste0("P", 1:6), each = 8), 2)
  d$NPX[d$Assay == "A1" & d$Subject != "P1"] <- NA
  result <- npx_ttest(d, pair_id = "Subject")
  expect_equal(nrow(result$result), 7)
  expect_true("OID10001" %in% result$excluded_assays)
})

test_that("unsupported model arguments cannot be silently ignored", {
  d <- runtime_fixture()
  expect_error(npx_analyze(d, pair_id = "Subject"), "not supported for ANOVA")
  expect_error(npx_analyze(d[d$Condition != "Followup", ], covariates = "Age"), "only for ANOVA")
})

test_that("stage diagnostics distinguish insufficient input and runtime errors", {
  insufficient <- structure(list(message = "Too few DEPs", call = NULL), class = c("npx_insufficient_input", "error", "condition"))
  expect_identical(npx_run_stage(stop(insufficient), "test")$status, "Skipped")
  expect_identical(npx_run_stage(stop("broken"), "test")$status, "Failed")
  expect_identical(npx_run_stage(list(files = character()), "test")$status, "Success")
})


test_that("test-only boxplot selection preserves full package defaults", {
  a <- list(deps = data.frame(Assay = c("A", "B", "C"), Adjusted_pval = c(0.03, 0.01, 0.02)))
  expect_null(formals(npx_pipeline)$boxplot_top_n)
  expect_null(npx_boxplot_selection(a))
  expect_identical(npx_boxplot_selection(a, 2), c("B", "C"))
  expect_equal(nrow(a$deps), 3)
  a$deps <- a$deps[FALSE, ]
  expect_identical(npx_boxplot_selection(a, 20), character())
})

test_that("pipeline retains diagnostics after failed posthoc and insufficient clusters", {
  d <- runtime_fixture()
  meta <- unique(d[c("SampleID", "Condition")])
  analysis <- npx_anova(d)
  analysis$deps <- analysis$deps[1:3, ]
  d$Condition <- NULL
  ppi_calls <- 0L
  local_mocked_bindings(
    npx_read_inputs = function(...) list(data = d, meta = meta, anno = data.frame()),
    npx_clean = function(...) list(data = d, check_log = NULL, check_log_clean = NULL),
    npx_analyze = function(...) analysis,
    npx_umap = function(...) NULL, npx_iqrplot = function(...) NULL,
    npx_distribution_boxplot = function(...) NULL,
    npx_heatmap = function(...) NULL, npx_dep_boxplot = function(...) NULL,
    npx_anova_posthoc = function(...) stop("posthoc unavailable"),
    npx_functional_analysis = function(...) stop("GO must not be called"),
    npx_ppi_network = function(analysis, ...) {
      ppi_calls <<- ppi_calls + 1L
      expect_equal(nrow(analysis$deps), 3)
      list(status = "Success", files = character())
    }
  )
  out <- tempfile()
  expect_warning(result <- npx_pipeline(d, meta, output_dir = out), "posthoc unavailable")
  expect_identical(result$cluster$status, "Skipped")
  expect_identical(result$cluster$code, "InsufficientInput")
  expect_identical(result$functional$code, "ClusteringUnavailable")
  expect_identical(result$ppi$status, "Success")
  expect_identical(ppi_calls, 1L)
  expect_true(file.exists(result$output_files$functional_status))
  expect_equal(nrow(read.csv(result$output_files$anova_posthoc)), 0)
  expect_setequal(list.files(out), c("metadata_used.csv", "QCreport_metadata.csv", "Results", "QC", "Tables", "log"))
  expect_identical(dirname(result$output_files$functional_status), file.path(out, "log"))
  expect_true(file.exists(result$output_files$session_info))
  expect_match(paste(readLines(result$output_files$run_log), collapse = "\n"), "posthoc unavailable")
  messages <- paste(readLines(result$output_files$run_log), collapse = "\n")
  expect_match(messages, "Assay boxplots started")
  expect_match(messages, "Assay boxplots failed")
  expect_match(messages, "GO functional analysis skipped")
  expect_match(messages, "PPI network analysis started")
  expect_match(messages, "PPI network analysis complete")
  expect_match(messages, "ANOVA analysis complete.*Assay boxplots=Failed; GO=Skipped; PPI=Success")

})

test_that("pipeline rejects missing model metadata before analysis and plots", {
  d <- runtime_fixture()
  meta <- unique(d[c("SampleID", "Condition")]); meta$Condition[1] <- " "
  d$Condition <- NULL
  local_mocked_bindings(
    npx_read_inputs = function(...) list(data = d, meta = meta, anno = data.frame()),
    npx_clean = function(...) list(data = d, check_log = NULL, check_log_clean = NULL)
  )
  expect_error(npx_pipeline(d, meta, output_dir = NULL), "Missing or blank model metadata")
})

test_that("actual NA QC entries do not propagate NA into sample status", {
  d <- runtime_fixture()
  d$SampleQC[1] <- NA_character_; d$SampleQC[2] <- "Warning"
  qc <- npx_sample_qc_status(d)
  expect_false(anyNA(qc$sample$QC))
  expect_identical(qc$sample$QC[qc$sample$SampleID == "S1"], "WARN")
})


test_that("print dispatches to the compact NPXplore summary", {
  result <- structure(list(analysis = list(method = "ttest_independent", levels = c("Healthy", "Disease")),
                           data = data.frame(Assay = "A")), class = "npxplore_result")
  text <- capture.output(print(result))
  expect_match(text[1], "NPXplore result")
  expect_lte(length(text), 5)
})
