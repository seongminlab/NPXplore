test_that("public QC preparation preserves original failures and metadata scope", {
  d <- runtime_fixture()
  meta <- unique(d[c("SampleID", "Condition")])
  meta <- meta[meta$SampleID != "S18", ]
  d$Block <- rep(c(1, 2), length.out = nrow(d))
  d$Condition <- NULL
  d$SampleQC[d$SampleID == "S1"] <- "FAIL"
  d$SampleQC[d$SampleID == "S2"] <- "WARN"
  inputs <- list(data = d, meta = meta)
  out <- tempfile()
  qc <- npx_qc_reports(d, meta = meta, output_dir = out)
  local_mocked_bindings(npx_clean = function(data) list(data = data, check_log = NULL, check_log_clean = NULL))
  prepared <- npx_apply_qc(inputs, qc, output_dir = out)
  expect_false("S1" %in% prepared$data$SampleID)
  expect_true("S2" %in% prepared$data$SampleID)
  expect_false("S18" %in% prepared$raw_data$SampleID)
  expect_true("S1" %in% prepared$qc_data$SampleID)
  expect_true(all(file.exists(prepared$files)))
  strict <- npx_apply_qc(inputs, qc, include_warn = FALSE)
  expect_false("S2" %in% strict$data$SampleID)
  qc$sample_qc <- qc$sample_qc[qc$sample_qc$SampleID != "S1", ]
  expect_error(npx_apply_qc(inputs, qc), "missing original sample")
})

test_that("public exports impute only the saved matrix and preserve statistical tables", {
  d <- runtime_fixture()
  d <- d[d$Condition != "Followup", ]
  analysis <- npx_ttest(d)
  d$NPX[1] <- NA_real_
  paths <- npx_export_tables(d, analysis, tempfile(), raw_data = d)
  expect_true(all(file.exists(paths)))
  expect_true(is.na(d$NPX[1]))
  matrix <- read.csv(paths[["npx"]], check.names = FALSE)
  expect_false(anyNA(matrix))
  expect_equal(nrow(read.csv(paths[["deps"]])), nrow(analysis$deps))
})

test_that("public plot export accepts result objects and long heatmap legends", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pheatmap")
  d <- runtime_fixture()
  d$Condition <- paste0(d$Condition, " very long sample group label requiring additional space")
  hm <- suppressWarnings(npx_heatmap(d))
  path <- tempfile(fileext = ".png")
  npx_save_plot(hm, path, width = 6, height = 5, dpi = 72)
  expect_gt(file.info(path)$size, 0)
  p <- ggplot2::ggplot(d, ggplot2::aes(NPX)) + ggplot2::geom_histogram(bins = 10)
  path <- tempfile(fileext = ".pdf")
  npx_save_plot(list(plot = p, width = 5, height = 4), path)
  expect_gt(file.info(path)$size, 0)
  expect_error(npx_save_plot(1, path), "plot must be")
})

 test_that("public preparation uses real Olink cleaning", {
   d <- runtime_fixture()
   meta <- unique(d[c("SampleID", "Condition")])
   d$Condition <- NULL
   prepared <- npx_apply_qc(list(data = d, meta = meta))
   expect_equal(nrow(prepared$data), nrow(d))
   expect_true("Condition" %in% names(prepared$data))
 })

test_that("NPX tables can be saved before analysis and retained during later exports", {
  d <- runtime_fixture()
  out <- tempfile()
  paths <- npx_export_tables(d, output_dir = out, raw_data = d)
  expect_setequal(names(paths), c("raw_data", "npx"))
  before <- readLines(paths[["npx"]])
  d$NPX <- d$NPX + 100
  analysis <- npx_anova(d)
  stats_paths <- npx_export_tables(d, analysis, output_dir = out, save_npx = FALSE)
  expect_false("npx" %in% names(stats_paths))
  expect_identical(readLines(paths[["npx"]]), before)
  expect_error(npx_export_tables(d, cluster_result = list()), "analysis is required")
  expect_error(npx_export_tables(d, save_npx = NA), "save_npx")
})

test_that("documented covariates work in matching ANOVA and posthoc exports", {
  d <- runtime_fixture()
  d <- tibble::as_tibble(d)
  d$Group <- d$Condition
  sample_index <- match(d$SampleID, unique(d$SampleID))
  d$SEX <- factor(ifelse(sample_index %% 2 == 0, "F", "M"))
  d$exposure <- (sample_index %% 5) + sample_index / 100
  analysis <- npx_anova(d, variable = "Group", covariates = c("SEX", "exposure"))
  posthoc <- npx_anova_posthoc(d, analysis, variable = "Group", covariates = c("SEX", "exposure"))
  expect_gt(nrow(analysis$deps), 0)
  expect_gt(nrow(posthoc), 0)
  clusters <- list(clusters = data.frame(Assay = analysis$deps$Assay,
                                        cluster = rep(1:2, length.out = nrow(analysis$deps))))
  paths <- npx_export_tables(d, analysis, output_dir = tempfile(), posthoc = posthoc,
                            cluster_result = clusters, save_npx = FALSE)
  expect_true(all(file.exists(paths)))
  expect_true("ClusterNo" %in% names(read.csv(paths[["clusters"]])))
})
