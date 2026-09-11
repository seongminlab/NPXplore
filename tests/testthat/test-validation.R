test_that("validation accepts minimal NPX data", {
  data <- data.frame(SampleID = "S1", SampleType = "SAMPLE", AssayType = "assay", Assay = "A", NPX = 1)
  meta <- data.frame(SampleID = "S1", Condition = "A")
  expect_true(npx_validate_inputs(data, meta))
})

test_that("validation catches duplicate metadata", {
  data <- data.frame(SampleID = "S1", SampleType = "SAMPLE", AssayType = "assay", Assay = "A", NPX = 1)
  meta <- data.frame(SampleID = c("S1", "S1"), Condition = "A")
  expect_error(npx_validate_inputs(data, meta), "one row")
})

test_that("NPX input preparation requires user-read data", {
  meta <- data.frame(SampleID = "S1", Condition = "A")
  expect_error(
    npx_read_inputs("not-an-npx-file", meta),
    "NPX data frame returned by OlinkAnalyze::read_npx"
  )
})

test_that("metadata input requires user-read data frame", {
  data <- data.frame(
    SampleID = "S1", SampleType = "SAMPLE", AssayType = "assay",
    Assay = "A", NPX = 1
  )
  expect_error(
    npx_read_inputs(data, "metadata.csv"),
    "Read the CSV first with read.csv"
  )
})

test_that("missing SampleID gives an actionable message", {
  data <- data.frame(SampleID = "S1", SampleType = "SAMPLE", AssayType = "assay", Assay = "A", NPX = 1)
  meta <- data.frame(PatientID = "S1", Condition = "A")
  expect_error(npx_validate_inputs(data, meta), "sample_id_col")
  expect_true(npx_validate_inputs(data, meta, sample_id_col = "PatientID"))
})

test_that("condition count selects two-condition analysis", {
  data <- data.frame(SampleID = rep(c("S1", "S2", "S3", "S4"), each = 1),
                     SampleType = "SAMPLE", AssayType = "assay", Assay = "A",
                     OlinkID = "OID12345", UniProt = "P1", Panel = "P",
                     NPX = c(1, 2, 3, 4), Condition = rep(c("A", "B"), 2))
  expect_warning(
    expect_error(npx_analyze(data, test = "bad"), "'arg' should be one of"),
    "control condition"
  )
})

test_that("analysis API uses variable and optional pair_id", {
  expect_identical(formals(npx_ttest)$variable, "Condition")
  expect_identical(formals(npx_wilcox)$variable, "Condition")
  expect_identical(formals(npx_anova)$variable, "Condition")
  expect_true(identical(npx_wilcon, npx_wilcon))
  expect_identical(formals(npx_analyze)$variable, "Condition")
  expect_null(formals(npx_analyze)$pair_id)
  expect_identical(formals(npx_pipeline)$variable, "Condition")
  expect_null(formals(npx_pipeline)$pair_id)
  expect_identical(formals(npx_pipeline)$include_warn, TRUE)
  expect_identical(formals(npx_pipeline)$include_fail, FALSE)
  expect_null(formals(npx_pipeline)$test)
  expect_identical(formals(npx_pipeline)$logFC, 0)
  expect_identical(formals(npx_pipeline)$`adj.p`, 0.05)
  expect_identical(formals(npx_pipeline)$run_functional_analysis, TRUE)
  expect_identical(formals(npx_heatmap)$scale_rows, FALSE)
  expect_silent(npx_prepare(
    data.frame(SampleID = "S1", SampleType = "SAMPLE", AssayType = "assay", Assay = "A", NPX = 1),
    data.frame(SampleID = "S1", Condition = "A"), include_warn = TRUE
  ))
  expect_silent(npx_prepare(
    data.frame(SampleID = "S1", SampleType = "SAMPLE", AssayType = "assay", Assay = "A", NPX = 1),
    data.frame(SampleID = "S1", Condition = "A"), include_warn = FALSE
  ))
})

test_that("analysis method and gene identifier defaults are configurable", {
  expect_null(formals(npx_analyze)$test)
  expect_identical(NPXplore:::npx_resolve_gene_ids(
    data.frame(Assay = c("GENE1", "GENE2"), UniProt = c("P1", "P2")),
    logFC = 0, adj.p = 0.05
  ), c("GENE1", "GENE2"))
  expect_identical(NPXplore:::npx_resolve_gene_ids(
    data.frame(UniProt = c("P1", "P2")), logFC = 0, adj.p = 0.05
  ), c("P1", "P2"))
})

test_that("DEP cutoffs use absolute Log2FC and adjusted p-value", {
  tab <- data.frame(Assay = c("A", "B", "C"), Log2FC = c(0.5, 1, -2),
                    Adjusted_pval = c(0.001, 0.04, 0.2))
  deps <- NPXplore:::npx_filter_deps(tab, logFC = 1, adj.p = 0.05)
  expect_identical(deps$Assay, "B")
  expect_error(NPXplore:::npx_validate_dep_cutoffs(logFC = -1), "non-negative")
  expect_error(NPXplore:::npx_validate_dep_cutoffs(adj.p = 0), "greater than 0")
})

test_that("visualization functions return plotting objects", {
  d <- tibble::tibble(
    SampleID = rep(paste0("S", 1:6), each = 3),
    Assay = rep(c("A", "B", "C"), 6),
    NPX = seq_len(18),
    Condition = rep(c("Healthy", "Disease"), each = 9)
  )
  individual <- npx_boxplot(d, assays = c("A", "B"))
  expect_length(individual$plots, 2)
  expect_true(all(vapply(individual$plots, inherits, logical(1), "ggplot")))
  expect_true(inherits(npx_boxplot(d, assays = c("A", "B"), mode = "combined")$plot, "ggplot"))
  dep_analysis <- list(
    method = "ttest_independent",
    deps = data.frame(Assay = "A", p.value = 0.01, Adjusted_pval = 0.02),
    result = data.frame(Assay = "A", p.value = 0.01, Adjusted_pval = 0.02)
  )
  dep_bp <- npx_dep_boxplot(d, dep_analysis, variable = "Condition")
  expect_true(inherits(dep_bp$plots[["A"]], "ggplot"))
  expect_identical(c(dep_bp$width, dep_bp$height), c(6, 9))
  um <- npx_umap(d)
  expect_true(inherits(um$plot, "ggplot"))
  expect_identical(c(um$width, um$height), c(6, 5))
  expect_true(inherits(npx_heatmap(d, assays = c("A", "B"))$matrix, "matrix"))
  distribution <- npx_distribution_boxplot(d)
  expect_true(inherits(distribution$plot, "ggplot"))
  expect_identical(c(distribution$width, distribution$height), c(16, 10))
  iqr <- npx_iqrplot(d)
  expect_true(inherits(iqr$plot, "ggplot"))
  expect_true(all(c("SampleID", "SampleMedian", "IQR_NPX", "N_Assay") %in% names(iqr$data)))
  expect_identical(c(iqr$width, iqr$height), c(10, 8))
})

test_that("control conditions are first in plots and heatmap annotations", {
  d <- tibble::tibble(
    SampleID = rep(paste0("S", 1:6), each = 3),
    Assay = rep(c("A", "B", "C"), 6),
    NPX = seq_len(18),
    Condition = rep(c("Disease", "Healthy"), each = 9)
  )
  bp <- npx_boxplot(d, assays = "A")
  expect_identical(levels(bp$data$Condition), c("Healthy", "Disease"))
  um <- npx_umap(d)
  expect_identical(levels(um$data$Condition), c("Healthy", "Disease"))
  hm <- npx_heatmap(d, assays = c("A", "B"))
  expect_identical(levels(hm$annotation$Condition), c("Healthy", "Disease"))
  expect_true(all(as.character(head(hm$annotation$Condition, 3)) == "Healthy"))
})

test_that("missing control preserves order and warns", {
  expect_warning(
    ord <- NPXplore:::npx_condition_levels(c("Group2", "Group1")),
    "control condition"
  )
  expect_identical(ord$plot, c("Group2", "Group1"))
  expect_identical(ord$test, c("Group2", "Group1"))
})

test_that("QC reports count all assay statuses for issue samples", {
  qc_data <- data.frame(
    SampleID = c(rep("S1", 4), "S2"),
    SampleType = "SAMPLE",
    AssayType = "assay",
    Assay = c("A", "B", "C", "D", "A"),
    Block = 1,
    AssayQC = "PASS",
    SampleQC = c("PASS", "WARN", "FAIL", "NA", "PASS"),
    Normalization = "NORMALIZED",
    Panel = "P",
    UniProt = paste0("P", 1:5),
    OlinkID = paste0("OID", 1:5),
    stringsAsFactors = FALSE
  )
  reports <- npx_qc_reports(qc_data)
  expect_identical(reports$warning_sample_statistics$SampleID, "S1")
  expect_identical(
    unname(as.integer(reports$warning_sample_statistics[1, c("PassAssay", "WarnedAssay", "FailedAssay", "ExcludedAssay")])),
    c(1L, 1L, 1L, 1L)
  )
})

test_that("sample QC uses all-evaluated-assays block failures", {
  qc_data <- data.frame(
    SampleID = c(rep("S_FAIL", 4), rep("S_WARN", 4)),
    SampleType = "SAMPLE",
    AssayType = "assay",
    Assay = rep(c("A", "B"), 4),
    Block = rep(c(1, 1, 2, 2), 2),
    AssayQC = "PASS",
    SampleQC = c("FAIL", "FAIL", "FAIL", "NA", "FAIL", "PASS", "WARN", "PASS"),
    Normalization = "NORMALIZED",
    Panel = "P",
    UniProt = paste0("P", 1:8),
    OlinkID = paste0("OID", 1:8),
    stringsAsFactors = FALSE
  )
  reports <- npx_qc_reports(qc_data)
  status <- reports$sample_qc[match(c("S_FAIL", "S_WARN"), reports$sample_qc$SampleID), ]
  expect_identical(as.character(status$QC), c("FAIL", "WARN"))
  expect_identical(as.integer(status$FailBlock), c(2L, 0L))
  expect_identical(as.integer(status$WarnBlock), c(0L, 2L))
  block_fail <- reports$block_qc |>
    dplyr::filter(SampleID == "S_FAIL") |>
    dplyr::pull(BlockQC)
  expect_identical(as.character(block_fail), c("FAIL", "FAIL"))

  qc_dir <- tempfile("npxplore-qc-")
  written <- NPXplore:::npx_write_qc(
    reports,
    qc_dir,
    meta = data.frame(SampleID = c("S_FAIL", "S_WARN"))
  )
  saved_sample_qc <- utils::read.csv(written$sample_qc, check.names = FALSE)
  expect_false("QCIssue" %in% names(saved_sample_qc))
})

test_that("pathway heatmaps use DEPs and show assay names", {
  skip_if_not_installed("pheatmap")
  data <- data.frame(
    SampleID = rep(c("S1", "S2", "S3", "S4"), each = 2),
    SampleType = "SAMPLE",
    AssayType = "assay",
    Assay = rep(c("Assay A", "Assay B"), 4),
    OlinkID = rep(c("OID1", "OID2"), 4),
    NPX = seq_len(8),
    Condition = rep(c("Control", "Disease"), each = 4),
    stringsAsFactors = FALSE
  )
  analysis <- list(
    deps = data.frame(OlinkID = c("OID1", "OID2"), Assay = c("Assay A", "Assay B"))
  )
  anno <- data.frame(
    OlinkID = c("OID1", "OID2"),
    Human_Pathway = c("Pathway 1", "Pathway 1")
  )
  out_dir <- tempfile("npxplore-pathway-")
  result <- npx_pathway_heatmaps(data, analysis, anno, output_dir = out_dir)
  expect_length(result$files, 1L)
  expect_true(file.exists(result$files[[1]]))
  expect_identical(result$mappings[["Pathway 1"]], c("Assay A", "Assay B"))
})

test_that("volcano plots use all results and save DEP outputs", {
  skip_if_not_installed("ggplot2")
  analysis <- list(
    levels = c("Control", "Disease"),
    result = data.frame(
      Assay = c("Assay A", "Assay B", "Assay C"),
      Log2FC = c(1, -1, 0),
      Adjusted_pval = c(0.01, 0.02, 0.5)
    )
  )
  out_dir <- tempfile("npxplore-volcano-")
  result <- npx_volcano(analysis, output_dir = out_dir)
  expect_true(all(c("Up (n=1)", "Down (n=1)", "NS (n=1)") %in% as.character(result$data$Significance)))
  expect_true(file.exists(result$files[["pdf"]]))
  expect_true(file.exists(result$files[["png"]]))
})

test_that("ANOVA cluster heatmap and lineplot save separate outputs", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("ggplot2")
  data <- data.frame(
    SampleID = rep(paste0("S", 1:4), each = 3),
    Assay = rep(c("Assay A", "Assay B", "Assay C"), 4),
    NPX = c(1, 1.2, 4, 2, 2.1, 3.8, 5, 5.2, 2, 6, 6.1, 1.8),
    Condition = rep(c("Control", "Disease"), each = 6)
  )
  cluster_result <- list(
    matrix = matrix(
      c(1, 2, 5, 6, 1.2, 2.1, 5.2, 6.1, 4, 3.8, 2, 1.8),
      nrow = 3,
      dimnames = list(c("Assay A", "Assay B", "Assay C"), paste0("S", 1:4))
    ),
    clusters = data.frame(
      Assay = c("Assay A", "Assay B", "Assay C"),
      cluster = factor(c(1, 1, 2))
    )
  )
  analysis <- list(method = "anova", deps = data.frame(Assay = c("Assay A", "Assay B", "Assay C")))
  out_dir <- tempfile("npxplore-clustering-")
  heatmap <- npx_cluster_heatmap(data, analysis, cluster_result = cluster_result, output_dir = out_dir)
  lineplot <- npx_cluster_lineplot(data, analysis, cluster_result = cluster_result, output_dir = out_dir)
  expect_true(file.exists(heatmap$files[["pdf"]]))
  expect_true(file.exists(heatmap$files[["png"]]))
  expect_true(file.exists(lineplot$files[["pdf"]]))
  expect_true(file.exists(lineplot$files[["png"]]))
  expect_identical(c(lineplot$width, lineplot$height), c(6, 8))
})

test_that("table output schemas are method-specific", {
  ttest_analysis <- list(
    method = "wilcox_independent",
    test_levels = c("Disease", "Healthy")
  )
  ttest_data <- data.frame(
    Assay = "A", OlinkID = "OID1", UniProt = "P1", Panel = "Panel1",
    estimate = 1, Disease = 2, Healthy = 1, p.value = 0.01,
    Adjusted_pval = 0.02, extra = "remove"
  )
  ttest_output <- NPXplore:::npx_statistical_table_for_output(ttest_data, ttest_analysis)
  expect_identical(
    names(ttest_output),
    c("Assay", "OlinkID", "UniProt", "Panel", "Log2FC", "Disease", "Healthy", "p.value", "Adjusted_pval")
  )

  anova_analysis <- list(method = "anova")
  anova_data <- data.frame(
    Assay = "A", OlinkID = "OID1", UniProt = "P1", term = "Condition",
    df = 1, sumsq = 2, meansq = 2, statistic = 3, p.value = 0.01,
    Adjusted_pval = 0.02, extra = "remove"
  )
  anova_output <- NPXplore:::npx_statistical_table_for_output(anova_data, anova_analysis)
  expect_identical(
    names(anova_output),
    c("Assay", "OlinkID", "UniProt", "term", "df", "sumsq", "meansq", "statistic", "p.value", "Adjusted_pval")
  )

  posthoc_output <- NPXplore:::npx_posthoc_table_for_output(
    data.frame(Assay = "A", OlinkID = "OID1", UniProt = "P1", term = "Condition",
               contrast = "Disease - Healthy", estimate = 1, Adjusted_pval = 0.02, extra = "remove")
  )
  expect_identical(
    names(posthoc_output),
    c("Assay", "OlinkID", "UniProt", "term", "contrast", "Log2FC", "Adjusted_pval")
  )
})
