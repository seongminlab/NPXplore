clustering_fixture <- function(n = 8L) {
  set.seed(819)
  mat <- matrix(stats::rnorm(n * 6), nrow = n)
  mat <- mat + rep(c(2, 12), length.out = n)
  rownames(mat) <- paste0("A", seq_len(n))
  colnames(mat) <- paste0("S", 1:6)
  data <- data.frame(Assay = rep(rownames(mat), times = 6),
                     SampleID = rep(colnames(mat), each = n),
                     NPX = as.vector(mat),
                     Group = rep(c("Control", "Disease", "Followup"), each = n * 2))
  list(data = data, analysis = list(method = "anova", deps = data.frame(Assay = rownames(mat))),
       matrix = mat)
}

test_that("NbClust handles modest DEP counts without requesting n centers", {
  skip_if_not_installed("NbClust")
  for (n in c(4L, 8L, 11L)) {
    input <- clustering_fixture(n)
    result <- npx_cluster_assays(input$data, input$analysis)
    expect_equal(nrow(result$clusters), n)
    expect_gte(result$optimal_k, 2L)
    expect_lte(result$optimal_k, n - 2L)
    expect_equal(result$matrix, input$matrix[rownames(result$matrix), colnames(result$matrix)])
    expect_identical(result$clustering_matrix, result$matrix)
    expect_false(result$scale_rows)
  }
})

test_that("too few assays or distinct profiles have actionable errors", {
  skip_if_not_installed("NbClust")
  for (n in 0:3) {
    input <- clustering_fixture(max(1, n))
    input$analysis$deps <- input$analysis$deps[seq_len(n), , drop = FALSE]
    expect_error(npx_cluster_assays(input$data, input$analysis), "At least four",
                  class = "npx_insufficient_input")
  }
  input <- clustering_fixture()
  input$data$NPX <- rep(seq_len(6), each = 8)
  expect_error(npx_cluster_assays(input$data, input$analysis), "distinct",
                class = "npx_insufficient_input")
})

test_that("optional assay scaling preserves NPX and standardizes assays", {
  skip_if_not_installed("NbClust")
  input <- clustering_fixture()
  result <- npx_cluster_assays(input$data, input$analysis, scale_rows = TRUE)
  expect_true(result$scale_rows)
  expect_equal(result$matrix, input$matrix[rownames(result$matrix), colnames(result$matrix)])
  expect_equal(unname(rowMeans(result$clustering_matrix)), rep(0, 8), tolerance = 1e-10)
  expect_equal(unname(apply(result$clustering_matrix, 1, stats::sd)), rep(1, 8))
  expect_identical(result$scaled_matrix, result$clustering_matrix)
})

test_that("multi-index cluster votes use the number row and deterministic ties", {
  skip_if_not_installed("NbClust")
  input <- clustering_fixture()
  local_mocked_bindings(NbClust = function(...) {
    list(Best.nc = rbind(Number_clusters = c(3, 2, 3, 2), Value_Index = c(20, 30, 10, 40)))
  }, .package = "NbClust")
  expect_equal(npx_cluster_assays(input$data, input$analysis)$optimal_k, 2)
  local_mocked_bindings(NbClust = function(...) {
    list(Best.nc = rbind(Number_clusters = c(3, 2, 3), Value_Index = c(20, 30, 10)))
  }, .package = "NbClust")
  expect_equal(npx_cluster_assays(input$data, input$analysis)$optimal_k, 3)
})

test_that("invalid cluster ranges and options fail before NbClust execution", {
  skip_if_not_installed("NbClust")
  input <- clustering_fixture()
  expect_error(npx_cluster_assays(input$data, input$analysis, list(min.nc = 7)), "range")
  expect_error(npx_cluster_assays(input$data, input$analysis, list(max.nc = NA)), "integer")
  expect_error(npx_cluster_assays(input$data, input$analysis, scale_rows = NA), "TRUE or FALSE")
})

test_that("heatmaps retain actual NPX and scale across samples per assay", {
  skip_if_not_installed("pheatmap")
  input <- clustering_fixture()
  raw <- npx_heatmap(input$data, variable = "Group", filename = tempfile(fileext = ".pdf"))
  expect_gt(max(raw$matrix), 5)
  expect_equal(raw$matrix, input$matrix[rownames(raw$matrix), colnames(raw$matrix)])
  scaled <- npx_heatmap(input$data, variable = "Group", scale_rows = TRUE,
                        filename = tempfile(fileext = ".pdf"))
  expect_equal(unname(rowMeans(scaled$matrix)), rep(0, 8), tolerance = 1e-10)
  expect_equal(unname(apply(scaled$matrix, 1, stats::sd)), rep(1, 8))
  one <- npx_heatmap(input$data, assays = "A2", variable = "Group",
                     filename = tempfile(fileext = ".pdf"))
  expect_equal(nrow(one$matrix), 1L)
  expect_gt(min(one$matrix), 5)
})

test_that("cluster heatmap uses custom grouping annotations and actual NPX", {
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("htmlwidgets")
  input <- clustering_fixture()
  clusters <- list(matrix = input$matrix,
                    clusters = data.frame(Assay = rownames(input$matrix), cluster = rep(1:2, 4)))
  local_mocked_bindings(saveWidget = function(...) invisible(NULL), .package = "htmlwidgets")
  result <- npx_cluster_heatmap(input$data, input$analysis, variable = "Group",
                                 cluster_result = clusters, output_dir = tempfile())
  expect_named(result$annotation_col, "Group")
  expect_equal(result$matrix, input$matrix[rownames(result$matrix), colnames(result$matrix)])
  expect_gt(max(result$matrix), 5)
  expect_true(all(file.exists(result$files[c("pdf", "png")])))
})

test_that("QC HTML exports survive unavailable self-contained conversion", {
  skip_if_not_installed("plotly")
  skip_if_not_installed("htmlwidgets")
  input <- clustering_fixture()
  saved <- list()
  local_mocked_bindings(ggsave = function(...) invisible(NULL), .package = "ggplot2")
  local_mocked_bindings(saveWidget = function(widget, file, selfcontained) {
    if (selfcontained) stop("pandoc unavailable")
    saved[[length(saved) + 1L]] <<- file
    invisible(NULL)
  }, .package = "htmlwidgets")
  expect_warning(iqr <- npx_iqrplot(input$data, variable = "Group", output_dir = tempfile()),
                   "Self-contained IQR HTML")
  expect_warning(distribution <- npx_distribution_boxplot(input$data, variable = "Group", output_dir = tempfile()),
                   "Self-contained distribution HTML")
  expect_equal(saved, list(iqr$files$html, distribution$files$html))
})
