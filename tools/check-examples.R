# Run from the checkout root; the demo data is excluded from built packages.
# Load the current source rather than a potentially stale installed NPXplore.
pkgload::load_all(".", quiet = TRUE)
paths <- c("example/NPXfile.parquet", "example/metadata_HC_vs_Cancer.csv",
           "example/metadata_human_dilute.csv")
stopifnot(all(file.exists(paths)))
data <- OlinkAnalyze::read_npx(paths[1])
for (variable in c("Condition", "Group")) {
  meta_path <- paths[if (variable == "Condition") 2L else 3L]
  inputs <- NPXplore::npx_read_inputs(data, read.csv(meta_path))
  prepared <- NPXplore::npx_apply_qc(inputs)
  result <- NPXplore::npx_analyze(prepared$data, variable = variable)
  stopifnot(nrow(prepared$data) > 0L, nrow(result$result) > 0L,
            all(c("Assay", "Adjusted_pval") %in% names(result$result)))
  message(variable, ": ", nrow(result$result), " statistical results; ",
          nrow(result$deps), " DEPs")
}
