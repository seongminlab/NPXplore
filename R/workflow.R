#' Read NPX and metadata files for step-by-step analysis
#' @param data_file NPX file path accepted by [OlinkAnalyze::read_npx()].
#' @param meta_file Metadata CSV path.
#' @param annotation_file Annotation XLSX path; defaults to the bundled annotation.
#' @param sample_id_col Sample identifier column in metadata.
#' @return A list with `data`, `meta`, and `anno`, accepted by [npx_apply_qc()].
#' @export
npx_read_files <- function(data_file, meta_file, annotation_file = npx_annotation_file(), sample_id_col = "SampleID") {
  if (!file.exists(data_file)) stop("NPX file does not exist: ", data_file)
  if (!file.exists(meta_file)) stop("Metadata file does not exist: ", meta_file)
  npx_read_inputs(OlinkAnalyze::read_npx(data_file),
    utils::read.csv(meta_file, check.names = FALSE), annotation_file, sample_id_col)
}

#' Apply original sample QC and prepare cleaned analysis data
#' @param inputs Result from [npx_read_files()] or [npx_read_inputs()].
#' @param qc Optional result from [npx_qc_reports()] on the original input data.
#' @param include_warn Retain samples with final WARN status.
#' @param include_fail Retain samples with final FAIL status. Cleaning still removes
#'   datapoints according to OlinkAnalyze's assay/data cleaning rules.
#' @param output_dir Optional root directory for metadata CSV exports.
#' @return A list with cleaned `data`, used `meta`, selected `analysis_meta`, original
#'   `qc`, `qc_meta`, metadata-matched `raw_data`, pre-exclusion `qc_data` for QC UMAP, check logs, and `files`.
#' @export
npx_apply_qc <- function(inputs, qc = NULL, include_warn = TRUE, include_fail = FALSE, output_dir = NULL) {
  if (!is.list(inputs) || !is.data.frame(inputs$data) || !is.data.frame(inputs$meta)) {
    stop("inputs must contain data and meta data frames from npx_read_files() or npx_read_inputs().")
  }
  for (value in list(include_warn, include_fail)) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop("QC inclusion flags must be TRUE or FALSE.")
  }
  npx_validate_inputs(inputs$data, inputs$meta)
  raw <- dplyr::filter(inputs$data, SampleID %in% inputs$meta$SampleID, SampleType == "SAMPLE")
  if (is.null(qc)) qc <- npx_qc_reports(raw)
  if (!is.list(qc) || !is.data.frame(qc$sample_qc) || !all(c("SampleID", "QC", "QCIssue", "FailBlock", "WarnBlock", "PassBlock") %in% names(qc$sample_qc))) {
    stop("qc must be a complete npx_qc_reports() result.")
  }
  original_ids <- unique(as.character(raw$SampleID[raw$AssayType == "assay"]))
  if (length(setdiff(original_ids, qc$sample_qc$SampleID))) stop("qc is missing original sample IDs; calculate QC before cleaning.")
  qc_meta <- npx_qc_metadata(qc, inputs$meta)
  metadata <- dplyr::select(inputs$meta, -dplyr::any_of(setdiff(names(qc_meta), "SampleID")))
  analysis_meta <- qc_meta |>
    dplyr::filter(.data$QC == "PASS" | (.data$QC == "WARN" & include_warn) | (.data$QC == "FAIL" & include_fail)) |>
    dplyr::inner_join(metadata, by = "SampleID")
  cleaned <- npx_clean(inputs$data)
  data <- npx_prepare(cleaned$data, analysis_meta, include_warn = TRUE, include_fail = TRUE)
  meta_used <- analysis_meta[match(unique(data$SampleID), analysis_meta$SampleID), , drop = FALSE]
  qc_data <- raw |>
    dplyr::filter(AssayType == "assay") |>
    dplyr::select("SampleID", "Assay", "NPX") |>
    dplyr::left_join(qc_meta, by = "SampleID")
  files <- character()
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files <- c(metadata_used = file.path(output_dir, "metadata_used.csv"),
               analysis_metadata = file.path(output_dir, "QCreport_metadata.csv"))
    utils::write.csv(dplyr::select(meta_used, -dplyr::any_of(c("FailBlock", "WarnBlock", "PassBlock", "QCIssue"))), files[["metadata_used"]], row.names = FALSE)
    utils::write.csv(dplyr::select(analysis_meta, -dplyr::any_of("QCIssue")), files[["analysis_metadata"]], row.names = FALSE)
  }
  list(data = data, meta = meta_used, analysis_meta = analysis_meta, qc = qc,
       qc_meta = qc_meta, raw_data = raw, qc_data = qc_data, check_log = cleaned$check_log,
       check_log_clean = cleaned$check_log_clean, files = files)
}

#' Export NPX, statistical, and DEP tables
#' @param data Prepared analysis data, such as `prepared$data` from [npx_apply_qc()].
#' @param analysis Optional result from [npx_ttest()], [npx_wilcox()], or [npx_anova()].
#'   Omit to export NPX data before statistical analysis.
#' @param save_npx Save the prepared NPX matrix. Set FALSE for later statistical exports.
#' @param output_dir Root output directory; tables are written under `Tables`.
#' @param raw_data Optional original NPX data. SAMPLE/assay rows are retained.
#'   Use `prepared$raw_data` to retain excluded samples from the original metadata.
#' @param posthoc Optional [npx_anova_posthoc()] result to export.
#' @param cluster_result Optional successful [npx_cluster_assays()] result.
#' @return Named paths of saved files. The exported NPX matrix uses assay median
#'   imputation; the supplied data and analysis objects are not modified.
#' @export
npx_export_tables <- function(data, analysis = NULL, output_dir = ".", raw_data = NULL, posthoc = NULL, cluster_result = NULL, save_npx = TRUE) {
  if (!is.null(analysis) && (!is.list(analysis) || !is.data.frame(analysis$result) || !is.data.frame(analysis$deps))) stop("analysis must contain result and deps tables.")
  if (!is.logical(save_npx) || length(save_npx) != 1L || is.na(save_npx)) stop("save_npx must be TRUE or FALSE.")
  if (!is.null(cluster_result) && is.null(analysis)) stop("analysis is required to export cluster assignments.")
  directory <- file.path(output_dir, "Tables")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  files <- character()
  save_table <- function(value, key, filename) {
    path <- file.path(directory, filename)
    utils::write.csv(value, path, row.names = FALSE)
    files[key] <<- path
  }
  if (!is.null(raw_data)) {
    save_table(dplyr::filter(raw_data, SampleType == "SAMPLE", AssayType == "assay"), "raw_data", "0.Raw_NPX_data.csv")
  }
  if (!is.null(analysis)) {
    save_table(npx_statistical_table_for_output(analysis$result, analysis), "analysis", "2.Statistical_analysis_results.csv")
    deps <- npx_statistical_table_for_output(analysis$deps, analysis)
    save_table(deps, "deps", "DEPs.csv")
    if (!is.null(analysis$excluded_assay_reasons)) save_table(analysis$excluded_assay_reasons, "excluded_assays", "Excluded_assays.csv")
  }
  if (save_npx) {
    mat <- npx_matrix(data)
    for (j in seq_len(ncol(mat))) {
      med <- stats::median(mat[, j], na.rm = TRUE)
      if (is.finite(med)) mat[is.na(mat[, j]), j] <- med
    }
    tab <- as.data.frame(t(mat))
    save_table(cbind(Assay = rownames(tab), tab), "npx", "1.NPX_table.csv")
  }
  if (!is.null(posthoc)) save_table(npx_posthoc_table_for_output(posthoc), "anova_posthoc", "3.ANOVA_posthoc_test.csv")
  if (!is.null(cluster_result)) {
    clusters <- npx_cluster_table(analysis$deps, cluster_result)
    save_table(dplyr::left_join(deps, clusters, by = "Assay"), "clusters", "DEPs_with_ClusterNo.csv")
  }
  files
}

#' Save an NPXplore plot without manual graphics-device management
#' @param plot A ggplot or a result from [npx_umap()] or [npx_heatmap()].
#' @param filename Output path ending in .pdf, .png, or .html.
#' @param width,height Optional dimensions in inches. Defaults to result dimensions
#'   or 8 by 6 inches for a bare ggplot. Heatmap export reserves label space.
#' @param dpi PNG resolution; defaults to 300.
#' @return The saved path, invisibly. Heatmaps support PDF/PNG; ggplots also HTML.
#' @export
npx_save_plot <- function(plot, filename, width = NULL, height = NULL, dpi = 300) {
  extension <- tolower(tools::file_ext(filename))
  if (!extension %in% c("pdf", "png", "html")) stop("filename must end in .pdf, .png, or .html.")
  object <- if (inherits(plot, "ggplot") || inherits(plot, "pheatmap")) plot else if (is.list(plot)) plot$plot else NULL
  if (is.null(object)) stop("plot must be a ggplot, npx_umap(), or npx_heatmap() result.")
  heatmap <- inherits(object, "pheatmap")
  if (is.null(width)) width <- if (heatmap) object$width else if (!inherits(plot, "ggplot")) plot$width else NULL
  if (is.null(height)) height <- if (heatmap) object$height else if (!inherits(plot, "ggplot")) plot$height else NULL
  if (is.null(width)) width <- if (heatmap) 12 else 8
  if (is.null(height)) height <- if (heatmap) 10 else 6
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (heatmap) {
    if (extension == "html") stop("Heatmap HTML is not supported by npx_save_plot(); use npx_cluster_heatmap() for cluster HTML.")
    npx_save_heatmap(object, filename, width, height, dpi)
  } else if (inherits(object, "ggplot")) {
    if (extension == "html") {
      if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) stop("plotly and htmlwidgets are required for HTML export.")
      widget <- plotly::ggplotly(object)
      tryCatch(htmlwidgets::saveWidget(widget, filename, selfcontained = TRUE), error = function(e) {
        npx_rethrow_time_limit(e)
        warning("Saving HTML with companion dependencies: ", conditionMessage(e), call. = FALSE)
        htmlwidgets::saveWidget(widget, filename, selfcontained = FALSE)
      })
    } else ggplot2::ggsave(filename, object, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
  } else stop("Unsupported plot object.")
  invisible(filename)
}
