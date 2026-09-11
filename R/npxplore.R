#' @importFrom dplyr .data
#' @importFrom dplyr any_of
NULL

utils::globalVariables(c("SampleID", "SampleType", "AssayType", "SampleQC",
                         "Block", "Assay", "AssayQC", "Normalization", "NPX",
                         "OlinkID", "n_levels"))

npx_validate_dep_cutoffs <- function(logFC = 0, adj.p = 0.05) {
  if (!is.numeric(logFC) || length(logFC) != 1L || is.na(logFC) || logFC < 0) {
    stop("logFC must be one non-negative number.")
  }
  if (!is.numeric(adj.p) || length(adj.p) != 1L || is.na(adj.p) || adj.p <= 0 || adj.p > 1) {
    stop("adj.p must be one number greater than 0 and less than or equal to 1.")
  }
  invisible(TRUE)
}

npx_condition_levels <- function(values, warn = TRUE) {
  values <- unique(trimws(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  keys <- gsub("[ -]+", "_", tolower(values))
  control_idx <- grepl("^(control|healthy|healthy_control|hc)$", keys)
  has_one_control <- sum(control_idx) == 1L
  if (!has_one_control && warn) {
    warning("A control condition (Control, Healthy, Healthy Control, or HC) was not uniquely identified; the current condition order will be used.", call. = FALSE)
  }
  plot_levels <- if (has_one_control) c(values[control_idx], values[!control_idx]) else values
  test_levels <- if (length(values) == 2L && has_one_control) rev(plot_levels) else plot_levels
  list(plot = plot_levels, test = test_levels,
       control = if (has_one_control) values[control_idx] else NULL,
       found = has_one_control)
}

npx_add_effect <- function(tab) {
  if ("Log2FC" %in% names(tab)) return(tab)
  if ("estimate" %in% names(tab)) {
    tab$Log2FC <- as.numeric(tab$estimate)
  } else if (all(c("estimate1", "estimate2") %in% names(tab))) {
    tab$Log2FC <- as.numeric(tab$estimate2) - as.numeric(tab$estimate1)
  }
  tab
}

npx_filter_deps <- function(tab, logFC = 0, adj.p = 0.05, require_effect = FALSE) {
  npx_validate_dep_cutoffs(logFC, adj.p)
  p_col <- intersect(c("Adjusted_pval", "adjusted_pval", "p.adj"), names(tab))[1]
  if (is.na(p_col)) stop("Analysis table has no adjusted p-value column.")
  tab <- npx_add_effect(tab)
  p <- as.numeric(tab[[p_col]])
  keep <- is.finite(p) & p < adj.p
  if ("Log2FC" %in% names(tab)) {
    keep <- keep & is.finite(as.numeric(tab$Log2FC)) & abs(as.numeric(tab$Log2FC)) >= logFC
  } else if (require_effect && logFC > 0) {
    warning("logFC was requested, but this analysis has no directional Log2FC column; only adj.p was applied.")
  }
  tab[keep, , drop = FALSE]
}

npx_complete_columns <- function(tab, columns) {
  for (column in setdiff(columns, names(tab))) tab[[column]] <- rep(NA, nrow(tab))
  tab[, columns, drop = FALSE]
}

npx_statistical_table_for_output <- function(tab, analysis) {
  if (!is.data.frame(tab)) stop("Statistical results must be a data frame.")
  method <- if (is.list(analysis) && !is.null(analysis$method)) as.character(analysis$method)[1] else ""
  if (method %in% c("ttest_independent", "ttest_paired", "wilcox_independent", "wilcox_paired")) {
    tab <- npx_add_effect(tab)
    condition_columns <- if (is.list(analysis) && !is.null(analysis$test_levels)) {
      as.character(analysis$test_levels)
    } else character()
    npx_complete_columns(
      tab,
      c("Assay", "OlinkID", "UniProt", "Panel", "Log2FC",
        condition_columns, "p.value", "Adjusted_pval")
    )
  } else if (identical(method, "anova")) {
    npx_complete_columns(
      tab,
      c("Assay", "OlinkID", "UniProt", "term", "df", "sumsq",
        "meansq", "statistic", "p.value", "Adjusted_pval")
    )
  } else {
    stop("Unsupported statistical analysis method: ", method)
  }
}

npx_posthoc_table_for_output <- function(tab) {
  if (!is.data.frame(tab)) stop("Post-hoc results must be a data frame.")
  npx_complete_columns(
    npx_add_effect(tab),
    c("Assay", "OlinkID", "UniProt", "term", "contrast", "Log2FC", "Adjusted_pval")
  )
}

#' Prepare already-read NPXplore inputs
#' @param data An NPX data frame returned by [OlinkAnalyze::read_npx()].
#' @param meta_file Metadata data frame returned by a user-side `read.csv()` call.
#' @param annotation_file Optional annotation XLSX file.
#' @param sample_id_col Metadata column containing sample identifiers. If the
#'   column is not named `SampleID`, provide its name here.
#' @return A list containing `data`, `meta`, and `anno`.
#' @export
npx_read_inputs <- function(data, meta_file, annotation_file = npx_annotation_file(), sample_id_col = "SampleID") {
  if (!is.data.frame(data)) {
    stop("data must be an NPX data frame returned by OlinkAnalyze::read_npx().")
  }
  if (!is.data.frame(meta_file)) {
    stop("meta_file must be a metadata data frame. Read the CSV first with read.csv().")
  }
  meta <- meta_file
  if (!file.exists(annotation_file)) stop("annotation_file does not exist: ", annotation_file)
  anno <- readxl::read_xlsx(annotation_file)
  npx_validate_inputs(data, meta, sample_id_col)
  if (sample_id_col != "SampleID") {
    if ("SampleID" %in% names(meta)) stop("Metadata contains both 'SampleID' and the specified sample_id_col.")
    names(meta)[names(meta) == sample_id_col] <- "SampleID"
  }
  list(data = data, meta = meta, anno = anno)
}

#' Locate NPXplore's built-in annotation
#' @return A file path.
#' @export
npx_annotation_file <- function() {
  path <- system.file("extdata", "HT_panel_pathway_info.xlsx", package = "NPXplore")
  if (!nzchar(path)) stop("Built-in annotation is not installed.")
  path
}

#' Validate NPX data and metadata
#' @param data An Olink NPX data frame.
#' @param meta A metadata data frame.
#' @param sample_id_col Metadata sample identifier column.
#' @export
npx_validate_inputs <- function(data, meta, sample_id_col = "SampleID") {
  if (!is.data.frame(data) || !is.data.frame(meta)) stop("data and meta must be data frames.")
  missing_data <- setdiff(c("SampleID", "SampleType", "AssayType", "Assay", "NPX"), names(data))
  if (length(missing_data)) stop("data is missing required columns: ", paste(missing_data, collapse = ", "))
  if (!sample_id_col %in% names(meta)) {
    if (identical(sample_id_col, "SampleID")) {
      stop("Metadata does not contain a 'SampleID' column. Please specify the column containing sample identifiers with sample_id_col, for example sample_id_col = 'PatientID'.")
    }
    stop("Metadata column not found: ", sample_id_col, ". Please provide the correct sample_id_col name.")
  }
  if (anyDuplicated(meta[[sample_id_col]])) stop("Metadata must contain one row per ", sample_id_col, ".")
  if (anyNA(meta[[sample_id_col]]) || any(!nzchar(trimws(as.character(meta[[sample_id_col]]))))) stop("Metadata sample identifiers cannot be missing or blank.")
  invisible(TRUE)
}

#' Summarise block-level and sample-level Olink QC status
#'
#' A block is FAIL only when every evaluated assay in that block is FAIL.
#' Assays with `SampleQC == "NA"` are excluded from the denominator. A sample
#' is FAIL when at least two blocks are FAIL; otherwise any WARN/FAIL record
#' makes the sample WARN. This preserves the original Olink issue flag while
#' adding a conservative whole-sample exclusion rule.
#' @param data An Olink NPX data frame.
#' @return A list with `block` and `sample` QC summaries.
#' @keywords internal
npx_sample_qc_status <- function(data) {
  required <- c("SampleID", "SampleType", "AssayType", "Assay", "Block", "SampleQC")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("data is missing required QC columns: ", paste(missing, collapse = ", "))
  rows <- data |>
    dplyr::filter(SampleType == "SAMPLE", AssayType == "assay") |>
    dplyr::mutate(
      .SampleQC = toupper(trimws(as.character(.data$SampleQC))),
      .SampleQC = dplyr::coalesce(dplyr::na_if(.data$.SampleQC, ""), "NA"),
      .SampleQC = dplyr::if_else(.data$.SampleQC == "WARNING", "WARN", .data$.SampleQC)
    )
  row_cols <- c(intersect(c("SampleID", "Panel", "Block", "Assay"), names(rows)), ".SampleQC")
  rows <- dplyr::distinct(rows[, row_cols, drop = FALSE])
  group_cols <- intersect(c("SampleID", "Panel", "Block"), names(rows))
  block <- rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      TotalAssay = dplyr::n_distinct(.data$Assay),
      NAAssay = sum(.data$.SampleQC == "NA"),
      EvaluatedAssay = .data$TotalAssay - .data$NAAssay,
      PassAssay = sum(.data$.SampleQC == "PASS"),
      WarnAssay = sum(.data$.SampleQC == "WARN"),
      FailAssay = sum(.data$.SampleQC == "FAIL"),
      BlockQC = dplyr::case_when(
        .data$EvaluatedAssay > 0 & .data$FailAssay == .data$EvaluatedAssay ~ "FAIL",
        .data$WarnAssay > 0 | .data$FailAssay > 0 ~ "WARN",
        TRUE ~ "PASS"
      ),
      .groups = "drop"
    )
  sample <- block |>
    dplyr::group_by(SampleID) |>
    dplyr::summarise(
      FailBlock = sum(.data$BlockQC == "FAIL"),
      WarnBlock = sum(.data$BlockQC == "WARN"),
      PassBlock = sum(.data$BlockQC == "PASS"),
      QCIssue = any(.data$BlockQC %in% c("WARN", "FAIL")),
      QC = dplyr::case_when(
        .data$FailBlock >= 2 ~ "FAIL",
        .data$QCIssue ~ "WARN",
        TRUE ~ "PASS"
      ),
      .groups = "drop"
    )
  list(block = block, sample = sample)
}

#' Run OlinkAnalyze quality checks and cleaning
#' @param data An Olink NPX data frame.
#' @return A list with `data`, `check_log`, `check_log_clean`, and
#'   `failed_samples`. `failed_samples` is reported for review; whole samples
#'   are not removed automatically because the reviewed metadata determines
#'   which samples are used for analysis.
#' @export
npx_clean <- function(data) {
  if (is.data.frame(data)) data <- tibble::as_tibble(data)
  check_log <- OlinkAnalyze::check_npx(data)
  data_clean <- OlinkAnalyze::clean_npx(data, check_log = check_log)
  qc_status <- npx_sample_qc_status(data)
  failed_samples <- qc_status$sample$SampleID[qc_status$sample$QC == "FAIL"]
  list(data = data_clean, check_log = check_log,
       check_log_clean = OlinkAnalyze::check_npx(data_clean),
       failed_samples = failed_samples,
       block_qc = qc_status$block,
       sample_qc = qc_status$sample)
}

#' Prepare sample assay data by joining metadata
#' @param data Cleaned Olink NPX data.
#' @param meta Metadata.
#' @param sample_id_col Metadata sample identifier column.
#' @param include_warn If `TRUE`, retain samples with Olink `WARN` status.
#'   Defaults to `TRUE`. If `FALSE`, WARN samples are removed. FAIL datapoints
#'   are handled by [npx_clean()], while the reviewed metadata determines which
#'   samples are included in the analysis.
#' @param include_fail If `TRUE`, retain samples with final sample-level
#'   `FAIL` status. Defaults to `FALSE`.
#' @return A tibble containing assay rows and metadata.
#' @export
npx_prepare <- function(data, meta, sample_id_col = "SampleID", include_warn = TRUE, include_fail = FALSE) {
  for (arg in list(include_warn = include_warn, include_fail = include_fail)) {
    if (!is.logical(arg) || length(arg) != 1L || is.na(arg)) {
      stop("include_warn and include_fail must be one non-missing logical value: TRUE or FALSE.")
    }
  }
  npx_validate_inputs(data, meta, sample_id_col)
  if (sample_id_col != "SampleID") names(meta)[names(meta) == sample_id_col] <- "SampleID"
  if ("QC" %in% names(meta)) {
    allowed <- c("PASS")
    if (include_warn) allowed <- c(allowed, "WARN")
    if (include_fail) allowed <- c(allowed, "FAIL")
    meta <- dplyr::filter(meta, toupper(trimws(as.character(.data$QC))) %in% allowed)
  }
  out <- dplyr::filter(data, SampleID %in% meta$SampleID, SampleType == "SAMPLE", AssayType == "assay")
  out <- dplyr::inner_join(out, tibble::as_tibble(meta), by = "SampleID")
  if (!nrow(out)) stop("No assay rows remain after metadata/QC filtering.")
  out
}

#' Create assay and sample QC summary tables
#' @param data NPX data. The function applies the Rmd sample and assay filters
#'   internally before creating the QC tables.
#' @param meta Optional metadata; QC is restricted to its SampleID values.
#' @param output_dir Optional output root for QC CSV files.
#' @return A list of QC report data frames, with saved file paths when requested.
#' @export
npx_qc_reports <- function(data, meta = NULL, output_dir = NULL) {
  if (!is.null(meta)) {
    if (!is.data.frame(meta) || !"SampleID" %in% names(meta)) stop("meta must contain SampleID.")
    data <- dplyr::filter(data, SampleID %in% meta$SampleID)
  }
  required <- c("SampleID", "SampleType", "AssayType", "Assay", "Block", "AssayQC", "SampleQC", "Normalization")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("data is missing required QC columns: ", paste(missing, collapse = ", "))
  ht <- dplyr::filter(data, SampleType == "SAMPLE")
  ht_assay <- dplyr::filter(ht, AssayType == "assay")
  blocks <- sort(unique(ht_assay$Block))
  blocks <- blocks[!is.na(blocks)]
  assay_qc <- dplyr::summarise(
    dplyr::group_by(ht_assay, Block),
    NumberOfAssay = dplyr::n_distinct(Assay),
    WarnedAssay = dplyr::n_distinct(Assay[toupper(trimws(as.character(AssayQC))) == "WARN"]),
    ExcludedAssay = dplyr::n_distinct(Assay[toupper(trimws(as.character(AssayQC))) == "NA"]),
    .groups = "drop"
  )
  if (length(blocks)) {
    assay_qc <- assay_qc[match(blocks, assay_qc$Block), , drop = FALSE]
  }
  excluded <- dplyr::distinct(dplyr::filter(ht, Normalization == "EXCLUDED"), dplyr::across(dplyr::any_of(c("AssayQC", "Panel", "Block", "UniProt", "Assay", "OlinkID"))))
  warned <- dplyr::distinct(dplyr::filter(ht, toupper(trimws(as.character(AssayQC))) == "WARN"), dplyr::across(dplyr::any_of(c("AssayQC", "Panel", "Block", "UniProt", "Assay", "OlinkID"))))
  qc_status <- npx_sample_qc_status(ht_assay)
  issue_ids <- qc_status$sample$SampleID[qc_status$sample$QCIssue]
  issue_ids <- issue_ids[!is.na(issue_ids) & nzchar(issue_ids)]
  # Find issue samples first, then count every assay status for those samples.
  # This matches the Rmd: Issue is used only to select SampleID, while the
  # subsequent PASS/WARN/FAIL/NA counts use the complete ht.assay table.
  issues <- dplyr::filter(ht_assay, SampleID %in% issue_ids)
  sample_qc <- dplyr::summarise(
    dplyr::group_by(issues, SampleID),
    PassAssay = dplyr::n_distinct(Assay[toupper(trimws(as.character(SampleQC))) == "PASS"]),
    WarnedAssay = dplyr::n_distinct(Assay[toupper(trimws(as.character(SampleQC))) == "WARN"]),
    FailedAssay = dplyr::n_distinct(Assay[toupper(trimws(as.character(SampleQC))) == "FAIL"]),
    ExcludedAssay = dplyr::n_distinct(Assay[toupper(trimws(as.character(SampleQC))) == "NA"]),
    .groups = "drop"
  )
  result <- list(assay = assay_qc, excluded_assay = excluded, warned_assay = warned,
       warning_sample_statistics = sample_qc,
       block_qc = qc_status$block,
       sample_qc = qc_status$sample)
  if (!is.null(output_dir)) result$files <- npx_write_qc(result, output_dir, meta)
  result
}

#' Prepare grouping levels and factor order for differential analysis
#' @keywords internal
npx_prepare_analysis <- function(data, variable = "Condition", min_levels = 2L) {
  data <- tibble::as_tibble(data)
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  condition_order <- npx_condition_levels(data[[variable]])
  if (length(condition_order$plot) < min_levels) {
    stop("At least ", min_levels, " condition levels are required.")
  }
  data[[variable]] <- factor(
    trimws(as.character(data[[variable]])),
    levels = condition_order$test
  )
  list(data = data, order = condition_order, levels = condition_order$plot)
}

#' Keep assays with usable observations in both two-condition groups
#' @keywords internal
npx_prepare_two_group_data <- function(data, variable) {
  excluded_assays <- character()
  if ("OlinkID" %in% names(data)) {
    valid <- data |>
      dplyr::filter(!is.na(.data$NPX), !is.na(.data[[variable]])) |>
      dplyr::group_by(.data$OlinkID, .data[[variable]]) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::count(.data$OlinkID, name = "n_levels") |>
      dplyr::filter(.data$n_levels == 2L) |>
      dplyr::pull(.data$OlinkID)
    excluded_assays <- setdiff(unique(as.character(data$OlinkID)), as.character(valid))
    data <- dplyr::filter(data, .data$OlinkID %in% valid)
  }
  if (!nrow(data)) stop("No assays have usable NPX observations in every condition.")
  list(data = data, excluded_assays = excluded_assays)
}

# Validate individual t-tests before calling the vectorized Olink implementation.
npx_ttest_usable <- function(data, variable, pair_id = NULL) {
  if (!is.null(pair_id) && !pair_id %in% names(data)) stop("Pairing column not found: ", pair_id)
  keys <- intersect(c("Assay", "OlinkID", "UniProt", "Panel"), names(data))
  groups <- dplyr::group_split(dplyr::group_by(data, dplyr::across(dplyr::all_of(keys))))
  levels <- levels(data[[variable]])
  reasons <- lapply(groups, function(g) {
    x <- g$NPX[g[[variable]] == levels[1] & !is.na(g[[variable]])]
    y <- g$NPX[g[[variable]] == levels[2] & !is.na(g[[variable]])]
    if (!is.null(pair_id)) {
      g <- g[!is.na(g[[pair_id]]), , drop = FALSE]
      ids <- as.character(g[[pair_id]])
      a <- g[[variable]] == levels[1]; b <- g[[variable]] == levels[2]
      a[is.na(a)] <- FALSE; b[is.na(b)] <- FALSE
      if (anyDuplicated(ids[a]) || anyDuplicated(ids[b])) {
        stop("Each pair must have only one observation per condition and assay.")
      }
      common <- intersect(ids[a], ids[b])
      x <- g$NPX[a][match(common, ids[a])]; y <- g$NPX[b][match(common, ids[b])]
    }
    reason <- tryCatch({
      stats::t.test(x, y, paired = !is.null(pair_id)); NULL
    }, error = function(e) conditionMessage(e))
    if (is.null(reason)) return(NULL)
    data.frame(OlinkID = as.character(g$OlinkID[1]), Reason = reason)
  })
  excluded <- dplyr::bind_rows(reasons)
  if (!nrow(excluded)) excluded <- data.frame(OlinkID = character(), Reason = character())
  kept <- dplyr::filter(data, !.data$OlinkID %in% excluded$OlinkID)
  if (!nrow(kept)) stop("No assays support a t-test: ", paste(unique(excluded$Reason), collapse = "; "))
  if (nrow(excluded)) message("Excluded ", nrow(excluded), " assay(s) that cannot support a t-test; see analysis$excluded_assay_reasons.")
  list(data = kept, excluded = excluded)
}

#' Assemble a consistent NPXplore differential-analysis result
#' @keywords internal
npx_finalize_analysis <- function(result, method, variable, condition_order,
                                  excluded_assays = character(), logFC = 0,
                                  adj.p = 0.05, require_effect = FALSE) {
  result <- npx_add_effect(result)
  deps <- npx_filter_deps(
    result,
    logFC = logFC,
    adj.p = adj.p,
    require_effect = require_effect
  )
  list(
    method = method,
    variable = variable,
    levels = condition_order$plot,
    test_levels = condition_order$test,
    control = condition_order$control,
    excluded_assays = excluded_assays,
    logFC = logFC,
    adj.p = adj.p,
    result = result,
    deps = deps
  )
}

#' Calculate two-condition differential proteins with an Olink t-test
#' @param data Prepared assay data.
#' @param variable Metadata column used as the grouping variable.
#' @param pair_id Optional pairing identifier for paired tests.
#' @param check_log Optional Olink QC log.
#' @param logFC Minimum absolute Log2 fold-change for DEPs.
#' @param adj.p Maximum adjusted p-value for DEPs.
#' @return A list containing the complete result in `result` and filtered DEPs in `deps`.
#' @export
npx_ttest <- function(data, variable = "Condition", pair_id = NULL,
                      check_log = NULL, logFC = 0, adj.p = 0.05) {
  npx_validate_dep_cutoffs(logFC, adj.p)
  prepared <- npx_prepare_analysis(data, variable, min_levels = 2L)
  if (length(prepared$levels) != 2L) {
    stop("npx_ttest() requires exactly two condition levels.")
  }
  two_group <- npx_prepare_two_group_data(prepared$data, variable)
  usable <- npx_ttest_usable(two_group$data, variable, pair_id)
  args <- list(df = usable$data, variable = variable, check_log = check_log)
  if (!is.null(pair_id)) args$pair_id <- pair_id
  result <- do.call(OlinkAnalyze::olink_ttest, args)
  method <- if (is.null(pair_id)) "ttest_independent" else "ttest_paired"
  final <- npx_finalize_analysis(result, method, variable, prepared$order,
                        unique(c(two_group$excluded_assays, usable$excluded$OlinkID)), logFC, adj.p,
                        require_effect = TRUE)
  missing_group <- data.frame(OlinkID = two_group$excluded_assays,
                              Reason = rep("No usable observations in one or both conditions", length(two_group$excluded_assays)))
  final$excluded_assay_reasons <- dplyr::bind_rows(missing_group, usable$excluded)
  final
}

#' Calculate two-condition differential proteins with an Olink Wilcoxon test
#' @param data Prepared assay data.
#' @param variable Metadata column used as the grouping variable.
#' @param pair_id Optional pairing identifier for paired tests.
#' @param check_log Optional Olink QC log.
#' @param logFC Minimum absolute Log2 fold-change for DEPs.
#' @param adj.p Maximum adjusted p-value for DEPs.
#' @return A list containing the complete result in `result` and filtered DEPs in `deps`.
#' @export
npx_wilcox <- function(data, variable = "Condition", pair_id = NULL,
                       check_log = NULL, logFC = 0, adj.p = 0.05) {
  npx_validate_dep_cutoffs(logFC, adj.p)
  prepared <- npx_prepare_analysis(data, variable, min_levels = 2L)
  if (length(prepared$levels) != 2L) {
    stop("npx_wilcox() requires exactly two condition levels.")
  }
  two_group <- npx_prepare_two_group_data(prepared$data, variable)
  args <- list(df = two_group$data, variable = variable, check_log = check_log)
  if (!is.null(pair_id)) args$pair_id <- pair_id
  result <- do.call(OlinkAnalyze::olink_wilcox, args)
  method <- if (is.null(pair_id)) "wilcox_independent" else "wilcox_paired"
  npx_finalize_analysis(result, method, variable, prepared$order,
                        two_group$excluded_assays, logFC, adj.p,
                        require_effect = TRUE)
}

#' Backward-compatible alias for [npx_wilcox()].
#' @param ... Arguments passed to [npx_wilcox()].
#' @export
npx_wilcon <- function(...) npx_wilcox(...)

#' Calculate multi-condition differential proteins with one-way ANOVA
#' @param data Prepared assay data.
#' @param variable Metadata column used as the grouping variable.
#' @param covariates Optional covariate column names for ANOVA.
#' @param check_log Optional Olink QC log.
#' @param logFC Minimum absolute Log2 fold-change for DEPs.
#' @param adj.p Maximum adjusted p-value for DEPs.
#' @return A list containing the complete result in `result` and filtered DEPs in `deps`.
#' @export
npx_anova <- function(data, variable = "Condition", covariates = NULL,
                      check_log = NULL, logFC = 0, adj.p = 0.05) {
  npx_validate_dep_cutoffs(logFC, adj.p)
  prepared <- npx_prepare_analysis(data, variable, min_levels = 3L)
  if (!is.null(covariates)) {
    missing_covariates <- setdiff(covariates, names(prepared$data))
    if (length(missing_covariates)) {
      stop("Covariate column(s) not found: ", paste(missing_covariates, collapse = ", "))
    }
  }
  result <- OlinkAnalyze::olink_anova(
    df = prepared$data,
    variable = variable,
    covariates = covariates,
    check_log = check_log
  )
  npx_finalize_analysis(result, "anova", variable, prepared$order,
                        logFC = logFC, adj.p = adj.p)
}

#' Run the appropriate differential analysis
#' @param data Prepared assay data.
#' @param variable Metadata column used as the grouping variable.
#' @param test Optional method: `ttest`, `wilcox`, or `anova`. If `NULL`,
#'   two groups use `npx_ttest()` and three or more groups use `npx_anova()`.
#' @param pair_id Optional pairing identifier for paired two-group tests.
#' @param covariates Optional covariate column names for ANOVA.
#' @param check_log Optional Olink QC log.
#' @param logFC Minimum absolute Log2 fold-change for DEPs.
#' @param adj.p Maximum adjusted p-value for DEPs.
#' @return A result list from `npx_ttest()`, `npx_wilcox()`, or `npx_anova()`.
#' @export
npx_analyze <- function(data, variable = "Condition", test = NULL, pair_id = NULL,
                        covariates = NULL, check_log = NULL, logFC = 0,
                        adj.p = 0.05) {
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  condition_order <- npx_condition_levels(data[[variable]])
  values <- condition_order$plot
  if (length(values) < 2L) stop("At least two condition levels are required.")
  if (!is.null(test)) test <- match.arg(test, c("ttest", "wilcox", "anova"))
  selected_test <- if (is.null(test)) if (length(values) == 2L) "ttest" else "anova" else test
  if (length(values) > 2L && selected_test %in% c("ttest", "wilcox")) {
    stop("ttest and wilcox require exactly two condition levels; use test = 'anova'.")
  }
  if (length(values) == 2L && selected_test == "anova") {
    stop("anova requires three or more condition levels.")
  }
  if (selected_test == "anova" && !is.null(pair_id)) stop("pair_id is not supported for ANOVA; repeated-measures analysis requires a different model.")
  if (selected_test != "anova" && length(covariates)) stop("covariates are supported only for ANOVA in NPXplore.")
  if (selected_test == "ttest") {
    npx_ttest(data, variable, pair_id, check_log, logFC, adj.p)
  } else if (selected_test == "wilcox") {
    npx_wilcox(data, variable, pair_id, check_log, logFC, adj.p)
  } else {
    npx_anova(data, variable, covariates, check_log, logFC, adj.p)
  }
}

#' Calculate ANOVA post-hoc contrasts
#' @param data Prepared assay data.
#' @param analysis An object returned by [npx_analyze()].
#' @param variable Metadata grouping column.
#' @param covariates Optional ANOVA covariates.
#' @param check_log Optional Olink QC log.
#' @param logFC Minimum absolute post-hoc Log2 fold-change.
#' @param adj.p Maximum adjusted post-hoc p-value.
#' @return A filtered post-hoc result table.
#' @export
npx_anova_posthoc <- function(data, analysis = NULL, variable = "Condition", covariates = NULL,
                             check_log = NULL, logFC = 0, adj.p = 0.05) {
  npx_validate_dep_cutoffs(logFC, adj.p)
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  ids <- NULL
  if (!is.null(analysis)) {
    deps <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$deps else analysis
    if (is.data.frame(deps) && "OlinkID" %in% names(deps)) ids <- unique(as.character(deps$OlinkID))
  }
  args <- list(df = data, olinkid_list = ids, variable = variable,
               covariates = covariates, effect = variable,
               check_log = check_log)
  if (is.null(ids)) args$olinkid_list <- NULL
  out <- do.call(OlinkAnalyze::olink_anova_posthoc, args)
  out <- npx_add_effect(out)
  if ("Threshold" %in% names(out)) out$Threshold <- NULL
  if (!nrow(out)) return(out)
  npx_filter_deps(out, logFC = logFC, adj.p = adj.p)
}

#' Cluster differential assays for multi-condition analysis
#' @param data Prepared assay data.
#' @param analysis An object returned by [npx_analyze()].
#' @param nbclust_args Arguments passed to [NbClust::NbClust()]. Cluster counts
#'   are capped below both the assay count and number of distinct profiles to
#'   allow the package's additional `k + 1` calculation. Multiple indices use
#'   majority vote; ties select the smaller cluster count.
#' @param seed Random seed used by clustering.
#' @param scale_rows If `TRUE`, standardize each assay across samples before
#'   clustering. Defaults to `FALSE`, using NPX values directly. This affects
#'   clustering only; the returned `matrix` retains the imputed NPX values.
#' @return A list containing the NPX `matrix`, `clustering_matrix`, optional
#'   `scaled_matrix`, cluster assignments, and NbClust result. Fewer than four
#'   variable assays or three distinct profiles raises an
#'   `npx_insufficient_input` error rather than assigning arbitrary clusters.
#' @export
npx_cluster_assays <- function(data, analysis, nbclust_args = list(index = "silhouette"), seed = 123,
                                scale_rows = FALSE) {
  if (!requireNamespace("NbClust", quietly = TRUE)) stop("Package 'NbClust' is required for npx_cluster_assays().")
  if (!is.logical(scale_rows) || length(scale_rows) != 1L || is.na(scale_rows)) stop("scale_rows must be TRUE or FALSE.")
  if (!is.list(nbclust_args)) stop("nbclust_args must be a list.")
  insufficient <- function(message) {
    stop(structure(list(message = message, call = NULL),
                   class = c("npx_insufficient_input", "error", "condition")))
  }
  if (is.list(analysis) && !is.null(analysis$method) && !identical(analysis$method, "anova")) {
    stop("Assay clustering is available only for ANOVA results.")
  }
  deps <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$deps else analysis
  if (!is.data.frame(deps) || !"Assay" %in% names(deps)) stop("analysis must contain a DEP table with an 'Assay' column.")
  assays <- unique(as.character(deps$Assay)); assays <- assays[!is.na(assays) & nzchar(assays)]
  if (length(assays) < 4L) insufficient("At least four DEP assays are required for cluster-count optimization; no clusters were assigned.")
  available <- unique(as.character(data$Assay[is.finite(data$NPX)]))
  if (length(intersect(assays, available)) < 4L) insufficient("At least four DEP assays with finite NPX values are required for clustering; no clusters were assigned.")
  mat <- t(npx_matrix(data, assays = assays))
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  for (i in seq_len(nrow(mat))) {
    med <- stats::median(mat[i, ], na.rm = TRUE)
    mat[i, is.na(mat[i, ])] <- med
  }
  mat <- mat[apply(mat, 1, function(x) all(is.finite(x)) && length(x) > 1L && stats::sd(x) > 0), , drop = FALSE]
  if (nrow(mat) < 4L) insufficient("At least four variable DEP assays are required for clustering; no clusters were assigned.")
  clustering_matrix <- if (scale_rows) t(scale(t(mat))) else mat
  distinct_profiles <- nrow(unique(as.data.frame(clustering_matrix)))
  if (distinct_profiles < 3L) insufficient("At least three distinct DEP assay profiles are required for cluster-count optimization; no clusters were assigned.")
  max_count <- min(nrow(mat) - 2L, distinct_profiles - 1L)
  defaults <- list(data = clustering_matrix, distance = "euclidean", min.nc = 2,
                   max.nc = min(10L, max_count), method = "kmeans",
                   index = "silhouette")
  args <- utils::modifyList(defaults, nbclust_args)
  args$data <- clustering_matrix
  valid_count <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x == as.integer(x) && x >= 2L
  if (!valid_count(args$min.nc) || !valid_count(args$max.nc)) stop("NbClust min.nc and max.nc must be finite integers of at least 2.")
  if (args$max.nc > max_count) {
    warning("NbClust max.nc reduced to ", max_count, " for the available assay profiles.", call. = FALSE)
    args$max.nc <- max_count
  }
  if (args$min.nc > args$max.nc) stop("NbClust cluster-count range is empty: min.nc exceeds the supported max.nc (", args$max.nc, ").")
  set.seed(seed)
  nb <- do.call(NbClust::NbClust, args)
  best <- nb$Best.nc
  candidates <- if (is.null(dim(best))) {
    if ("Number_clusters" %in% names(best)) best[["Number_clusters"]] else best[1]
  } else if ("Number_clusters" %in% rownames(best)) {
    best["Number_clusters", ]
  } else stop("NbClust did not return a Number_clusters row.")
  candidates <- as.numeric(candidates)
  candidates <- candidates[is.finite(candidates) & candidates == floor(candidates) &
                             candidates >= args$min.nc & candidates <= args$max.nc]
  if (!length(candidates)) stop("NbClust did not return a valid cluster count.")
  votes <- table(candidates)
  optimal_k <- min(as.integer(names(votes)[votes == max(votes)]))
  set.seed(seed)
  km <- stats::kmeans(clustering_matrix, centers = optimal_k, nstart = 50)
  clusters <- data.frame(Assay = rownames(mat), cluster = as.integer(km$cluster), stringsAsFactors = FALSE)
  clusters$cluster <- factor(clusters$cluster, levels = sort(unique(clusters$cluster)))
  list(matrix = mat, clustering_matrix = clustering_matrix,
       scaled_matrix = if (scale_rows) clustering_matrix else NULL, scale_rows = scale_rows,
       clusters = clusters, optimal_k = optimal_k, nbclust = nb, kmeans = km)
}

npx_safe_filename <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))

npx_qc_metadata <- function(qc, meta) {
  ids <- unique(as.character(meta$SampleID))
  status <- qc$sample_qc
  out <- data.frame(SampleID = ids, stringsAsFactors = FALSE) |>
    dplyr::left_join(status, by = "SampleID")
  out$QC[is.na(out$QC)] <- "PASS"
  out$QCIssue[is.na(out$QCIssue)] <- FALSE
  out$FailBlock[is.na(out$FailBlock)] <- 0L
  out$WarnBlock[is.na(out$WarnBlock)] <- 0L
  out$PassBlock[is.na(out$PassBlock)] <- 0L
  out
}

npx_write_qc <- function(qc, output_dir, meta = NULL) {
  dir <- file.path(output_dir, "QC"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(assay = "1.Report_AssayQC.csv", excluded_assay = "2.Report_Excluded_Assay.csv",
             warned_assay = "3.Report_Warned_Assay.csv", block_qc = "4.Report_BlockQC.csv",
             warning_sample_statistics = "5.Report_Warning_Sample_Statistics.csv")
  out <- list()
  for (nm in names(files)) { path <- file.path(dir, files[[nm]]); utils::write.csv(qc[[nm]], path, row.names = FALSE); out[[nm]] <- path }
  sample_file <- file.path(dir, "4.Report_SampleQC.csv")
  meta_qc <- if (is.null(meta)) qc$sample_qc else npx_qc_metadata(qc, meta)
  sample_qc_output <- dplyr::select(meta_qc, -dplyr::any_of("QCIssue"))
  utils::write.csv(sample_qc_output, sample_file, row.names = FALSE)
  out$sample_qc <- sample_file
  out
}

npx_try <- function(expr, label) {
  tryCatch(expr, error = function(e) { warning(label, " skipped: ", conditionMessage(e), call. = FALSE); NULL })
}

npx_try_timeout <- function(expr, label, seconds = Inf) {
  if (is.infinite(seconds)) return(npx_try(expr, label))
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
  setTimeLimit(elapsed = seconds, transient = TRUE)
  npx_try(expr, label)
}

npx_rethrow_time_limit <- function(e) {
  limits <- c(gettext("reached elapsed time limit", domain = "R"),
              gettext("reached CPU time limit", domain = "R"))
  if (any(endsWith(conditionMessage(e), limits))) stop(e)
}

npx_cluster_table <- function(deps, cluster_result) {
  table <- cluster_result$clusters
  if (!is.data.frame(table) || !all(c("Assay", "cluster") %in% names(table))) {
    stop("cluster_result must contain a clusters table with Assay and cluster columns.")
  }
  table <- data.frame(Assay = trimws(as.character(table$Assay)),
                      ClusterNo = trimws(as.character(table$cluster)), stringsAsFactors = FALSE)
  table <- unique(table[!is.na(table$Assay) & nzchar(table$Assay) &
                          !is.na(table$ClusterNo) & nzchar(table$ClusterNo), , drop = FALSE])
  table <- table[table$Assay %in% as.character(deps$Assay), , drop = FALSE]
  if (anyDuplicated(table$Assay)) stop("An Assay has conflicting cluster assignments.")
  if (!nrow(table)) stop("No DEP assays match a valid cluster assignment.")
  ids <- unique(table$ClusterNo)
  ids <- if (all(grepl("^[0-9]+$", ids))) as.character(sort(as.integer(ids))) else sort(ids)
  table <- table[order(match(table$ClusterNo, ids), table$Assay), , drop = FALSE]
  attr(table, "unassigned_assays") <- setdiff(unique(as.character(deps$Assay)), table$Assay)
  table
}

npx_stage_message <- function(stage, value) {
  status <- if (is.null(value$status)) "NotRun" else value$status
  label <- switch(status, Success = "complete", Partial = "finished with partial results",
                  Failed = "failed", Error = "failed", Skipped = "skipped", NotRun = "not run", status)
  reason <- if (is.null(value$reason) || !nzchar(value$reason)) "" else paste0(" Reason: ", value$reason)
  message("[NPXplore] ", stage, " ", label, ".", reason)
  invisible(value)
}

npx_run_stage <- function(expr, stage, seconds = Inf) {
  message("[NPXplore] ", stage, " started.")
  if (is.finite(seconds)) {
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = TRUE), add = TRUE)
    setTimeLimit(elapsed = seconds, transient = TRUE)
  }
  result <- tryCatch({
    value <- expr
    if (is.null(value)) stop("Stage returned no result.")
    if (is.null(value$status)) value$status <- "Success"
    value
  }, error = function(e) {
    insufficient <- inherits(e, "npx_insufficient_input")
    timed_out <- any(endsWith(conditionMessage(e), c(
      gettext("reached elapsed time limit", domain = "R"),
      gettext("reached CPU time limit", domain = "R"))))
    code <- if (insufficient) "InsufficientInput" else if (timed_out) "TimeLimit" else "RuntimeError"
    status <- if (insufficient) "Skipped" else "Failed"
    list(status = status, code = code, reason = conditionMessage(e),
         files = character(), results = list())
  })
  npx_stage_message(stage, result)
  result
}

npx_blocked_stage <- function(reason, code = "ClusteringUnavailable") {
  list(status = "Skipped", code = code, reason = reason, files = character(), results = list())
}

npx_stage_summary <- function(stages) {
  dplyr::bind_rows(lapply(names(stages), function(name) {
    value <- stages[[name]]
    data.frame(Stage = name, Status = if (is.null(value$status)) "NotRun" else value$status,
               Code = if (is.null(value$code)) "" else value$code,
               Reason = if (is.null(value$reason)) "" else value$reason, stringsAsFactors = FALSE)
  }))
}

npx_boxplot_selection <- function(analysis, top_n = NULL) {
  if (is.null(top_n)) return(NULL)
  deps <- analysis$deps
  if (!nrow(deps)) return(character())
  deps <- deps[order(deps$Adjusted_pval, as.character(deps$Assay), na.last = TRUE), , drop = FALSE]
  utils::head(unique(as.character(deps$Assay)), top_n)
}

npx_save_pipeline_html <- function(plot, file) {
  if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) return(NULL)
  npx_try({
    widget <- plotly::ggplotly(plot, tooltip = "text")
    tryCatch(htmlwidgets::saveWidget(widget, file = file, selfcontained = TRUE), error = function(e) {
      warning("Self-contained HTML unavailable; saving companion dependencies: ", conditionMessage(e), call. = FALSE)
      htmlwidgets::saveWidget(widget, file = file, selfcontained = FALSE)
    })
    file
  }, "Interactive HTML export")
}

#' Run the complete NPXplore pipeline
#' @param data_file NPX data frame already read with [OlinkAnalyze::read_npx()].
#'   NPXplore does not read NPX files internally.
#' @param meta_file Metadata data frame returned by a user-side `read.csv()` call.
#' @param annotation_file Optional annotation XLSX.
#' @param sample_id_col Metadata column containing sample identifiers. Defaults
#'   to `SampleID`.
#' @param variable Metadata column used as the grouping variable.
#' @param test Analysis method. Defaults to automatic selection from the
#'   number of levels in `variable`.
#' @param pair_id Optional pairing identifier.
#' @param covariates Optional ANOVA covariates.
#' @param include_warn If `TRUE`, retain WARN samples. Defaults to `TRUE`.
#'   If `FALSE`, WARN samples are removed.
#' @param include_fail If `TRUE`, retain samples with final sample-level
#'   `FAIL` status. Defaults to `FALSE`.
#' @param output_dir Directory in which to save the metadata used for analysis.
#' @param run_functional_analysis If `TRUE`, run GO and PPI analysis. Defaults
#'   to `TRUE`.
#' @param enrichment_timeout Best-effort R elapsed time limit for each complete
#'   GO/PPI stage, including downloads and figure export. Defaults to `Inf`
#'   (no stage limit); first-time STRING downloads can take several minutes.
#'   This is not a guaranteed HTTP deadline: external libraries may handle
#'   time-limit errors internally.
#' @param logFC Minimum absolute Log2 fold-change for DEPs. Defaults to `0`.
#' @param adj.p Maximum adjusted p-value for DEPs. Defaults to `0.05`.
#' @param boxplot_top_n Optional positive integer limiting individual boxplots to
#'   the top DEPs by adjusted p-value. NULL (default) plots all assays. This
#'   changes only boxplot exports, not statistical tests, DEP cutoffs, GO or PPI.
#' @param cluster_scale_rows If `TRUE`, use per-assay z-scores only for ANOVA
#'   clustering. Defaults to `FALSE` to retain raw NPX distances. Statistical
#'   testing, fold changes, and raw NPX heatmaps are not transformed.
#' @return A structured pipeline result. Files are organized under `Results`,
#'   `QC`, `Tables`, and `log`, with the two metadata CSV files at the root.
#'   The `output_files` list includes diagnostic, execution-log and session paths.
#' @export
npx_pipeline <- function(data_file, meta_file, annotation_file = npx_annotation_file(), sample_id_col = "SampleID", variable = "Condition", test = NULL, pair_id = NULL, covariates = NULL, include_warn = TRUE, include_fail = FALSE, output_dir = "NPXplore_output", logFC = 0, adj.p = 0.05, run_functional_analysis = TRUE, enrichment_timeout = Inf, cluster_scale_rows = FALSE, boxplot_top_n = NULL) {
  npx_validate_dep_cutoffs(logFC, adj.p)
  if (!is.logical(include_warn) || length(include_warn) != 1L || is.na(include_warn) ||
      !is.logical(include_fail) || length(include_fail) != 1L || is.na(include_fail)) {
    stop("include_warn and include_fail must be one non-missing logical value: TRUE or FALSE.")
  }
  if (!is.logical(run_functional_analysis) || length(run_functional_analysis) != 1L || is.na(run_functional_analysis)) {
    stop("run_functional_analysis must be one non-missing logical value: TRUE or FALSE.")
  }
  if (!is.logical(cluster_scale_rows) || length(cluster_scale_rows) != 1L || is.na(cluster_scale_rows)) {
    stop("cluster_scale_rows must be TRUE or FALSE.")
  }
  if (!is.numeric(enrichment_timeout) || length(enrichment_timeout) != 1L || is.na(enrichment_timeout) || enrichment_timeout <= 0) {
    stop("enrichment_timeout must be one positive number of seconds.")
  }
  if (!is.null(boxplot_top_n) && (!is.numeric(boxplot_top_n) || length(boxplot_top_n) != 1L ||
      !is.finite(boxplot_top_n) || boxplot_top_n < 1 || boxplot_top_n != as.integer(boxplot_top_n))) {
    stop("boxplot_top_n must be NULL or one positive integer.")
  }
  log_connection <- NULL
  log_file <- NULL
  if (!is.null(output_dir)) {
    log_dir <- file.path(output_dir, "log")
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(log_dir, "NPXplore_run.log")
    log_connection <- file(log_file, open = "wt")
    on.exit(close(log_connection), add = TRUE)
    writeLines(c(paste("NPXplore", utils::packageVersion("NPXplore")),
                 paste("Started:", Sys.time()), paste("Variable:", variable),
                 paste("logFC:", logFC, "adj.p:", adj.p)), log_connection)
  }
  write_log <- function(kind, condition) {
    if (!is.null(log_connection)) {
      writeLines(paste0("[", kind, "] ", conditionMessage(condition)), log_connection)
      flush(log_connection)
    }
  }
  withCallingHandlers({
  inputs <- npx_read_inputs(data_file, meta_file, annotation_file, sample_id_col)
  # QC reports intentionally use the original data, matching the Rmd. The
  # metadata/SAMPLE/assay collection is kept separate from cleaned analysis
  # data so QC counts are not changed by clean_npx().
  raw_ht <- dplyr::filter(inputs$data, SampleID %in% inputs$meta$SampleID, SampleType == "SAMPLE")
  raw_ht_assay <- dplyr::filter(raw_ht, AssayType == "assay")
  qc <- npx_qc_reports(raw_ht)
  qc_meta <- npx_qc_metadata(qc, inputs$meta)
  qc_columns <- setdiff(names(qc_meta), "SampleID")
  analysis_input_meta <- dplyr::select(inputs$meta, -dplyr::any_of(qc_columns))
  analysis_meta <- qc_meta |>
    dplyr::filter(
      .data$QC == "PASS" |
        (.data$QC == "WARN" & include_warn) |
        (.data$QC == "FAIL" & include_fail)
    ) |>
    dplyr::inner_join(analysis_input_meta, by = "SampleID")
  cleaned <- npx_clean(inputs$data)
  prepared <- npx_prepare(cleaned$data, analysis_meta, include_warn = TRUE, include_fail = TRUE)
  required_metadata <- unique(c(variable, pair_id, covariates))
  missing_columns <- setdiff(required_metadata, names(prepared))
  if (length(missing_columns)) stop("Model metadata column(s) not found: ", paste(missing_columns, collapse = ", "))
  incomplete <- vapply(prepared[required_metadata], function(x) any(is.na(x) | !nzchar(trimws(as.character(x)))), logical(1))
  if (any(incomplete)) stop("Missing or blank model metadata in: ", paste(required_metadata[incomplete], collapse = ", "), ". Complete or remove those samples before running the pipeline.")
  meta_used <- analysis_meta[match(unique(prepared$SampleID), analysis_meta$SampleID), , drop = FALSE]
  metadata_file <- NULL
  analysis_metadata_file <- NULL
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    analysis_metadata_file <- file.path(output_dir, "QCreport_metadata.csv")
    qc_report_metadata_output <- dplyr::select(
      analysis_meta,
      -dplyr::any_of("QCIssue")
    )
    utils::write.csv(qc_report_metadata_output, analysis_metadata_file, row.names = FALSE)
    metadata_file <- file.path(output_dir, "metadata_used.csv")
    metadata_used_output <- dplyr::select(
      meta_used,
      -dplyr::any_of(c("FailBlock", "WarnBlock", "PassBlock", "QCIssue"))
    )
    utils::write.csv(metadata_used_output, metadata_file, row.names = FALSE)
  }
  message("[NPXplore] Running statistical analysis for ", dplyr::n_distinct(prepared$SampleID), " samples.")
  analysis <- npx_analyze(prepared, variable, test, pair_id, covariates, cleaned$check_log_clean, logFC, adj.p)
  functional <- NULL; ppi <- NULL; cluster_result <- NULL; bp <- list(status = "NotRun")
  deps_file <- NULL; output_files <- list(metadata_used = metadata_file, analysis_metadata = analysis_metadata_file, run_log = log_file)
  if (!is.null(output_dir)) {
    output_files$qc <- npx_write_qc(qc, output_dir, inputs$meta)
    tables_dir <- file.path(output_dir, "Tables")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    raw_file <- file.path(tables_dir, "0.Raw_NPX_data.csv")
    utils::write.csv(raw_ht_assay, raw_file, row.names = FALSE)
    output_files$raw_data <- raw_file
    statistical_output <- npx_statistical_table_for_output(analysis$result, analysis)
    deps_output <- npx_statistical_table_for_output(analysis$deps, analysis)
    deps_file <- file.path(tables_dir, "DEPs.csv")
    utils::write.csv(deps_output, deps_file, row.names = FALSE)
    output_files$deps <- deps_file
    result_file <- file.path(tables_dir, "2.Statistical_analysis_results.csv")
    utils::write.csv(statistical_output, result_file, row.names = FALSE)
    output_files$analysis <- result_file
    if (!is.null(analysis$excluded_assay_reasons)) {
      excluded_file <- file.path(tables_dir, "Excluded_assays.csv")
      utils::write.csv(analysis$excluded_assay_reasons, excluded_file, row.names = FALSE)
      output_files$excluded_assays <- excluded_file
    }
    npx_file <- file.path(tables_dir, "1.NPX_table.csv")
    npx_mat <- npx_matrix(prepared)
    for (j in seq_len(ncol(npx_mat))) {
      med <- stats::median(npx_mat[, j], na.rm = TRUE)
      if (is.finite(med)) npx_mat[is.na(npx_mat[, j]), j] <- med
    }
    npx_tab <- as.data.frame(t(npx_mat)); npx_tab <- cbind(Assay = rownames(npx_tab), npx_tab)
    utils::write.csv(npx_tab, npx_file, row.names = FALSE); output_files$npx <- npx_file
    figures_dir <- file.path(output_dir, "Results")
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      message("[NPXplore] Saving UMAP and sample QC figures.")
      u <- npx_try(npx_umap(prepared, variable = variable), "UMAP")
      if (!is.null(u)) {
        umap_dir <- file.path(figures_dir, "UMAP")
        dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
        f <- file.path(umap_dir, "UMAP.pdf")
        ggplot2::ggsave(f, u$plot, width = u$width, height = u$height, units = "in")
        png_file <- file.path(umap_dir, "UMAP.png")
        ggplot2::ggsave(
          png_file, u$plot, width = u$width, height = u$height,
          units = "in", dpi = 300, bg = "white"
        )
        output_files$umap <- c(pdf = f, png = png_file)
        if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
          html_file <- file.path(umap_dir, "UMAP.html")
          output_files$umap_html <- npx_save_pipeline_html(u$plot, html_file)
        }
      }
      qc_umap_data <- dplyr::left_join(
        raw_ht_assay[, c("SampleID", "Assay", "NPX"), drop = FALSE],
        qc_meta,
        by = "SampleID"
      )
      qc_umap <- npx_try(npx_umap(qc_umap_data, variable = "QC"), "QC UMAP")
      if (!is.null(qc_umap)) {
        f <- file.path(output_dir, "QC", "UMAP_QC_status.pdf")
        dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
        ggplot2::ggsave(f, qc_umap$plot, width = qc_umap$width, height = qc_umap$height, units = "in")
        png_file <- file.path(output_dir, "QC", "UMAP_QC_status.png")
        ggplot2::ggsave(png_file, qc_umap$plot, width = qc_umap$width, height = qc_umap$height,
                        units = "in", dpi = 600, bg = "white")
        output_files$qc_umap <- c(pdf = f, png = png_file)
        if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
          html_file <- file.path(output_dir, "QC", "UMAP_QC_status.html")
          output_files$qc_umap_html <- npx_save_pipeline_html(qc_umap$plot, html_file)
        } else {
          warning("Packages 'plotly' and 'htmlwidgets' are required for the interactive QC UMAP HTML output; PDF and PNG were saved.", call. = FALSE)
        }
      }
      distribution_dir <- file.path(figures_dir, "Distribution")
      iqr <- npx_try(npx_iqrplot(prepared, variable = variable, output_dir = distribution_dir), "IQR plot")
      if (!is.null(iqr)) output_files$iqr <- unlist(iqr$files, use.names = TRUE)
      distribution <- npx_try(npx_distribution_boxplot(prepared, variable = variable, output_dir = distribution_dir), "sample distribution boxplot")
      if (!is.null(distribution)) {
        output_files$distribution_boxplot <- unlist(distribution$files, use.names = TRUE)
      }
    }
    message("[NPXplore] Saving the full NPX heatmap; large panels can take several minutes.")
    heatmap_dir <- file.path(figures_dir, "Heatmap")
    dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
    heatmap_pdf <- file.path(heatmap_dir, "Heatmap_NPX.pdf")
    all_hm <- npx_try(npx_heatmap(prepared, variable = variable, scale_rows = FALSE, filename = heatmap_pdf), "analysis heatmap")
    if (!is.null(all_hm)) {
      heatmap_png <- file.path(heatmap_dir, "Heatmap_NPX.png")
      npx_save_plot(all_hm, heatmap_png, width = 12, height = 10, dpi = 600)
      output_files$analysis_heatmap <- c(pdf = heatmap_pdf, png = heatmap_png)
    }
    dep_assays <- unique(as.character(analysis$deps$Assay)); dep_assays <- dep_assays[!is.na(dep_assays) & nzchar(dep_assays)]
    dep_figures_dir <- file.path(figures_dir, "DEPs")
    dir.create(dep_figures_dir, recursive = TRUE, showWarnings = FALSE)
    box_dir <- file.path(figures_dir, "DEPs", "Boxplot")
    message("[NPXplore] Saving individual boxplots: ", if (is.null(boxplot_top_n)) "all assays" else paste0("top ", boxplot_top_n, " DEPs"), ".")
    bp <- npx_run_stage(
      npx_dep_boxplot(
        prepared,
        analysis = analysis,
        meta = analysis_meta,
        variable = variable,
        assays = npx_boxplot_selection(analysis, boxplot_top_n),
        output_dir = box_dir
      ),
      "Assay boxplots"
    )
    if (identical(bp$status, "Success")) {
      output_files$dep_boxplot <- unlist(bp$files, use.names = TRUE)
    }
    if (length(dep_assays)) {
      message("[NPXplore] Saving DEP heatmaps and analysis figures.")
      dep_heatmap_file <- file.path(dep_figures_dir, "DEPs_heatmap.pdf")
      hm <- npx_try(npx_heatmap(prepared, variable = variable, assays = dep_assays, scale_rows = FALSE, filename = dep_heatmap_file), "DEP heatmap")
      if (!is.null(hm)) output_files$heatmap <- dep_heatmap_file
      if (length(analysis$levels) == 2L) {
        volcano_dir <- file.path(figures_dir, "DEPs", "VolcanoPlot")
        vc <- npx_try(
          npx_volcano(analysis, logFC = logFC, adj.p = adj.p, output_dir = volcano_dir),
          "volcano plot"
        )
        if (!is.null(vc) && !is.null(vc$files)) output_files$volcano <- vc$files
      }
      if ("Human_Pathway" %in% names(inputs$anno)) {
        pathway_dir <- file.path(figures_dir, "DEPs", "Heatmap_Pathway")
        pathway_heatmaps <- npx_try(
          npx_pathway_heatmaps(
            prepared,
            analysis = analysis,
            anno = inputs$anno,
            variable = variable,
            output_dir = pathway_dir
          ),
          "DEP pathway heatmaps"
        )
        if (!is.null(pathway_heatmaps)) {
          output_files$pathway_heatmaps <- pathway_heatmaps$files
        }
      }
    }
    cluster_result <- NULL
    if (length(analysis$levels) > 2L) {
      posthoc <- npx_try(npx_anova_posthoc(prepared, analysis, variable, covariates, cleaned$check_log_clean, logFC, adj.p), "ANOVA posthoc")
      if (is.null(posthoc)) posthoc <- data.frame()
      posthoc_file <- file.path(tables_dir, "3.ANOVA_posthoc_test.csv")
      utils::write.csv(npx_posthoc_table_for_output(posthoc), posthoc_file, row.names = FALSE)
      output_files$anova_posthoc <- posthoc_file
      cl <- npx_run_stage(npx_cluster_assays(prepared, analysis, scale_rows = cluster_scale_rows), "Clustering")
      cluster_result <- cl
      if (identical(cl$status, "Success")) {
        cluster_file <- file.path(tables_dir, "DEPs_with_ClusterNo.csv")
        cluster_export <- cl$clusters |>
          dplyr::transmute(Assay = as.character(.data$Assay), ClusterNo = as.integer(.data$cluster))
        clustered_deps <- dplyr::left_join(deps_output, cluster_export, by = "Assay")
        utils::write.csv(clustered_deps, cluster_file, row.names = FALSE)
        output_files$clusters <- cluster_file
        cluster_dir <- file.path(figures_dir, "DEPs", "Clustering")
        cluster_hm <- npx_try(
          npx_cluster_heatmap(prepared, analysis, variable = variable, cluster_result = cl, output_dir = cluster_dir),
          "cluster heatmap"
        )
        if (!is.null(cluster_hm)) output_files$cluster_heatmap <- cluster_hm$files
        cluster_line <- npx_try(
          npx_cluster_lineplot(prepared, analysis, variable = variable, cluster_result = cl, output_dir = cluster_dir),
          "cluster lineplot"
        )
        if (!is.null(cluster_line)) output_files$cluster_lineplot <- cluster_line$files
      }
    }
    cluster_blocked <- identical(analysis$method, "anova") && !identical(cluster_result$status, "Success")
    functional <- if (run_functional_analysis && cluster_blocked) npx_blocked_stage(
      paste0("GO not attempted because clustering is unavailable: ", cluster_result$reason)
    ) else if (run_functional_analysis) npx_run_stage(
      npx_functional_analysis(
        analysis = analysis,
        data = prepared,
        cluster_result = cluster_result,
        output_dir = output_dir,
        databases = npx_go_databases,
        padj_cutoff = 0.05,
        top_n = 10
      ),
      "GO functional analysis",
      enrichment_timeout
    ) else NULL
    if (run_functional_analysis && cluster_blocked) npx_stage_message("GO functional analysis", functional)
    if (!is.null(functional)) {
      output_files$gene_ontology <- functional$files
    }
    ppi <- if (run_functional_analysis) npx_run_stage(
      npx_ppi_network(
        analysis = analysis,
        data = prepared,
        anno = inputs$anno,
        cluster_result = cluster_result,
        output_dir = output_dir,
        score_threshold = 700,
        min_proteins = 3,
        width = 10,
        height = 8,
        dpi = 600,
        seed = 123
      ),
      "PPI network analysis",
      enrichment_timeout
    ) else NULL
    if (!is.null(ppi)) {
      output_files$ppi_network <- ppi$files
      output_files$ppi_summary <- ppi$summary_file
      output_files$string_mapping <- ppi$mapping_file
      output_files$ppi_pathway_mapping <- ppi$pathway_mapping_file
    }
  }
  if (!run_functional_analysis || is.null(output_dir)) {
    reason <- if (!run_functional_analysis) "run_functional_analysis is FALSE." else "output_dir is NULL; functional output stages were not run."
    functional <- ppi <- list(status = "NotRun", code = "Disabled", reason = reason, files = character(), results = list())
    npx_stage_message("GO functional analysis", functional)
    npx_stage_message("PPI network analysis", ppi)
  }
  diagnostics <- npx_stage_summary(list(Clustering = cluster_result, GO = functional, PPI = ppi))
  if (!is.null(output_dir)) {
    status_file <- file.path(output_dir, "log", "Functional_analysis_status.csv")
    utils::write.csv(diagnostics, status_file, row.names = FALSE)
    output_files$functional_status <- status_file
    session_file <- file.path(output_dir, "log", "sessionInfo.txt")
    writeLines(utils::capture.output(utils::sessionInfo()), session_file)
    output_files$session_info <- session_file
    writeLines(paste("Completed:", Sys.time()), log_connection)
  }
  method_label <- switch(analysis$method, anova = "ANOVA", ttest = "T-test", wilcox = "Wilcoxon", analysis$method)
  message("[NPXplore] ", method_label, " analysis complete. Stage status: Assay boxplots=", bp$status,
          "; GO=", functional$status, "; PPI=", ppi$status, ".")
  structure(list(data = prepared, meta = meta_used, anno = inputs$anno, check_log = cleaned$check_log, check_log_clean = cleaned$check_log_clean, qc = qc, analysis = analysis, cluster = cluster_result, functional = functional, ppi = ppi, diagnostics = diagnostics, output_files = output_files), class = "npxplore_result")
  }, message = function(c) write_log("message", c),
     warning = function(c) write_log("warning", c),
     error = function(c) write_log("error", c))
}

#' Print an NPXplore result
#' @param x An `npxplore_result` object.
#' @param ... Unused additional arguments.
#' @export
print.npxplore_result <- function(x, ...) { cat("NPXplore result\n  method: ", x$analysis$method, "\n  conditions: ", paste(x$analysis$levels, collapse = ", "), "\n  assay rows: ", nrow(x$data), "\n", sep = ""); invisible(x) }
