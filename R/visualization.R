#' Make an assay matrix from prepared NPX data
#' @param data Prepared Olink assay data.
#' @param assays Optional assay names to retain.
#' @return A numeric matrix with samples in rows.
#' @keywords internal
npx_matrix <- function(data, assays = NULL) {
  required <- c("SampleID", "Assay", "NPX")
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("data is missing required columns: ", paste(missing, collapse = ", "))
  if (!is.null(assays)) data <- dplyr::filter(data, Assay %in% assays)
  if (!nrow(data)) stop("No assays remain for visualization.")
  wide <- data |>
    dplyr::select(SampleID, Assay, NPX) |>
    dplyr::group_by(SampleID, Assay) |>
    dplyr::summarise(NPX = if (all(is.na(NPX))) NA_real_ else stats::median(NPX, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Assay, values_from = NPX)
  ids <- wide$SampleID
  matrix <- as.matrix(wide[, setdiff(names(wide), "SampleID"), drop = FALSE])
  rownames(matrix) <- ids
  matrix <- matrix[, colSums(!is.na(matrix)) > 0, drop = FALSE]
  if (!ncol(matrix)) stop("No assays contain usable NPX values.")
  matrix
}

npx_umap_dimensions <- function(umap_data, base_width = 8, min_height = 4, max_height = 10) {
  x_range <- diff(range(umap_data$UMAP1, finite = TRUE))
  y_range <- diff(range(umap_data$UMAP2, finite = TRUE))
  if (!is.finite(x_range) || x_range <= 0) x_range <- 1
  if (!is.finite(y_range) || y_range <= 0) y_range <- 1
  height <- base_width * y_range / x_range
  height <- max(min_height, min(max_height, height))
  width <- height * x_range / y_range
  list(width = width, height = height)
}

#' Calculate and plot a UMAP embedding
#' @param data Prepared assay data.
#' @param variable Metadata column used for point colour.
#' @param assays Optional assays to use.
#' @param seed Random seed.
#' @param n_neighbors Number of UMAP neighbours.
#' @return A list with `data`, a ggplot object in `plot`, and the output
#'   dimensions in `width` and `height` (6 x 5 inches for analysis UMAP and
#'   8 x 8 inches for QC UMAP).
#' @export
npx_umap <- function(data, variable = "Condition", assays = NULL, seed = 123, n_neighbors = 15) {
  if (!requireNamespace("umap", quietly = TRUE)) stop("Package 'umap' is required for npx_umap().")
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  mat <- npx_matrix(data, assays)
  if (ncol(mat) < 2L) stop("UMAP requires at least 2 variable assays.")
  for (j in seq_len(ncol(mat))) {
    missing <- is.na(mat[, j])
    mat[missing, j] <- stats::median(mat[, j], na.rm = TRUE)
  }
  mat <- mat[, apply(mat, 2, stats::sd, na.rm = TRUE) > 0, drop = FALSE]
  if (nrow(mat) < 3L || ncol(mat) < 2L) stop("UMAP requires at least 3 samples and 2 variable assays.")
  cfg <- umap::umap.defaults
  cfg$n_components <- 2L
  cfg$n_neighbors <- max(2L, min(as.integer(n_neighbors), nrow(mat) - 1L))
  cfg$random_state <- as.integer(seed)
  set.seed(seed)
  embedding <- umap::umap(mat, config = cfg)$layout
  out <- data.frame(UMAP1 = embedding[, 1], UMAP2 = embedding[, 2], SampleID = rownames(mat), stringsAsFactors = FALSE)
  labels <- unique(data[, c("SampleID", variable), drop = FALSE])
  out <- dplyr::left_join(out, labels, by = "SampleID")
  condition_order <- if (identical(variable, "QC")) {
    qc_levels <- unique(trimws(as.character(out[[variable]])))
    qc_levels <- qc_levels[!is.na(qc_levels) & nzchar(qc_levels)]
    preferred <- c("PASS", "WARN", "FAIL")
    list(plot = c(intersect(preferred, qc_levels), setdiff(qc_levels, preferred)))
  } else npx_condition_levels(out[[variable]])
  out[[variable]] <- factor(trimws(as.character(out[[variable]])), levels = condition_order$plot)
  out$hover_text <- paste0(
    "SampleID: ", out$SampleID,
    "<br>", variable, ": ", out[[variable]],
    "<br>UMAP1: ", round(out$UMAP1, 3),
    "<br>UMAP2: ", round(out$UMAP2, 3)
  )
  is_qc <- identical(variable, "QC")
  out <- out |>
    dplyr::filter(
      !is.na(.data[[variable]]),
      is.finite(.data$UMAP1),
      is.finite(.data$UMAP2)
    )
  palette <- if (is_qc) {
    c(PASS = "royalblue2", WARN = "goldenrod2", FAIL = "firebrick2")
  } else {
    c("#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
      "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52")
  }
  palette <- stats::setNames(rep(palette, length.out = length(condition_order$plot)), condition_order$plot)
  if (is_qc) {
    plot <- ggplot2::ggplot(
      out,
      ggplot2::aes(x = .data[["UMAP1"]], y = .data[["UMAP2"]],
                   fill = .data[[variable]], label = .data[["SampleID"]],
                   text = .data[["hover_text"]])
    ) +
      ggplot2::geom_point(size = 3, shape = 21, colour = "black", stroke = 0.6, alpha = 0.9) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE, na.value = "grey70") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 3),
        axis.line = ggplot2::element_blank()
      ) +
      ggplot2::labs(title = "UMAP: QC status", fill = "QC")
  } else {
    plot <- ggplot2::ggplot(
      out,
      ggplot2::aes(x = .data[["UMAP1"]], y = .data[["UMAP2"]],
                   fill = .data[[variable]], label = .data[["SampleID"]],
                   text = .data[["hover_text"]])
    ) +
      ggplot2::geom_point(shape = 21, size = 3.5, colour = "black", stroke = 0.8, alpha = 0.9) +
      ggplot2::scale_fill_manual(values = palette, drop = FALSE, na.value = "grey70") +
      ggplot2::theme_classic() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 3),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(colour = "black", linewidth = 1.0),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text = ggplot2::element_text(size = 13, colour = "black"),
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
        legend.title = ggplot2::element_text(size = 14, face = "bold"),
        legend.text = ggplot2::element_text(size = 12),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA)
      ) +
      ggplot2::labs(title = "UMAP", x = "UMAP1", y = "UMAP2", fill = variable)
  }
  dimensions <- if (is_qc) list(width = 8, height = 8) else list(width = 6, height = 5)
  list(data = out, plot = plot, width = dimensions$width, height = dimensions$height)
}

#' Plot NPX distributions for selected assays
#' @param data Prepared assay data.
#' @param variable Metadata column used for the x-axis.
#' @param assays Assays to plot. Defaults to the first 12 assays.
#' @param mode Plot mode: `"individual"` (default) returns one plot per assay;
#'   `"combined"` returns one faceted plot containing all selected assays.
#' @param output_dir Optional directory. In `"individual"` mode, one PDF per
#'   assay is saved here.
#' @return A list with plotting data, a named `plots` list, and `plot` for
#'   combined mode. Individual plots are also available in `plots`.
#' @export
npx_boxplot <- function(data, variable = "Condition", assays = NULL,
                        mode = "individual", output_dir = NULL) {
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  mode <- match.arg(mode, c("individual", "combined"))
  if (is.null(assays)) assays <- unique(as.character(data$Assay))[seq_len(min(12L, length(unique(data$Assay))))]
  plot_data <- dplyr::filter(data, Assay %in% assays, !is.na(NPX))
  if (!nrow(plot_data)) stop("No observations remain for boxplot.")
  condition_order <- npx_condition_levels(plot_data[[variable]])
  plot_data[[variable]] <- factor(trimws(as.character(plot_data[[variable]])), levels = condition_order$plot)
  make_plot <- function(one_assay) {
    ggplot2::ggplot(dplyr::filter(plot_data, Assay == one_assay), ggplot2::aes(x = .data[[variable]], y = NPX, fill = .data[[variable]])) +
      ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
      ggplot2::theme_classic() + ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = one_assay, x = variable, y = "NPX")
  }
  plots <- stats::setNames(lapply(assays, make_plot), assays)
  if (!is.null(output_dir)) {
    if (mode != "individual") stop("output_dir is supported only in individual mode.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    for (assay in names(plots)) {
      file_name <- paste0(gsub("[^A-Za-z0-9_.-]+", "_", assay), ".pdf")
      ggplot2::ggsave(file.path(output_dir, file_name), plots[[assay]], width = 6, height = 5, units = "in")
    }
  }
  combined <- NULL
  if (mode == "combined") {
    combined <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[variable]], y = NPX, fill = .data[[variable]])) +
      ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 0.8) +
      ggplot2::facet_wrap(~Assay, scales = "free_y") + ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none") + ggplot2::labs(x = variable, y = "NPX")
  }
  list(data = plot_data, plots = plots, plot = combined, output_dir = output_dir)
}

#' Plot per-assay boxplots annotated with differential-analysis statistics
#'
#' The function selects the assay-specific p-value and adjusted p-value from
#' an analysis object returned by [npx_ttest()], [npx_wilcox()], or
#' [npx_anova()]. Condition levels are ordered with Control/Healthy/HC first,
#' matching the other NPXplore figures.
#' @param data NPX assay data. It may already contain `variable` or receive it
#'   from `meta` by `SampleID`.
#' @param analysis A result returned by `npx_ttest()`, `npx_wilcox()`, or
#'   `npx_anova()`, or a corresponding result data frame.
#' @param meta Optional metadata data frame containing `SampleID` and
#'   `variable`.
#' @param variable Metadata column used for the x-axis.
#' @param assays Assays to plot. Defaults to all assays in `data`, matching
#'   the supplied example code. Use this argument to restrict output to a
#'   selected DEP set.
#' @param output_dir Optional output directory. PDF and 600 dpi PNG files are
#'   written there when supplied.
#' @return A list with `data`, named `plots`, `files`, `method`, `width`, and
#'   `height`.
#' @export
npx_dep_boxplot <- function(data, analysis, meta = NULL, variable = "Condition",
                            assays = NULL, output_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for npx_dep_boxplot().")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for npx_dep_boxplot().")
  }
  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop("Package 'ggtext' is required for npx_dep_boxplot().")
  }
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (!"SampleID" %in% names(data)) stop("data must contain a 'SampleID' column.")

  stat_table <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$result else analysis
  if (!is.data.frame(stat_table)) stop("analysis must be an analysis result list or data frame.")
  method <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$method else NULL
  if (is.null(method) || !length(method) || is.na(method)) {
    method <- if ("estimate" %in% names(stat_table)) "ttest" else "anova"
  }
  method_label <- if (grepl("wilcox", method, ignore.case = TRUE)) {
    "Wilcoxon"
  } else if (grepl("anova", method, ignore.case = TRUE)) {
    "ANOVA"
  } else {
    "t-test"
  }

  if (!is.null(meta)) {
    if (!is.data.frame(meta) || !all(c("SampleID", variable) %in% names(meta))) {
      stop("meta must contain 'SampleID' and the grouping column '", variable, "'.")
    }
    data <- dplyr::select(data, -dplyr::any_of(variable)) |>
      dplyr::left_join(
        dplyr::distinct(meta[, c("SampleID", variable), drop = FALSE]),
        by = "SampleID"
      )
  }
  if (!variable %in% names(data)) {
    stop("Grouping column not found in data or meta: ", variable)
  }
  if (!"Assay" %in% names(data)) stop("data must contain an 'Assay' column.")
  if (!all(c("p.value", "Adjusted_pval") %in% names(stat_table))) {
    stop("analysis must contain 'p.value' and 'Adjusted_pval' columns.")
  }

  if (is.null(assays)) {
    assays <- unique(as.character(data$Assay))
  }
  assays <- unique(as.character(assays))
  assays <- assays[!is.na(assays) & nzchar(assays)]
  if (!length(assays)) stop("No assays were selected for DEP boxplots.")

  plot_data <- data |>
    dplyr::mutate(
      .SampleID = as.character(.data$SampleID),
      .Condition = trimws(as.character(.data[[variable]])),
      .Assay = as.character(.data$Assay)
    ) |>
    dplyr::filter(
      .data$.Assay %in% assays,
      !is.na(.data$NPX),
      is.finite(.data$NPX),
      !is.na(.data$.Condition),
      nzchar(.data$.Condition)
    )
  if (!nrow(plot_data)) stop("No observations remain for DEP boxplots.")

  condition_levels <- npx_condition_levels(plot_data$.Condition)$plot
  plot_data$.Condition <- factor(plot_data$.Condition, levels = condition_levels)
  plot_data$Condition_num <- as.numeric(plot_data$.Condition)
  plot_data <- dplyr::filter(plot_data, !is.na(.data$Condition_num))

  plotly_palette <- c(
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
  )
  condition_colors <- stats::setNames(
    rep(plotly_palette, length.out = length(condition_levels)),
    condition_levels
  )

  stat_assay_col <- if ("Assay" %in% names(stat_table)) "Assay" else NULL
  if (is.null(stat_assay_col)) stop("analysis must contain an 'Assay' column.")
  stat_table$.Assay <- as.character(stat_table[[stat_assay_col]])

  format_stat <- function(x) {
    if (!length(x) || is.na(x) || !is.finite(as.numeric(x))) return("NA")
    format.pval(as.numeric(x), digits = 3, eps = 1e-300)
  }

  make_plot <- function(assay) {
    d <- dplyr::filter(plot_data, .data$.Assay == assay)
    test_result <- stat_table |>
      dplyr::filter(.data$.Assay == assay) |>
      dplyr::slice(1)
    range_df <- d |>
      dplyr::group_by(.data$Condition_num) |>
      dplyr::summarise(
        ymin = min(.data$NPX, na.rm = TRUE),
        ymax = max(.data$NPX, na.rm = TRUE),
        .groups = "drop"
      )
    p_text <- if (nrow(test_result) == 1L) {
      paste0(
        "<i>P</i> value = ", format_stat(test_result$p.value),
        "   |   Adjusted <i>P</i> = ", format_stat(test_result$Adjusted_pval)
      )
    } else {
      paste0("No ", method_label, " result")
    }
    ggplot2::ggplot(
      d,
      ggplot2::aes(x = .data$Condition_num, y = .data$NPX)
    ) +
      ggplot2::geom_errorbar(
        data = range_df,
        ggplot2::aes(x = .data$Condition_num, ymin = .data$ymin, ymax = .data$ymax),
        inherit.aes = FALSE, width = 0.10, colour = "black", linewidth = 0.7
      ) +
      ggplot2::geom_boxplot(
        ggplot2::aes(group = .data$.Condition, fill = .data$.Condition),
        width = 0.55, outlier.shape = NA, colour = "black", linewidth = 0.7,
        alpha = 1
      ) +
      ggplot2::geom_jitter(
        width = 0.08, size = 1.35, alpha = 0.8, colour = "black", shape = 16
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(condition_levels), labels = condition_levels,
        limits = c(0.5, length(condition_levels) + 0.5)
      ) +
      ggplot2::scale_fill_manual(values = condition_colors[condition_levels], drop = FALSE) +
      ggplot2::scale_y_continuous(
        breaks = scales::breaks_pretty(n = 5),
        labels = scales::label_number(accuracy = 1),
        expand = ggplot2::expansion(mult = c(0.05, 0.12))
      ) +
      ggplot2::labs(title = assay, subtitle = p_text, x = variable, y = "NPX") +
      ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.5),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(colour = "black", linewidth = 1),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 11, colour = "black"),
        axis.text.y = ggplot2::element_text(size = 13, colour = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5, size = 18, face = "bold"),
        plot.subtitle = ggtext::element_markdown(hjust = 0.5, size = 11),
        legend.position = "none"
      )
  }

  plots <- stats::setNames(lapply(assays, make_plot), assays)
  files <- list(pdf = character(), png = character())
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    safe_name <- function(x) {
      x <- gsub("[^A-Za-z0-9_-]+", "_", x)
      x <- gsub("_+", "_", x)
      gsub("^_|_$", "", x)
    }
    files$pdf <- stats::setNames(
      vapply(names(plots), function(assay) {
        file <- file.path(output_dir, paste0(safe_name(assay), "_boxplot.pdf"))
        ggplot2::ggsave(file, plots[[assay]], width = 6, height = 9,
                        units = "in", device = grDevices::pdf)
        file
      }, character(1L)), names(plots))
    files$png <- stats::setNames(
      vapply(names(plots), function(assay) {
        file <- file.path(output_dir, paste0(safe_name(assay), "_boxplot.png"))
        ggplot2::ggsave(file, plots[[assay]], width = 6, height = 9,
                        units = "in", dpi = 600, bg = "white")
        file
      }, character(1L)), names(plots))
  }
  list(
    data = plot_data,
    plots = plots,
    files = files,
    method = method,
    width = 6,
    height = 9
  )
}

#' Plot NPX distribution by sample and condition
#' @param data Prepared assay data.
#' @param variable Metadata column containing sample groups.
#' @param output_dir Optional directory for PDF, PNG, and HTML outputs.
#' @return A list containing the imputed long-format data, static plot,
#'   optional interactive plot, and output file paths.
#' @export
npx_distribution_boxplot <- function(data, variable = "Condition", output_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for npx_distribution_boxplot().")
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  required <- c("SampleID", variable)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("data is missing required columns: ", paste(missing, collapse = ", "))

  mat <- npx_matrix(data)
  for (j in seq_len(ncol(mat))) {
    med <- stats::median(mat[, j], na.rm = TRUE)
    if (is.finite(med)) mat[is.na(mat[, j]), j] <- med
  }

  meta_condition <- data |>
    dplyr::select(dplyr::all_of(c("SampleID", variable))) |>
    dplyr::mutate(
      SampleID = as.character(.data$SampleID),
      .condition = trimws(as.character(.data[[variable]]))
    ) |>
    dplyr::distinct(.data$SampleID, .keep_all = TRUE)
  condition_order <- npx_condition_levels(meta_condition$.condition)$plot
  sample_order <- meta_condition$SampleID[meta_condition$SampleID %in% rownames(mat)]

  npx_long <- as.data.frame(mat, check.names = FALSE) |>
    tibble::rownames_to_column("SampleID") |>
    tidyr::pivot_longer(
      cols = -SampleID,
      names_to = "Assay",
      values_to = "NPX"
    ) |>
    dplyr::left_join(meta_condition[, c("SampleID", ".condition"), drop = FALSE], by = "SampleID") |>
    dplyr::filter(
      is.finite(.data$NPX),
      !is.na(.data$.condition),
      .data$.condition %in% condition_order
    ) |>
    dplyr::mutate(
      SampleID = factor(.data$SampleID, levels = sample_order),
      Condition = factor(.data$.condition, levels = condition_order)
    )

  plotly_palette <- c(
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
  )
  condition_colors <- stats::setNames(
    rep(plotly_palette, length.out = length(condition_order)),
    condition_order
  )

  static_plot <- ggplot2::ggplot(
    npx_long,
    ggplot2::aes(x = .data[["SampleID"]], y = .data[["NPX"]], fill = .data[["Condition"]])
  ) +
    ggplot2::geom_boxplot(
      width = 0.7, outlier.shape = 16, outlier.size = 1.2,
      outlier.alpha = 0.9, colour = "black"
    ) +
    ggplot2::scale_x_discrete(limits = sample_order) +
    ggplot2::scale_fill_manual(values = condition_colors, drop = FALSE) +
    ggplot2::labs(
      title = "NPX Distribution by Sample",
      x = "Sample ID", y = "NPX", fill = variable
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 1),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 13, colour = "black"),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )

  interactive_plot <- NULL
  files <- list(pdf = NULL, png = NULL, html = NULL)
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files$pdf <- file.path(output_dir, "boxplot.pdf")
    files$png <- file.path(output_dir, "boxplot.png")
    ggplot2::ggsave(files$pdf, static_plot, width = 16, height = 10, units = "in")
    ggplot2::ggsave(files$png, static_plot, width = 16, height = 10, units = "in", dpi = 600, bg = "white")
  }

  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    interactive_plot <- plotly::plot_ly(
      data = npx_long,
      x = ~SampleID, y = ~NPX, color = ~Condition,
      colors = condition_colors, type = "box", boxpoints = "outliers",
      marker = list(size = 3, color = "black"), line = list(width = 1),
      hovertemplate = paste0(
        paste0("SampleID: %{x}<br>NPX: %{y}<br>", variable, ": %{fullData.name}<br><extra></extra>")
      )
    ) |>
      plotly::layout(
        title = "NPX Distribution by Sample",
        xaxis = list(title = "Sample ID", categoryorder = "array", categoryarray = sample_order, tickangle = 90, tickfont = list(size = 8)),
        yaxis = list(title = "NPX"),
        boxmode = "overlay",
        legend = list(title = list(text = variable))
      )
    if (!is.null(output_dir)) {
      files$html <- file.path(output_dir, "boxplot.html")
      html_saved <- tryCatch({
        htmlwidgets::saveWidget(interactive_plot, file = files$html, selfcontained = TRUE)
        TRUE
      }, error = function(e) {
        warning("Self-contained distribution HTML could not be saved; saving HTML with external dependencies instead: ", conditionMessage(e), call. = FALSE)
        tryCatch({
          htmlwidgets::saveWidget(interactive_plot, file = files$html, selfcontained = FALSE)
          TRUE
        }, error = function(e2) {
          warning("Distribution HTML output was skipped: ", conditionMessage(e2), call. = FALSE)
          FALSE
        })
      })
      if (!html_saved) files$html <- NULL
    }
  }

  list(data = npx_long, plot = static_plot, interactive = interactive_plot, files = files,
       width = 16, height = 10)
}

#' Plot an NPX heatmap
#' @param data Prepared assay data.
#' @param variable Metadata column used for column annotation.
#' @param assays Optional assays to retain. Defaults to all assays.
#' @param scale_rows If `TRUE`, z-score each assay.
#' @param filename Optional PDF/PNG path. Export space expands to fit annotation text.
#' @param show_rownames Show assay names on heatmap rows.
#' @param show_colnames Show sample names on heatmap columns.
#' @param fontsize_row Row-label font size.
#' @param cellheight Optional row cell height in points.
#' @param width Minimum heatmap export width in inches; expands for long annotation labels.
#' @param height Minimum heatmap export height in inches.
#' @param main Plot title.
#' @return A list with `matrix`, `annotation`, and the pheatmap result.
#' @export
npx_heatmap <- function(data, variable = "Condition", assays = NULL, scale_rows = FALSE,
                        filename = NULL, show_rownames = FALSE, show_colnames = TRUE,
                        fontsize_row = 8, cellheight = NA, width = 12, height = 10,
                        main = "NPX Heatmap by Assay and Condition") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Package 'pheatmap' is required for npx_heatmap().")
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  mat <- npx_matrix(data, assays)
  for (j in seq_len(ncol(mat))) {
    med <- stats::median(mat[, j], na.rm = TRUE)
    if (is.finite(med)) mat[is.na(mat[, j]), j] <- med
  }
  mat[!is.finite(mat)] <- NA_real_
  annotation <- data |>
    dplyr::select(dplyr::all_of(c("SampleID", variable))) |>
    dplyr::mutate(
      SampleID = as.character(.data$SampleID),
      .condition = trimws(as.character(.data[[variable]]))
    ) |>
    dplyr::distinct(.data$SampleID, .keep_all = TRUE) |>
    dplyr::filter(
      !is.na(.data$.condition),
      nzchar(.data$.condition),
      .data$SampleID %in% rownames(mat)
    )
  condition_order <- npx_condition_levels(annotation$.condition)$plot
  annotation <- annotation |>
    dplyr::mutate(.condition = factor(.data$.condition, levels = condition_order))
  annotation <- as.data.frame(annotation)
  rownames(annotation) <- annotation$SampleID
  annotation_col <- data.frame(annotation[[".condition"]], row.names = rownames(annotation))
  names(annotation_col) <- variable
  sample_order <- rownames(annotation_col)[order(annotation_col[[variable]])]
  mat <- mat[sample_order, , drop = FALSE]
  annotation_col <- annotation_col[sample_order, , drop = FALSE]
  if (scale_rows) mat <- scale(mat)
  mat[!is.finite(mat)] <- NA_real_
  mat <- t(mat)
  condition_counts <- table(annotation_col[[variable]])
  gaps_col <- cumsum(condition_counts)
  if (length(gaps_col) > 1L) gaps_col <- gaps_col[-length(gaps_col)] else gaps_col <- NULL
  plotly_palette <- c(
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
  )
  condition_colors <- stats::setNames(
    rep(plotly_palette, length.out = length(condition_order)),
    condition_order
  )
  heatmap_colors <- grDevices::colorRampPalette(c("royalblue2", "white", "brown3"))(100)
  heatmap_range <- range(mat, finite = TRUE)
  if (any(!is.finite(heatmap_range))) stop("No finite NPX values remain for heatmap plotting.")
  if (diff(heatmap_range) == 0) heatmap_range <- heatmap_range + c(-0.5, 0.5)
  heatmap_breaks <- seq(heatmap_range[1], heatmap_range[2], length.out = 101)
  hm <- pheatmap::pheatmap(
    mat = mat,
    color = heatmap_colors,
    breaks = heatmap_breaks,
    annotation_col = annotation_col,
    annotation_colors = stats::setNames(list(condition_colors), variable),
    cluster_rows = nrow(mat) > 1L,
    cluster_cols = FALSE,
    gaps_col = gaps_col,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    angle_col = 90,
    fontsize_col = 6,
    border_color = NA,
    na_col = "grey80",
    main = main,
    filename = NA_character_,
    width = width,
    height = height,
    fontsize_row = fontsize_row,
    cellheight = cellheight,
    silent = !is.null(filename)
  )
  if (!is.null(filename)) hm <- npx_save_heatmap(hm, filename, width, height)
  list(matrix = mat, annotation = annotation_col, plot = hm)
}

# Reserve the actual drawn annotation text width, including the legend key.
# pheatmap's estimated legend column can be narrower than its text grobs.
npx_heatmap_layout <- function(gt, width, height) {
  grDevices::pdf(NULL, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::pushViewport(gt$vp)
  widths <- grid::convertWidth(gt$widths, "inches", valueOnly = TRUE)
  heights <- grid::convertHeight(gt$heights, "inches", valueOnly = TRUE)
  for (name in c("annotation_legend", "col_annotation_names")) {
    index <- which(gt$layout$name == name)
    if (!length(index)) next
    grob <- gt$grobs[[index]]
    texts <- if (inherits(grob, "text")) list(grob) else as.list(grob$children)
    right <- vapply(texts, function(text) {
      if (!inherits(text, "text")) return(0)
      max(grid::convertWidth(text$x + grid::grobWidth(text), "inches", valueOnly = TRUE))
    }, numeric(1))
    column <- gt$layout$l[index]
    widths[column] <- max(widths[column], right + 0.15)
  }
  # Long labels must not leave a zero/negative-width heatmap body.
  matrix_column <- gt$layout$l[gt$layout$name == "matrix"]
  widths[matrix_column] <- max(widths[matrix_column], min(3, width / 2))
  gt$widths <- grid::unit(widths, "inches")
  gt$heights <- grid::unit(heights, "inches")
  grid::popViewport()
  list(gtable = gt, width = max(width, sum(widths)) + 0.2,
       height = max(height, sum(heights)) + 0.2)
}

npx_save_heatmap <- function(hm, filename, width, height, dpi = 300) {
  layout <- npx_heatmap_layout(hm$gtable, width, height)
  extension <- tolower(tools::file_ext(filename))
  if (identical(extension, "pdf")) {
    grDevices::pdf(filename, width = layout$width, height = layout$height)
  } else {
    device <- switch(extension, png = grDevices::png, jpeg = grDevices::jpeg,
      jpg = grDevices::jpeg, tiff = grDevices::tiff, bmp = grDevices::bmp,
      stop("Unsupported heatmap file extension: ", extension))
    device(filename, width = layout$width, height = layout$height, units = "in", res = dpi)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(layout$gtable)
  hm$gtable <- layout$gtable
  hm$width <- layout$width
  hm$height <- layout$height
  hm
}

npx_pathway_mapping <- function(analysis, anno) {
  if (!is.data.frame(anno)) stop("anno must be an annotation data frame.")
  required_anno <- c("OlinkID", "Human_Pathway")
  missing_anno <- setdiff(required_anno, names(anno))
  if (length(missing_anno)) {
    stop("anno is missing required column(s): ", paste(missing_anno, collapse = ", "))
  }
  if (!is.list(analysis) || !is.data.frame(analysis$deps)) {
    stop("analysis must contain a data.frame named 'deps'.")
  }
  if (!"OlinkID" %in% names(analysis$deps) || !"Assay" %in% names(analysis$deps)) {
    stop("analysis$deps must contain 'OlinkID' and 'Assay' columns.")
  }

  deps <- dplyr::select(analysis$deps, OlinkID, Assay) |>
    dplyr::mutate(
      OlinkID = as.character(.data$OlinkID),
      Assay = trimws(as.character(.data$Assay))
    ) |>
    dplyr::filter(!is.na(.data$OlinkID), !is.na(.data$Assay), nzchar(.data$Assay)) |>
    dplyr::distinct()
  pathway_map <- dplyr::select(anno, dplyr::all_of(required_anno)) |>
    dplyr::mutate(
      OlinkID = as.character(.data$OlinkID),
      Human_Pathway = trimws(as.character(.data[["Human_Pathway"]]))
    ) |>
    dplyr::filter(
      !is.na(.data[["OlinkID"]]),
      !is.na(.data[["Human_Pathway"]]),
      nzchar(.data[["Human_Pathway"]])
    ) |>
    dplyr::distinct() |>
    dplyr::inner_join(deps, by = "OlinkID")

  pathway_map
}

#' Save pathway-specific heatmaps for significant DEPs
#' @param data Prepared NPX data used for the differential analysis.
#' @param analysis Result from [npx_analyze()], [npx_ttest()], [npx_wilcox()],
#'   or [npx_anova()].
#' @param anno Annotation data frame containing `OlinkID` and `Human_Pathway`.
#' @param variable Metadata column used for column annotation.
#' @param output_dir Directory for pathway heatmaps.
#' @return A list containing generated files and pathway/assay mappings.
#' @export
npx_pathway_heatmaps <- function(data, analysis, anno, variable = "Condition",
                                 output_dir = "Results/DEPs/Heatmap_Pathway") {
  pathway_map <- npx_pathway_mapping(analysis, anno)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  files <- character()
  mappings <- list()
  for (pathway in unique(pathway_map$Human_Pathway)) {
    assays <- unique(pathway_map$Assay[pathway_map$Human_Pathway == pathway])
    assays <- assays[!is.na(assays) & nzchar(assays)]
    if (!length(assays)) next

    assay_matrix <- tryCatch(npx_matrix(data, assays = assays), error = function(e) NULL)
    if (is.null(assay_matrix)) next
    assays <- intersect(assays, colnames(assay_matrix))
    if (!length(assays)) next

    nr <- length(assays)
    nc <- nrow(assay_matrix)
    cell_h <- if (nr <= 100L) 15 else if (nr <= 500L) 9 else 7
    font_r <- if (nr <= 100L) 8 else if (nr <= 500L) 6 else 5
    safe_name <- npx_safe_filename(pathway)
    path_file <- file.path(output_dir, paste0("Heatmap_DEPs_", safe_name, ".pdf"))
    npx_heatmap(
      data,
      variable = variable,
      assays = assays,
      scale_rows = FALSE,
      filename = path_file,
      show_rownames = TRUE,
      show_colnames = FALSE,
      fontsize_row = font_r,
      cellheight = cell_h,
      width = max(12, 5 + nc * 0.12),
      height = max(8, (nr * cell_h + 350) / 72),
      main = paste0(nr, " DEPs: ", pathway)
    )
    if (file.exists(path_file)) {
      files <- c(files, path_file)
      mappings[[pathway]] <- assays
    }
  }
  list(files = files, mappings = mappings)
}

#' Save a cluster heatmap for an ANOVA DEP result
#' @param data Prepared NPX data.
#' @param analysis An object returned by [npx_analyze()], [npx_anova()], or a
#'   DEP result table with an `Assay` column.
#' @param variable Metadata grouping column used for sample annotations.
#' @param cluster_result Optional result from [npx_cluster_assays()].
#' @param output_dir Directory for PDF, PNG, and HTML outputs.
#' @return A list containing the cluster result, plot matrix, and output files.
#' @export
npx_cluster_heatmap <- function(data, analysis, variable = "Condition",
                                cluster_result = NULL,
                                output_dir = "Results/DEPs/Clustering") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) stop("Package 'pheatmap' is required for npx_cluster_heatmap().")
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  cl <- if (is.null(cluster_result)) npx_cluster_assays(data, analysis) else cluster_result
  if (!is.list(cl) || !is.data.frame(cl$clusters) || !all(c("Assay", "cluster") %in% names(cl$clusters))) {
    stop("cluster_result must be a result from npx_cluster_assays().")
  }

  cluster_table <- cl$clusters |>
    dplyr::transmute(
      Assay = as.character(.data[["Assay"]]),
      cluster = as.character(.data[["cluster"]])
    ) |>
    dplyr::filter(!is.na(.data$Assay), nzchar(.data$Assay), !is.na(.data$cluster), nzchar(.data$cluster)) |>
    dplyr::distinct(.data$Assay, .keep_all = TRUE)
  if (!nrow(cluster_table)) stop("No assay clusters are available for heatmap plotting.")

  cluster_values <- unique(cluster_table$cluster)
  cluster_levels <- if (all(grepl("^[0-9]+$", cluster_values))) {
    as.character(sort(as.integer(cluster_values)))
  } else sort(cluster_values)
  cluster_table$cluster <- factor(cluster_table$cluster, levels = cluster_levels)
  assay_matrix <- cl$matrix[cluster_table$Assay, , drop = FALSE]
  assay_matrix <- assay_matrix[rownames(assay_matrix) %in% cluster_table$Assay, , drop = FALSE]
  if (!nrow(assay_matrix) || !ncol(assay_matrix)) stop("No usable assay matrix is available for heatmap plotting.")

  annotation_col <- data |>
    dplyr::select(dplyr::all_of(c("SampleID", variable))) |>
    dplyr::mutate(
      SampleID = as.character(.data$SampleID),
      .condition = trimws(as.character(.data[[variable]]))
    ) |>
    dplyr::distinct(.data$SampleID, .keep_all = TRUE) |>
    dplyr::filter(!is.na(.data$.condition), nzchar(.data$.condition), .data$SampleID %in% colnames(assay_matrix))
  if (!nrow(annotation_col)) stop("No sample annotations are available for heatmap plotting.")
  condition_order <- npx_condition_levels(annotation_col$.condition)$plot
  annotation_col$.condition <- factor(annotation_col$.condition, levels = condition_order)
  annotation_col <- as.data.frame(annotation_col)
  rownames(annotation_col) <- annotation_col$SampleID
  annotation_col <- data.frame(annotation_col$.condition, row.names = rownames(annotation_col))
  names(annotation_col) <- variable
  sample_order <- rownames(annotation_col)[order(annotation_col[[variable]])]
  assay_matrix <- assay_matrix[, sample_order, drop = FALSE]
  annotation_col <- annotation_col[sample_order, , drop = FALSE]

  annotation_row <- data.frame(cluster = cluster_table$cluster, row.names = cluster_table$Assay)
  annotation_row <- annotation_row[rownames(assay_matrix), , drop = FALSE]
  row_order <- unlist(lapply(cluster_levels, function(clv) rownames(annotation_row)[annotation_row$cluster == clv]), use.names = FALSE)
  assay_matrix <- assay_matrix[row_order, , drop = FALSE]
  annotation_row <- annotation_row[row_order, , drop = FALSE]
  cluster_counts <- table(annotation_row$cluster)
  gaps_row <- if (length(cluster_counts) > 1L) cumsum(cluster_counts)[-length(cluster_counts)] else NULL
  condition_counts <- table(annotation_col[[variable]])
  gaps_col <- if (length(condition_counts) > 1L) cumsum(condition_counts)[-length(condition_counts)] else NULL

  plotly_palette <- c("#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
                      "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52")
  condition_colors <- stats::setNames(rep(plotly_palette, length.out = length(condition_order)), condition_order)
  cluster_colors <- stats::setNames(grDevices::hcl.colors(length(cluster_levels), palette = "Dark 3"), cluster_levels)
  heatmap_colors <- grDevices::colorRampPalette(c("royalblue2", "white", "brown3"))(100)
  heatmap_range <- range(assay_matrix, finite = TRUE)
  if (any(!is.finite(heatmap_range))) stop("No finite NPX values remain for cluster heatmap plotting.")
  if (diff(heatmap_range) == 0) heatmap_range <- heatmap_range + c(-0.5, 0.5)
  heatmap_breaks <- seq(heatmap_range[1], heatmap_range[2], length.out = 101)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_file <- file.path(output_dir, "DEPs_cluster_heatmap.pdf")
  png_file <- file.path(output_dir, "DEPs_cluster_heatmap.png")
  hm <- pheatmap::pheatmap(
      mat = assay_matrix,
      color = heatmap_colors,
      breaks = heatmap_breaks,
      annotation_col = annotation_col,
      annotation_row = annotation_row,
      annotation_colors = c(stats::setNames(list(condition_colors), variable), list(cluster = cluster_colors)),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      treeheight_row = 0,
      treeheight_col = 0,
      gaps_row = gaps_row,
      gaps_col = gaps_col,
      show_rownames = FALSE,
      show_colnames = TRUE,
      annotation_names_row = TRUE,
      annotation_names_col = TRUE,
      angle_col = 90,
      fontsize_col = 6,
      border_color = NA,
      na_col = "grey80",
      main = "DEPs Heatmap by Assay Cluster",
      filename = NA_character_,
      silent = TRUE
    )
  saved_heatmap <- npx_save_heatmap(hm, pdf_file, width = 6, height = 10)
  npx_save_heatmap(hm, png_file, width = 6, height = 10, dpi = 600)

  html_file <- file.path(output_dir, "DEPs_cluster_heatmap.html")
  interactive <- NULL
  files <- c(pdf = pdf_file, png = png_file)
  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    interactive <- plotly::plot_ly(
      x = colnames(assay_matrix),
      y = rownames(assay_matrix),
      z = assay_matrix,
      type = "heatmap",
      colors = c("royalblue2", "white", "brown3"),
      zmin = heatmap_range[1],
      zmax = heatmap_range[2],
      hovertemplate = "Assay: %{y}<br>SampleID: %{x}<br>NPX: %{z:.3f}<extra></extra>"
    ) |>
      plotly::layout(
        title = "DEPs Heatmap by Assay Cluster",
        xaxis = list(tickangle = 90),
        margin = list(l = 120, r = 40, t = 80, b = 120),
        shapes = list(list(type = "rect", xref = "paper", yref = "paper", x0 = 0, x1 = 1, y0 = 0, y1 = 1,
                           fillcolor = "rgba(0,0,0,0)", line = list(color = "black", width = 3)))
      ) |>
      plotly::config(responsive = TRUE, displayModeBar = TRUE)
    html_saved <- tryCatch({
      htmlwidgets::saveWidget(interactive, file = html_file, selfcontained = TRUE)
      TRUE
    }, error = function(e) {
      warning("Self-contained cluster heatmap HTML could not be saved; saving HTML with external dependencies instead: ", conditionMessage(e), call. = FALSE)
      tryCatch({ htmlwidgets::saveWidget(interactive, file = html_file, selfcontained = FALSE); TRUE }, error = function(e2) FALSE)
    })
    if (html_saved) files <- c(files, html = html_file)
  }
  list(cluster = cl, matrix = assay_matrix, annotation_col = annotation_col,
       annotation_row = annotation_row, plot = interactive, files = files,
       width = saved_heatmap$width, height = saved_heatmap$height)
}

#' Save cluster expression-pattern lineplots for an ANOVA DEP result
#' @param data Prepared NPX data.
#' @param analysis An object returned by [npx_analyze()], [npx_anova()], or a
#'   DEP result table with an `Assay` column.
#' @param variable Metadata grouping column.
#' @param cluster_result Optional result from [npx_cluster_assays()].
#' @param output_dir Directory for PDF, PNG, and HTML outputs.
#' @return A list containing the summary data, plot, and output files.
#' @export
npx_cluster_lineplot <- function(data, analysis, variable = "Condition",
                                 cluster_result = NULL,
                                 output_dir = "Results/DEPs/Clustering") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for npx_cluster_lineplot().")
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  cl <- if (is.null(cluster_result)) npx_cluster_assays(data, analysis) else cluster_result
  if (!is.list(cl) || !is.data.frame(cl$clusters) || !all(c("Assay", "cluster") %in% names(cl$clusters))) {
    stop("cluster_result must be a result from npx_cluster_assays().")
  }
  clusters <- cl$clusters |>
    dplyr::transmute(Assay = as.character(.data$Assay), cluster = as.character(.data$cluster)) |>
    dplyr::filter(!is.na(.data$Assay), nzchar(.data$Assay), !is.na(.data$cluster), nzchar(.data$cluster))
  cluster_values <- unique(clusters$cluster)
  cluster_levels <- if (all(grepl("^[0-9]+$", cluster_values))) as.character(sort(as.integer(cluster_values))) else sort(cluster_values)
  clusters$cluster <- factor(clusters$cluster, levels = cluster_levels)
  line_data <- data |>
    dplyr::select(SampleID, Assay, NPX, dplyr::all_of(variable)) |>
    dplyr::mutate(SampleID = as.character(.data$SampleID), Assay = as.character(.data$Assay),
                  .condition = trimws(as.character(.data[[variable]]))) |>
    dplyr::filter(!is.na(.data$.condition), nzchar(.data$.condition)) |>
    dplyr::inner_join(clusters, by = "Assay")
  if (!nrow(line_data)) stop("No data remain for cluster lineplot.")
  condition_order <- npx_condition_levels(line_data$.condition)$plot
  line_summary <- line_data |>
    dplyr::group_by(.data$cluster, .data$.condition) |>
    dplyr::summarise(
      Mean_NPX = mean(.data$NPX, na.rm = TRUE),
      N_NPX = sum(is.finite(.data$NPX)),
      SE_NPX = ifelse(.data$N_NPX > 1, stats::sd(.data$NPX, na.rm = TRUE) / sqrt(.data$N_NPX), 0),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      Condition = factor(.data$.condition, levels = condition_order),
      cluster = factor(as.character(.data$cluster), levels = cluster_levels)
    )
  plotly_palette <- c("#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
                      "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52")
  condition_colors <- stats::setNames(rep(plotly_palette, length.out = length(condition_order)), condition_order)
  p_lineplot <- ggplot2::ggplot(line_summary, ggplot2::aes(x = .data$Condition, y = .data$Mean_NPX, group = 1)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$Mean_NPX - .data$SE_NPX, ymax = .data$Mean_NPX + .data$SE_NPX), fill = "grey70", alpha = 0.35) +
    ggplot2::geom_line(linewidth = 1.2, colour = "black") +
    ggplot2::geom_point(ggplot2::aes(fill = .data$Condition), shape = 21, colour = "black", size = 3.5, stroke = 0.8) +
    ggplot2::facet_wrap(~ cluster, ncol = 1, scales = "free_y", labeller = ggplot2::as_labeller(function(x) paste0("Cluster ", x))) +
    ggplot2::scale_fill_manual(values = condition_colors, drop = FALSE) +
    ggplot2::labs(title = "Expression Pattern by Assay Cluster", x = variable, y = "Mean NPX", fill = variable) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20, face = "bold"),
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(size = 16, face = "bold", colour = "black", hjust = 0),
      axis.line = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 15, face = "bold"),
      axis.text = ggplot2::element_text(size = 12, colour = "black"), panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 2),
      panel.spacing = grid::unit(0.8, "lines"), legend.title = ggplot2::element_text(size = 13, face = "bold"), legend.text = ggplot2::element_text(size = 11)
    )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lineplot_height <- max(8, 3.5 * length(cluster_levels))
  pdf_file <- file.path(output_dir, "DEPs_Cluster_Lineplot.pdf")
  png_file <- file.path(output_dir, "DEPs_Cluster_Lineplot.png")
  html_file <- file.path(output_dir, "DEPs_Cluster_Lineplot.html")
  ggplot2::ggsave(pdf_file, p_lineplot, width = 6, height = lineplot_height, units = "in", device = grDevices::pdf)
  ggplot2::ggsave(png_file, p_lineplot, width = 6, height = lineplot_height, units = "in", dpi = 600, bg = "white")
  files <- c(pdf = pdf_file, png = png_file)
  interactive <- NULL
  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    interactive <- plotly::ggplotly(p_lineplot, tooltip = c("x", "y"), width = 600, height = max(800, 350 * length(cluster_levels)))
    html_saved <- tryCatch({ htmlwidgets::saveWidget(interactive, file = html_file, selfcontained = TRUE); TRUE }, error = function(e) {
      warning("Self-contained cluster lineplot HTML could not be saved; saving HTML with external dependencies instead: ", conditionMessage(e), call. = FALSE)
      tryCatch({ htmlwidgets::saveWidget(interactive, file = html_file, selfcontained = FALSE); TRUE }, error = function(e2) FALSE)
    })
    if (html_saved) files <- c(files, html = html_file)
  }
  list(data = line_summary, plot = p_lineplot, interactive = interactive, files = files,
       width = 6, height = lineplot_height, cluster = cl)
}

#' Plot differential-analysis results as a volcano plot
#' @param analysis An object returned by [npx_analyze()], or its result table.
#' @param logFC Minimum absolute Log2 fold-change for DEPs.
#' @param adj.p Maximum adjusted p-value for DEPs.
#' @param output_dir Optional directory for PDF, PNG, and HTML outputs.
#' @return A list with plotting data, a ggplot object, and optional output files.
#' @export
npx_volcano <- function(analysis, logFC = 0, adj.p = 0.05, output_dir = NULL) {
  tab <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$result else analysis
  if (!is.data.frame(tab)) stop("analysis must be an npx_analyze() result or data frame.")
  npx_validate_dep_cutoffs(logFC, adj.p)
  p_col <- intersect(c("Adjusted_pval", "adjusted_pval", "p.adj"), names(tab))[1]
  if (is.na(p_col)) stop("Analysis table has no adjusted p-value column.")
  tab <- npx_add_effect(tab)
  if (!"Log2FC" %in% names(tab)) stop("Analysis table has no Log2FC or estimate column.")

  effect <- as.numeric(tab$Log2FC)
  adjusted_p <- as.numeric(tab[[p_col]])
  finite_p <- adjusted_p[is.finite(adjusted_p) & adjusted_p > 0]
  p_for_plot <- adjusted_p
  p_for_plot[!is.finite(p_for_plot) | p_for_plot <= 0] <- .Machine$double.xmin
  plot_data <- tab
  plot_data$Log2FC <- effect
  plot_data$neg_log10_p <- -log10(p_for_plot)
  plot_data$Significance <- dplyr::case_when(
    adjusted_p < adj.p & effect > 0 & abs(effect) >= logFC ~ "Up",
    adjusted_p < adj.p & effect < 0 & abs(effect) >= logFC ~ "Down",
    TRUE ~ "NS"
  )
  plot_data <- plot_data[is.finite(plot_data$Log2FC) & is.finite(plot_data$neg_log10_p), , drop = FALSE]

  sig_count <- table(factor(plot_data$Significance, levels = c("Up", "Down", "NS")))
  sig_labels <- c(
    Up = paste0("Up (n=", unname(sig_count[["Up"]]), ")"),
    Down = paste0("Down (n=", unname(sig_count[["Down"]]), ")"),
    NS = paste0("NS (n=", unname(sig_count[["NS"]]), ")")
  )
  plot_data$Significance <- factor(
    unname(sig_labels[plot_data$Significance]),
    levels = unname(sig_labels)
  )

  label_col <- if ("Assay" %in% names(plot_data)) "Assay" else if ("UniProt" %in% names(plot_data)) "UniProt" else NULL
  plot_data$hover_text <- paste0(
    if (is.null(label_col)) "" else paste0(label_col, ": ", plot_data[[label_col]], "<br>"),
    "Log<sub>2</sub> fold change: ", round(plot_data$Log2FC, 3),
    "<br>Adjusted <i>P</i>: ", signif(as.numeric(plot_data[[p_col]]), 3)
  )
  condition_levels <- if (is.list(analysis) && !is.data.frame(analysis) && !is.null(analysis$levels)) {
    intersect(as.character(analysis$levels), names(plot_data))
  } else character()
  if (length(condition_levels)) {
    mean_text <- apply(plot_data[, condition_levels, drop = FALSE], 1, function(x) {
      paste0(condition_levels, ": ", round(as.numeric(x), 3), collapse = "<br>")
    })
    plot_data$hover_text <- paste0(plot_data$hover_text, "<br>", mean_text)
  }

  x_limit <- max(abs(plot_data$Log2FC), na.rm = TRUE)
  if (!is.finite(x_limit) || x_limit <= 0) x_limit <- 1
  x_limit <- x_limit + 0.2
  y_limit <- max(plot_data$neg_log10_p, na.rm = TRUE)
  if (!is.finite(y_limit) || y_limit <= 0) y_limit <- 1
  y_limit <- y_limit + 1
  sig_colors <- c(Up = "brown3", Down = "royalblue2", NS = "#666666")
  names(sig_colors) <- unname(sig_labels)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[["Log2FC"]],
      y = .data[["neg_log10_p"]],
      fill = .data[["Significance"]],
      text = .data[["hover_text"]]
    )
  ) +
    ggplot2::geom_point(shape = 21, size = 3, color = "black", stroke = 0.7, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = sig_colors, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, y_limit)) +
    ggplot2::labs(
      title = "Volcano Plot",
      x = expression(Log[2] ~ "fold change"),
      y = expression(-log[10] ~ "(Adjusted " * italic(P) * ")"),
      fill = "Significance"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 3),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 1.2),
      axis.title = ggplot2::element_text(size = 18, face = "bold"),
      axis.text = ggplot2::element_text(size = 14, color = "black"),
      plot.title = ggplot2::element_text(size = 22, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12)
    )

  files <- NULL
  interactive <- NULL
  if (!is.null(output_dir)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for volcano output.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_file <- file.path(output_dir, "Volcano_plot.pdf")
    png_file <- file.path(output_dir, "Volcano_plot.png")
    html_file <- file.path(output_dir, "Volcano_plot.html")
    ggplot2::ggsave(pdf_file, plot, width = 8, height = 10, units = "in", device = grDevices::pdf)
    ggplot2::ggsave(png_file, plot, width = 8, height = 10, units = "in", dpi = 600, bg = "white")
    files <- c(pdf = pdf_file, png = png_file)
    if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
      interactive <- plotly::ggplotly(plot, tooltip = "text", width = 800, height = 1000) |>
        plotly::layout(
          margin = list(t = 100, r = 40, b = 130, l = 130),
          shapes = list(list(
            type = "rect", xref = "paper", yref = "paper",
            x0 = 0, x1 = 1, y0 = 0, y1 = 1,
            fillcolor = "rgba(0,0,0,0)", line = list(color = "black", width = 3)
          ))
        ) |>
        plotly::config(responsive = TRUE, displayModeBar = TRUE)
      html_saved <- tryCatch({
        htmlwidgets::saveWidget(interactive, file = html_file, selfcontained = TRUE)
        TRUE
      }, error = function(e) {
        warning(
          "Self-contained volcano HTML could not be saved; saving HTML with external dependencies instead: ",
          conditionMessage(e),
          call. = FALSE
        )
        tryCatch({
          htmlwidgets::saveWidget(interactive, file = html_file, selfcontained = FALSE)
          TRUE
        }, error = function(e2) {
          warning("Volcano HTML output was skipped: ", conditionMessage(e2), call. = FALSE)
          FALSE
        })
      })
      if (html_saved) files <- c(files, html = html_file)
    } else {
      warning("Packages 'plotly' and 'htmlwidgets' are required for interactive volcano HTML output; PDF and PNG were saved.", call. = FALSE)
    }
  }
  list(data = plot_data, plot = plot, interactive = interactive, files = files, width = 8, height = 10,
       logFC = logFC, adj.p = adj.p, adjusted_p_column = p_col)
}

#' Plot sample median NPX against sample NPX IQR
#' @param data Prepared assay data.
#' @param variable Metadata column containing sample groups.
#' @param output_dir Optional directory for PDF, PNG, and HTML outputs.
#' @param filename Optional PDF path retained for compatibility. If supplied,
#'   the static plot is also saved to this path.
#' @return A list containing sample statistics, static plot, optional
#'   interactive plot, and output file paths.
#' @export
npx_iqrplot <- function(data, variable = "Condition", output_dir = NULL, filename = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required for npx_iqrplot().")
  if (!variable %in% names(data)) stop("Grouping column not found: ", variable)
  required <- c("SampleID", "NPX", variable)
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("data is missing required columns: ", paste(missing, collapse = ", "))

  meta_condition <- data |>
    dplyr::select(dplyr::all_of(c("SampleID", variable))) |>
    dplyr::mutate(
      SampleID = as.character(.data$SampleID),
      .condition = trimws(as.character(.data[[variable]]))
    ) |>
    dplyr::distinct(.data$SampleID, .keep_all = TRUE)
  condition_order <- npx_condition_levels(meta_condition$.condition)$plot
  sample_stats <- data |>
    dplyr::mutate(SampleID = as.character(.data$SampleID)) |>
    dplyr::group_by(.data$SampleID) |>
    dplyr::summarise(
      SampleMedian = stats::median(.data$NPX, na.rm = TRUE),
      IQR_NPX = stats::IQR(.data$NPX, na.rm = TRUE),
      N_Assay = sum(is.finite(.data$NPX)),
      .groups = "drop"
    ) |>
    dplyr::left_join(meta_condition[, c("SampleID", ".condition"), drop = FALSE], by = "SampleID") |>
    dplyr::filter(
      is.finite(.data$SampleMedian),
      is.finite(.data$IQR_NPX),
      .data$.condition %in% condition_order
    ) |>
    dplyr::mutate(
      Condition = factor(.data$.condition, levels = condition_order),
      hover_text = paste0(
        "SampleID: ", .data$SampleID,
        "<br>", variable, ": ", .data$Condition,
        "<br>Sample median NPX: ", round(.data$SampleMedian, 3),
        "<br>IQR of NPX: ", round(.data$IQR_NPX, 3),
        "<br>Number of assays: ", .data$N_Assay
      )
    )

  plotly_palette <- c(
    "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
    "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52"
  )
  condition_colors <- stats::setNames(
    rep(plotly_palette, length.out = length(condition_order)),
    condition_order
  )
  plot <- ggplot2::ggplot(
    sample_stats,
    ggplot2::aes(
      x = .data[["SampleMedian"]], y = .data[["IQR_NPX"]],
      fill = .data[["Condition"]], text = .data[["hover_text"]]
    )
  ) +
    ggplot2::geom_point(shape = 21, size = 3.5, colour = "black", stroke = 0.8, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = condition_colors, drop = FALSE) +
    ggplot2::labs(
      title = "Sample Median NPX vs IQR of NPX",
      x = "Sample median of NPX", y = "IQR of NPX", fill = variable
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.5),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 1),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text = ggplot2::element_text(size = 13, colour = "black"),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      legend.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 12)
    )

  interactive_plot <- NULL
  files <- list(pdf = NULL, png = NULL, html = NULL)
  if (!is.null(filename)) {
    ggplot2::ggsave(filename, plot, width = 10, height = 8, units = "in")
    files$pdf <- filename
  }
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    files$pdf <- file.path(output_dir, "IQRplot.pdf")
    files$png <- file.path(output_dir, "IQRplot.png")
    ggplot2::ggsave(files$pdf, plot, width = 10, height = 8, units = "in")
    ggplot2::ggsave(files$png, plot, width = 10, height = 8, units = "in", dpi = 600, bg = "white")
  }
  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    interactive_plot <- plotly::ggplotly(plot, tooltip = "text")
    if (!is.null(output_dir)) {
      files$html <- file.path(output_dir, "IQRplot.html")
      html_saved <- tryCatch({
        htmlwidgets::saveWidget(interactive_plot, file = files$html, selfcontained = TRUE)
        TRUE
      }, error = function(e) {
        warning("Self-contained IQR HTML could not be saved; saving HTML with external dependencies instead: ", conditionMessage(e), call. = FALSE)
        tryCatch({
          htmlwidgets::saveWidget(interactive_plot, file = files$html, selfcontained = FALSE)
          TRUE
        }, error = function(e2) {
          warning("IQR HTML output was skipped: ", conditionMessage(e2), call. = FALSE)
          FALSE
        })
      })
      if (!html_saved) files$html <- NULL
    }
  }
  list(data = sample_stats, plot = plot, interactive = interactive_plot, files = files,
       width = 10, height = 8)
}
