npx_go_databases <- c(
  Biological_Process = "GO_Biological_Process_2026",
  Molecular_Function = "GO_Molecular_Function_2026",
  Cellular_Component = "GO_Cellular_Component_2026",
  KEGG_pathway = "KEGG_2026"
)

npx_go_invalid_response <- function(...) {
  stop(structure(list(message = paste0(...), call = NULL),
                 class = c("npx_go_invalid_response", "error", "condition")))
}

#' Run GO/KEGG enrichment through Enrichr
#' @param genes Character vector of gene symbols, or a data frame containing
#'   `Assay` and/or `UniProt`. `Assay` is preferred and `UniProt` is the fallback.
#' @param databases Enrichr library names.
#' @param analysis Optional analysis result. Its filtered `deps` table is used.
#' @param logFC Minimum absolute Log2 fold-change for filtering a data frame.
#' @param adj.p Maximum adjusted p-value for filtering a data frame.
#' @return A named list of complete Enrichr result tables.
#' @export
npx_go_enrichment <- function(genes = NULL, databases = unname(npx_go_databases), analysis = NULL, logFC = 0, adj.p = 0.05) {
  if (!requireNamespace("enrichR", quietly = TRUE)) stop("Package 'enrichR' is required for npx_go_enrichment().")
  if (!is.null(analysis)) genes <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$deps else analysis
  if (is.data.frame(genes)) genes <- npx_resolve_gene_ids(genes, logFC = logFC, adj.p = adj.p)
  genes <- unique(trimws(as.character(genes))); genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) stop("genes must contain at least one non-empty gene symbol.")
  if (!length(databases) || any(is.na(databases)) || any(!nzchar(databases))) {
    stop("databases must contain at least one non-empty Enrichr database name.")
  }
  # enrichR 3.4 initializes these in .onAttach(), which :: does not run.
  # Supply missing defaults locally and preserve the user's site/proxy settings.
  defaults <- list(
    enrichR.base.address = paste0(getOption("enrichR.sites.base.address", "https://maayanlab.cloud/"), "Enrichr/"),
    enrichR.live = TRUE, enrichR.quiet = FALSE
  )
  option_names <- c(names(defaults), "stringsAsFactors")
  old_options <- stats::setNames(lapply(option_names, getOption), option_names)
  on.exit(options(old_options), add = TRUE)
  options(defaults[vapply(old_options[names(defaults)], is.null, logical(1L))])
  message(
    "[Enrichr] submitting ", length(genes), " Assay/gene identifiers to: ",
    paste(databases, collapse = ", ")
  )
  result <- tryCatch(
    enrichR::enrichr(genes, databases),
    error = function(e) {
      npx_rethrow_time_limit(e)
      stop(
        "Enrichr request failed for ", length(genes),
        " Assay/gene identifiers and database(s) [",
        paste(databases, collapse = ", "), "]: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (!is.list(result)) {
    npx_go_invalid_response(
      "Enrichr returned an invalid result for database(s) [",
      paste(databases, collapse = ", "), "]. Expected a named list of data frames."
    )
  }
  missing <- setdiff(databases, names(result))
  if (length(missing)) {
    npx_go_invalid_response(
      "Enrichr response is missing database(s): ", paste(missing, collapse = ", "),
      ". Check library names with enrichR::listEnrichrDbs()."
    )
  }
  required <- c("Term", "P.value", "Adjusted.P.value", "Combined.Score", "Genes")
  for (database in databases) {
    table <- result[[database]]
    if (!is.data.frame(table) || !all(required %in% names(table))) {
      npx_go_invalid_response("Enrichr returned invalid columns for database '", database,
           "'. Expected: ", paste(required, collapse = ", "))
    }
  }
  result[databases]
}

npx_resolve_gene_ids <- function(data, logFC = 0, adj.p = 0.05) {
  if (!is.data.frame(data)) stop("data must be a data frame.")
  if (all(c("Adjusted_pval", "Log2FC") %in% names(data))) data <- npx_filter_deps(data, logFC = logFC, adj.p = adj.p)
  else if ("Adjusted_pval" %in% names(data)) data <- npx_filter_deps(data, logFC = 0, adj.p = adj.p)
  id_col <- if ("Assay" %in% names(data)) "Assay" else if ("UniProt" %in% names(data)) "UniProt" else NULL
  if (is.null(id_col)) stop("No gene identifier found. Expected 'Assay' or 'UniProt'.")
  data[[id_col]]
}

npx_enrichr_table <- function(result, genes) {
  if (is.null(result) || !is.data.frame(result) || !nrow(result)) return(NULL)
  if (!"Genes" %in% names(result)) result$Genes <- ""
  if (!"Adjusted.P.value" %in% names(result)) result$Adjusted.P.value <- NA_real_
  if (!"P.value" %in% names(result)) result$P.value <- NA_real_
  if (!"Combined.Score" %in% names(result)) result$Combined.Score <- NA_real_
  result |>
    dplyr::mutate(
      Count = lengths(strsplit(as.character(.data$Genes), ";")),
      GeneRatio = .data$Count / length(genes),
      minusLog10Padj = -log10(pmax(as.numeric(.data$Adjusted.P.value), .Machine$double.xmin))
    ) |>
    dplyr::arrange(.data$Adjusted.P.value, dplyr::desc(.data$Combined.Score))
}

npx_go_dotplot <- function(result, title, top_n = 10, p_cutoff = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || is.null(result) || !nrow(result)) return(NULL)
  plot_data <- result |>
    dplyr::filter(!is.na(.data$Adjusted.P.value), .data$Adjusted.P.value <= p_cutoff,
                  is.finite(.data$Combined.Score), is.finite(.data$GeneRatio)) |>
    dplyr::arrange(dplyr::desc(.data$Combined.Score)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::arrange(dplyr::desc(.data$GeneRatio))
  if (!nrow(plot_data)) return(NULL)
  plot_data$Term <- vapply(
    as.character(plot_data$Term),
    function(term) paste(strwrap(term, width = 55), collapse = "\n"),
    character(1L)
  )
  plot_data$Term <- factor(plot_data$Term, levels = rev(unique(plot_data$Term)))
  ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$GeneRatio, y = .data$Term, size = .data$Count, color = .data$minusLog10Padj,
    text = paste0("Term: ", .data$Term, "<br>GeneRatio: ", round(.data$GeneRatio, 3),
                  "<br>Count: ", .data$Count, "<br>Adjusted P: ", signif(.data$Adjusted.P.value, 3))
  )) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.15))) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::scale_size(range = c(3, 10)) +
    ggplot2::labs(title = title, x = "GeneRatio", y = NULL, size = "Count", color = expression(-log[10]("Adjusted P"))) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black", lineheight = 0.9),
      axis.text.x = ggplot2::element_text(size = 10, colour = "black", angle = 45, hjust = 1, vjust = 1),
      axis.title.x = ggplot2::element_text(size = 13, face = "bold"), panel.grid.minor = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 11), legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(t = 15, r = 25, b = 20, l = 15)
    )
}

npx_go_save_plot <- function(plot, file_prefix, output_dir, width = 9, height = 6) {
  if (is.null(plot)) return(character())
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character(); errors <- character()
  for (format in c("pdf", "png", "html")) {
    file <- file.path(output_dir, paste0(file_prefix, ".", format))
    saved <- tryCatch({
      if (format == "pdf") ggplot2::ggsave(file, plot, width = width, height = height, units = "in", device = grDevices::pdf, limitsize = FALSE)
      else if (format == "png") ggplot2::ggsave(file, plot, width = width, height = height, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
      else {
        if (!requireNamespace("plotly", quietly = TRUE) || !requireNamespace("htmlwidgets", quietly = TRUE)) stop("Packages 'plotly' and 'htmlwidgets' are required for GO HTML output.")
        widget <- plotly::ggplotly(plot, tooltip = "text", width = width * 100, height = height * 100)
        tryCatch(htmlwidgets::saveWidget(widget, file, selfcontained = TRUE), error = function(e) {
          npx_rethrow_time_limit(e)
          warning("Self-contained GO HTML could not be saved; saving HTML with external dependencies instead: ", conditionMessage(e), call. = FALSE)
          htmlwidgets::saveWidget(widget, file, selfcontained = FALSE)
        })
      }
      TRUE
    }, error = function(e) { npx_rethrow_time_limit(e); e })
    if (inherits(saved, "error")) errors <- c(errors, paste0(format, ": ", conditionMessage(saved)))
    else files <- c(files, stats::setNames(file, format))
  }
  attr(files, "errors") <- errors
  files
}

npx_go_write_workbook <- function(results, file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) stop("Package 'openxlsx' is required for GO Excel output.")
  workbook <- openxlsx::createWorkbook(); if (!length(results)) results <- list(Result = NULL)
  for (nm in names(results)) {
    sheet <- substr(gsub("[^A-Za-z0-9_]", "_", nm), 1, 31); if (!nzchar(sheet)) sheet <- "Result"
    openxlsx::addWorksheet(workbook, sheet)
    result <- results[[nm]]
    if (is.null(result) || !is.data.frame(result) || !nrow(result)) {
      openxlsx::writeData(workbook, sheet, data.frame(Message = "No enrichment terms returned")); next
    }
    result <- result |>
      dplyr::select(-dplyr::any_of(c("Old.P.value", "Old.Adjusted.P.value"))) |>
      dplyr::relocate(dplyr::any_of(c("Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Count", "GeneRatio", "minusLog10Padj", "Genes")), .before = 1)
    openxlsx::writeData(workbook, sheet, result); openxlsx::addFilter(workbook, sheet, rows = 1, cols = seq_len(ncol(result))); openxlsx::freezePane(workbook, sheet, firstRow = TRUE)
    openxlsx::setColWidths(workbook, sheet, cols = seq_len(ncol(result)), widths = "auto")
    if ("Term" %in% names(result)) openxlsx::setColWidths(workbook, sheet, which(names(result) == "Term"), widths = 45)
    if ("Genes" %in% names(result)) openxlsx::setColWidths(workbook, sheet, which(names(result) == "Genes"), widths = 60)
  }
  openxlsx::saveWorkbook(workbook, file, overwrite = TRUE); file
}

npx_go_diagnostics <- function(label, database, name, genes, status = "Success", code = "Completed", reason = "", stage = "Query", terms = NA_integer_, significant = NA_integer_) {
  data.frame(Set = label, Database = database, Label = name, Stage = stage,
             Status = status, Code = code, GeneCount = length(genes),
             TermCount = terms, SignificantTermCount = significant,
             Reason = reason, stringsAsFactors = FALSE)
}

npx_go_status <- function(diagnostics) {
  statuses <- diagnostics$Status
  if (!length(statuses) || all(statuses == "Skipped")) return("Skipped")
  statuses <- statuses[statuses != "Skipped"]
  if (all(statuses == "Failed")) return("Failed")
  if (any(statuses %in% c("Failed", "Partial"))) return("Partial")
  "Success"
}

npx_go_analyze_set <- function(dep_data, label, output_dir, databases = npx_go_databases, padj_cutoff = 0.05, top_n = 10, diagnostics_dir = output_dir) {
  if (!is.data.frame(dep_data) || !"Assay" %in% names(dep_data)) stop("GO data must contain an Assay column.")
  if (!length(databases) || anyNA(databases) || any(!nzchar(databases))) stop("databases must contain non-empty Enrichr database names.")
  genes <- unique(trimws(as.character(dep_data$Assay))); genes <- genes[!is.na(genes) & nzchar(genes)]
  labels <- names(databases)
  if (is.null(labels)) labels <- rep("", length(databases))
  unnamed <- is.na(labels) | !nzchar(labels)
  labels[unnamed] <- databases[unnamed]
  names(databases) <- make.unique(labels, sep = "_")
  safe_label <- npx_safe_filename(label)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character(); results <- list(); diagnostics <- list()
  for (nm in names(databases)) {
    database <- databases[[nm]]
    diagnostic <- npx_go_diagnostics(label, database, nm, genes)
    if (!length(genes)) {
      diagnostic$Status <- "Skipped"; diagnostic$Code <- "InputEmpty"
      diagnostic$Stage <- "Input"; diagnostic$Reason <- "No non-empty Assay identifiers in this DEP set."
      diagnostics[[nm]] <- diagnostic
      next
    }
    message("[GeneOntology] ", label, " / ", database, ": querying ", length(genes), " Assay identifiers.")
    result <- tryCatch({
      raw <- npx_go_enrichment(genes = genes, databases = database)
      table <- raw[[database]]
      required <- c("Term", "P.value", "Adjusted.P.value", "Combined.Score", "Genes")
      if (!is.data.frame(table) || !all(required %in% names(table))) npx_go_invalid_response("Invalid Enrichr table for ", database)
      if (nrow(table) && (!is.numeric(table$Adjusted.P.value) || !is.numeric(table$Combined.Score))) npx_go_invalid_response("Enrichr adjusted p-values and Combined.Score must be numeric for ", database)
      table
    }, error = function(e) { npx_rethrow_time_limit(e); e })
    if (inherits(result, "error")) {
      diagnostic$Status <- "Failed"
      diagnostic$Code <- if (inherits(result, "npx_go_invalid_response")) "InvalidResponse" else "QueryError"
      diagnostic$Reason <- conditionMessage(result)
      warning("GO ", label, " / ", database, " [", diagnostic$Code, "]: ", diagnostic$Reason, call. = FALSE)
      diagnostics[[nm]] <- diagnostic
      next
    }
    diagnostic$TermCount <- nrow(result)
    diagnostic$SignificantTermCount <- sum(is.finite(result$Adjusted.P.value) & result$Adjusted.P.value <= padj_cutoff)
    results[nm] <- list(npx_enrichr_table(result, genes))
    if (!nrow(result)) {
      diagnostic$Code <- "EmptyResult"; diagnostic$Reason <- "Database returned no enrichment terms."
    } else if (!diagnostic$SignificantTermCount) {
      diagnostic$Code <- "NoSignificantTerms"
      diagnostic$Reason <- paste0("No terms pass adjusted p-value <= ", padj_cutoff, "; full results retained.")
    } else {
      diagnostic$Stage <- "Plot"
      plot <- tryCatch(npx_go_dotplot(results[[nm]], paste(label, gsub("_", " ", nm)), top_n = top_n, p_cutoff = padj_cutoff),
                       error = function(e) { npx_rethrow_time_limit(e); e })
      if (inherits(plot, "error") || is.null(plot)) {
        diagnostic$Status <- "Partial"; diagnostic$Code <- "PlotError"
        diagnostic$Reason <- if (inherits(plot, "error")) conditionMessage(plot) else "Significant terms exist but no plot could be constructed."
      } else {
        diagnostic$Stage <- "Export"
        saved <- tryCatch(npx_go_save_plot(plot, paste0(safe_label, "_GO_", npx_safe_filename(nm)), output_dir),
                          error = function(e) { npx_rethrow_time_limit(e); e })
        if (inherits(saved, "error")) {
          diagnostic$Status <- "Partial"; diagnostic$Code <- "OutputError"; diagnostic$Reason <- conditionMessage(saved)
        } else {
          files <- c(files, saved)
          if (length(attr(saved, "errors"))) {
            diagnostic$Status <- "Partial"; diagnostic$Code <- "OutputError"
            diagnostic$Reason <- paste(attr(saved, "errors"), collapse = "; ")
          }
        }
      }
    }
    message("[GeneOntology] ", label, " / ", database, ": ", diagnostic$Code, " (", diagnostic$TermCount, " terms; ", diagnostic$SignificantTermCount, " significant).")
    diagnostics[[nm]] <- diagnostic
  }
  diagnostics <- do.call(rbind, diagnostics); rownames(diagnostics) <- NULL
  if (length(results)) {
    excel_file <- file.path(output_dir, paste0("GeneOntology_", safe_label, ".xlsx"))
    saved <- tryCatch(npx_go_write_workbook(results, excel_file), error = function(e) { npx_rethrow_time_limit(e); e })
    if (inherits(saved, "error")) {
      diagnostics <- rbind(diagnostics, npx_go_diagnostics(label, NA_character_, "Workbook", genes,
        "Partial", "OutputError", conditionMessage(saved), "Export"))
    } else files <- c(files, excel = saved)
  }
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics_file <- file.path(diagnostics_dir, paste0("GeneOntology_", safe_label, "_diagnostics.csv"))
  saved <- tryCatch({ utils::write.csv(diagnostics, diagnostics_file, row.names = FALSE); TRUE },
                    error = function(e) { npx_rethrow_time_limit(e); e })
  if (inherits(saved, "error")) {
    diagnostics <- rbind(diagnostics, npx_go_diagnostics(label, NA_character_, "Diagnostics", genes,
      "Partial", "OutputError", conditionMessage(saved), "Export"))
    warning("GO diagnostics output failed: ", conditionMessage(saved), call. = FALSE)
  } else files <- c(files, diagnostics = diagnostics_file)
  list(status = npx_go_status(diagnostics), reason = paste(unique(diagnostics$Reason[nzchar(diagnostics$Reason)]), collapse = "; "),
       files = files, results = results, genes = genes, diagnostics = diagnostics,
       diagnostics_file = if (!inherits(saved, "error")) diagnostics_file else NULL)
}

npx_go_run_set <- function(dep_data, label, output_dir, ...) {
  tryCatch(npx_go_analyze_set(dep_data, label, output_dir, ...), error = function(e) {
    npx_rethrow_time_limit(e)
    warning("GO set failed: ", conditionMessage(e), call. = FALSE)
    list(status = "Failed", reason = conditionMessage(e), files = character(), results = list(),
         diagnostics = npx_go_diagnostics(label, NA_character_, NA_character_, character(),
           "Failed", "InputError", conditionMessage(e), "Input"))
  })
}

#' Run functional analysis for two-group or ANOVA DEP results
#' @param analysis Result from [npx_analyze()], [npx_ttest()], [npx_wilcox()], or [npx_anova()].
#' @param data Prepared NPX data, required for automatic ANOVA clustering.
#' @param cluster_result Optional result from [npx_cluster_assays()].
#' @param output_dir Pipeline root directory. GO outputs are saved under
#'   `Results/GeneOntology` and diagnostics under `log`.
#' @param databases Enrichr databases.
#' @param padj_cutoff Adjusted p-value cutoff used for plot terms.
#' @param top_n Number of terms shown in each dotplot.
#' @return A list containing GO results and output files.
#' @export
npx_functional_analysis <- function(analysis, data = NULL, cluster_result = NULL, output_dir = ".", databases = npx_go_databases, padj_cutoff = 0.05, top_n = 10) {
  if (!is.list(analysis) || !is.data.frame(analysis$deps)) stop("analysis must contain a data.frame named 'deps'.")
  log_dir <- file.path(output_dir, "log"); dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  go_dir <- file.path(output_dir, "Results", "GeneOntology"); dir.create(go_dir, recursive = TRUE, showWarnings = FALSE); files <- character(); results <- list()
  if (length(analysis$levels) == 2L) {
    up <- analysis$deps[is.finite(as.numeric(analysis$deps$Log2FC)) & analysis$deps$Log2FC > 0, , drop = FALSE]
    down <- analysis$deps[is.finite(as.numeric(analysis$deps$Log2FC)) & analysis$deps$Log2FC < 0, , drop = FALSE]
    for (direction in c(UP = "up_regulated", DOWN = "down_regulated")) {
      dep_data <- if (direction == "up_regulated") up else down
      set_result <- npx_go_run_set(dep_data, direction, go_dir, databases, padj_cutoff, top_n, diagnostics_dir = log_dir)
      if (!is.null(set_result)) { files <- c(files, set_result$files); results[[direction]] <- set_result }
    }
  } else if (identical(analysis$method, "anova")) {
    cluster_table <- tryCatch({
      if (is.null(cluster_result)) {
        if (is.null(data)) stop("data is required to calculate ANOVA cluster-specific GO analysis.")
        cluster_result <- npx_cluster_assays(data, analysis)
      }
      npx_cluster_table(analysis$deps, cluster_result)
    }, error = function(e) { npx_rethrow_time_limit(e); e })
    if (inherits(cluster_table, "error")) {
      status <- if (inherits(cluster_table, "npx_insufficient_input") || !nrow(analysis$deps)) "Skipped" else "Failed"
      code <- if (!nrow(analysis$deps)) "InputEmpty" else if (inherits(cluster_table, "npx_insufficient_input")) "InsufficientInput" else "ClusteringError"
      results$Clustering <- list(status = status, reason = conditionMessage(cluster_table), files = character(), results = list(),
        diagnostics = npx_go_diagnostics("Clustering", NA_character_, NA_character_, character(), status, code, conditionMessage(cluster_table), "Clustering"))
      cluster_table <- data.frame(Assay = character(), ClusterNo = character())
    }
    unassigned <- attr(cluster_table, "unassigned_assays")
    if (length(unassigned)) results$Unassigned <- list(status = "Skipped", reason = "DEP assays lack a cluster assignment.", files = character(), results = list(),
      diagnostics = npx_go_diagnostics("Unassigned", NA_character_, NA_character_, unassigned, "Skipped", "UnassignedCluster", paste(unassigned, collapse = "; "), "Clustering"))
    for (cluster_id in unique(cluster_table$ClusterNo)) {
      cluster_data <- dplyr::inner_join(analysis$deps, cluster_table[cluster_table$ClusterNo == cluster_id, , drop = FALSE], by = "Assay")
      cluster_dir <- file.path(go_dir, paste0("Cluster_", npx_safe_filename(cluster_id)))
      set_result <- npx_go_run_set(cluster_data, paste0("Cluster_", cluster_id), cluster_dir, databases, padj_cutoff, top_n, diagnostics_dir = log_dir)
      if (!is.null(set_result)) { files <- c(files, set_result$files); results[[paste0("Cluster_", cluster_id)]] <- set_result }
    }
  } else stop("Functional GO analysis requires a two-group test or ANOVA result.")
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics")); rownames(diagnostics) <- NULL
  diagnostic_file <- file.path(log_dir, "GO_diagnostics.csv")
  saved <- tryCatch({ utils::write.csv(diagnostics, diagnostic_file, row.names = FALSE); TRUE },
                    error = function(e) { npx_rethrow_time_limit(e); e })
  if (inherits(saved, "error")) {
    diagnostics <- rbind(diagnostics, npx_go_diagnostics("GO", NA_character_, "Diagnostics", character(),
      "Partial", "OutputError", conditionMessage(saved), "Export"))
    warning("GO diagnostics output failed: ", conditionMessage(saved), call. = FALSE)
  } else files <- c(files, diagnostics = diagnostic_file)
  list(status = npx_go_status(diagnostics), reason = paste(unique(diagnostics$Reason[nzchar(diagnostics$Reason)]), collapse = "; "),
       files = files, results = results, cluster = cluster_result, diagnostics = diagnostics,
       diagnostics_file = if (!inherits(saved, "error")) diagnostic_file else NULL)
}

#' Query STRING protein-protein interactions
#' @param genes Gene symbols or identifiers accepted by STRINGdb.
#' @param species NCBI taxonomy identifier; defaults to human.
#' @param score_threshold Minimum combined STRING score (0-1000).
#' @param analysis Optional analysis result whose filtered `deps` table supplies identifiers.
#' @param logFC Minimum absolute Log2 fold-change for filtering.
#' @param adj.p Maximum adjusted p-value for filtering.
#' @return A list containing mapped proteins and interaction edges.
#' @export
npx_string_ppi <- function(genes = NULL, species = 9606, score_threshold = 400, analysis = NULL, logFC = 0, adj.p = 0.05) {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) stop("Package 'STRINGdb' is required for npx_string_ppi().")
  if (!is.numeric(score_threshold) || length(score_threshold) != 1L || !is.finite(score_threshold) || score_threshold < 0 || score_threshold > 1000) stop("score_threshold must be between 0 and 1000.")
  if (!is.null(analysis)) genes <- if (is.list(analysis) && !is.data.frame(analysis)) analysis$deps else analysis
  if (is.data.frame(genes)) genes <- npx_resolve_gene_ids(genes, logFC = logFC, adj.p = adj.p)
  genes <- unique(as.character(genes)); genes <- genes[!is.na(genes) & nzchar(genes)]
  if (!length(genes)) stop("genes must contain at least one non-empty identifier.")
  db <- STRINGdb::STRINGdb$new(version = "12.0", species = as.integer(species), score_threshold = as.integer(score_threshold))
  mapped <- db$map(data.frame(gene = genes, stringsAsFactors = FALSE), "gene", removeUnmappedRows = TRUE)
  if (!nrow(mapped)) return(list(mapped = mapped, interactions = data.frame()))
  interactions <- db$get_interactions(mapped$STRING_id); interactions <- interactions[interactions$combined_score >= score_threshold, , drop = FALSE]
  list(mapped = mapped, interactions = interactions)
}

npx_string_aliases <- function(string_db) {
  aliases <- as.data.frame(npx_ppi_query(string_db$get_aliases()), stringsAsFactors = FALSE)
  string_id_col <- grep("STRING", names(aliases), ignore.case = TRUE, value = TRUE)[1]
  alias_col <- grep("alias", names(aliases), ignore.case = TRUE, value = TRUE)[1]
  source_col <- grep("source", names(aliases), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(string_id_col) || is.na(alias_col)) {
    npx_ppi_stop("STRING alias table must contain STRING_id and alias columns.", "ValidationError")
  }
  data.frame(
    STRING_id = as.character(aliases[[string_id_col]]),
    alias = as.character(aliases[[alias_col]]),
    source = if (!is.na(source_col)) as.character(aliases[[source_col]]) else NA_character_,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(!is.na(.data$STRING_id), nzchar(.data$STRING_id), !is.na(.data$alias), nzchar(.data$alias))
}

npx_string_mapping <- function(dep_data, string_db) {
  input <- dep_data |>
    dplyr::transmute(
      Assay = as.character(.data$Assay),
      UniProt = as.character(.data$UniProt)
    ) |>
    dplyr::filter(!is.na(.data$Assay), nzchar(.data$Assay)) |>
    dplyr::distinct(.data$Assay, .keep_all = TRUE)
  aliases <- npx_string_aliases(string_db)
  gene_mapping <- aliases |>
    dplyr::filter(.data$alias %in% input$Assay) |>
    dplyr::mutate(mapping_priority = dplyr::case_when(
      grepl("HGNC", .data$source, ignore.case = TRUE) ~ 1,
      grepl("symbol", .data$source, ignore.case = TRUE) ~ 2,
      TRUE ~ 10
    )) |>
    dplyr::arrange(.data$alias, .data$mapping_priority) |>
    dplyr::group_by(.data$alias) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::transmute(Assay = .data$alias, STRING_id_gene = .data$STRING_id, STRING_source_gene = .data$source)
  uniprot_mapping <- input |>
    dplyr::left_join(gene_mapping, by = "Assay") |>
    dplyr::filter(is.na(.data$STRING_id_gene), !is.na(.data$UniProt), nzchar(.data$UniProt)) |>
    dplyr::select("UniProt") |>
    dplyr::distinct() |>
    dplyr::inner_join(
      aliases |>
        dplyr::filter(.data$alias %in% input$UniProt) |>
        dplyr::mutate(mapping_priority = dplyr::case_when(
          grepl("UniProt", .data$source, ignore.case = TRUE) ~ 1,
          grepl("Swiss", .data$source, ignore.case = TRUE) ~ 2,
          TRUE ~ 10
        )) |>
        dplyr::arrange(.data$alias, .data$mapping_priority) |>
        dplyr::group_by(.data$alias) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup(),
      by = c("UniProt" = "alias")
    ) |>
    dplyr::transmute(UniProt = .data$UniProt, STRING_id_uniprot = .data$STRING_id, STRING_source_uniprot = .data$source)
  input |>
    dplyr::left_join(gene_mapping, by = "Assay") |>
    dplyr::left_join(uniprot_mapping, by = "UniProt") |>
    dplyr::mutate(
      STRING_id = dplyr::coalesce(.data$STRING_id_gene, .data$STRING_id_uniprot),
      STRING_mapping_source = dplyr::case_when(
        !is.na(.data$STRING_id_gene) ~ "Gene_symbol",
        !is.na(.data$STRING_id_uniprot) ~ "UniProt",
        TRUE ~ "Unmapped"
      ),
      STRING_alias_source = dplyr::coalesce(.data$STRING_source_gene, .data$STRING_source_uniprot)
    ) |>
    dplyr::distinct(.data$Assay, .keep_all = TRUE)
}

# Preserve the failing stage without losing the original service message.
npx_ppi_stop <- function(message, category) {
  stop(structure(list(message = message, call = NULL, category = category),
    class = c("npx_ppi_error", "error", "condition")))
}

npx_ppi_query <- function(expr) {
  tryCatch(expr, error = function(e) {
    npx_rethrow_time_limit(e)
    npx_ppi_stop(conditionMessage(e), "QueryError")
  })
}

npx_ppi_node_colors <- function(values, color_mode, max_value = NULL, cluster_colors = NULL) {
  if (identical(color_mode, "cluster")) {
    colors <- unname(cluster_colors[as.character(values)])
    colors[is.na(colors)] <- "#BDBDBD"
    return(colors)
  }
  values <- as.numeric(values)
  scale <- if (identical(color_mode, "log2fc")) {
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_value, max_value))
  } else {
    ggplot2::scale_color_gradientn(colors = c("#2C7BB6", "#FFFFBF", "#D7191C"), limits = c(0, max_value))
  }
  scale$map(values)
}

npx_ppi_graph <- function(dep_data, string_db, score_threshold, min_proteins, color_column, mapping = NULL) {
  if (!all(c("Assay", "UniProt", color_column) %in% names(dep_data))) {
    npx_ppi_stop("PPI data must contain Assay, UniProt, and the selected node color column.", "ValidationError")
  }
  if (!nrow(dep_data)) return(list(status = "Skipped", category = "InputEmpty", reason = "No input proteins", mapped = data.frame()))
  if (length(unique(dep_data$Assay[!is.na(dep_data$Assay) & nzchar(dep_data$Assay)])) < min_proteins) return(list(status = "Skipped", category = "InsufficientInput", reason = "Too few input proteins", mapped = data.frame()))
  if (is.null(mapping)) mapping <- npx_string_mapping(dep_data, string_db)
  color_values <- if (identical(color_column, "ClusterNo")) as.character else as.numeric
  mapped <- dep_data |>
    dplyr::mutate(Assay = as.character(.data$Assay), UniProt = as.character(.data$UniProt)) |>
    dplyr::left_join(mapping, by = c("Assay", "UniProt")) |>
    dplyr::filter(!is.na(.data$STRING_id), nzchar(.data$STRING_id)) |>
    dplyr::mutate(color_value = color_values(.data[[color_column]])) |>
    dplyr::distinct(.data$Assay, .keep_all = TRUE)
  if (length(unique(mapped$STRING_id)) < min_proteins) return(list(status = "Skipped", category = "InsufficientMapped", reason = "Too few mapped proteins", mapped = mapped))
  interactions <- npx_ppi_query(string_db$get_interactions(unique(mapped$STRING_id)))
  if (is.null(interactions) || !nrow(interactions)) return(list(status = "Skipped", category = "NoInteractions", reason = "No STRING interactions", mapped = mapped))
  if (!all(c("from", "to", "combined_score") %in% names(interactions))) npx_ppi_stop("STRING interactions must contain from, to, and combined_score columns.", "ValidationError")
  interactions <- interactions |>
    dplyr::filter(.data$combined_score >= score_threshold) |>
    dplyr::rename(from_STRING = "from", to_STRING = "to")
  if (!nrow(interactions)) return(list(status = "Skipped", category = "NoHighConfidenceInteractions", reason = "No high-confidence interactions", mapped = mapped))
  id_map <- mapped |>
    dplyr::select("STRING_id", "Assay", "color_value") |>
    dplyr::distinct(.data$STRING_id, .keep_all = TRUE)
  edges <- interactions |>
    dplyr::left_join(id_map |>
      dplyr::select("STRING_id", from = "Assay"), by = c("from_STRING" = "STRING_id")) |>
    dplyr::left_join(id_map |>
      dplyr::select("STRING_id", to = "Assay"), by = c("to_STRING" = "STRING_id")) |>
    dplyr::filter(!is.na(.data$from), !is.na(.data$to), .data$from != .data$to) |>
    dplyr::mutate(gene1 = pmin(.data$from, .data$to), gene2 = pmax(.data$from, .data$to)) |>
    dplyr::group_by(.data$gene1, .data$gene2) |>
    dplyr::slice_max(.data$combined_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(from = .data$gene1, to = .data$gene2, combined_score = .data$combined_score)
  if (!nrow(edges)) return(list(status = "Skipped", category = "NoInteractions", reason = "No mapped PPI edges", mapped = mapped))
  connected <- unique(c(edges$from, edges$to))
  nodes <- mapped |>
    dplyr::filter(.data$Assay %in% connected) |>
    dplyr::transmute(name = .data$Assay, color_value = .data$color_value) |>
    dplyr::distinct(.data$name, .keep_all = TRUE)
  graph <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(combined_score = "max", "ignore"))
  igraph::V(graph)$degree <- igraph::degree(graph)
  igraph::V(graph)$betweenness <- igraph::betweenness(graph, normalized = TRUE)
  igraph::V(graph)$harmonic <- igraph::harmonic_centrality(graph, normalized = TRUE)
  list(status = "Success", category = "Success", reason = "", graph = graph, edges = edges, nodes = nodes, mapped = mapped)
}

npx_ppi_save_network <- function(network, label, output_dir, color_mode, score_threshold, width, height, dpi, seed) {
  if (!identical(network$status, "Success")) return(network)
  files <- character()
  output_errors <- character()
  output_status <- c(pdf = NA, png = NA, html = NA, edges = NA, nodes = NA)
  save_output <- function(format, path, expr) {
    tryCatch({
      force(expr)
      files[format] <<- path
      output_status[format] <<- TRUE
      TRUE
    }, error = function(e) {
      npx_rethrow_time_limit(e)
      output_status[format] <<- FALSE
      output_errors <<- c(output_errors, paste0(format, ": ", conditionMessage(e)))
      FALSE
    })
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  safe <- npx_safe_filename(label)
  graph <- network$graph
  node_stats <- data.frame(
    Assay = igraph::V(graph)$name,
    ColorValue = igraph::V(graph)$color_value,
    Degree = igraph::V(graph)$degree,
    Betweenness = igraph::V(graph)$betweenness,
    Harmonic_Centrality = igraph::V(graph)$harmonic,
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(.data$Degree))
  table_dir <- file.path(output_dir, "STRING_tables"); dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  edge_file <- file.path(table_dir, paste0(safe, "_STRING_edges.csv")); node_file <- file.path(table_dir, paste0(safe, "_STRING_nodes.csv"))
  numeric_column <- switch(color_mode, log2fc = "Log2FC", adjusted_p = "Adjusted_pval", cluster = "ClusterNo")
  if (numeric_column %in% names(network$mapped)) {
    node_stats <- dplyr::left_join(node_stats, network$mapped |>
      dplyr::select(dplyr::all_of(c("Assay", numeric_column))) |>
      dplyr::distinct(.data$Assay, .keep_all = TRUE), by = "Assay")
  }
  cluster_colors <- network$cluster_colors
  max_value <- NULL
  if (!identical(color_mode, "cluster")) {
    max_value <- max(abs(node_stats$ColorValue), na.rm = TRUE)
    if (!is.finite(max_value) || max_value == 0) max_value <- 1
  }
  node_stats$NodeColor <- npx_ppi_node_colors(node_stats$ColorValue, color_mode, max_value, cluster_colors)
  save_output("edges", edge_file, utils::write.csv(network$edges, edge_file, row.names = FALSE))
  save_output("nodes", node_file, utils::write.csv(node_stats, node_file, row.names = FALSE))
  set.seed(seed)
  layout_matrix <- tryCatch(igraph::layout_with_fr(graph), error = function(e) { npx_rethrow_time_limit(e); e })
  if (inherits(layout_matrix, "error")) {
    output_errors <- c(output_errors, paste0("layout: ", conditionMessage(layout_matrix)))
    output_status[c("pdf", "png", "html")] <- FALSE
    return(utils::modifyList(network, list(files = files, node_statistics = node_stats,
      output_status = output_status, output_errors = output_errors,
      status = if (length(files)) "Partial" else "Error", category = "OutputError",
      reason = paste(output_errors, collapse = "; "))))
  }
  layout_data <- data.frame(
    Assay = igraph::V(graph)$name,
    x = layout_matrix[, 1],
    y = layout_matrix[, 2],
    stringsAsFactors = FALSE
  )
  node_stats <- dplyr::left_join(node_stats, layout_data, by = "Assay")
  save_output("nodes", node_file, utils::write.csv(node_stats, node_file, row.names = FALSE))
  p <- tryCatch({
  p <- ggraph::ggraph(graph, layout = "manual", x = layout_matrix[, 1], y = layout_matrix[, 2]) +
    ggraph::geom_edge_link(ggplot2::aes(width = .data$combined_score, alpha = .data$combined_score), colour = "grey65") +
    ggraph::scale_edge_width(range = c(0.4, 2.5), name = "STRING score") +
    ggraph::scale_edge_alpha(range = c(0.3, 0.85), guide = "none") +
    ggraph::geom_node_point(ggplot2::aes(fill = .data$color_value, size = .data$degree), shape = 21, color = "grey20", stroke = 0.35) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name), size = 3.5, vjust = -0.8, check_overlap = TRUE) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_size_continuous(range = c(4, 11), name = "Degree") +
    ggplot2::labs(title = label, subtitle = paste0("STRING PPI | score >= ", score_threshold, " | ", igraph::vcount(graph), " proteins | ", igraph::ecount(graph), " interactions")) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0.5), plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5), plot.margin = ggplot2::margin(25, 25, 25, 25))
  if (identical(color_mode, "log2fc")) {
    p <- p + ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_value, max_value), name = expression(Log[2] * FC))
  } else if (identical(color_mode, "cluster")) {
    p <- p + ggplot2::scale_fill_manual(values = cluster_colors, limits = names(cluster_colors),
      drop = FALSE, na.value = "#BDBDBD", name = "Cluster No.")
  } else {
    p <- p + ggplot2::scale_fill_gradientn(colors = c("#2C7BB6", "#FFFFBF", "#D7191C"), limits = c(0, max_value), name = expression(-log[10]("Adjusted P")))
  }
  p
  }, error = function(e) {
    npx_rethrow_time_limit(e)
    output_errors <<- c(output_errors, paste0("plot: ", conditionMessage(e)))
    NULL
  })
  pdf_file <- file.path(output_dir, paste0(safe, "_STRING.pdf")); png_file <- file.path(output_dir, paste0(safe, "_STRING.png"))
  if (!is.null(p)) {
    save_output("pdf", pdf_file, ggplot2::ggsave(pdf_file, p, width = width, height = height, units = "in", device = grDevices::pdf, limitsize = FALSE))
    save_output("png", png_file, ggplot2::ggsave(png_file, p, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE))
  } else output_status[c("pdf", "png")] <- FALSE
  if (requireNamespace("visNetwork", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    html_file <- file.path(output_dir, paste0(safe, "_STRING.html"))
    save_output("html", html_file, {
    color_label <- switch(color_mode, log2fc = "Log2FC", adjusted_p = "-log10 adjusted p-value", cluster = "Cluster No.")
    title_text <- if (identical(color_mode, "log2fc")) {
      paste0("Assay: ", node_stats$Assay, "<br>Log2FC: ", signif(node_stats$Log2FC, 4), "<br>Node color: ", node_stats$NodeColor, "<br>Degree: ", node_stats$Degree)
    } else if (identical(color_mode, "cluster")) {
      paste0("Assay: ", node_stats$Assay, "<br>Cluster No.: ",
        ifelse(is.na(node_stats$ClusterNo), "Unassigned", node_stats$ClusterNo),
        "<br>Node color: ", node_stats$NodeColor, "<br>Degree: ", node_stats$Degree)
    } else {
      paste0("Assay: ", node_stats$Assay, "<br>Adjusted_pval: ", signif(node_stats$Adjusted_pval, 4), "<br>", color_label, ": ", signif(node_stats$ColorValue, 4), "<br>Node color: ", node_stats$NodeColor, "<br>Degree: ", node_stats$Degree)
    }
    node_data <- node_stats |>
      dplyr::mutate(
        id = .data$Assay,
        label = .data$Assay,
        value = .data$Degree,
        title = title_text,
        color.background = .data$NodeColor,
        color.border = "#333333",
        color.highlight.background = .data$NodeColor,
        color.hover.background = .data$NodeColor,
        borderWidth = 1,
        x = .data$x * 100,
        y = -.data$y * 100,
        fixed = TRUE
      ) |>
      dplyr::select("id", "label", "title", "value", "color.background", "color.border",
                    "color.highlight.background", "color.hover.background", "borderWidth", "x", "y", "fixed")
    edge_data <- network$edges |>
      dplyr::transmute(from = .data$from, to = .data$to, title = paste0("STRING score: ", .data$combined_score), width = .data$combined_score / 200)
    widget <- visNetwork::visNetwork(node_data, edge_data) |>
      visNetwork::visNodes(shape = "dot", scaling = list(min = 8, max = 35), borderWidth = 1, borderWidthSelected = 3) |>
      visNetwork::visEdges(color = list(color = "grey70", highlight = "black", hover = "black"), smooth = FALSE, selectionWidth = 3) |>
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE), nodesIdSelection = TRUE) |>
      visNetwork::visInteraction(hover = TRUE, navigationButtons = TRUE, multiselect = TRUE) |>
      visNetwork::visPhysics(enabled = FALSE)
    if (identical(color_mode, "cluster")) {
      legend_nodes <- data.frame(label = paste("Cluster", names(cluster_colors)),
        shape = "dot", color = unname(cluster_colors), size = 15, stringsAsFactors = FALSE)
      if (anyNA(node_stats$ClusterNo)) legend_nodes <- rbind(legend_nodes,
        data.frame(label = "Unassigned", shape = "dot", color = "#BDBDBD", size = 15))
    } else {
      legend_values <- if (identical(color_mode, "log2fc")) c(-max_value, 0, max_value) else c(0, max_value / 2, max_value)
      legend_nodes <- data.frame(label = format(signif(legend_values, 3), trim = TRUE),
        shape = "dot", color = npx_ppi_node_colors(legend_values, color_mode, max_value),
        size = 15, stringsAsFactors = FALSE)
    }
    widget <- visNetwork::visLegend(widget, useGroups = FALSE, addNodes = legend_nodes,
                                    position = "right", width = 0.18, main = color_label)
    tryCatch(visNetwork::visSave(widget, html_file, selfcontained = TRUE, background = "white"), error = function(e) {
      npx_rethrow_time_limit(e)
      visNetwork::visSave(widget, html_file, selfcontained = FALSE, background = "white")
    })
    })
  }
  network$files <- files
  network$plot <- p
  network$node_statistics <- node_stats
  network$output_status <- output_status
  network$output_errors <- output_errors
  if (length(output_errors)) {
    network$status <- if (length(files)) "Partial" else "Error"
    network$category <- "OutputError"
    network$reason <- paste(output_errors, collapse = "; ")
  }
  network
}

#' Draw STRING PPI networks for two-group or ANOVA results
#' @param analysis Result from an NPXplore statistical analysis.
#' @param data Prepared NPX data, used for automatic ANOVA clustering.
#' @param anno Optional annotation data with `OlinkID` and `Human_Pathway`.
#'   Adds pathway networks using the same membership as [npx_pathway_heatmaps()].
#' @param cluster_result Optional result from [npx_cluster_assays()].
#' @param output_dir Pipeline output directory.
#' @param score_threshold Minimum STRING interaction score.
#' @param min_proteins Minimum mapped proteins per network.
#' @param width,height Figure dimensions in inches.
#' @param dpi PNG resolution.
#' @param seed Network layout seed.
#' @return A list containing per-network `results`, successful output `files`,
#'   `status`, `reason`, and the `summary` and `mapping` tables with their file paths.
#'   `PPI_summary.csv` distinguishes empty or insufficient inputs, mapping shortages,
#'   absent interactions, query errors, validation errors, and output errors.
#'   PDF/PNG/HTML fields are TRUE for saved formats, FALSE for failures, and NA
#'   when unattempted (including unavailable optional HTML packages).
#'   `STRING_ID_mapping.csv` retains unmapped inputs; `NotQueried` indicates that
#'   STRING mapping was not completed. A network may return `Partial` when some
#'   output formats failed while others were saved.
#' @export
npx_ppi_network <- function(analysis, data = NULL, anno = NULL, cluster_result = NULL, output_dir = ".", score_threshold = 700, min_proteins = 3, width = 10, height = 8, dpi = 600, seed = 123) {
  if (!is.list(analysis) || !is.data.frame(analysis$deps)) stop("analysis must contain a data.frame named 'deps'.")
  if (!is.numeric(score_threshold) || length(score_threshold) != 1L || !is.finite(score_threshold) || score_threshold < 0 || score_threshold > 1000) stop("score_threshold must be between 0 and 1000.")
  if (!is.numeric(min_proteins) || length(min_proteins) != 1L || !is.finite(min_proteins) || min_proteins < 1 || min_proteins != as.integer(min_proteins)) stop("min_proteins must be a positive integer.")
  root <- file.path(output_dir, "Results", "STRING_PPI")
  log_root <- file.path(output_dir, "log")
  results <- list(); files <- character(); summary_rows <- list()
  pathway_map <- data.frame(OlinkID = character(), Human_Pathway = character(), Assay = character())
  mapping <- data.frame(Assay = character(), UniProt = character(), STRING_id = character(), STRING_mapping_source = character())
  record <- function(label, saved, dep_data) {
    results[[label]] <<- saved
    if (!is.null(saved$files)) files <<- c(files, saved$files)
    output <- function(format) {
      if (!is.null(saved$output_status)) return(unname(saved$output_status[format]))
      if (format %in% names(saved$files)) return(TRUE)
      NA
    }
    summary_rows[[label]] <<- data.frame(
      Set = label, Status = saved$status,
      Category = if (is.null(saved$category)) saved$status else saved$category,
      Reason = if (is.null(saved$reason)) "" else saved$reason,
      InputCount = if ("Assay" %in% names(dep_data)) length(unique(dep_data$Assay[!is.na(dep_data$Assay) & nzchar(dep_data$Assay)])) else nrow(dep_data),
      MappedCount = if (!is.null(saved$mapped) && "STRING_id" %in% names(saved$mapped)) {
        ids <- saved$mapped$STRING_id
        length(unique(ids[!is.na(ids) & nzchar(ids)]))
      } else if ("Assay" %in% names(dep_data) && any(mapping$STRING_mapping_source != "NotQueried")) {
        ids <- mapping$STRING_id[mapping$Assay %in% dep_data$Assay]
        length(unique(ids[!is.na(ids) & nzchar(ids)]))
      } else NA_integer_,
      NodeCount = if (is.null(saved$nodes)) 0L else nrow(saved$nodes),
      EdgeCount = if (is.null(saved$edges)) 0L else nrow(saved$edges),
      PDF = output("pdf"), PNG = output("png"), HTML = output("html"),
      stringsAsFactors = FALSE
    )
    message("PPI: ", label, " - ", saved$status, " [", summary_rows[[label]]$Category, "]", if (nzchar(summary_rows[[label]]$Reason)) paste0(": ", summary_rows[[label]]$Reason))
  }
  failed <- function(e, category = "ProcessingError") {
    npx_rethrow_time_limit(e)
    list(status = "Error", category = if (is.null(e$category)) category else e$category, reason = conditionMessage(e))
  }
  finish <- function() {
    summary <- dplyr::bind_rows(summary_rows)
    errors <- character()
    write_diagnostic <- function(value, filename) {
      path <- file.path(log_root, filename)
      tryCatch({
        dir.create(log_root, recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(value, path, row.names = FALSE)
        path
      }, error = function(e) {
        npx_rethrow_time_limit(e)
        errors <<- c(errors, paste0(filename, ": ", conditionMessage(e)))
        NA_character_
      })
    }
    mapping_file <- write_diagnostic(mapping, "STRING_ID_mapping.csv")
    summary_file <- write_diagnostic(summary, "PPI_summary.csv")
    pathway_file <- write_diagnostic(pathway_map, "PPI_pathway_mapping.csv")
    status <- if (any(summary$Status %in% c("Error", "Partial"))) {
      if (any(summary$Status %in% c("Success", "Partial"))) "Partial" else "Failed"
    } else if (any(summary$Status == "Success")) "Success" else "Skipped"
    if (length(errors)) {
      status <- if (length(files)) "Partial" else "Failed"
      warning("PPI diagnostics could not be fully saved: ", paste(errors, collapse = "; "), call. = FALSE)
    }
    list(results = results, files = files, output_dir = root, summary = summary, mapping = mapping,
      summary_file = summary_file, mapping_file = mapping_file,
      pathway_mapping = pathway_map, pathway_mapping_file = pathway_file, status = status,
      reason = paste(c(unique(summary$Reason[nzchar(summary$Reason)]), errors), collapse = "; "))
  }
  if (!nrow(analysis$deps)) {
    record("All", list(status = "Skipped", category = "InputEmpty", reason = "No DEPs passed the analysis thresholds", mapped = data.frame()), analysis$deps)
    return(finish())
  }
  two_group <- length(analysis$levels) == 2L
  required <- c("Assay", "UniProt", if (two_group) "Log2FC" else "Adjusted_pval")
  missing <- setdiff(required, names(analysis$deps))
  if (length(missing)) {
    record("All", list(status = "Error", category = "ValidationError", reason = paste("PPI DEP data is missing:", paste(missing, collapse = ", "))), analysis$deps)
    return(finish())
  }
  mapping <- analysis$deps |>
    dplyr::transmute(Assay = as.character(.data$Assay), UniProt = as.character(.data$UniProt), STRING_id = NA_character_, STRING_mapping_source = "NotQueried") |>
    dplyr::distinct()
  groups <- list()
  cluster_colors <- NULL
  cluster_table <- NULL
  cluster_block <- NULL
  add_group <- function(label, dep_data, directory, color_mode, color_column, export_label = label) {
    groups[[label]] <<- list(data = dep_data, directory = directory, color_mode = color_mode,
                            color_column = color_column, export_label = export_label)
  }
  dep_data <- analysis$deps
  if (two_group) {
    add_group("All_DEPs", dep_data, root, "log2fc", "Log2FC", "ALL_DEPs")
  } else if (identical(analysis$method, "anova")) {
    dep_data$MinusLog10AdjustedP <- -log10(pmax(as.numeric(dep_data$Adjusted_pval), .Machine$double.xmin))
    # The all-DEP network does not depend on successful clustering.
    add_group("All_DEPs", dep_data, root, "adjusted_p", "MinusLog10AdjustedP", "ALL_DEPs")
    cluster_error <- tryCatch({
      if (is.null(cluster_result)) {
        if (is.null(data)) npx_ppi_stop("data is required for automatic ANOVA clustering.", "ValidationError")
        cluster_result <- npx_cluster_assays(data, analysis)
      }
      if (!is.null(cluster_result$status) && !identical(cluster_result$status, "Success")) {
        npx_ppi_stop(paste("Clustering unavailable:", cluster_result$reason), "MissingCluster")
      }
      cluster_table <- npx_cluster_table(dep_data, cluster_result)
      NULL
    }, error = function(e) { npx_rethrow_time_limit(e); e })
    if (inherits(cluster_error, "error")) {
      cluster_block <- if (inherits(cluster_error, "npx_insufficient_input") ||
          identical(cluster_error$category, "MissingCluster")) {
        list(status = "Skipped", category = "MissingCluster", reason = conditionMessage(cluster_error))
      } else {
        list(status = "Error", category = "ValidationError", reason = conditionMessage(cluster_error))
      }
      record("Clusters", cluster_block, dep_data)
    } else {
      unassigned <- attr(cluster_table, "unassigned_assays")
      if (length(unassigned)) record("Unassigned", list(status = "Skipped", category = "MissingCluster",
        reason = "DEPs without a valid cluster assignment; retained in all-DEP and pathway networks (grey)."),
        dep_data[dep_data$Assay %in% unassigned, , drop = FALSE])
      dep_data <- dplyr::left_join(dplyr::select(dep_data, -dplyr::any_of("ClusterNo")), cluster_table, by = "Assay")
      cluster_ids <- unique(cluster_table$ClusterNo)
      cluster_colors <- stats::setNames(grDevices::hcl.colors(length(cluster_ids), "Dark 3"), cluster_ids)
      for (cluster_id in cluster_ids) {
        add_group(paste0("Cluster_", cluster_id),
          dep_data[!is.na(dep_data$ClusterNo) & dep_data$ClusterNo == cluster_id, , drop = FALSE],
          file.path(root, paste0("Cluster_", npx_safe_filename(cluster_id))), "adjusted_p", "MinusLog10AdjustedP",
          paste0("Cluster", cluster_id, "_DEPs"))
      }
    }
  } else {
    record("All", list(status = "Error", category = "ValidationError",
      reason = "PPI network analysis requires a two-group test or ANOVA result."), dep_data)
    return(finish())
  }
  if (!is.null(anno)) {
    pathway_error <- tryCatch({
      pathway_map <- npx_pathway_mapping(analysis, anno)
      if (!is.null(cluster_table)) pathway_map <- dplyr::left_join(pathway_map, cluster_table, by = "Assay")
      NULL
    }, error = function(e) { npx_rethrow_time_limit(e); e })
    if (inherits(pathway_error, "error")) {
      record("Pathways", failed(pathway_error, "ValidationError"), dep_data)
    } else if (!nrow(pathway_map)) {
      record("Pathways", list(status = "Skipped", category = "NoPathwayAnnotation",
        reason = "No DEP assays match a non-empty Human_Pathway annotation."), dep_data)
    } else if (!two_group && is.null(cluster_table)) {
      if (is.null(cluster_block)) {
        cluster_block <- list(status = "Skipped", category = "MissingCluster",
          reason = "Pathway cluster colors require valid ANOVA cluster assignments; all-DEP PPI is still attempted.")
      }
      record("Pathways", cluster_block, dep_data)
    } else {
      pathways <- unique(pathway_map$Human_Pathway)
      safe_pathways <- make.unique(vapply(pathways, npx_safe_filename, character(1)), sep = "_")
      for (i in seq_along(pathways)) {
        pathway <- pathways[[i]]
        assays <- pathway_map$Assay[pathway_map$Human_Pathway == pathway]
        add_group(paste0("Pathway: ", pathway), dep_data[dep_data$Assay %in% assays, , drop = FALSE],
          file.path(root, "Pathways", safe_pathways[[i]]),
          if (two_group) "log2fc" else "cluster", if (two_group) "Log2FC" else "ClusterNo",
          paste0(safe_pathways[[i]], "_DEPs"))
      }
    }
  }
  string_db <- NULL; mapping_error <- NULL; initialized <- FALSE
  initialize_string <- function() {
    initialized <<- TRUE
    tryCatch({
      for (pkg in c("STRINGdb", "ggplot2", "ggraph", "igraph")) if (!requireNamespace(pkg, quietly = TRUE)) npx_ppi_stop(paste0("Package '", pkg, "' is required for npx_ppi_network()."), "DependencyError")
      string_db <<- npx_ppi_query(STRINGdb::STRINGdb$new(version = "12.0", species = 9606, score_threshold = as.integer(score_threshold), input_directory = ""))
      message("PPI: mapping ", nrow(analysis$deps), " input proteins to STRING identifiers...")
      mapping <<- npx_string_mapping(analysis$deps, string_db)
    }, error = function(e) { mapping_error <<- failed(e) })
  }
  for (label in names(groups)) {
    group <- groups[[label]]
    dep_data <- group$data
    message("PPI: ", label, " (", nrow(dep_data), " input proteins)")
    input_count <- length(unique(dep_data$Assay[!is.na(dep_data$Assay) & nzchar(dep_data$Assay)]))
    saved <- tryCatch({
      if (input_count >= min_proteins && !initialized) initialize_string()
      if (input_count >= min_proteins && !is.null(mapping_error)) {
        mapping_error
      } else {
        network <- npx_ppi_graph(dep_data, string_db, score_threshold, min_proteins, group$color_column, mapping = mapping)
        network$cluster_colors <- cluster_colors
        npx_ppi_save_network(network, group$export_label, group$directory, group$color_mode, score_threshold, width, height, dpi, seed)
      }
    }, error = function(e) failed(e))
    if (identical(saved$status, "Error")) warning("PPI ", label, " failed: ", saved$reason, call. = FALSE)
    record(label, saved, dep_data)
  }
  finish()
}
