#!/usr/bin/env Rscript
# Offline verification: no figures, enrichment services, or pipeline stage overrides.
# With an output root, compare every exported core table to fresh default results.
source(file.path(npx_check_script_directory <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))), "regenerate-examples.R"))

npx_example_compare_csv <- function(expected, actual, label = basename(actual)) {
  if (!file.exists(actual)) stop(label, ": missing export: ", actual)
  read <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  a <- read(expected); b <- read(actual)
  if (!identical(names(a), names(b)) || nrow(a) != nrow(b)) stop(label, ": exported columns or row count differ from current default results.")
  difference <- all.equal(a, b, tolerance = 1e-10, check.attributes = FALSE)
  if (!isTRUE(difference)) stop(label, ": export differs from current defaults: ", paste(difference, collapse = "; "))
  invisible(TRUE)
}

npx_check_examples <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) && args[1] %in% c("-h", "--help")) {
    cat("Usage: Rscript tools/check-examples.R [output-root]\n",
        "Without output-root, validate current default QC and statistical results only.\n", sep = "")
    return(invisible(NULL))
  }
  if (length(args) > 1L) stop("Too many arguments; use --help.")
  root <- npx_example_root()
  output_root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else NULL
  pkgload::load_all(root, quiet = TRUE, export_all = FALSE)
  pinned <- npx_example_manifest(root)
  data <- OlinkAnalyze::read_npx(file.path(root, "example/NPXfile.parquet"))
  specs <- npx_example_specs()
  temp <- tempfile("npx-example-check-"); dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  counts <- list()
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    meta <- utils::read.csv(file.path(root, spec$Metadata), check.names = FALSE, stringsAsFactors = FALSE)
    inputs <- NPXplore::npx_read_inputs(data, meta)
    expected_dir <- file.path(temp, spec$Directory)
    prepared <- NPXplore::npx_apply_qc(inputs, output_dir = expected_dir)
    analysis <- NPXplore::npx_analyze(prepared$data, variable = spec$Variable, check_log = prepared$check_log_clean)
    stopifnot(nrow(prepared$data) > 0L, nrow(analysis$result) > 0L,
              all(c("Assay", "Adjusted_pval") %in% names(analysis$result)),
              !any(prepared$meta$QC == "FAIL"))
    counts[[i]] <- data.frame(Example = spec$Directory, Variable = spec$Variable,
                              Samples = nrow(prepared$meta), Assays = nrow(analysis$result), DEPs = nrow(analysis$deps))
    if (!is.null(output_root)) {
      actual <- file.path(output_root, spec$Directory)
      if (!dir.exists(actual)) stop("Missing example directory: ", actual)
      posthoc <- if (length(analysis$levels) > 2L) {
        NPXplore::npx_anova_posthoc(prepared$data, analysis, variable = spec$Variable,
                                  check_log = prepared$check_log_clean)
      } else NULL
      exported <- NPXplore::npx_export_tables(prepared$data, analysis, output_dir = expected_dir,
                                              raw_data = prepared$raw_data, posthoc = posthoc)
      qc_files <- getFromNamespace("npx_write_qc", "NPXplore")(prepared$qc, expected_dir, inputs$meta)
      files <- c(prepared$files, exported, unlist(qc_files, use.names = FALSE))
      for (path in files) {
        relative <- substring(path, nchar(expected_dir) + 2L)
        npx_example_compare_csv(path, file.path(actual, relative), paste(spec$Directory, relative))
      }
      log_dir <- file.path(actual, "log")
      source_file <- file.path(log_dir, "source_manifest.csv")
      if (!file.exists(source_file)) stop("Missing source manifest: regenerate examples with tools/regenerate-examples.R.")
      saved <- utils::read.csv(source_file, stringsAsFactors = FALSE)
      if (!identical(saved, pinned)) stop("Source or input SHA256 differs from the recorded snapshot: ", spec$Directory)
      expected_inputs <- file.path(expected_dir, "input_manifest.csv")
      utils::write.csv(pinned[pinned$Kind == "Input", , drop = FALSE], expected_inputs, row.names = FALSE)
      npx_example_compare_csv(expected_inputs, file.path(log_dir, "input_manifest.csv"))
      defaults <- npx_example_defaults()
      expected_defaults <- file.path(expected_dir, "pipeline_defaults.csv")
      utils::write.csv(data.frame(Argument = names(defaults), Default = unname(defaults)), expected_defaults, row.names = FALSE)
      npx_example_compare_csv(expected_defaults, file.path(log_dir, "pipeline_defaults.csv"))
      expected_membership <- file.path(expected_dir, "sample_membership.csv")
      utils::write.csv(npx_example_membership(meta, spec$Variable, list(qc = prepared$qc, meta = prepared$meta)), expected_membership, row.names = FALSE)
      npx_example_compare_csv(expected_membership, file.path(log_dir, "sample_membership.csv"))
      run <- utils::read.csv(file.path(log_dir, "run_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
      expected <- list(Example = spec$Directory, Variable = spec$Variable,
                       PackageVersion = as.character(utils::packageVersion("NPXplore")), RunStatus = "Complete",
                       SourceSHA256 = npx_example_manifest_hash(pinned, "Source"), InputSHA256 = npx_example_manifest_hash(pinned, "Input"),
                       MetadataSamples = nrow(meta), Samples = nrow(prepared$meta), PreparedAssays = length(unique(prepared$data$Assay)),
                       StatisticalResults = nrow(analysis$result), DEPs = nrow(analysis$deps))
      expected <- c(expected, setNames(as.list(defaults), paste0("Default_", names(defaults))))
      if (nrow(run) != 1L || !all(names(expected) %in% names(run))) stop("Incomplete run manifest: ", spec$Directory)
      for (key in names(expected)) if (!identical(as.character(run[[key]]), as.character(expected[[key]]))) stop("Stale run manifest field ", key, ": ", spec$Directory)
      if (!grepl("^[0-9a-f]{40}$", run$SourceGitCommit) || is.na(run$SourceDirty) || is.na(run$Seed)) stop("Invalid Git/seed provenance: ", spec$Directory)
      stages <- utils::read.csv(file.path(log_dir, "Functional_analysis_status.csv"), stringsAsFactors = FALSE)
      for (stage in c("Clustering", "GO", "PPI")) {
        if (!identical(as.character(run[[paste0(stage, "Status")]]), stages$Status[match(stage, stages$Stage)])) stop("Stage status differs from run manifest: ", stage)
      }
      if (!file.exists(file.path(log_dir, "sessionInfo.txt"))) stop("Missing R session provenance.")
    }
    message(spec$Variable, ": ", nrow(prepared$meta), " samples; ", nrow(analysis$result), " statistical results; ", nrow(analysis$deps), " DEPs", if (!is.null(output_root)) "; exports verified" else "")
  }
  npx_example_assert_manifest(root, pinned)
  invisible(do.call(rbind, counts))
}

if (sys.nframe() == 0L) npx_check_examples()
