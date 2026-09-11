test_that("STRING retrieval errors remain distinct from empty interactions", {
  skip_if_not_installed("igraph")
  dep_data <- data.frame(Assay = c("A", "B", "C"), UniProt = c("P1", "P2", "P3"), Log2FC = 1)
  mapping <- data.frame(Assay = dep_data$Assay, UniProt = dep_data$UniProt, STRING_id = paste0("s", 1:3))
  failing_db <- list(get_interactions = function(ids) stop("HTTP 503 service unavailable"))
  expect_error(
    NPXplore:::npx_ppi_graph(dep_data, failing_db, 700, 3, "Log2FC", mapping),
    "HTTP 503 service unavailable"
  )
  empty_db <- list(get_interactions = function(ids) data.frame())
  empty <- NPXplore:::npx_ppi_graph(dep_data, empty_db, 700, 3, "Log2FC", mapping)
  expect_identical(empty$status, "Skipped")
  expect_identical(empty$reason, "No STRING interactions")
})

test_that("two-group PPI builds one all-DEP network with log2FC colors", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep_data <- data.frame(
    Assay = LETTERS[1:9], UniProt = paste0("P", 1:9), OlinkID = paste0("OID", 1:9),
    Log2FC = rep(c(1, 1, -1), each = 3), Human_Pathway = "existing annotation"
  )
  anno <- data.frame(OlinkID = dep_data$OlinkID, Human_Pathway = rep(c("A", "B", "C"), each = 3))
  alias_calls <- 0L
  constructor_threshold <- NULL
  fake_db <- list(
    get_aliases = function() {
      alias_calls <<- alias_calls + 1L
      data.frame(STRING_id = paste0("s", 1:9), alias = dep_data$Assay, source = "HGNC")
    },
    get_interactions = function(ids) {
      data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900))
    }
  )
  local_mocked_bindings(
    STRINGdb = list(new = function(version, species, score_threshold, input_directory) {
      constructor_threshold <<- score_threshold
      fake_db
    }),
    .package = "STRINGdb"
  )
  local_mocked_bindings(
    npx_ppi_save_network = function(network, label, output_dir, color_mode, ...) {
      if (identical(network$status, "Success")) network$files <- paste0(label, ".pdf")
      network$color_mode <- color_mode
      network
    },
    .package = "NPXplore"
  )
  result <- suppressMessages(npx_ppi_network(
    list(deps = dep_data, levels = c("control", "case"), method = "ttest"),
    anno = anno, score_threshold = 750, output_dir = tempfile()
  ))
  expect_identical(alias_calls, 1L)
  expect_identical(constructor_threshold, 750L)
  expect_identical(names(result$results), c("All_DEPs", "Pathway: A", "Pathway: B", "Pathway: C"))
  expect_identical(result$results$All_DEPs$status, "Success")
  expect_identical(result$results$All_DEPs$color_mode, "log2fc")
  expect_equal(result$summary$InputCount, c(9, 3, 3, 3))
  expect_identical(unname(result$files), c("ALL_DEPs.pdf", "A_DEPs.pdf", "B_DEPs.pdf", "C_DEPs.pdf"))
  expect_true(all(vapply(result$results, function(x) identical(x$color_mode, "log2fc"), logical(1))))
  expect_setequal(result$pathway_mapping$Human_Pathway, c("A", "B", "C"))
})

test_that("Wilcoxon PPI also builds one all-DEP network", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep_data <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Log2FC = c(-1, 0.5, 2))
  fake_db <- list(
    get_aliases = function() data.frame(STRING_id = paste0("s", 1:3), alias = dep_data$Assay, source = "HGNC"),
    get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900))
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(npx_ppi_save_network = function(network, label, ...) {
    network$label <- label
    network
  }, .package = "NPXplore")
  result <- suppressMessages(npx_ppi_network(list(deps = dep_data, levels = c("A", "B"), method = "wilcox"), output_dir = tempfile()))
  expect_identical(names(result$results), "All_DEPs")
  expect_identical(result$results$All_DEPs$label, "ALL_DEPs")
})

test_that("small PPI inputs skip STRING retrieval", {
  dep_data <- data.frame(Assay = "A", UniProt = "P1", Log2FC = 1)
  unused_db <- list(get_aliases = function() stop("should not download aliases"))
  result <- NPXplore:::npx_ppi_graph(dep_data, unused_db, 700, 3, "Log2FC")
  expect_identical(result$status, "Skipped")
  expect_identical(result$reason, "Too few input proteins")
})

test_that("PPI minimum protein count uses unique STRING proteins", {
  dep_data <- data.frame(Assay = c("A", "B", "C"), UniProt = c("P1", "P2", "P3"), Log2FC = 1)
  mapping <- data.frame(Assay = dep_data$Assay, UniProt = dep_data$UniProt, STRING_id = c("s1", "s1", "s2"))
  unused_db <- list(get_interactions = function(ids) stop("should not query a two-protein network"))
  result <- NPXplore:::npx_ppi_graph(dep_data, unused_db, 700, 3, "Log2FC", mapping)
  expect_identical(result$status, "Skipped")
  expect_identical(result$reason, "Too few mapped proteins")
})

test_that("ANOVA PPI processes clusters with adjusted p-value colors", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep_data <- data.frame(Assay = LETTERS[1:6], UniProt = paste0("P", 1:6), Adjusted_pval = c(0, 0.001, 0.01, 0.02, 0.03, 0.04))
  alias_calls <- 0L
  fake_db <- list(
    get_aliases = function() {
      alias_calls <<- alias_calls + 1L
      data.frame(STRING_id = paste0("s", 1:6), alias = dep_data$Assay)
    },
    get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900))
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(
    npx_ppi_save_network = function(network, label, output_dir, color_mode, ...) {
      network$files <- file.path(output_dir, paste0(label, ".pdf"))
      network$color_mode <- color_mode
      network
    },
    .package = "NPXplore"
  )
  result <- suppressMessages(npx_ppi_network(
    list(deps = dep_data, levels = c("A", "B", "C"), method = "anova"),
    cluster_result = list(clusters = data.frame(Assay = dep_data$Assay, cluster = rep(1:2, each = 3))), output_dir = tempfile()
  ))
  expect_identical(alias_calls, 1L)
  expect_identical(names(result$results), c("All_DEPs", "Cluster_1", "Cluster_2"))
  expect_true(all(vapply(result$results, function(x) identical(x$status, "Success"), logical(1))))
  expect_identical(result$results$Cluster_1$color_mode, "adjusted_p")
  expect_true(all(is.finite(result$results$Cluster_1$nodes$color_value)))
  expect_equal(result$results$Cluster_1$nodes$color_value[2:3], c(3, 2))
  expect_match(result$files[[1]], "Results/STRING_PPI/ALL_DEPs.pdf", fixed = TRUE)
  expect_true(all(grepl("Results/STRING_PPI/Cluster_", result$files[-1], fixed = TRUE)))
  expect_identical(result$results$All_DEPs$color_mode, "adjusted_p")
})

test_that("PPI reports empty, insufficient, and invalid inputs without downloading", {
  local_mocked_bindings(npx_string_mapping = function(...) stop("must not download"), .package = "NPXplore")
  out <- tempfile()
  empty <- suppressMessages(npx_ppi_network(list(deps = data.frame()), output_dir = out))
  expect_identical(empty$status, "Skipped")
  expect_identical(empty$summary$Category, "InputEmpty")
  expect_true(file.exists(empty$summary_file))
  expect_true(file.exists(empty$mapping_file))
  dep <- data.frame(Assay = "A", UniProt = "P1", OlinkID = "OID1", Log2FC = 1)
  small <- suppressMessages(npx_ppi_network(list(deps = dep, levels = c("A", "B")),
    anno = data.frame(OlinkID = "OID1", Human_Pathway = "Pathway"), output_dir = tempfile()))
  expect_identical(small$results$All_DEPs$category, "InsufficientInput")
  expect_identical(small$mapping$STRING_mapping_source, "NotQueried")
  expect_true(all(is.na(small$summary$PDF)))
  invalid <- suppressMessages(npx_ppi_network(list(deps = data.frame(Assay = "A"), levels = c("A", "B")), output_dir = tempfile()))
  expect_identical(invalid$status, "Failed")
  expect_identical(invalid$summary$Category, "ValidationError")
})

test_that("STRING mapping query failure is saved as a query failure", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), OlinkID = paste0("OID", 1:3), Log2FC = 1)
  local_mocked_bindings(STRINGdb = list(new = function(...) list(get_aliases = function() stop("HTTP 503 aliases unavailable"))), .package = "STRINGdb")
  expect_warning(result <- suppressMessages(npx_ppi_network(list(deps = dep, levels = c("A", "B")),
    output_dir = tempfile())), "HTTP 503")
  expect_identical(result$results$All_DEPs$category, "QueryError")
  expect_match(result$results$All_DEPs$reason, "HTTP 503 aliases unavailable")
  expect_identical(result$status, "Failed")
  expect_true(file.exists(result$summary_file))
  expect_equal(nrow(result$mapping), 3)
})

test_that("PPI distinguishes mapping shortage, absent edges, and score filtering", {
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Log2FC = 1)
  mapping <- data.frame(Assay = dep$Assay, UniProt = dep$UniProt, STRING_id = paste0("s", 1:3))
  missing <- mapping; missing$STRING_id[2:3] <- NA_character_
  result <- NPXplore:::npx_ppi_graph(dep, NULL, 700, 3, "Log2FC", missing)
  expect_identical(result$category, "InsufficientMapped")
  empty_db <- list(get_interactions = function(ids) data.frame())
  expect_identical(NPXplore:::npx_ppi_graph(dep, empty_db, 700, 3, "Log2FC", mapping)$category, "NoInteractions")
  weak_db <- list(get_interactions = function(ids) data.frame(from = "s1", to = "s2", combined_score = 400))
  expect_identical(NPXplore:::npx_ppi_graph(dep, weak_db, 700, 3, "Log2FC", mapping)$category, "NoHighConfidenceInteractions")
  bad_db <- list(get_interactions = function(ids) data.frame(wrong_column = 1))
  expect_error(NPXplore:::npx_ppi_graph(dep, bad_db, 700, 3, "Log2FC", mapping), class = "npx_ppi_error")
})

test_that("PPI mapping export retains unmapped proteins", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), OlinkID = paste0("OID", 1:3), Log2FC = 1)
  fake_db <- list(get_aliases = function() data.frame(STRING_id = "s1", alias = "A", source = "HGNC"))
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  result <- suppressMessages(npx_ppi_network(list(deps = dep, levels = c("A", "B")),
    anno = data.frame(OlinkID = dep$OlinkID, Human_Pathway = "Pathway"), output_dir = tempfile()))
  exported <- utils::read.csv(result$mapping_file)
  expect_equal(nrow(exported), 3)
  expect_equal(sum(exported$STRING_mapping_source == "Unmapped"), 2)
  expect_identical(result$results$All_DEPs$category, "InsufficientMapped")
  expect_equal(result$summary$MappedCount[result$summary$Set == "All_DEPs"], 1)
})

test_that("PPI interactive node colors preserve continuous numeric values", {
  skip_if_not_installed("ggplot2")
  colors <- NPXplore:::npx_ppi_node_colors(c(-2, -1, 0, 1, 2), "log2fc", 2)
  expect_equal(length(unique(colors)), 5)
  expect_identical(toupper(colors[c(1, 3, 5)]), c("#2166AC", "#FFFFFF", "#B2182B"))
  colors <- NPXplore:::npx_ppi_node_colors(c(0, 1, 2, 3, 4), "adjusted_p", 4)
  expect_equal(length(unique(colors)), 5)
  expect_identical(toupper(colors[c(1, 3, 5)]), c("#2C7BB6", "#FFFFBF", "#D7191C"))
})

test_that("a failed PPI output format preserves tables and other formats", {
  for (pkg in c("igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Log2FC = c(1, 2, 3))
  mapping <- data.frame(Assay = dep$Assay, UniProt = dep$UniProt, STRING_id = paste0("s", 1:3))
  db <- list(get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900)))
  graph <- NPXplore:::npx_ppi_graph(dep, db, 700, 3, "Log2FC", mapping)
  local_mocked_bindings(ggsave = function(filename, ...) {
    if (grepl("\\.pdf$", filename)) stop("PDF device failed")
    file.create(filename)
  }, .package = "ggplot2")
  result <- NPXplore:::npx_ppi_save_network(graph, "Test", tempfile(), "log2fc", 700, 4, 4, 72, 1)
  expect_identical(result$status, "Partial")
  expect_identical(result$category, "OutputError")
  expect_match(result$reason, "PDF device failed")
  expect_false(result$output_status[["pdf"]])
  expect_true(result$output_status[["png"]])
  expect_true(all(file.exists(result$files[c("edges", "nodes", "png")])))
  exported_nodes <- utils::read.csv(result$files[["nodes"]])
  expect_true("Log2FC" %in% names(exported_nodes))
  expect_true("NodeColor" %in% names(exported_nodes))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", exported_nodes$NodeColor)))
})

test_that("PPI HTML serializes log2FC node colors and fixed coordinates", {
  for (pkg in c("igraph", "ggraph", "ggplot2", "visNetwork", "htmlwidgets")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Log2FC = c(-2, 0, 2))
  mapping <- data.frame(Assay = dep$Assay, UniProt = dep$UniProt, STRING_id = paste0("s", 1:3))
  db <- list(get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900)))
  graph <- NPXplore:::npx_ppi_graph(dep, db, 700, 3, "Log2FC", mapping)
  local_mocked_bindings(ggsave = function(filename, ...) file.create(filename), .package = "ggplot2")
  result <- NPXplore:::npx_ppi_save_network(graph, "MixedLog2FC", tempfile(), "log2fc", 700, 4, 4, 72, 1)
  html <- paste(readLines(result$files[["html"]], warn = FALSE), collapse = "\n")
  exported_nodes <- utils::read.csv(result$files[["nodes"]])
  expect_true(result$output_status[["html"]])
  expect_true(all(c("#2166AC", "#FFFFFF", "#B2182B") %in% toupper(exported_nodes$NodeColor)))
  expect_true(all(vapply(exported_nodes$NodeColor, grepl, logical(1), x = html, fixed = TRUE)))
  expect_match(html, "Log2FC")
  expect_match(html, "\"fixed\":\\[true")
})

test_that("PPI HTML serializes adjusted p-value node colors for ANOVA networks", {
  for (pkg in c("igraph", "ggraph", "ggplot2", "visNetwork", "htmlwidgets")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Adjusted_pval = c(0.05, 0.01, 0.0001))
  dep$MinusLog10AdjustedP <- -log10(dep$Adjusted_pval)
  mapping <- data.frame(Assay = dep$Assay, UniProt = dep$UniProt, STRING_id = paste0("s", 1:3))
  db <- list(get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900)))
  graph <- NPXplore:::npx_ppi_graph(dep, db, 700, 3, "MinusLog10AdjustedP", mapping)
  local_mocked_bindings(ggsave = function(filename, ...) file.create(filename), .package = "ggplot2")
  result <- NPXplore:::npx_ppi_save_network(graph, "Cluster_1", tempfile(), "adjusted_p", 700, 4, 4, 72, 1)
  html <- paste(readLines(result$files[["html"]], warn = FALSE), collapse = "\n")
  exported_nodes <- utils::read.csv(result$files[["nodes"]])
  expect_true(result$output_status[["html"]])
  expect_true("Adjusted_pval" %in% names(exported_nodes))
  expect_true(all(vapply(exported_nodes$NodeColor, grepl, logical(1), x = html, fixed = TRUE)))
  expect_match(html, "Adjusted_pval")
  expect_match(html, "-log10 adjusted p-value")
})

test_that("time limit errors escape STRING diagnostics", {
  expect_error(NPXplore:::npx_ppi_query(stop("reached elapsed time limit")), "reached elapsed time limit")
  tryCatch(NPXplore:::npx_ppi_query(stop("reached elapsed time limit")), error = function(e) expect_false(inherits(e, "npx_ppi_error")))
})

test_that("layout failure preserves the quantitative network tables", {
  for (pkg in c("igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Log2FC = c(-2, 0, 2))
  mapping <- data.frame(Assay = dep$Assay, UniProt = dep$UniProt, STRING_id = paste0("s", 1:3))
  db <- list(get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = c(800, 900)))
  graph <- NPXplore:::npx_ppi_graph(dep, db, 700, 3, "Log2FC", mapping)
  local_mocked_bindings(layout_with_fr = function(...) stop("layout failed"), .package = "igraph")
  result <- NPXplore:::npx_ppi_save_network(graph, "All_DEPs", tempfile(), "log2fc", 700, 4, 4, 72, 1)
  expect_identical(result$status, "Partial")
  expect_true(all(file.exists(result$files[c("nodes", "edges")])))
  expect_match(result$reason, "layout failed")
})

test_that("ANOVA pathway colors preserve cluster identity across overlapping pathways", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2", "visNetwork")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:6], OlinkID = paste0("OID", 1:6),
    UniProt = paste0("P", 1:6), Adjusted_pval = rep(0.01, 6))
  anno <- data.frame(OlinkID = paste0("OID", c(1, 2, 3, 6, 2, 3, 4, 5, 6, 99)),
    Human_Pathway = c(rep("Path X", 4), rep("Path Y", 4), "Small", "Non DEP"))
  clusters <- list(clusters = data.frame(Assay = LETTERS[1:5], cluster = c(1, 2, 1, 3, 2)))
  fake_db <- list(
    get_aliases = function() data.frame(STRING_id = paste0("s", 1:6), alias = dep$Assay, source = "HGNC"),
    get_interactions = function(ids) data.frame(from = head(ids, -1), to = tail(ids, -1), combined_score = 800)
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(ggsave = function(filename, ...) file.create(filename), .package = "ggplot2")
  out <- tempfile()
  result <- suppressMessages(npx_ppi_network(list(deps = dep, method = "anova", levels = c("A", "B", "C")),
    anno = anno, cluster_result = clusters, output_dir = out))
  x <- result$results[["Pathway: Path X"]]
  y <- result$results[["Pathway: Path Y"]]
  expect_identical(x$status, "Success")
  expect_identical(y$status, "Success")
  expect_setequal(x$nodes$name, c("A", "B", "C", "F"))
  expect_setequal(y$nodes$name, c("B", "C", "D", "E"))
  expect_identical(x$cluster_colors, y$cluster_colors)
  nx <- x$node_statistics; ny <- y$node_statistics
  expect_identical(nx$NodeColor[nx$Assay == "B"], ny$NodeColor[ny$Assay == "B"])
  expect_identical(nx$NodeColor[nx$Assay == "A"], nx$NodeColor[nx$Assay == "C"])
  expect_false(identical(nx$NodeColor[nx$Assay == "A"], nx$NodeColor[nx$Assay == "B"]))
  expect_identical(nx$NodeColor[nx$Assay == "F"], "#BDBDBD")
  expect_true(is.na(nx$ClusterNo[nx$Assay == "F"]))
  built <- ggplot2::ggplot_build(x$plot)$data[[2]]
  expect_setequal(toupper(built$fill), toupper(nx$NodeColor))
  html <- paste(readLines(x$files[["html"]], warn = FALSE), collapse = "\n")
  expect_match(html, "Cluster No.", fixed = TRUE)
  expect_match(html, "Unassigned")
  expect_true(all(vapply(nx$NodeColor, grepl, logical(1), x = html, fixed = TRUE)))
  expect_match(x$files[["pdf"]], "Results/STRING_PPI/Pathways/Path_X/Path_X_DEPs_STRING.pdf", fixed = TRUE)
  expect_identical(result$results[["Pathway: Small"]]$category, "InsufficientInput")
  expect_false("Pathway: Non DEP" %in% names(result$results))
  expect_setequal(result$pathway_mapping$Assay, dep$Assay)
  expect_true(file.exists(result$pathway_mapping_file))
  expect_equal(result$summary$InputCount[result$summary$Set == "All_DEPs"], 6)
})

test_that("ANOVA all-DEP PPI survives unavailable clustering", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], OlinkID = paste0("OID", 1:3),
    UniProt = paste0("P", 1:3), Adjusted_pval = c(0.01, 0.02, 0.03))
  fake_db <- list(
    get_aliases = function() data.frame(STRING_id = paste0("s", 1:3), alias = dep$Assay),
    get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = 800)
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(npx_ppi_save_network = function(network, ...) network, .package = "NPXplore")
  result <- suppressMessages(npx_ppi_network(list(deps = dep, method = "anova", levels = c("A", "B", "C")),
    anno = data.frame(OlinkID = dep$OlinkID, Human_Pathway = "Test"),
    cluster_result = list(status = "Skipped", reason = "Too few DEPs"), output_dir = tempfile()))
  expect_identical(result$results$All_DEPs$status, "Success")
  expect_equal(result$results$All_DEPs$nodes$color_value, -log10(dep$Adjusted_pval))
  expect_identical(result$results$Clusters$category, "MissingCluster")
  expect_identical(result$results$Pathways$category, "MissingCluster")
})

test_that("malformed ANOVA cluster table is a validation error without discarding all-DEP PPI", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], OlinkID = paste0("OID", 1:3),
    UniProt = paste0("P", 1:3), Adjusted_pval = c(0.01, 0.02, 0.03))
  fake_db <- list(
    get_aliases = function() data.frame(STRING_id = paste0("s", 1:3), alias = dep$Assay),
    get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = 800)
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(npx_ppi_save_network = function(network, ...) network, .package = "NPXplore")
  result <- suppressMessages(npx_ppi_network(list(deps = dep, method = "anova", levels = c("A", "B", "C")),
    anno = data.frame(OlinkID = dep$OlinkID, Human_Pathway = "Test"),
    cluster_result = list(clusters = data.frame(Assay = dep$Assay, wrong = 1:3)), output_dir = tempfile()))
  expect_identical(result$results$All_DEPs$status, "Success")
  expect_identical(result$results$Clusters$status, "Error")
  expect_identical(result$results$Clusters$category, "ValidationError")
  expect_identical(result$results$Pathways$status, "Error")
  expect_identical(result$results$Pathways$category, "ValidationError")
  expect_identical(result$status, "Partial")
})

test_that("conflicting ANOVA cluster assignments are validation errors without discarding all-DEP PPI", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], OlinkID = paste0("OID", 1:3),
    UniProt = paste0("P", 1:3), Adjusted_pval = c(0.01, 0.02, 0.03))
  fake_db <- list(
    get_aliases = function() data.frame(STRING_id = paste0("s", 1:3), alias = dep$Assay),
    get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = 800)
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(npx_ppi_save_network = function(network, ...) network, .package = "NPXplore")
  result <- suppressMessages(npx_ppi_network(list(deps = dep, method = "anova", levels = c("A", "B", "C")),
    anno = data.frame(OlinkID = dep$OlinkID, Human_Pathway = "Test"),
    cluster_result = list(clusters = data.frame(Assay = c("A", "A", "B", "C"), cluster = c(1, 2, 1, 2))),
    output_dir = tempfile()))
  expect_identical(result$results$All_DEPs$status, "Success")
  expect_identical(result$results$Clusters$status, "Error")
  expect_identical(result$results$Clusters$category, "ValidationError")
  expect_match(result$results$Clusters$reason, "conflicting cluster assignments")
  expect_identical(result$results$Pathways$category, "ValidationError")
  expect_identical(result$status, "Partial")
})

test_that("invalid pathway annotation does not discard all-DEP PPI", {
  for (pkg in c("STRINGdb", "igraph", "ggraph", "ggplot2")) skip_if_not_installed(pkg)
  dep <- data.frame(Assay = LETTERS[1:3], UniProt = paste0("P", 1:3), Log2FC = c(-1, 1, 2))
  fake_db <- list(
    get_aliases = function() data.frame(STRING_id = paste0("s", 1:3), alias = dep$Assay),
    get_interactions = function(ids) data.frame(from = ids[1:2], to = ids[2:3], combined_score = 800)
  )
  local_mocked_bindings(STRINGdb = list(new = function(...) fake_db), .package = "STRINGdb")
  local_mocked_bindings(npx_ppi_save_network = function(network, ...) network, .package = "NPXplore")
  result <- suppressMessages(npx_ppi_network(list(deps = dep, levels = c("A", "B")),
    anno = data.frame(wrong = "column"), output_dir = tempfile()))
  expect_identical(result$results$All_DEPs$status, "Success")
  expect_identical(result$results$Pathways$category, "ValidationError")
  expect_identical(result$status, "Partial")
})
