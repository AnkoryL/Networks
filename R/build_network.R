#' Build interaction network and save results
#'
#' @param filtered_data Data frame of filtered GO annotations (GO_ID, ENSEMBL, EVIDENCE, SYMBOL)
#' @param gene_df Data frame of input genes (ENSEMBL)
#' @param duplicated_symbols Character vector of duplicated symbols (optional)
#' @param label_type One of "SYMBOL", "ENSEMBL", "both", "SYMBOLalias"
#' @param layout_name Name of igraph layout function (e.g., "layout_nicely")
#' @param threshold Minimum number of shared GO terms to draw edge
#' @param output_folder_path Path to save results
#' @param connection_type String to annotate connection type
#'
#' @return List with elements: graph (igraph object), final_table (data.frame), plot_path (string)
#' @export
build_network <- function(filtered_data, gene_df, duplicated_symbols = NULL,
                          label_type = "SYMBOL", layout_name = "layout_nicely",
                          threshold = 5, output_folder_path, connection_type = "common_go_term") {

  if (!requireNamespace("igraph", quietly = TRUE)) stop("igraph package is required")
  graph_from_data_frame <- igraph::graph_from_data_frame
  graph_from_adjacency_matrix <- igraph::graph_from_adjacency_matrix
  neighbors <- igraph::neighbors
  degree <- igraph::degree
  delete_edges <- igraph::delete_edges
  delete_vertices <- igraph::delete_vertices

  gene_go_terms <- stats::aggregate(GO_ID ~ ENSEMBL, filtered_data, toString)

  go_gene_graph <- graph_from_data_frame(
    d = filtered_data[, c("GO_ID", "ENSEMBL")],
    directed = FALSE
  )

  gene_matrix <- matrix(0, nrow = nrow(gene_go_terms), ncol = nrow(gene_go_terms))
  rownames(gene_matrix) <- gene_go_terms$ENSEMBL
  colnames(gene_matrix) <- gene_go_terms$ENSEMBL

  edge_info <- list()
  for (i in 1:(nrow(gene_go_terms)-1)) {
    for (j in (i+1):nrow(gene_go_terms)) {
      gi <- rownames(gene_matrix)[i]
      gj <- rownames(gene_matrix)[j]

      terms_i <- names(neighbors(go_gene_graph, gi))
      terms_j <- names(neighbors(go_gene_graph, gj))

      common_terms <- intersect(terms_i, terms_j)
      if (length(common_terms) == 0) next

      evi_i <- filtered_data$EVIDENCE[filtered_data$ENSEMBL == gi & filtered_data$GO_ID %in% common_terms]
      evi_j <- filtered_data$EVIDENCE[filtered_data$ENSEMBL == gj & filtered_data$GO_ID %in% common_terms]

      edge_info[[paste0(gi, "--", gj)]] <- paste0(common_terms, "(", evi_i, "—", evi_j, ")", collapse = "; ")

      gene_matrix[i, j] <- length(common_terms)
      gene_matrix[j, i] <- length(common_terms)
    }
  }

  gene_gene_graph <- graph_from_adjacency_matrix(gene_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  igraph::E(gene_gene_graph)$tooltip <- paste0(names(edge_info), "_", unlist(edge_info))

  igraph::V(gene_gene_graph)$type <- "gene"
  igraph::V(gene_gene_graph)$tooltip <- paste0(igraph::V(gene_gene_graph)$name, "_GO_terms:", gene_go_terms$GO_ID[match(igraph::V(gene_gene_graph)$name, gene_go_terms$ENSEMBL)])

  gene_gene_graph <- delete_edges(gene_gene_graph, igraph::E(gene_gene_graph)[igraph::E(gene_gene_graph)$weight < threshold])
  gene_gene_graph <- delete_vertices(gene_gene_graph, which(degree(gene_gene_graph) == 0))

  ## Label assignment
  if (label_type == "ENSEMBL") {
    igraph::V(gene_gene_graph)$label <- igraph::V(gene_gene_graph)$name
  } else if (label_type == "SYMBOL") {
    igraph::V(gene_gene_graph)$label <- filtered_data$SYMBOL[match(igraph::V(gene_gene_graph)$name, filtered_data$ENSEMBL)]
  } else if (label_type == "both") {
    sym <- filtered_data$SYMBOL[match(igraph::V(gene_gene_graph)$name, filtered_data$ENSEMBL)]
    igraph::V(gene_gene_graph)$label <- paste0(sym, " (", igraph::V(gene_gene_graph)$name, ")")
  } else if (label_type == "SYMBOLalias") {
    if (is.null(duplicated_symbols)) stop("SYMBOLalias requires duplicated_symbols")
    symbol_by_ensembl <- filtered_data$SYMBOL[match(igraph::V(gene_gene_graph)$name, filtered_data$ENSEMBL)]
    is_dup <- symbol_by_ensembl %in% duplicated_symbols
    igraph::V(gene_gene_graph)$label <- ifelse(is_dup, paste0(symbol_by_ensembl, " (", igraph::V(gene_gene_graph)$name, ")"), symbol_by_ensembl)
  }

  igraph::V(gene_gene_graph)$color <- "lightblue"

  layout_fun <- getExportedValue("igraph", layout_name)
  coords <- layout_fun(gene_gene_graph)

  plot_path <- file.path(output_folder_path, "gene_network_plot.png")
  grDevices::png(plot_path, width = 4800, height = 3240, res = 300)
  plot(
    gene_gene_graph,
    vertex.size = 8,
    vertex.label.cex = 0.6,
    vertex.label.dist = 1.5,
    edge.width = igraph::E(gene_gene_graph)$weight / 10,
    edge.color = "gray",
    edge.label.cex = 0.5,
    edge.label.color = "darkgray",
    edge.label.font = 1,
    edge.label.dist = 0.01,
    layout = coords
  )
  grDevices::dev.off()

  ## Final table
  edge_df <- data.frame(edge = names(edge_info), values = unlist(edge_info), stringsAsFactors = FALSE)
  edge_df <- cbind(edge_df, do.call(rbind, strsplit(edge_df$edge, "--")))
  names(edge_df)[3:4] <- c("gene1", "gene2")
  edge_df$connection_type <- connection_type
  edge_df$edge <- NULL
  edge_df <- edge_df[, c("gene1", "gene2", "connection_type", "values")]

  out_path_table <- file.path(output_folder_path, "gene_interaction_output_table_base.csv")
  utils::write.table(edge_df, file = out_path_table, sep = "\t", quote = FALSE, row.names = FALSE)

  return(list(
    graph = gene_gene_graph,
    final_table = edge_df,
    plot_path = plot_path
  ))
}
