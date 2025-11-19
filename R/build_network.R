#' Build a GO-based gene–gene network
#'
#' @name build_go_network
#'
#' @param filtered_data Data frame with GO annotations after filtering.
#' @param gene_df Data frame with gene identifiers.
#' @param output_folder_path Directory where plots and tables will be saved.
#' @param duplicated_symbols Optional character vector of duplicated gene
#' @param label_type Label mode for graph vertices
#' @param layout_name Name of igraph layout function to use.
#' @param threshold Minimum number of shared GO terms for an edge to be kept.
#' @param connection_type Name of the connection type used in output tables.
#'
#' @return A list with:
#'   \item{graph}{igraph object}
#'   \item{final_table}{Edge table with shared GO information}
#'   \item{plot_path}{Path to saved plot}
#'   \item{connectivity}{Fraction of vertices in the largest component}
#'
#' @export
build_go_network <- function(
    filtered_data,
    gene_df,
    output_folder_path,
    duplicated_symbols = NULL,
    label_type,
    layout_name,
    threshold,
    connection_type) {

  if (is.null(filtered_data) || nrow(filtered_data) < 2) {
    stop(error_messages$no_genes_after_filtering)
  }

  intersection <- intersect(filtered_data$Ensembl_gene_id, unique(gene_df$Ensembl_gene_id))
  if (length(intersection) == 0) stop(error_messages$no_genes_overlap)

  gene_go_terms <- stats::aggregate(GO_ID ~ Ensembl_gene_id, filtered_data, toString)

  go_gene_graph <- igraph::graph_from_data_frame(
    d = filtered_data[, c("GO_ID", "Ensembl_gene_id")],
    directed = FALSE
  )

  all_vertices <- igraph::V(go_gene_graph)$name
  gene_list <- gene_go_terms$Ensembl_gene_id
  valid_genes <- gene_list[gene_list %in% all_vertices]

  if (length(valid_genes) < 2) stop("Not enough valid genes to build network.")

  edge_info <- list()

  gene_matrix <- matrix(0, nrow = length(valid_genes), ncol = length(valid_genes),
                        dimnames = list(valid_genes, valid_genes))

  for (i in 1:(length(valid_genes) - 1)) {
    for (j in (i + 1):length(valid_genes)) {
      gi <- valid_genes[i]
      gj <- valid_genes[j]

      terms_i <- names(igraph::neighbors(go_gene_graph, gi))
      terms_j <- names(igraph::neighbors(go_gene_graph, gj))
      common_terms <- intersect(terms_i, terms_j)
      if (length(common_terms) == 0) next

      values_vec <- sapply(common_terms, function(go) {
        evi_i <- filtered_data$EVIDENCE[filtered_data$Ensembl_gene_id == gi & filtered_data$GO_ID == go]
        evi_j <- filtered_data$EVIDENCE[filtered_data$Ensembl_gene_id == gj & filtered_data$GO_ID == go]
        paste0(go, "(", evi_i, "—", evi_j, ")")
      })

      edge_info[[paste0(gi, "--", gj)]] <- paste(values_vec, collapse = "; ")
      gene_matrix[i, j] <- length(common_terms)
      gene_matrix[j, i] <- length(common_terms)
    }
  }

  gene_gene_graph <- igraph::graph_from_adjacency_matrix(gene_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  gene_gene_graph <- igraph::delete_edges(gene_gene_graph, igraph::E(gene_gene_graph)[igraph::E(gene_gene_graph)$weight < threshold])
  gene_gene_graph <- igraph::delete_vertices(gene_gene_graph, which(igraph::degree(gene_gene_graph) == 0))
  vertices <- igraph::V(gene_gene_graph)$name


  g_components <- igraph::components(gene_gene_graph)
  max_component_size <- max(g_components$csize)
  total_vertices <- igraph::vcount(gene_gene_graph)
  connectivity_ratio <- max_component_size / total_vertices


  if (label_type == "SYMBOL") {
    igraph::V(gene_gene_graph)$label <- filtered_data$SYMBOL[match(vertices, filtered_data$Ensembl_gene_id)]
  } else if (label_type == "ENSEMBL") {
    igraph::V(gene_gene_graph)$label <- vertices
  } else if (label_type == "both") {
    sym <- filtered_data$SYMBOL[match(vertices, filtered_data$Ensembl_gene_id)]
    igraph::V(gene_gene_graph)$label <- paste0(sym, " (", vertices, ")")
  } else if (label_type == "SYMBOLalias") {
    if (is.null(duplicated_symbols)) stop(error_messages$no_duplicated_symbol)
    symbol_by_ensembl <- filtered_data$SYMBOL[match(vertices, filtered_data$Ensembl_gene_id)]
    is_dup <- symbol_by_ensembl %in% duplicated_symbols
    igraph::V(gene_gene_graph)$label <- ifelse(is_dup, paste0(symbol_by_ensembl, " (", vertices, ")"), symbol_by_ensembl)
  }

  igraph::V(gene_gene_graph)$color <- "lightblue"

  layout_fun <- getExportedValue("igraph", layout_name)
  coords <- layout_fun(gene_gene_graph)

  plot_path <- file.path(output_folder_path, "gene_network_plot.png")
  grDevices::png(plot_path, width = 4800, height = 3240, res = 300)
  plot(gene_gene_graph, vertex.size = 8, vertex.label.cex = 0.6,
       vertex.label.dist = 1.5, edge.width = igraph::E(gene_gene_graph)$weight / 10,
       edge.color = "gray", edge.label.cex = 0.5, edge.label.color = "darkgray",
       edge.label.font = 1, edge.label.dist = 0.01, layout = coords)
  grDevices::dev.off()

  if (length(edge_info) == 0) stop("No gene interactions found.")

  edge_info_filtered <- edge_info[sapply(edge_info, function(x) {
    length(strsplit(x, "; ")[[1]]) >= threshold
  })]

  if (length(edge_info_filtered) == 0) stop("No gene interactions found above threshold.")

  edge_df <- data.frame(
    gene1 = sapply(strsplit(names(edge_info), "--"), `[`, 1),
    gene2 = sapply(strsplit(names(edge_info), "--"), `[`, 2),
    connection_type_info = paste0("GO",";",connection_type),
    supplemental_info <- paste("undirected"),
    values = unlist(edge_info),
    stringsAsFactors = FALSE
  )

  if (label_type == "SYMBOL") {
    edge_df$gene1 <- filtered_data$SYMBOL[match(edge_df$gene1, filtered_data$Ensembl_gene_id)]
    edge_df$gene2 <- filtered_data$SYMBOL[match(edge_df$gene2, filtered_data$Ensembl_gene_id)]
  } else if (label_type == "both") {
    sym1 <- filtered_data$SYMBOL[match(edge_df$gene1, filtered_data$Ensembl_gene_id)]
    sym2 <- filtered_data$SYMBOL[match(edge_df$gene2, filtered_data$Ensembl_gene_id)]
    edge_df$gene1 <- paste0(sym1, " (", edge_df$gene1, ")")
    edge_df$gene2 <- paste0(sym2, " (", edge_df$gene2, ")")
  } else if (label_type == "SYMBOLalias") {
    symbol_dict <- filtered_data %>%
      dplyr::filter(!is.na(Ensembl_gene_id), !is.na(SYMBOL)) %>%
      dplyr::distinct(Ensembl_gene_id, SYMBOL) %>%
      dplyr::mutate(is_dup = SYMBOL %in% duplicated_symbols,
                    display_label = ifelse(is_dup, paste0(SYMBOL, " (", ENSEMBL, ")"), SYMBOL))
    edge_df$gene1 <- symbol_dict$display_label[match(edge_df$gene1, symbol_dict$Ensembl_gene_id)]
    edge_df$gene2 <- symbol_dict$display_label[match(edge_df$gene2, symbol_dict$Ensembl_gene_id)]
  }

  return(list(
    graph = gene_gene_graph,
    final_table = edge_df,
    plot_path = plot_path,
    connectivity = connectivity_ratio
  ))
}

#' Build a KEGG-based gene–gene network
#'
#' @name build_kegg_network
#'
#' @param filtered_data Data frame with KEGG annotation after filtering.
#' @param gene_df Data frame with the original gene list.
#' @param output_folder_path Directory where plots and tables will be saved.
#' @param label_type Vertex label style
#' @param layout_name Name of igraph layout function.
#' @param connection_type Text identifier of the connection type for output tables.
#' @param threshold Minimum number of shared pathways required for edge retention.
#'
#' @return A list with:
#'   \item{graph}{igraph object}
#'   \item{final_table}{Edge table with shared pathway information}
#'   \item{plot_path}{Path to PNG plot}
#'   \item{connectivity}{Fraction of vertices in the largest connected component}
#'
#' @export
build_kegg_network <- function(
    filtered_data,
    gene_df,
    output_folder_path,
    label_type,
    layout_name,
    connection_type,
    threshold) {

  if (is.null(filtered_data) || nrow(filtered_data) < 2) stop("No genes available for network.")


  intersection <- intersect(filtered_data$Gene_symbol, unique(gene_df$Gene_symbol))
  if (length(intersection) == 0) stop("No genes overlap between annotation and input list.")

  gene_pathways <- stats::aggregate(pathway_id ~ Gene_symbol, filtered_data, toString)

  gene_list <- unique(  filtered_data$Gene_symbol)
  gene_matrix <- matrix(0, nrow = length(gene_list), ncol = length(gene_list),
                        dimnames = list(gene_list, gene_list))
  edge_info <- list()

  for (i in 1:(length(gene_list) - 1)) {
    for (j in (i + 1):length(gene_list)) {
      gi <- gene_list[i]
      gj <- gene_list[j]

      pathways_i <- filtered_data$pathway_id[  filtered_data$Gene_symbol == gi]
      pathways_j <- filtered_data$pathway_id[  filtered_data$Gene_symbol == gj]
      common_paths <- intersect(pathways_i, pathways_j)
      if (length(common_paths) == 0) next

      edge_info[[paste0(gi, "--", gj)]] <- paste(common_paths, collapse = "; ")
      gene_matrix[i, j] <- length(common_paths)
      gene_matrix[j, i] <- length(common_paths)
    }
  }

  gene_gene_graph <- igraph::graph_from_adjacency_matrix(gene_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  gene_gene_graph <- igraph::delete_edges(gene_gene_graph, igraph::E(gene_gene_graph)[igraph::E(gene_gene_graph)$weight < threshold])
  gene_gene_graph <- igraph::delete_vertices(gene_gene_graph, which(igraph::degree(gene_gene_graph) == 0))
  vertices <- igraph::V(gene_gene_graph)$name

  g_components <- igraph::components(gene_gene_graph)
  max_component_size <- max(g_components$csize)
  total_vertices <- igraph::vcount(gene_gene_graph)
  connectivity_ratio <- max_component_size / total_vertices


  if (label_type == "SYMBOL") {
    igraph::V(gene_gene_graph)$label <- vertices
  } else if (label_type == "both") {
    igraph::V(gene_gene_graph)$label <- paste0(vertices, " (", vertices, ")")
  }

  igraph::V(gene_gene_graph)$color <- "lightblue"

  layout_fun <- getExportedValue("igraph", layout_name)
  coords <- layout_fun(gene_gene_graph)

  plot_path <- file.path(output_folder_path, "gene_network_plot.png")
  grDevices::png(plot_path, width = 4800, height = 3240, res = 300)
  plot(gene_gene_graph, vertex.size = 8, vertex.label.cex = 0.6,
       vertex.label.dist = 1.5, edge.width = igraph::E(gene_gene_graph)$weight,
       edge.color = "gray", layout = coords)
  grDevices::dev.off()

  if (length(edge_info) == 0) stop("No gene interactions found.")

  edge_info_filtered <- edge_info[sapply(edge_info, function(x) {
    length(strsplit(x, "; ")[[1]]) >= threshold
  })]

  if (length(edge_info_filtered) == 0) stop("No gene interactions found above threshold.")

  edge_df <- data.frame(
    gene1 = sapply(strsplit(names(edge_info_filtered), "--"), `[`, 1),
    gene2 = sapply(strsplit(names(edge_info_filtered), "--"), `[`, 2),
    connection_type_info = paste0("KEGG",";",connection_type),
    supplemental_info <- paste("undirected"),
    values = unlist(edge_info_filtered),
    stringsAsFactors = FALSE
  )

  return(list(
    graph = gene_gene_graph,
    final_table = edge_df,
    plot_path = plot_path,
    connectivity = connectivity_ratio
  ))
}
