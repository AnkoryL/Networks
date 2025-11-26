#' Build KEGG relation graph for a gene list
#'
#' @name relation_in_kegg_pathway
#'
#' @param genes_list_path Path to the input file with a gene list.
#' @param species_prefix Character string defining which internal
#' @param output_folder_path Directory for saving results. Created if absent.
#' @param label_type Label mode for graph vertices. `"SYMBOL"`.
#' @param layout_name Name of the igraph layout function to use.
#' @param use_ortholog = "none"
#' @param subtype_filter Vector of subtype names to exclude from relations.
#'
#' @export
relation_in_kegg_pathway <- function(
    genes_list_path,
    species_prefix,
    output_folder_path = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Networks_results"),
    use_ortholog = "none",
    label_type = "SYMBOL",
    layout_name = "layout_nicely",
    subtype_filter = "compound") {

  message("[1/10] Starting pipeline...")

  connection_type <- as.character(match.call()[[1]])

    message("[2/10] Reading gene list from file...")

  gene_df <- prepare_genes(genes_list_path, species_prefix, output_folder_path)

  if (use_ortholog != "none") {
    ortholog_df <- load_ortholog(from_species = species_prefix, to_species = use_ortholog)
    gene_df <- map_genes_to_orthologs(gene_df, ortholog_df, output_folder_path)

    species_prefix <- use_ortholog
  }

  message("[3/10] Logging pipeline parameters...")
  output_folder_path <- normalizePath(output_folder_path, winslash = "/", mustWork = FALSE)
  if (!dir.exists(output_folder_path)) dir.create(output_folder_path, recursive = TRUE)

  log_file <- file.path(output_folder_path, "log_parameters.txt")

  log_lines <- c(
    sprintf("Species: %s", species_prefix),
    sprintf("Genes list path: %s", genes_list_path),
    sprintf("Output folder: %s", output_folder_path),
    sprintf("Connection type: %s", connection_type),
    sprintf("Label type: %s", label_type),
    sprintf("Layout: %s", layout_name)
  )
  writeLines(log_lines, con = log_file)


  message("[4/10] Loading internal KEGG annotation data...")

  kegg_df<- load_annotation(source_type="KEGG", species_prefix=species_prefix,connection_type="relation_in_kegg_pathway")

  message("[5/10] Annotating genes with relation in KEGG pathways...")
  message("[6/10] Checking coverage of annotation...")
  annotated_genes <- unique(c(kegg_df$from, kegg_df$to)[c(kegg_df$from, kegg_df$to) %in% gene_df$Gene_symbol])
  not_mapped <- setdiff(gene_df$Gene_symbol, annotated_genes)

  out_not_mapped <- file.path(
    output_folder_path,
    paste0(species_prefix, "_not_mapped_table.txt")
  )

  write.table(
    not_mapped,
    file = out_not_mapped,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message("Annotated genes: ", length(annotated_genes))
  message("Not mapped genes: ", length(not_mapped))

  write(sprintf("%d  genes(s) not mapped", length(not_mapped)),
        file = log_file, append = TRUE)


    message("[7/10] Filtering...")
  filtered_df <- kegg_df[
    kegg_df$from %in% gene_df$Gene_symbol &
    kegg_df$to %in% gene_df$Gene_symbol &
    kegg_df$from != kegg_df$to,]


    if (!is.null(subtype_filter)) {
      filtered_df <- filtered_df[!filtered_df$subtype_name %in% subtype_filter, ]
  }

  write.table(filtered_df, file = file.path(output_folder_path,"KEGG_relation_filtered.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

  message("[8/10] Building network...")

  gene_graph <- igraph::graph_from_data_frame(filtered_df[, c("from", "to")], directed = TRUE)

  vertices <- igraph::V(gene_graph)$name
  if (label_type == "SYMBOL") {
    igraph::V(gene_graph)$label <- vertices
  } else if (label_type == "both") {
    igraph::V(gene_graph)$label <- paste0(vertices, " (", vertices, ")")
  }
  igraph::V(gene_graph)$color <- "lightblue"

  layout_fun <- getExportedValue("igraph", layout_name)
  coords <- layout_fun(gene_graph)

  plot_path <- file.path(output_folder_path, "gene_relation_network_plot.png")
  grDevices::png(plot_path, width = 4800, height = 3240, res = 300)
  plot(gene_graph, vertex.size = 8, vertex.label.cex = 0.6,
       vertex.label.dist = 1.5, edge.color = "gray", layout = coords)
  grDevices::dev.off()

  g_components <- igraph::components(gene_graph)
  connectivity_ratio <- max(g_components$csize) / igraph::vcount(gene_graph)

  edge_df <- filtered_df
  edge_df$connection_type_info <- paste(
    "KEGG",
    "relation_in_KEGG_pathway",
    sep = ";"
  )
  edge_df$supplemental_info <- paste(
    edge_df$pathway_id,
    edge_df$relation_type,
    edge_df$subtype_name,
    sep = ";"
  )

  edge_df <- edge_df[, c("from", "to", "connection_type_info", "supplemental_info")]

  final_table <- edge_df
  connectivity <- connectivity_ratio
  message("[9/10] Saving results...")

  out_path_table <- file.path(output_folder_path, paste0(species_prefix, "_", connection_type, "_gene_interaction_output_table_for_cytoscape.csv"))
  write.table(final_table, file = out_path_table, sep = "\t", quote = FALSE, row.names = FALSE)

  write(sprintf("Network connectivity: %.3f", connectivity),
        file = log_file, append = TRUE)

  message("[10/10] Pipeline completed successfully!")
  message(" → Output folder: ", output_folder_path)

}
