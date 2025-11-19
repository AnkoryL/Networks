#' Run the KEGG common-pathway pipeline
#'
#' @name common_kegg_pathway
#' @description Run the pipeline for common KEGG pathways for a given gene list.
#' @param genes_list_path Path to gene list
#' @param species_prefix Species prefix, e.g., "human"
#' @param output_folder_path Output folder
#' @param min_common_pathways Minimum shared pathways for edge inclusion
#' @param use_ortholog = "none"
#' @return List with network graph and summary tables
#' @export
common_kegg_pathway <- function(
    genes_list_path,
    species_prefix,
    output_folder_path = "./results",
    use_ortholog = "none",
    label_type = "SYMBOL",
    layout_name = "layout_nicely",
    threshold = 5
) {
  message("[1/10] Starting pipeline...")

  connection_type <- as.character(match.call()[[1]])

  message("[2/10] Reading gene list from file...")

  # gene_df <-read_input_file(genes_list_path)
  gene_df <- prepare_genes(genes_list_path, species_prefix, output_folder_path)

  if (use_ortholog != "none") {
    ortholog_df <- load_ortholog(from_species = species_prefix, to_species = use_ortholog)
    gene_df <- map_genes_to_orthologs(gene_df, ortholog_df, output_folder_path)

    species_prefix <- use_ortholog
  }

  message("[3/10] Logging pipeline parameters...")
    log_file <- file.path(output_folder_path, "log_parameters.txt")
    dir.create(output_folder_path, showWarnings = FALSE, recursive = TRUE)

    log_lines <- c(
      sprintf("Species: %s", species_prefix),
      sprintf("Genes list path: %s", genes_list_path),
      sprintf("Output folder: %s", normalizePath(output_folder_path, winslash = "/", mustWork = FALSE)),
      sprintf("Connection type: %s", connection_type),
      sprintf("Use ortholog: %s", use_ortholog),
      sprintf("Label type: %s", label_type),
      sprintf("Threshold: %s", threshold),
      sprintf("Layout: %s", layout_name)
    )
    writeLines(log_lines, con = log_file)

    message("[4/10] Loading internal KEGG annotation data...")

    kegg_df<- load_annotation(source_type="KEGG", species_prefix=species_prefix,connection_type="common_kegg_pathway")

    message("[5/10] Annotating genes with KEGG common pathways...")
    message("[6/10] Checking coverage of annotation...")
    annotated <- annotate_with_kegg_common_pathways(gene_df, species_prefix, kegg_df, output_folder_path)

    annotated_genes <- annotated$annotated_genes
    not_mapped <- annotated$not_mapped

       write(sprintf("%d  genes(s) not mapped", not_mapped),
          file = log_file, append = TRUE)

    write.table(annotated_genes, file = file.path(output_folder_path,"annotated_genes_kegg.txt"), sep ="\t", row.names = FALSE, quote = FALSE)

    message("[7/10] Filtering: not applicable.")

    message("[8/10] Building KEGG network...")
    network <- build_kegg_network(
      filtered_data = annotated_genes,
      gene_df = gene_df,
      label_type = label_type,
      layout_name = layout_name,
      threshold = threshold,
      output_folder_path = output_folder_path,
      connection_type = connection_type
    )

    message("[9/10] Saving results...")
    out_path_table <- file.path(output_folder_path, paste0(species_prefix, "_", connection_type, "_gene_interaction_output_table_for_cytoscape.csv"))
    write.table(network$final_table, file = out_path_table, sep = "\t", quote = FALSE, row.names = FALSE)

    write(sprintf("Network connectivity: %.3f", network$connectivity),
          file = log_file, append = TRUE)

    message("[10/10] Pipeline completed successfully!")
    message(" → Output folder: ", output_folder_path)

}

