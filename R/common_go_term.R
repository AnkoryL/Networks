#' Run the full GO-based gene interaction pipeline using internal GO annotations
#'
#' @param genes_list_path Path to gene list (CSV or TXT).
#' @param species_prefix Short prefix (e.g., "human", "mouse").
#' @param go_type GO category: BP, CC, MF, IP.
#' @param output_folder_path Output directory (default: "./results").
#' @param config_evidence_path Optional path to config_evidence.txt. If missing, uses internal.
#' @param level_from GO tree level (default: 5).
#' @param level_to GO tree level (default: 20).
#' @param label_type Gene label type (default: "SYMBOL").
#' @param layout_name Graph layout (default: "layout_nicely").
#' @param threshold Threshold for edge filtering (default: 5).
#' @param use_ortholog = "none"
#'
#' @return Data frame of filtered annotations (invisible).
#' @export
#'
common_go_term <- function(
    genes_list_path,
    species_prefix,
    go_type,
    output_folder_path = file.path(Sys.getenv("USERPROFILE"), "Desktop", "Networks_results"),
    use_ortholog = "none",
    config_evidence_path = default_evidence_experimental_only(),
    level_from = 5,
    level_to = 20,
    label_type = "SYMBOL",
    layout_name = "layout_nicely",
    threshold = 5
) {
  message("[1/10] Starting pipeline...")

    connection_type <- as.character(match.call()[[1]])

    message("[2/10] Reading gene list...")
    # gene_df <- read_input_file(genes_list_path)
    gene_df <- prepare_genes(genes_list_path, species_prefix, output_folder_path)

    if (use_ortholog != "none") {
      ortholog_df <- load_ortholog(from_species = species_prefix, to_species = use_ortholog)
      gene_df <- map_genes_to_orthologs(gene_df, ortholog_df, output_folder_path)

      species_prefix <- use_ortholog
    }


  message("[3/10] Checking/creating output folder and logging parameters...")
  log_file <- file.path(output_folder_path, "log_parameters.txt")
  dir.create(output_folder_path, showWarnings = FALSE, recursive = TRUE)
  log_lines <- c(
    sprintf("Species: %s", species_prefix),
    sprintf("Genes list path: %s", genes_list_path),
    sprintf("GO type: %s", go_type),
    sprintf("Output folder: %s", output_folder_path),
    sprintf("Connection type: %s", connection_type),
    sprintf("Use ortholog: %s", use_ortholog),
    sprintf("GO level range: %s to %s", level_from, level_to),
    sprintf("Label type: %s", label_type),
    sprintf("Threshold: %s", threshold),
    sprintf("Layout: %s", layout_name)
  )
  writeLines(log_lines, con = log_file)


  message("[4/10] Loading internal GO annotation data...")
  go_df<- load_annotation(source_type="GO", species_prefix=species_prefix, go_type = go_type)

  message("[5/10] Annotating genes with GO terms...")
  message("[6/10] Checking annotation coverage...")

  annotated <- annotate_with_go_terms(gene_df, go_df, output_folder_path, species_prefix)

  annotated_genes <- annotated$annotated_genes
  not_mapped <- annotated$not_mapped


  write(sprintf("%d  genes(s) not mapped", not_mapped),
        file = log_file, append = TRUE)


  message("[7/10] Filtering GO annotations...")
  filtered_data <- filter_go_annotations(annotated_genes, config_evidence_path, level_from, level_to)
  if (is.null(filtered_data) || nrow(filtered_data) < 2) stop(error_messages$no_genes_after_filtering)

  message("[8/10] Building network...")
  network <- build_go_network(
    filtered_data = filtered_data,
    gene_df = gene_df,
    label_type = label_type,
    layout_name = layout_name,
    threshold = threshold,
    output_folder_path = output_folder_path,
    connection_type = connection_type
  )

  message("[9/10] Save results...")
  out_path_table <- file.path(output_folder_path, paste0(species_prefix, "_", connection_type, "_gene_interaction_output_table_for_cytoscape.csv"))
  write.table(network$final_table, file = out_path_table, sep = "\t", quote = FALSE, row.names = FALSE)

  write(sprintf("Network connectivity: %.3f", network$connectivity),
        file = log_file, append = TRUE)

  message("[10/10] Pipeline completed successfully!")
  message(" → Output folder: ", output_folder_path)

  return(invisible(filtered_data))
}

