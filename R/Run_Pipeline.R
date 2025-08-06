#' Run the full GO-based gene interaction pipeline
#'
#' @param genes_list_path Path to gene list (CSV or TXT).
#' @param species_prefix Short prefix (e.g., "human", "mouse").
#' @param go_type GO category: BP, CC, MF, IP.
#' @param config_path Optional path to main_config.txt for loading defaults.
#' @param config_evidence_path Optional path to config_evidence.txt. If missing, uses internal.
#' @param output_folder_path Output directory (default: "./results").
#' @param connection_type Type of connection (default: "common_go_term").
#' @param level_from GO tree level (default: 5).
#' @param level_to GO tree level (default: 20).
#' @param label_type Gene label type (default: "SYMBOL").
#' @param layout_name Graph layout (default: "layout_nicely").
#' @param threshold Threshold for edge filtering (default: 5).
#'
#' @return NULL
#' @export
run_pipeline <- function(
    genes_list_path,
    species_prefix,
    go_type,
    output_folder_path = "./results",
    config_path = NULL,
    config_evidence_path = default_evidence_experimental_only(),
    connection_type = "common_go_term",
    level_from = 5,
    level_to = 20,
    label_type = "SYMBOL",
    layout_name = "layout_nicely",
    threshold = 5
) {
  message("[1/10] Starting pipeline...")

  if (!is.null(config_path)) {
    if (!file.exists(config_path)) {
      stop("main_config.txt not found at: ", config_path)
    }
    message("[2/10] Reading configuration from file...")
    config <- read_config_txt(config_path)

    if (missing(genes_list_path))    genes_list_path    <- config$genes_list_path
    if (missing(species_prefix))     species_prefix     <- config$species_prefix
    if (missing(go_type))            go_type            <- config$go_type
    if (missing(output_folder_path)) output_folder_path <- config$output_folder_path
    if (missing(connection_type))    connection_type    <- config$connection_type
    if (missing(level_from))         level_from         <- config$level_from
    if (missing(level_to))           level_to           <- config$level_to
    if (missing(label_type))         label_type         <- config$label_type
    if (missing(layout_name))        layout_name        <- config$layout
    if (missing(threshold))          threshold          <- config$threshold
    if (is.null(config_evidence_path) || config_evidence_path == "") {
      config_evidence_path <- default_evidence_experimental_only()
    }
  } else {
    message("[2/10] Using provided arguments...")
  }

  go_file_path <- get_go_file_by_prefix(species_prefix, go_type)

  message("[3/10] Logging parameters to output folder...")
  log_lines <- c(
    sprintf("   Species: %s", species_prefix),
    sprintf("   GO file path: %s", go_file_path),
    sprintf("   Genes list path: %s", genes_list_path),
    sprintf("   Config evidence path: %s", config_evidence_path),
    sprintf("   Output folder path: %s", normalizePath(output_folder_path, winslash = "/", mustWork = FALSE)),
    sprintf("   Connection type: %s", connection_type),
    sprintf("   GO level range: %s to %s", level_from, level_to),
    sprintf("   Label type: %s", label_type),
    sprintf("   Threshold weight: %s", threshold),
    sprintf("   Layout option: %s", layout_name)
  )

  log_file <- file.path(output_folder_path, "log_parameters.txt")
  writeLines(log_lines, con = log_file)



  message("[4/10] Reading gene list...")

  gene_df <- prepare_genes(genes_list_path, species_prefix, output_folder_path)

  message("[5/10] Reading GO annotation file...")

  go_df <- read_go_annotations(go_file_path, species_prefix)

  message("[6/10] Annotating genes with GO terms...")

  annotated_genes <- annotate_with_go_terms(gene_df, go_df, output_folder_path, species_prefix)

  message("[7/10] Checking annotation coverage...")
  missing <- annotated_genes %>%
    dplyr::filter(is.na(SYMBOL) | is.na(GENE_ID) | is.na(GO_ID)) %>%
    dplyr::pull(ENSEMBL) %>%
    unique()

  utils::write.table(missing, file = file.path(output_folder_path, "term_gene_table_not_mapped.txt"),
                     sep = "\t", row.names = FALSE, quote = FALSE)

  message("[8/10] Filtering GO terms by evidence and level...")
  filtered_data <- filter_go_annotations(annotated_genes, config_evidence_path, level_from, level_to)

  message("[9/10] Building network...")
  network <- build_network(
    filtered_data,
    gene_df,
    duplicated_symbols = if (exists("duplicated_symbols")) duplicated_symbols else NULL,
    label_type = label_type,
    layout_name = layout_name,
    threshold = threshold,
    output_folder_path = output_folder_path,
    connection_type = connection_type
  )

  final_path <- file.path(output_folder_path, "gene_interaction_output_table_for_cytoscape.csv")

  if (label_type == "SYMBOL") {
    sym1 <- filtered_data$SYMBOL[match(network$final_table$gene1, filtered_data$ENSEMBL)]
    sym2 <- filtered_data$SYMBOL[match(network$final_table$gene2, filtered_data$ENSEMBL)]
    network$final_table$gene1 <- sym1
    network$final_table$gene2 <- sym2
  } else if (label_type == "both") {
    sym1 <- filtered_data$SYMBOL[match(network$final_table$gene1, filtered_data$ENSEMBL)]
    sym2 <- filtered_data$SYMBOL[match(network$final_table$gene2, filtered_data$ENSEMBL)]
    network$final_table$gene1 <- paste0(sym1, " (", network$final_table$gene1, ")")
    network$final_table$gene2 <- paste0(sym2, " (", network$final_table$gene2, ")")
  } else if (label_type == "SYMBOLalias") {
    symbol_dict <- filtered_data %>%
      dplyr::filter(!is.na(ENSEMBL), !is.na(SYMBOL)) %>%
      dplyr::distinct(ENSEMBL, SYMBOL) %>%
      dplyr::mutate(
        is_dup = SYMBOL %in% duplicated_symbols,
        display_label = ifelse(is_dup, paste0(SYMBOL, " (", ENSEMBL, ")"), SYMBOL)
      )
    network$final_table$gene1 <- symbol_dict$display_label[match(network$final_table$gene1, symbol_dict$ENSEMBL)]
    network$final_table$gene2 <- symbol_dict$display_label[match(network$final_table$gene2, symbol_dict$ENSEMBL)]
  }

  utils::write.table(network$final_table, file = final_path, sep = "\t", quote = FALSE, row.names = FALSE)

  message("[10/10] Results saved.")
  message("Pipeline completed successfully!\n → Output folder: ", output_folder_path)
  invisible()
}
