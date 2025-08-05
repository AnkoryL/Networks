#' Run the full GO-based gene interaction pipeline
#'
#' @param config_path Path to the main_config.txt
#'
#' @return NULL
#' @export
run_pipeline <- function(config_path) {
  message("[1/10] Starting pipeline...")

  if (!file.exists(config_path)) stop("main_config.txt not found at given path.")

  message("[2/10] Reading configuration...")
  config <- read_config_txt(config_path)

  species_prefix        <- config$special_prefix
  genes_list_path       <- config$genes_list_path
  go_file_path          <- config$go_file_path
  config_evidence_path  <- config$config_evidence_path
  output_folder_path    <- config$output_folder_path
  connection_type       <- config$connection_type
  level_from            <- config$level_from
  level_to              <- config$level_to
  label_type            <- config$label_type
  layout_name           <- config$layout
  threshold             <- config$threshold

  message("[3/10] Reading gene list...")
  genes_input <- read_input_file(genes_list_path)

  if ("Ensembl_gene_id" %in% colnames(genes_input)) {
    genes <- genes_input$Ensembl_gene_id
    gene_df <- data.frame(ENSEMBL = genes, stringsAsFactors = FALSE)
  } else if ("Gene_symbol" %in% colnames(genes_input)) {
    message("Gene symbols detected. Annotating to Ensembl...")
    symbol_result <- annotate_symbols(genes_input, species_prefix, output_folder_path)
    gene_df <- symbol_result$gene_df
    duplicated_symbols <- symbol_result$duplicated_symbols
  } else {
    stop("Input must contain either Ensembl_gene_id or Gene_symbol column.")
  }

  message("[4/10] Reading GO annotation file...")
  go_lines <- readLines(go_file_path)
  go_df <- do.call(rbind, lapply(go_lines, parse_go_line))
  go_df <- go_df[-1, ]

  entrez_keys <- unique(go_df$GENE_ID)
  annotation_db <- get_annotation_db(species_prefix)

  suppressMessages(go_annots <- AnnotationDbi::select(
    annotation_db,
    keys = entrez_keys,
    columns = c("SYMBOL", "ENSEMBL"),
    keytype = "ENTREZID"
  ))

  go_annotated <- dplyr::left_join(go_df, go_annots, by = c("GENE_ID" = "ENTREZID"))

  message("[5/10] Annotating genes with GO terms...")
  annotated_genes <- dplyr::left_join(gene_df, go_annotated, by = "ENSEMBL")

  out_term_path <- file.path(output_folder_path, paste0(species_prefix, "_term_gene_table.csv"))
  utils::write.table(annotated_genes, file = out_term_path, sep = "\t", row.names = FALSE, quote = FALSE)

  message("[6/10] Checking annotation coverage...")
  missing <- annotated_genes %>%
    dplyr::filter(is.na(SYMBOL) | is.na(GENE_ID) | is.na(GO_ID)) %>%
    dplyr::pull(ENSEMBL) %>%
    unique()

  utils::write.table(missing, file = file.path(output_folder_path, "term_gene_table_not_mapped.txt"),
                     sep = "\t", row.names = FALSE, quote = FALSE)

  message("[7/10] Filtering GO terms by evidence and level...")
  filtered_data <- filter_go_annotations(annotated_genes, config_evidence_path, level_from, level_to)

  message("[8/10] Building network...")
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

  message("[9/10] Results saved.")
  message("[10/10] Pipeline completed successfully!\n → Output folder: ", output_folder_path)
  invisible()
}
