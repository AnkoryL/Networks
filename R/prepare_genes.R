#' Prepare gene list from input file
#'
#' @param genes_list_path Path to gene list file.
#' @param species_prefix Expected species.
#' @param output_folder_path Path for outputs (e.g., duplicated symbols).
#'
#' @return A data frame with Ensembl gene IDs.
#' @export
prepare_genes <- function(genes_list_path, species_prefix, output_folder_path) {
  genes_raw <- read_input_file(genes_list_path)

  if ("Ensembl_gene_id" %in% colnames(genes_raw)) {
    message("      Detected Ensembl_gene_id column — using it directly.")
    gene_df <- data.frame(ENSEMBL = unlist(genes_raw$Ensembl_gene_id), stringsAsFactors = FALSE)

  } else if ("Gene_symbol" %in% colnames(genes_raw)) {
    message("      Mapping Gene_symbol to Ensembl...")

    annotation_db <- switch(
      species_prefix,
      "human"     = org.Hs.eg.db,
      "mouse"     = org.Mm.eg.db,
      "macaque"   = org.Mmu.eg.db,
      "zebrafish" = org.Dr.eg.db,
      stop("Unsupported species prefix.")
    )

    suppressMessages(annotations <- AnnotationDbi::select(
      annotation_db,
      keys = AnnotationDbi::keys(annotation_db),
      columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
    ))

    data <- data.frame(Gene_symbol = genes_raw$Gene_symbol)

    mapped <- dplyr::left_join(
      data,
      annotations[, c("ENSEMBL", "SYMBOL", "ENTREZID")],
      by = c("Gene_symbol" = "SYMBOL")
    )

    duplicated_map <- mapped %>%
      filter(!is.na(ENSEMBL)) %>%
      group_by(Gene_symbol) %>%
      summarise(
        n_ensembl = dplyr::n_distinct(ENSEMBL),
        ensembl_ids = paste(unique(ENSEMBL), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_ensembl > 1)

    if (nrow(duplicated_map) > 0) {
      write.table(
        duplicated_map,
        file = file.path(output_folder_path, "duplicated_symbol_mappings.txt"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }

    gene_df <- mapped %>%
      filter(!is.na(ENSEMBL)) %>%
      dplyr::select(ENSEMBL) %>%
      dplyr::distinct()

  } else {
    stop("Input must contain either 'Ensembl_gene_id' or 'Gene_symbol' column.")
  }

  inferred <- infer_species(gene_df$ENSEMBL[1])
  if (is.na(inferred)) stop("Could not infer species from Ensembl ID.")
  if (inferred != species_prefix) stop(sprintf(
    "Species mismatch: expected '%s', got '%s'", species_prefix, inferred
  ))

  message(sprintf("      Species confirmed: %s", inferred))

  gene_df
}
