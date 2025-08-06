#' Prepare gene list from input file
#'
#' @param genes_list_path Path to gene list file.
#' @param species_prefix Expected species.
#' @param output_folder_path Path for outputs (e.g., duplicated symbols).
#'
#' @return A data frame with Ensembl gene IDs.
#' @export
prepare_genes <- function(genes_list_path, species_prefix, output_folder_path) {

    genes_list <- read_input_file(genes_list_path)
  # genes_list_vector <- unlist(genes_list[, 1])  # Assuming first column has Ensembl IDs

  if ("Ensembl_gene_id" %in% colnames(genes_list)) {
    message("      Detected Ensembl_gene_id column — using it directly.")
    genes_list <- unlist(genes_list$Ensembl_gene_id)
    gene_df <- data.frame(ENSEMBL = genes_list, row.names = NULL, stringsAsFactors = FALSE)

    infer_species <- function(ensembl_id) {
      dplyr::case_when(
        grepl("^ENSG",    ensembl_id, ignore.case = TRUE) ~ "human",
        grepl("^ENSMUSG", ensembl_id, ignore.case = TRUE) ~ "mouse",
        grepl("^ENSMMUG", ensembl_id, ignore.case = TRUE) ~ "macaque",
        grepl("^ENSDARG", ensembl_id, ignore.case = TRUE) ~ "zebrafish",
        TRUE ~ NA_character_
      )
    }

    inferred_species <- infer_species(gene_df[1,])


    if (is.na(inferred_species)) {
      stop("   Species could not be automatically inferred from input. Please check the matching of your Ensembl IDs, the entered species, and the GO annotation.")
    } else {
      message(sprintf("      Automatically defined species: %s \n", inferred_species))
    }

    if (species_prefix!=inferred_species) {
      stop(sprintf(
        "   The entered species ('%s') differs from the automatically determined one ('%s')",
        species_prefix, inferred_species
      ))
    }


  } else if ("Gene_symbol" %in% colnames(genes_list)) {
    message("      Ensembl_gene_id column not found — attempting to annotate Gene_symbol to Ensembl.")

    #label_type <- "SYMBOLboth"
    label_type <- "SYMBOLalias"

    annotation_db <- switch(
      species_prefix,
      "human"     = org.Hs.eg.db::org.Hs.eg.db,
      "mouse"     = org.Mm.eg.db::org.Mm.eg.db,
      "macaque"   = org.Mmu.eg.db::org.Mmu.eg.db,
      "zebrafish" = org.Dr.eg.db::org.Dr.eg.db,
      stop("Unsupported species prefix.")
    )

    suppressMessages(annotations <- AnnotationDbi::select(
      annotation_db,
      keys = AnnotationDbi::keys(annotation_db),
      columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
    ))

    data <- data.frame(Gene_symbol = genes_list$Gene_symbol)

    annotationdbi_data_full <- dplyr::left_join(
      data,
      annotations[, c("ENSEMBL", "SYMBOL", "ENTREZID")],
      by = c("Gene_symbol" = "SYMBOL")
    )

    # Check Ensembl duplicates
    duplicated_symbol_mappings <- annotationdbi_data_full %>%
      dplyr::filter(!is.na(ENSEMBL)) %>%
      dplyr::group_by(Gene_symbol) %>%
      dplyr::summarise(
        n_ensembl = dplyr::n_distinct(ENSEMBL),
        ensembl_ids = paste(unique(ENSEMBL), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_ensembl > 1)

    duplicated_symbols <- duplicated_symbol_mappings$Gene_symbol

    output_duplicated_file_path <- file.path(output_folder_path, "duplicated_symbol_mappings.txt")
    write.table(
      duplicated_symbol_mappings,
      file = output_duplicated_file_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    ## Delete Ensembl duplicates
    #annotationdbi_data_full <- annotationdbi_data_full %>%
    # filter(!is.na(ENSEMBL)) %>%
    # group_by(Gene_symbol) %>%
    # slice(1) %>%
    # ungroup()


    genes_list <- annotationdbi_data_full %>%
      dplyr::filter(!is.na(ENSEMBL)) %>%
      dplyr::rename(Ensembl_gene_id = ENSEMBL)

    message("      Gene symbol -> Ensembl      completed.")

    genes_list <- unlist(genes_list$Ensembl_gene_id)
    gene_df <- data.frame(ENSEMBL = genes_list, row.names = NULL, stringsAsFactors = FALSE)

    infer_species <- function(ensembl_id) {
      dplyr::case_when(
        grepl("^ENSG",    ensembl_id, ignore.case = TRUE) ~ "human",
        grepl("^ENSMUSG", ensembl_id, ignore.case = TRUE) ~ "mouse",
        grepl("^ENSMMUG", ensembl_id, ignore.case = TRUE) ~ "macaque",
        grepl("^ENSDARG", ensembl_id, ignore.case = TRUE) ~ "zebrafish",
        TRUE ~ NA_character_
      )

    }

    inferred_species <- infer_species(gene_df[1,])

    if (is.na(inferred_species)) {
      stop("   Species could not be automatically inferred from input. Please check the matching of your Ensembl IDs, the entered species, and the GO annotation.")
    } else {
      message(sprintf("      Automatically defined species: %s\n", inferred_species))
    }

    if (species_prefix!=inferred_species) {
      stop(sprintf(
        "\n   The entered species ('%s') differs from the automatically determined one ('%s').\n",
        species_prefix, inferred_species
      ))
    }

  } else {
    stop("\n   Input must contain either 'Ensembl_gene_id' or 'Gene_symbol' column.\n")
  }
    return(gene_df)
}

