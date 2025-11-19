#' Prepare gene list from input file
#'
#' @param genes_list_path Path to gene list file.
#' @param species_prefix Expected species.
#' @param output_folder_path Path for outputs (e.g., duplicated symbols).
#'
#' @return A data frame with Ensembl gene IDs.
#' @export
#'

prepare_genes <- function(genes_list_path, species_prefix, output_folder_path) {

  genes_list <- read_input_file(genes_list_path)

  infer_species <- function(ensembl_id) {
    dplyr::case_when(
      grepl("^ENSG",    ensembl_id, ignore.case = TRUE) ~ "human",
      grepl("^ENSMUSG", ensembl_id, ignore.case = TRUE) ~ "mouse",
      grepl("^ENSMMUG", ensembl_id, ignore.case = TRUE) ~ "macaque",
      grepl("^ENSDARG", ensembl_id, ignore.case = TRUE) ~ "zebrafish",
      TRUE ~ NA_character_
    )
  }

  annotation_db <- switch(
    species_prefix,
    "human"     = org.Hs.eg.db::org.Hs.eg.db,
    "mouse"     = org.Mm.eg.db::org.Mm.eg.db,
    "macaque"   = org.Mmu.eg.db::org.Mmu.eg.db,
    "zebrafish" = org.Dr.eg.db::org.Dr.eg.db,
    stop(error_messages$unsupported_species_prefix)
  )

  suppressMessages(
    annotations <- AnnotationDbi::select(
      annotation_db,
      keys = AnnotationDbi::keys(annotation_db),
      columns = c("ENSEMBL", "SYMBOL")
    )
  )

  if ("Ensembl_gene_id" %in% colnames(genes_list) &&
      "Gene_symbol" %in% colnames(genes_list)) {

    message("      Detected both Ensembl_gene_id and Gene_symbol — using input as is.")

    df <- genes_list[, c("Ensembl_gene_id", "Gene_symbol")]

    inferred <- infer_species(df$Ensembl_gene_id[1])
    if (is.na(inferred)) stop(error_messages$species_not_automatically_inferred)
    if (species_prefix != inferred)
      stop(sprintf(error_messages$species_diffes, species_prefix, inferred))

    return(df)
  }

  if ("Ensembl_gene_id" %in% colnames(genes_list)) {

    message("      Detected Ensembl_gene_id column — will enrich with Gene_symbol.")

    vec <- unlist(genes_list$Ensembl_gene_id)

    inferred <- infer_species(vec[1])
    if (is.na(inferred)) stop(error_messages$species_not_automatically_inferred)
    if (species_prefix != inferred)
      stop(sprintf(error_messages$species_diffes, species_prefix, inferred))

    df <- data.frame(Ensembl_gene_id = vec, stringsAsFactors = FALSE)

    df <- dplyr::left_join(
      df,
      annotations[, c("ENSEMBL", "SYMBOL")],
      by = c("Ensembl_gene_id" = "ENSEMBL")
    )

    names(df)[names(df) == "SYMBOL"] <- "Gene_symbol"

    out_gene_list_path <- file.path(output_folder_path, "double_index_gene_list.txt")
    write.table(
      df,
      file = out_gene_list_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    return(df)
  }

  if ("Gene_symbol" %in% colnames(genes_list)) {

    message("      No Ensembl IDs found — mapping Gene_symbol → Ensembl.")

    symbols <- genes_list$Gene_symbol

    annotation_full <- dplyr::left_join(
      data.frame(Gene_symbol = symbols),
      annotations[, c("ENSEMBL", "SYMBOL")],
      by = c("Gene_symbol" = "SYMBOL")
    )

    duplicated_symbol_mappings <- annotation_full %>%
      dplyr::filter(!is.na(ENSEMBL)) %>%
      dplyr::group_by(Gene_symbol) %>%
      dplyr::summarise(
        n_ensembl = dplyr::n_distinct(ENSEMBL),
        ensembl_ids = paste(unique(ENSEMBL), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_ensembl > 1)

    write.table(
      duplicated_symbol_mappings,
      file = file.path(output_folder_path, "duplicated_symbol_mappings.txt"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    annotation_full <- annotation_full %>%
      dplyr::filter(!is.na(ENSEMBL))

    inferred <- infer_species(annotation_full$ENSEMBL[1])
    if (is.na(inferred)) stop(error_messages$species_not_automatically_inferred)
    if (species_prefix != inferred)
      stop(sprintf(error_messages$species_diffes, species_prefix, inferred))

    df <- annotation_full %>%
      dplyr::rename(Ensembl_gene_id = ENSEMBL) %>%
      dplyr::select(Ensembl_gene_id, Gene_symbol)

    out_gene_list_path <- file.path(output_folder_path, "double_index_gene_list.txt")
    write.table(
      df,
      file = out_gene_list_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    message("      Gene_symbol → Ensembl completed.")

    return(df)
  }

  stop(error_messages$unsupported_genes_input)
}

