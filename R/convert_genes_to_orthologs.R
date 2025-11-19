#' Convert a gene list into orthologs using internal ortholog tables
#'
#' @param genes_list_path Path to a file with input genes.
#'   Must contain either `Ensembl_gene_id` or `Gene_symbol`.
#' @param from_species Species of the input genes ("human", "mouse", "macaque", "zebrafish").
#' @param to_species Species to convert into ("human", "mouse", "macaque").
#' @param output_folder_path Folder where results will be saved.
#'
#' @details
#' This function loads internal ortholog tables, filters out rows without valid
#' original/ortholog Ensembl IDs, keeps only the rows corresponding to the input
#' gene list, and constructs a ready-to-use table of ortholog Ensembl IDs and symbols.
#'
#' Resulting tables are saved as TXT files in `output_folder_path`.
#'
#' @export
convert_genes_to_orthologs <- function(
    genes_list_path,
    from_species,
    to_species,
    output_folder_path = "./results_orthologs"
) {
  message("[1/4] Reading input genes...")
  dir.create(output_folder_path, showWarnings = FALSE, recursive = TRUE)

  gene_df <- prepare_genes(
    genes_list_path = genes_list_path,
    species_prefix = from_species,
    output_folder_path = output_folder_path
  )

  message("[2/4] Loading ortholog table...")
  ortholog_df <- load_ortholog(from_species = from_species, to_species = to_species)

  message("[3/4] Mapping genes to orthologs...")
  mapping <- map_genes_to_orthologs(gene_df, ortholog_df,output_folder_path)

  message("[4/4] Saving results...")

  out_map <- file.path(
    output_folder_path,
    paste0(from_species, "_to_", to_species, "_ortholog_list.txt")
  )
  write.table(
    mapping,
    file = out_map,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  openxlsx::write.xlsx(
    mapping,
    file = file.path(
      output_folder_path,
      paste0(from_species, "_to_", to_species, "_ortholog_list.xlsx")
    )
  )

  message("Conversion completed.")

}
