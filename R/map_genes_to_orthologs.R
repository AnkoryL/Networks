#' Map input genes to orthologs using a precompiled ortholog table
#' @param gene_df Data frame of input genes. Must contain a column `Ensembl_gene_id`.
#' @param ortholog_df Data frame of orthologs. Must contain at least two columns:
#'
#' @return A list with two elements:
#'   - `mapped_genes` — data frame with columns `Ensembl_gene_id` and `Gene_symbol` ready for downstream analysis
#'   - `orthologs` — the filtered ortholog table with original and ortholog columns preserved

#' @export
map_genes_to_orthologs <- function(gene_df, ortholog_df, output_folder_path) {

  homolog_ensembl_col <- grep("_homolog_ensembl_gene$", colnames(ortholog_df), value = TRUE)
  homolog_symbol_col <- grep("_homolog_associated_gene_name$", colnames(ortholog_df), value = TRUE)

  if (length(homolog_ensembl_col) != 1 || length(homolog_symbol_col) != 1) {
    stop("not found orthologs")
  }

  orthologs <- ortholog_df[
    ortholog_df$ensembl_gene_id != "" &
      ortholog_df[[homolog_ensembl_col]] != "",
  ]

  orthologs <- orthologs[orthologs$ensembl_gene_id %in% gene_df$Ensembl_gene_id, ]

  out_file <- file.path(output_folder_path, "orthologs_list_filtered_by_input_genes.txt")
  write.table(
    orthologs,
    file = out_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  openxlsx::write.xlsx(
    orthologs,
    file = file.path(output_folder_path, "orthologs_list_filtered_by_input_genes.xlsx")
    )


  mapped_for_analysis <- data.frame(
    Ensembl_gene_id = orthologs[[homolog_ensembl_col]],
    Gene_symbol = ifelse(orthologs[[homolog_symbol_col]] == "", NA,
                         orthologs[[homolog_symbol_col]]),
    stringsAsFactors = FALSE
  )

  return(mapped_for_analysis)
}
