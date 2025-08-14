#' Annotate genes with GO terms from parsed cluego file
#' @export
annotate_with_go_terms <- function(gene_df, go_df, output_folder_path, species_prefix) {
  annotated <- dplyr::left_join(gene_df, go_df, by = "ENSEMBL")
  out_term_path <- file.path(output_folder_path, paste0(species_prefix, "_term_gene_table.csv"))
  utils::write.table(annotated, file = out_term_path, sep = "\t", row.names = FALSE, quote = FALSE)
  annotated
}
