#' Annotate genes with GO terms from parsed cluego file
#' @export
annotate_with_go_terms <- function(gene_df, go_df, output_path) {
  annotated <- dplyr::left_join(gene_df, go_df, by = "ENSEMBL")
  utils::write.table(annotated, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  annotated
}
