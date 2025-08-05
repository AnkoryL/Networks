#' Save unmapped genes to file
#' @export
check_annotation_coverage <- function(annotated_genes, genes_list, output_path) {
  missing <- annotated_genes %>%
    dplyr::filter(is.na(SYMBOL) | is.na(GENE_ID) | is.na(GO_ID)) %>%
    dplyr::pull(ENSEMBL) %>%
    na.omit() %>%
    unique()

  utils::write.table(missing, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

  list(
    total = length(genes_list),
    unmapped = length(missing),
    unmapped_ids = missing
  )
}
