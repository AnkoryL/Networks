#' Save unmapped genes to file
#' @export
check_annotation_coverage <- function(annotated_genes, genes_list, output_path) {
  genes_list <- unique(genes_list)
  annotated_filtered <- annotated_genes %>%
    dplyr::filter(ENSEMBL %in% genes_list)

  missing <- annotated_filtered %>%
    dplyr::filter(is.na(SYMBOL) | is.na(GENE_ID) | is.na(GO_ID)) %>%
    dplyr::pull(ENSEMBL) %>%
    na.omit() %>%
    unique()

  if (length(missing) > 0) {
    utils::write.table(missing, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    file.create(output_path)
  }

  list(
    total = length(genes_list),
    unmapped = length(missing),
    unmapped_ids = missing
  )
}
