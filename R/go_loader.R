#' Read and annotate GO data
#'
#' @param go_file_path Path to GO annotation file.
#' @param species_prefix Species prefix for mapping.
#'
#' @return A data frame with GO + SYMBOL + ENSEMBL.
#' @export
read_go_annotations <- function(go_file_path, species_prefix) {
  go_data <- readLines(go_file_path)
  parsed <- do.call(rbind, lapply(go_data, parse_go_line))
  parsed <- parsed[-1, ]  # skip header or dummy row

  unique_keys <- unique(parsed$GENE_ID)

  annotation_db <- switch(
    species_prefix,
    "human"     = org.Hs.eg.db,
    "mouse"     = org.Mm.eg.db,
    "macaque"   = org.Mmu.eg.db,
    "zebrafish" = org.Dr.eg.db
  )

  suppressMessages(go_annotations <- AnnotationDbi::select(
    annotation_db,
    keys = unique_keys,
    columns = c("SYMBOL", "ENSEMBL"),
    keytype = "ENTREZID"
  ))

  dplyr::left_join(parsed, go_annotations, by = c("GENE_ID" = "ENTREZID"))
}
