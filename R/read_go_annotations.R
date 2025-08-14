#' Read and annotate GO data
#'
#' @param go_file_path Path to GO annotation file.
#' @param species_prefix Species prefix for mapping.
#'
#' @return A data frame with GO + SYMBOL + ENSEMBL.
#' @export
read_go_annotations <- function(go_file_path, species_prefix) {
  go_data <- readLines(go_file_path)
  input_go_data <- do.call(rbind, lapply(go_data, parse_go_line))
  input_go_data <- input_go_data[-1, ]

  unique_keys <- unique(input_go_data$GENE_ID)

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
    keytype = "ENTREZID")
  )

  annotationdbi_go_data <- input_go_data %>%
    dplyr::left_join(go_annotations, by = c("GENE_ID" = "ENTREZID"),
              relationship = "many-to-many")

}
