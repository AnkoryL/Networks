#' Filter GO annotations by level and evidence
#' @export
#'
filter_go_annotations <- function(annotated_genes, config_evidence_path, level_from, level_to) {
  evidence_conf <- read.delim(config_evidence_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  valid_evidence <- evidence_conf$Evidence[evidence_conf$Boolean == "TRUE"]

  df <- annotated_genes %>%
    dplyr::mutate(Level = strsplit(as.character(Level), ",")) %>%
    dplyr::mutate(Level = lapply(Level, as.numeric)) %>%
    dplyr::filter(sapply(Level, function(x) any(x >= level_from & x <= level_to))) %>%
    dplyr::mutate(Level = sapply(Level, function(x) paste(x, collapse = ","))) %>%
    dplyr::filter(EVIDENCE %in% valid_evidence) %>%
    dplyr::select(GO_ID, ENSEMBL, EVIDENCE, SYMBOL)

  df
}
