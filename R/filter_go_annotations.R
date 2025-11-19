#' Filter GO annotations by level and evidence, returning a cleaned table for network building
#'
#' @param annotated_genes Data.frame with columns: GO_ID, ENSEMBL, EVIDENCE, SYMBOL, Level (comma-separated string)
#' @param config_evidence_path Path to evidence configuration file
#' @param level_from Minimum GO level to include
#' @param level_to Maximum GO level to include
#'
#' @return Data.frame with columns GO_ID, Ensembl_gene_id , EVIDENCE, SYMBOL filtered by level and evidence
filter_go_annotations <- function(annotated_genes, config_evidence_path, level_from, level_to) {
  config_evidence <- read.delim(config_evidence_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  valid_evidence <- config_evidence$Evidence[config_evidence$Boolean == "TRUE"]

  data <- na.omit(annotated_genes)

  data$Level <- strsplit(as.character(data$Level), ",")
  data$Level <- lapply(data$Level, function(x) as.numeric(trimws(x)))

  valid_level <- function(x) any(x >= level_from & x <= level_to)
  data <- data[sapply(data$Level, valid_level), ]

  data$Level <- sapply(data$Level, function(x) paste(x, collapse = ","))

  data <- data[data$EVIDENCE %in% valid_evidence, ]

  filtered_data <- data[, c("GO_ID", "Ensembl_gene_id", "EVIDENCE", "SYMBOL")]

  return(filtered_data)
}
