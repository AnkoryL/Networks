#' Parse a single GO annotation line
#'
#' @param line A character string from the GO annotation file.
#'
#' @return A data frame with columns: GO_ID, Level, Category, GENE_ID, EVIDENCE.
#' @keywords internal
#' @export
parse_go_line <- function(line) {
  parts <- unlist(strsplit(line, "\t"))

  go_id <- parts[1]
  level <- parts[2]
  category <- parts[3]
  genes_with_evidence <- strsplit(parts[4], "\\|")[[1]]

  data.frame(
    GO_ID = rep(go_id, length(genes_with_evidence)),
    Level = rep(level, length(genes_with_evidence)),
    Category = rep(category, length(genes_with_evidence)),
    GENE_ID = gsub("^[^:]+:", "", genes_with_evidence),
    EVIDENCE = gsub(":.*$", "", genes_with_evidence),
    stringsAsFactors = FALSE
  )
}
