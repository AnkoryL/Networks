#' Read a gene list file in various formats
#'
#' @param path Path to a gene list file (csv, txt, tsv, xls, xlsx).
#'
#' @return A data frame with gene identifiers.
#' @export
read_input_file <- function(path) {
  file_ext <- tools::file_ext(path)

  if (file_ext %in% c("csv")) {
    message("      Gene-list is CSV file!")
    return(read.csv(path, sep = "\t", stringsAsFactors = FALSE))
  } else if (file_ext %in% c("txt", "tsv")) {
    message("      Gene-list is TXT/TSV file!")
    return(read.delim(path, stringsAsFactors = FALSE))
  } else if (file_ext %in% c("xls", "xlsx")) {
    message("      Gene-list is Excel file!")
    return(openxlsx::read.xlsx(path))
  } else {
    stop(error_messages$unsupported_genes_file_format)
  }
}
