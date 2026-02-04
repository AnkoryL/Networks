#' Read a gene list file in various formats, auto-detecting extension
#'
#' @param path Path to the gene list file (can be name only; function will detect extension).
#' @param folder Folder where the file is located. Default: current working directory.
#'
#' @return A data frame with gene identifiers.
#' @export
read_input_file <- function(path, folder = ".") {

  if (dirname(path) != ".") {
    folder <- dirname(path)
    path_name <- basename(path)
  } else {
    path_name <- path
  }

  if (!grepl("\\.", path_name)) {
    exts <- c("csv", "txt", "tsv", "xls", "xlsx")
    files_found <- list.files(
      folder,
      pattern = paste0("^", path_name, "\\.(", paste(exts, collapse = "|"), ")$"),
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(files_found) == 0) {
      stop(
        "File '", path_name,
        "' with supported extension not found in folder ", folder
      )
    }

    if (length(files_found) > 1) {
      message("      Multiple files found, using: ", basename(files_found[1]))
    }

    path <- files_found[1]
  }

  file_ext <- tolower(tools::file_ext(path))

  if (file_ext == "csv") {
    message("      Gene-list is CSV file!")
    read.csv(path, sep = "\t", stringsAsFactors = FALSE)
  } else if (file_ext %in% c("txt", "tsv")) {
    message("      Gene-list is TXT/TSV file!")
    read.delim(path, stringsAsFactors = FALSE)
  } else if (file_ext %in% c("xls", "xlsx")) {
    message("      Gene-list is Excel file!")
    openxlsx::read.xlsx(path)
  } else {
    stop(error_messages$unsupported_genes_file_format)
  }
}
