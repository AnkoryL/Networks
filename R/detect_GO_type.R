#' Detect GO type from file name
#'
#' @param go_file_path Path to GO file.
#'
#' @return GO category prefix: BP, MF, CC, IP.
#' @export
detect_go_type <- function(go_file_path) {
  type <- dplyr::case_when(
    grepl("BiologicalProcess", go_file_path)     ~ "BP",
    grepl("CellularComponent", go_file_path)     ~ "CC",
    grepl("ImmuneSystemProcess", go_file_path)   ~ "IP",
    grepl("MolecularFunction", go_file_path)     ~ "MF"
  )
  if (is.na(type)) {
    stop("GO file name must contain a recognizable type (e.g., BiologicalProcess).")
  }
  type
}
