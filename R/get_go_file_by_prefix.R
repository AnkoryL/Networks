#' Locate internal GO annotation file
#'
#' @param species_prefix Short species code: human, mouse, zebrafish, macaque
#' @param go_type GO category prefix: BP, MF, CC, IP
#'
#' @return Full system path to matching GO annotation file
#' @export
get_go_file_by_prefix <- function(species_prefix, go_type) {
  type_map <- c(
    BP = "BiologicalProcess",
    MF = "MolecularFunction",
    CC = "CellularComponent",
    IP = "ImmuneSystemProcess"
  )

  species_map <- c(
    human = "HomoSapiens",
    mouse = "MusMusculus",
    macaque = "Macacamulatta",
    zebrafish = "Daniorerio"
  )

  if (!go_type %in% names(type_map)) {
    stop(error_messages$unsupported_go_type)
  }
  if (!species_prefix %in% names(species_map)) {
    stop(error_messages$unsupported_species_prefix)
  }

  full_type <- type_map[[go_type]]
  full_species <- species_map[[species_prefix]]

  annotation_dir <- system.file("extdata/annotations", package = "Networks")
  if (annotation_dir == "") {
    stop(error_messages$not_found_annotation_directory)
  }

  all_files <- list.files(annotation_dir, full.names = TRUE)

  pattern <- paste0("(?i)", full_species, ".*", full_type, ".*\\.txt(\\.gz)?$")
  matched <- grep(pattern, all_files, value = TRUE)

  if (length(matched) == 0) {
    stop(sprintf(
      error_messages$not_found_annotation_file,
      species_prefix, full_species, full_type
    ))
  } else if (length(matched) > 1) {
    warning(sprintf("Multiple matching annotation files found. Using: %s", basename(matched[1])))
  }

  return(matched[1])
}
