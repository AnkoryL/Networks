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
    stop("Unknown go_type. Use one of: BP, MF, CC, IP")
  }
  if (!species_prefix %in% names(species_map)) {
    stop("Unknown species_prefix. Use one of: human, mouse, macaque, zebrafish")
  }

  full_type <- type_map[[go_type]]
  full_species <- species_map[[species_prefix]]

  annotation_dir <- system.file("extdata/annotations", package = "Networks")
  if (annotation_dir == "") {
    stop("Annotation directory not found in the package.")
  }

  all_files <- list.files(annotation_dir, full.names = TRUE)

  pattern <- paste0("(?i)", full_species, ".*", full_type, ".*\\.txt(\\.gz)?$")
  matched <- grep(pattern, all_files, value = TRUE)

  if (length(matched) == 0) {
    stop(sprintf(
      "No annotation file found for species = '%s' (%s), type = '%s'",
      species_prefix, full_species, full_type
    ))
  } else if (length(matched) > 1) {
    warning(sprintf("Multiple matching annotation files found. Using: %s", basename(matched[1])))
  }

  return(matched[1])
}
