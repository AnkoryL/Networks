#' Load and validate configuration for gene network analysis
#'
#' @param config_path Path to a main_config.txt file.
#'
#' @return A named list of validated configuration parameters.
#' @export
load_and_validate_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(paste0("File not found: ", config_path))
  }

  config <- read_config_txt(config_path)

  required_keys <- c(
    "genes_list_path", "species_prefix",  "go_type", "go_file_path", "config_evidence_path",
    "output_folder_path", "connection_type", "level_from", "level_to",
    "label_type", "threshold", "layout"
  )

  missing <- setdiff(required_keys, names(config))
  if (length(missing) > 0) {
    stop(error_messages$missing_config_parameters, paste(missing, collapse = ", "))
  }

  if (!(config$label_type %in% c("SYMBOL", "ENSEMBL", "both"))) {
    stop(error_messages$invalid_label_type)
  }

  message("      √ Config loaded successfully.")

  config
}
