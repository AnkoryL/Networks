#' Read main configuration from a tab-delimited file
#'
#' @param config_path Path to a main_config.txt file.
#'
#' @return A named list of configuration parameters.
#' @export
read_config_txt <- function(config_path) {
  config <- read.delim(config_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  config <- config[!is.na(config$parameter) & config$parameter != "", ]
  config_list <- setNames(as.list(config$value), config$parameter)

  numeric_fields <- c("level_from", "level_to", "threshold")
  for (field in numeric_fields) {
    if (field %in% names(config_list)) {
      config_list[[field]] <- as.numeric(config_list[[field]])
    }
  }

  return(config_list)
}
