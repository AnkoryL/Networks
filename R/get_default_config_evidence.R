#' Get default path to config_evidence.txt
#' @return Path to internal config_evidence.txt
#' @export
default_evidence_experimental_only <- function() {
  system.file("extdata/templates/config_evidence.txt", package = "Networks")
}
