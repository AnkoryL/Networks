#' Create working directory with example config files
#'
#' @param path Path to create the config project folder. Defaults to "my_project".
#' @param copy_example_configs Logical; if TRUE, copies example config files into the project folder. Default is FALSE.
#' @param create_empty_gene_list Logical; if TRUE, creates an empty gene list file inside the project folder. Default is FALSE.
#' @export
create_project_folder <- function(path = "my_project", copy_example_configs = FALSE, create_empty_gene_list = FALSE) {
  # if (dir.exists(path)) stop("Folder already exists. Choose another name or delete the existing one.")

  dir.create(path, recursive = TRUE)
  dir.create(file.path(path, "results"))

  # Copy configs
  main_config_template <- system.file("extdata/templates", "main_config.txt", package = "Networks")
  evidence_config_template <- system.file("extdata/templates", "config_evidence.txt", package = "Networks")

  if (!file.exists(main_config_template) || !file.exists(evidence_config_template)) {
    stop("Internal config templates not found in the package.")
  }

  if (copy_example_configs == TRUE) {
  file.copy(main_config_template, file.path(path, "main_config.txt"))
  file.copy(evidence_config_template, file.path(path, "config_evidence.txt"))
  }

  # (Optional) create empty gene list
  if (create_empty_gene_list == TRUE) {
  gene_list_path <- file.path(path, "my_genes.csv")
  readr::write_csv(data.frame(Ensembl_gene_id = character()), gene_list_path)
  }

  message("Project created at: ", normalizePath(path))

}
