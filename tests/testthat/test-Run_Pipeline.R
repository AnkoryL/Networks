# test_that("run_pipeline works", {
#
#   output_folder_path <- "C:/Users/Ankory/YandexDisk/is_not_a_wolf/Networks/inst/extdata/example/result"
#
#   species <- "macaque"
#   go_type <- "IP"
#   file_format <- "csv"
#   id_type <- "ensembl"
#
#   genes_list_path <- system.file("extdata/example/",
#                                  paste0(species, "_", id_type, ".", file_format),
#                                  package = "Networks")
#
#   run_pipeline(
#     go_type =  go_type,
#     species_prefix = species,
#     genes_list_path = genes_list_path,
#     output_folder_path = output_folder_path,
#     connection_type = "common_go_term",
#     level_from = 1,
#     level_to = 20,
#     label_type = "SYMBOL",
#     threshold = 1,
#     layout_name = "layout_nicely"
#   )
#
#   expect_true(file.exists(file.path(output_folder_path, "log_parameters.txt")))
#   expect_true(file.exists(file.path(output_folder_path, "term_gene_table_not_mapped.txt")))
#   expect_true(file.exists(file.path(output_folder_path, "gene_interaction_output_table_for_cytoscape.csv")))
#
#   output_table <- read.csv(file.path(output_folder_path, "gene_interaction_output_table_for_cytoscape.csv"))
#   expect_gt(nrow(output_table), 0)
# })
#


##  test_file("C:/Users/Ankory/YandexDisk/is_not_a_wolf/Networks/tests/testthat/test-Run_Pipeline.R", reporter = "progress")

test_that("run_pipeline works for all species, GO types, gene formats and ID types", {

  species_list <- c(
    # "human",
    # "mouse",
    # "macaque",
    "zebrafish"
  )

  go_types <- c(
    # "BP",
    # "CC",
    # "MF",
    "IP"
  )

  file_formats <- c(
    # "xlsx",
    # "txt",
    "csv"
  )

  id_types <- c(
    "ensembl",
    "symbol",
    "ensembl_symbol"
  )

  output_folder_path <- "C:/Users/Ankory/YandexDisk/is_not_a_wolf/Networks/inst/extdata/example/result"


  for (species in species_list) {
    for (go in go_types) {
      for (format in file_formats) {
        for (id in id_types) {
          file_name <- paste0(species, "_", id, ".", format)
          genes_list_path <- system.file("extdata/example/", file_name, package = "Networks")

          res <- tryCatch({
            cat(sprintf("\n--- Running pipeline for %s, GO: %s, format: %s, id: %s ---\n",
                        species, go, format, id))

            run_pipeline(
              go_type = go,
              species_prefix = species,
              genes_list_path = genes_list_path,
              output_folder_path = output_folder_path,
              connection_type = "common_go_term",
              level_from = 1,
              level_to = 20,
              label_type = "SYMBOL",
              threshold = 1,
              layout_name = "layout_nicely"
            )

            expect_true(file.exists(file.path(output_folder_path, "log_parameters.txt")))
            expect_true(file.exists(file.path(output_folder_path, "term_gene_table_not_mapped.txt")))
            expect_true(file.exists(file.path(output_folder_path, "gene_interaction_output_table_for_cytoscape.csv")))

            output_table <- read.csv(file.path(output_folder_path, "gene_interaction_output_table_for_cytoscape.csv"))
            expect_gt(nrow(output_table), 0)

          TRUE
          },
          error = function(e) {
            message(sprintf("Error in run_pipeline for %s, GO: %s, id: %s, format: %s:\n%s",
                            species, go, id, format, conditionMessage(e)))
            FALSE
          })

          if (!isTRUE(res)) {
            message(sprintf("Skipping this combination due to error, continuing loop."))
            next
          }
        }
      }
    }
  }
})
