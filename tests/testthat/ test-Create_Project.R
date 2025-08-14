test_that("create_project_folder creates expected structure", {
  tmp <- withr::local_tempdir()

  test_main <- file.path(tmp, "main_config.txt")
  test_evidence <- file.path(tmp, "config_evidence.txt")
  writeLines("fake main config", test_main)
  writeLines("fake evidence config", test_evidence)

  with_mocked_bindings(
    system.file = function(...) {
      args <- list(...)
      if (args[[2]] == "main_config.txt") return(test_main)
      if (args[[2]] == "config_evidence.txt") return(test_evidence)
      stop("Unexpected system.file call")
    },
    {
      project_path <- file.path(tmp, "test_project")
      create_project_folder(
        path = project_path,
        copy_example_configs = TRUE,
        create_empty_gene_list = TRUE
      )

      expect_true(dir.exists(project_path))
      expect_true(dir.exists(file.path(project_path, "results")))
      expect_true(file.exists(file.path(project_path, "main_config.txt")))
      expect_true(file.exists(file.path(project_path, "config_evidence.txt")))
      expect_true(file.exists(file.path(project_path, "my_genes.csv")))

      gene_list <- readr::read_csv(file.path(project_path, "my_genes.csv"), show_col_types = FALSE)
      expect_named(gene_list, "Ensembl_gene_id")
      expect_equal(nrow(gene_list), 0)
    }
  )
})
