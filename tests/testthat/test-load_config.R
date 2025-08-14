test_that("load_and_validate_config errors if required keys are missing", {
  tmp_file <- tempfile(fileext = ".txt")
  write.table(data.frame(parameter = c("genes_list_path", "species_prefix"),
                         value = c("some_path", "human")),
              file = tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    load_and_validate_config(tmp_file),
    regexp = "Missing required config parameters"
  )

  unlink(tmp_file)
})

test_that("load_and_validate_config errors if label_type invalid", {
  tmp_file <- tempfile(fileext = ".txt")
  params <- data.frame(
    parameter = c(
      "genes_list_path", "species_prefix",  "go_type", "go_file_path", "config_evidence_path",
      "output_folder_path", "connection_type", "level_from", "level_to",
      "label_type", "threshold", "layout"
    ),
    value = c(
      "genes.txt", "human", "BP", "go_file.txt", "config_evidence.txt",
      "out_folder", "common_go_term", "10", "20",
      "INVALID_LABEL", "5", "layout_nicely"
    )
  )
  write.table(params, file = tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    load_and_validate_config(tmp_file),
    regexp = "Invalid label_type"
  )

  unlink(tmp_file)
})

test_that("load_and_validate_config loads valid config correctly", {
  tmp_file <- tempfile(fileext = ".txt")
  params <- data.frame(
    parameter = c(
      "genes_list_path", "species_prefix",  "go_type", "go_file_path", "config_evidence_path",
      "output_folder_path", "connection_type", "level_from", "level_to",
      "label_type", "threshold", "layout"
    ),
    value = c(
      "genes.txt", "human", "BP", "go_file.txt", "config_evidence.txt",
      "out_folder", "common_go_term", "10", "20",
      "SYMBOL", "5", "layout_nicely"
    )
  )
  write.table(params, file = tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_message(
    res <- load_and_validate_config(tmp_file),
    regexp = "Config loaded successfully"
  )

  message("Returned config parameters:")
  for (nm in names(res)) {
    message(sprintf("  %s: %s", nm, res[[nm]]))
  }

  expect_type(res, "list")
  expect_true(all(names(res) %in% params$parameter))

  unlink(tmp_file)
})
