test_that("get_go_file_by_prefix returns correct file path and handles errors", {
  path <- get_go_file_by_prefix("human", "BP")
  expect_true(file.exists(path))
  expect_match(basename(path), "HomoSapiens.*BiologicalProcess.*\\.txt(\\.gz)?$", perl = TRUE)

  expect_error(
    get_go_file_by_prefix("human", "UNKNOWN"),
    "Unknown go_type"
  )

  expect_error(
    get_go_file_by_prefix("species", "BP"),
    "Unknown species_prefix"
  )

 })
