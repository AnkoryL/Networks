test_that("project config is created", {
  temp_dir <- tempfile("test_config_")
  create_project_config(temp_dir, open_main_config = FALSE)
  expect_true(file.exists(file.path(temp_dir, "main_config.txt")))
  expect_true(file.exists(file.path(temp_dir, "config_evidence.txt")))
  expect_true(dir.exists(file.path(temp_dir, "input")))
  expect_true(dir.exists(file.path(temp_dir, "results")))
})
