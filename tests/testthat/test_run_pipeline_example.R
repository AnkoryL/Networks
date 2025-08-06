test_that("Pipeline runs on example config", {
  example_path <- system.file("extdata/example", package = "yourPackage")
  config_path <- file.path(example_path, "example_main_config.txt")
  expect_error(run_pipeline(config_path), NA)
})
