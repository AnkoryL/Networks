testthat::test_that("read_config_txt reads config from templates correctly", {
  config_path <- system.file("extdata/templates/main_config.txt", package = "Networks")
  testthat::expect_true(file.exists(config_path))

  cfg <- read_config_txt(config_path)

  testthat::expect_type(cfg, "list")
  testthat::expect_true(all(names(cfg) != ""))

  for (field in c("level_from", "level_to", "threshold")) {
    if (field %in% names(cfg)) {
      testthat::expect_type(cfg[[field]], "double")
    }
  }

  if ("species" %in% names(cfg)) {
    testthat::expect_type(cfg$species, "character")
  }
})
