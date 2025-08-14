testthat::test_that("read_input_file reads CSV, TXT/TSV, and Excel correctly", {
  # CSV
  csv_path <- system.file("extdata/example/test_genes.csv", package = "Networks")
  testthat::expect_true(file.exists(csv_path))
  df_csv <- read_input_file(csv_path)
  testthat::expect_s3_class(df_csv, "data.frame")
  testthat::expect_true(nrow(df_csv) > 0)

  # TXT
  txt_path <- system.file("extdata/example/test_genes.txt", package = "Networks")
  testthat::expect_true(file.exists(txt_path))
  df_txt <- read_input_file(txt_path)
  testthat::expect_s3_class(df_txt, "data.frame")
  testthat::expect_true(nrow(df_txt) > 0)

  # Excel
  xlsx_path <- system.file("extdata/example/test_genes.xlsx", package = "Networks")
  testthat::expect_true(file.exists(xlsx_path))
  df_xlsx <- read_input_file(xlsx_path)
  testthat::expect_s3_class(df_xlsx, "data.frame")
  testthat::expect_true(nrow(df_xlsx) > 0)

  # Unsupported
  tmp <- tempfile(fileext = ".pdf")
  writeLines("fake content", tmp)
  testthat::expect_error(read_input_file(tmp), "Unsupported file format")
})
