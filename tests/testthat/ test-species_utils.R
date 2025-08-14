testthat::test_that("infer_species correctly identifies species from Ensembl IDs", {
  expect_identical(infer_species("ENSG00000139618"), "human")
  expect_identical(infer_species("ENSMUSG00000064372"), "mouse")
  expect_identical(infer_species("ENSMMUG00000012345"), "macaque")
  expect_identical(infer_species("ENSDARG00000009876"), "zebrafish")

  expect_identical(infer_species("ensg00000139618"), "human")
  expect_identical(infer_species("ensmusg00000064372"), "mouse")

  expect_true(is.na(infer_species("XYZ00000000001")))

  expect_true(is.na(infer_species("")))
})
