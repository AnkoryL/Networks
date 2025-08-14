testthat::test_that("read_go_annotations reads and annotates GO data correctly", {
  go_file_path <- system.file("extdata/annotations/HomoSapiens_GO_BiologicalProcess-EBI-UniProt-GOA-ACAP-ARAP_04.04.2025_00h00.txt.gz", package = "Networks")

  testthat::expect_true(file.exists(go_file_path))

  result <- read_go_annotations(go_file_path, "human")

  testthat::expect_s3_class(result, "data.frame")
  testthat::expect_true(all(c("GO_ID", "Level", "Category", "GENE_ID", "EVIDENCE") %in% colnames(result)))
  testthat::expect_gt(nrow(result), 0)
})

