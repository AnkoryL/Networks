
test_that("annotate_with_go_terms correctly joins data and writes output file", {
  tmp_dir <- tempdir()
  species <- "human"

  gene_df <- data.frame(
    ENSEMBL = c("ENSG0001", "ENSG0002"),
    stringsAsFactors = FALSE
  )

  go_df <- data.frame(
    ENSEMBL = c("ENSG0001", "ENSG0003"),
    GO_ID = c("GO:00001", "GO:00002"),
    stringsAsFactors = FALSE
  )

  result <- annotate_with_go_terms(gene_df, go_df, tmp_dir, species)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("ENSEMBL", "GO_ID") %in% colnames(result)))
  expect_equal(nrow(result), nrow(gene_df)) # Join left, все строки gene_df должны остаться

  expected_file <- file.path(tmp_dir, paste0(species, "_term_gene_table.csv"))
  expect_true(file.exists(expected_file))

  file_content <- read.delim(expected_file, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(result, file_content)
})
