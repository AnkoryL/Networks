
test_that("prepare_genes works with Ensembl_gene_id column", {
  tmp_dir <- tempdir()
  tmp_file <- tempfile(fileext = ".tsv")

  df <- data.frame(Ensembl_gene_id = c("ENSG000001", "ENSG000002"), stringsAsFactors = FALSE)
  write.table(df, tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  res <- prepare_genes(tmp_file, "human", tmp_dir)

  expect_true(is.data.frame(res))
  expect_equal(colnames(res), "ENSEMBL")
  expect_equal(nrow(res), 2)
  expect_equal(res$ENSEMBL, df$Ensembl_gene_id)
})

test_that("prepare_genes works with Gene_symbol column and outputs duplicates file", {
  tmp_dir <- tempdir()
  tmp_file <- tempfile(fileext = ".tsv")

  df <- data.frame(Gene_symbol = c("BRCA1", "TP53", "BRCA1"), stringsAsFactors = FALSE)
  write.table(df, tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  res <- prepare_genes(tmp_file, "human", tmp_dir)

  expect_true(is.data.frame(res))
  expect_equal(colnames(res), "ENSEMBL")
  expect_true(all(grepl("^ENS", res$ENSEMBL)))

  dup_file <- file.path(tmp_dir, "duplicated_symbol_mappings.txt")
  expect_true(file.exists(dup_file))

  dup_content <- read.delim(dup_file, sep = "\t", stringsAsFactors = FALSE)
  expect_true("Gene_symbol" %in% colnames(dup_content))
})

test_that("prepare_genes errors on missing required column", {
  tmp_file <- tempfile(fileext = ".tsv")
  df <- data.frame(Other_col = c("a", "b"), stringsAsFactors = FALSE)
  write.table(df, tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    prepare_genes(tmp_file, "human", tempdir()),
    regexp = "Input must contain either 'Ensembl_gene_id' or 'Gene_symbol' column"
  )
})

test_that("prepare_genes errors on species mismatch", {
  tmp_file <- tempfile(fileext = ".tsv")
  df <- data.frame(Ensembl_gene_id = c("ENSMUSG0000001"), stringsAsFactors = FALSE)
  write.table(df, tmp_file, sep = "\t", quote = FALSE, row.names = FALSE)

  expect_error(
    prepare_genes(tmp_file, "human", tempdir()),
    regexp = "differs from the automatically determined one"
  )
})
