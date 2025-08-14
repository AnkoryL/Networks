test_that("check_annotation_coverage correctly detects unmapped genes", {
  tmp_file <- tempfile(fileext = ".txt")

  annotated <- data.frame(
    ENSEMBL = c("ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004"),
    SYMBOL  = c("GENE1", NA, "GENE3", "GENE4"),
    GENE_ID = c("1", NA, "3", "4"),
    GO_ID   = c("GO:00001", NA, "GO:00003", NA),
    stringsAsFactors = FALSE
  )

  gene_list <- c("ENSG0001", "ENSG0002", "ENSG0003", "ENSG0004")

  res <- check_annotation_coverage(annotated, gene_list, tmp_file)

  file_contents <- readLines(tmp_file)

  expect_type(res, "list")
  expect_named(res, c("total", "unmapped", "unmapped_ids"))
  expect_equal(res$total, 4)
  expect_equal(sort(res$unmapped_ids), sort(c("ENSG0002", "ENSG0004")))
  expect_equal(res$unmapped, 2)

  expect_true(file.exists(tmp_file))
  expect_equal(sort(file_contents), sort(c("ENSG0002", "ENSG0004")))
})


test_that("check_annotation_coverage returns zero unmapped when all genes are annotated", {
  tmp_file <- tempfile(fileext = ".txt")

  annotated <- data.frame(
    ENSEMBL = c("ENSG0001", "ENSG0002"),
    SYMBOL  = c("GENE1", "GENE2"),
    GENE_ID = c("1", "2"),
    GO_ID   = c("GO:00001", "GO:00002"),
    stringsAsFactors = FALSE
  )

  gene_list <- c("ENSG0001", "ENSG0002")

  res <- check_annotation_coverage(annotated, gene_list, tmp_file)

  file_contents <- readLines(tmp_file)

  expect_equal(res$total, 2)
  expect_equal(res$unmapped, 0)
  expect_equal(length(res$unmapped_ids), 0)

  expect_true(file.exists(tmp_file))
  expect_equal(length(file_contents), 0) # empty file
})
