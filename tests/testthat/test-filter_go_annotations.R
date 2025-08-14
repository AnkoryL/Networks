test_that("filter_go_annotations filters by level range and valid evidence", {
  tmp_conf <- tempfile(fileext = ".tsv")

  evidence_conf <- data.frame(
    Evidence = c("EXP", "IDA", "IEA"),
    Boolean = c("TRUE", "FALSE", "TRUE"),
    stringsAsFactors = FALSE
  )
  write.table(evidence_conf, tmp_conf, sep = "\t", row.names = FALSE, quote = FALSE)

  annotated <- data.frame(
    GO_ID    = c("GO:00001", "GO:00002", "GO:00003"),
    ENSEMBL  = c("ENSG0001", "ENSG0002", "ENSG0003"),
    EVIDENCE = c("EXP", "IDA", "IEA"),
    SYMBOL   = c("GENE1", "GENE2", "GENE3"),
    Level    = c("3,5", "2,7", "4"), # levels as comma-separated strings
    stringsAsFactors = FALSE
  )

  res <- filter_go_annotations(annotated, tmp_conf, level_from = 1, level_to = 20)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("GO_ID", "ENSEMBL", "EVIDENCE", "SYMBOL"))

  expect_equal(nrow(res), 2)
  expect_true(all(res$EVIDENCE %in% c("EXP", "IEA")))
  expect_true(all(res$GO_ID %in% c("GO:00001", "GO:00003")))
})


test_that("filter_go_annotations returns error when no records match", {
  tmp_conf <- tempfile(fileext = ".tsv")

  evidence_conf <- data.frame(
    Evidence = c("EXP"),
    Boolean = c("FALSE"),
    stringsAsFactors = FALSE
  )
  write.table(evidence_conf, tmp_conf, sep = "\t", row.names = FALSE, quote = FALSE)

  annotated <- data.frame(
    GO_ID    = c("GO:00001", "GO:00002"),
    ENSEMBL  = c("ENSG0001", "ENSG0002"),
    EVIDENCE = c("EXP", "EXP"),
    SYMBOL   = c("GENE1", "GENE2"),
    Level    = c("3", "4"),
    stringsAsFactors = FALSE
  )

  expect_error(
    filter_go_annotations(annotated, tmp_conf, level_from = 1, level_to = 20),
    regexp = "No entries passed the filtering criteria"
  )
})
