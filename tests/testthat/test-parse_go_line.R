test_that("parse_go_line parses GO annotation line correctly", {
  line <- "GO:0008150\t1\tBiological Process\tEXP:GENE1|IDA:GENE2|IEA:GENE3"

  result <- parse_go_line(line)

  expected <- data.frame(
    GO_ID = rep("GO:0008150", 3),
    Level = rep("1", 3),
    Category = rep("Biological Process", 3),
    GENE_ID = c("GENE1", "GENE2", "GENE3"),
    EVIDENCE = c("EXP", "IDA", "IEA"),
    stringsAsFactors = FALSE
  )

  expect_equal(result, expected)
})
