test_that("build_network works for multiple label and connection types", {
  tmpdir <- tempdir()

  filtered_data <- data.frame(
    GO_ID = c("GO:0001", "GO:0001", "GO:0002", "GO:0002", "GO:0003", "GO:0003"),
    ENSEMBL = c("ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004", "ENSG000005", "ENSG000006"),
    EVIDENCE = c("EXP", "EXP", "IDA", "IDA", "TAS", "TAS"),
    SYMBOL = c("GENE1", "GENE2", "GENE1", "GENE3", "GENE2", "GENE4"),
    stringsAsFactors = FALSE
  )

  gene_df <- data.frame(
    ENSEMBL = c("ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004"),
    stringsAsFactors = FALSE
  )

  duplicated_symbols_stub <- data.frame(
    ENSEMBL = c("ENSG000001", "ENSG000002"),
    SYMBOL = c("GENE1", "GENE1"),
    stringsAsFactors = FALSE
  )

  label_types <- c("SYMBOL", "SYMBOLalias", "both")
  connection_types <- c("common_go_term")

  for (lt in label_types) {
    for (ct in connection_types) {
      res <- build_network(
        filtered_data = filtered_data,
        gene_df = gene_df,
        duplicated_symbols = if (lt == "SYMBOLalias") duplicated_symbols_stub else NULL,
        label_type = lt,
        layout_name = "layout_nicely",
        threshold = 1,
        output_folder_path = tmpdir,
        connection_type = ct
      )

      expect_s3_class(res$graph, "igraph")
      expect_true(file.exists(res$plot_path))
      expect_true(file.exists(file.path(tmpdir, "gene_interaction_output_table_base.csv")))
      expect_gt(nrow(res$final_table), 0)
      expect_true(all(c("gene1", "gene2", "connection_type", "values") %in% colnames(res$final_table)))
    }
  }
})
