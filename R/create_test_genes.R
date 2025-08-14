#' Create genes files for multiple organisms
#'
#' This function generates gene annotation tables for a predefined set of organisms,
#' including Ensembl IDs and gene symbols, and saves them as TXT, CSV, and XLSX files
#' in the specified output directory.
#'
#' @param output_folder_path Character. Path to the directory where the annotation files
#'        will be saved. The directory will be created if it does not exist.
#'
#' @return NULL. Files are written to disk as a side effect.
#'
#' @details
#' The function supports human, mouse, macaque, and zebrafish organisms.
#' For each organism, three files are generated:
#' \itemize{
#'   \item{Ensembl-only gene IDs}
#'   \item{Gene symbol-only IDs}
#'   \item{Merged table with Ensembl and gene symbols}
#' }
#'
#' @export
#'

# create_test_genes("C:/Users/Ankory/YandexDisk/is_not_a_wolf/Networks/inst/extdata/example/")
create_test_genes <- function(output_folder_path) {
  if (!dir.exists(output_folder_path)) {
    dir.create(output_folder_path, recursive = TRUE)
  }

  create_path <- function(filename) {
    file.path(output_folder_path, filename)
  }

  organisms_data <- list(
    human = list(
      ensembl = c(
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000157764",
        "ENSG00000146648",
        "ENSG00000139618",
        "ENSG00000073282",
        "ENSG00000136997"
      ),
      symbols = c(
        "TP53",
        "BRCA1",
        "BRAF",
        "EGFR",
        "BRCA2",
        "TP63",
        "MYC")
    ),
    mouse = list(
      ensembl = c(
        "ENSMMUG00000014601",
        "ENSMUSG00000017146",
        "ENSMUSG00000002413",
        "ENSMUSG00000020122",
        "ENSMUSG00000041147",
        "ENSMUSG00000022510",
        "ENSMUSG00000022346"
      ),
      symbols = c(
        "Trp53",
        "Brca1",
        "Braf",
        "Egfr",
        "Brca2",
        "Trp63",
        "Myc")
    ),
    macaque = list(
      ensembl = c(
        "ENSMMUG00000008639",
        "ENSMMUG00000001329",
        "ENSMMUG00000042793",
        "ENSMMUG00000022394",
        "ENSMMUG00000007197",
        "ENSMMUG00000016466",
        "ENSMMUG00000014601"
      ),
      symbols = c(
        "TP53",
        "BRCA1",
        "BRAF",
        "EGFR",
        "BRCA2",
        "TP63",
        "MYC")
    ),
    zebrafish = list(
      ensembl = c(
        "ENSDARG00000035559",
        "ENSDARG00000098242",
        "ENSDARG00000017661",
        "ENSDARG00000013847",
        "ENSDARG00000079015",
        "ENSDARG00000044356",
        "ENSDARG00000101637"
      ),
      symbols = c(
        "tp53",
        "brcc3",
        "braf",
        "egfra",
        "brca2",
        "tp63",
        "ccnd1")
    )
  )

  save_data <- function(data, suffix) {
    write.table(data, create_path(paste0(suffix, ".txt")),  sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data, create_path(paste0(suffix, ".csv")), sep = "\t", row.names = FALSE, quote = FALSE)
    openxlsx::write.xlsx(data, create_path(paste0(suffix, ".xlsx")))
  }

  create_data <- function(organism, data) {
    merged_data <- data.frame(
      Ensembl_gene_id = data$ensembl,
      Gene_symbol = data$symbols,
      stringsAsFactors = FALSE)

    ensembl_only <- data.frame(
      Ensembl_gene_id = data$ensembl,
      stringsAsFactors = FALSE)

    symbol_only <- data.frame(
      Gene_symbol = data$symbols,
      stringsAsFactors = FALSE)

    save_data(ensembl_only, paste0(organism, "_ensembl"))
    save_data(symbol_only, paste0(organism, "_symbol"))
    save_data(merged_data, paste0(organism, "_ensembl_symbol"))
  }

  for (org in names(organisms_data)) {
    message(paste0("Processing organism: ", org, "\n"))
    create_data(org, organisms_data[[org]])
  }

  message("All files have been saved successfully!\n")
}
