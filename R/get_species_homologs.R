#' Retrieve homologs for all supported species
#'
#' Queries Ensembl via biomaRt and retrieves homologs from human, mouse, macaque
#' and zebrafish for the input species. Saves CSV (and optionally RDS) files.
#'
#' @param species_input Character. One of: "human", "mouse", "macaque", "zebrafish".
#' @param output_dir Character. Output directory.
#' @param save_rda Logical. Whether to save .rda files.
#' @return A named list of data frames with homolog tables.
#' @export
get_species_homologs <- function(species_input,
                                 output_dir = "./results",
                                 save_csv = FALSE,
                                 save_rda = FALSE) {

  species_map <- list(
    human     = list(dataset = "hsapiens_gene_ensembl",   prefix = "hsapiens"),
    mouse     = list(dataset = "mmusculus_gene_ensembl",  prefix = "mmusculus"),
    macaque   = list(dataset = "mmulatta_gene_ensembl",   prefix = "mmulatta"),
    zebrafish = list(dataset = "drerio_gene_ensembl",     prefix = "drerio")
  )

  if (!species_input %in% names(species_map))
    stop("Unsupported species: ", species_input)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  dataset_name <- species_map[[species_input]]$dataset

  ensembl <- biomaRt::useMart("ensembl")
  dataset_en <- biomaRt::useDataset(dataset_name, mart = ensembl)


  all_targets <- species_map
  all_targets[[species_input]] <- NULL

  out <- list()

  for (tgt in names(all_targets)) {

    prefix <- all_targets[[tgt]]$prefix

    attrs <- c(
      "ensembl_gene_id",
      "external_gene_name",
      paste0(prefix, "_homolog_ensembl_gene"),
      paste0(prefix, "_homolog_associated_gene_name"),
      paste0(prefix, "_homolog_orthology_type"),
      "description"
    )

    res <- biomaRt::getBM(attributes = attrs, mart = dataset_en)

    if (save_csv) {
    csv_name <- file.path(output_dir, paste0(dataset_name, "_to_", tgt, "_orthologs.csv"))
    write.table(res, csv_name, sep = "\t", quote = FALSE, row.names = FALSE)
}
    if (save_rda) {
      rda_name <- file.path(output_dir, paste0(dataset_name, "_to_", tgt, "_orthologs.rda"))
      saveRDS(res, rda_name)
    }

    out[[tgt]] <- res
  }


}
