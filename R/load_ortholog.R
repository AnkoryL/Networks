#' Load ortholog mapping
#'
#' @param from_species Species of input genes ("human", "mouse", "macacaque").
#' @param to_species Target species for mapping ("human", "mouse", "macacaque").
#'
#' @return Data frame with ortholog mapping, e.g., Ensembl IDs.
#' @export
load_ortholog <- function(from_species, to_species) {

  key <- paste0(from_species, "_to_", to_species)

  ortholog_path <- switch(
    key,

    "human_to_macaque"   = system.file("extdata/orthologs/hsapiens_gene_ensembl_to_macaque_orthologs.rda", package = "Networks"),
    "human_to_mouse"     = system.file("extdata/orthologs/hsapiens_gene_ensembl_to_mouse_orthologs.rda", package = "Networks"),
    "human_to_zebrafish" = system.file("extdata/orthologs/hsapiens_gene_ensembl_to_zebrafish_orthologs.rda", package = "Networks"),


    "macaque_to_human"     = system.file("extdata/orthologs/mmulatta_gene_ensembl_to_human_orthologs.rda", package = "Networks"),
    "macaque_to_mouse"     = system.file("extdata/orthologs/mmulatta_gene_ensembl_to_mouse_orthologs.rda", package = "Networks"),
    "macaque_to_zebrafish" = system.file("extdata/orthologs/mmulatta_gene_ensembl_to_zebrafish_orthologs.rda", package = "Networks"),


    "mouse_to_human"     = system.file("extdata/orthologs/mmusculus_gene_ensembl_to_human_orthologs.rda", package = "Networks"),
    "mouse_to_macaque"   = system.file("extdata/orthologs/mmusculus_gene_ensembl_to_macaque_orthologs.rda", package = "Networks"),
    "mouse_to_zebrafish" = system.file("extdata/orthologs/mmusculus_gene_ensembl_to_zebrafish_orthologs.rda", package = "Networks"),


    "zebrafish_to_human"   = system.file("extdata/orthologs/drerio_gene_ensembl_to_human_orthologs.rda", package = "Networks"),
    "zebrafish_to_mouse"   = system.file("extdata/orthologs/drerio_gene_ensembl_to_mouse_orthologs.rda", package = "Networks"),
    "zebrafish_to_macaque" = system.file("extdata/orthologs/drerio_gene_ensembl_to_macaque_orthologs.rda", package = "Networks"),

    stop(sprintf("No ortholog mapping available for %s to %s", from_species, to_species))
  )

  env <- new.env()
  load(ortholog_path, envir = env)

  ortholog_df <- env[[ls(env)[1]]]

  return(ortholog_df)
}
