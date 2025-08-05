#' Infer species from Ensembl gene ID
#'
#' @param ensembl_id A character string with Ensembl ID (e.g. ENSG...).
#'
#' @return A species prefix (human, mouse, etc.) or NA if not matched.
#' @export
infer_species <- function(ensembl_id) {
  dplyr::case_when(
    grepl("^ENSG",    ensembl_id, ignore.case = TRUE) ~ "human",
    grepl("^ENSMUSG", ensembl_id, ignore.case = TRUE) ~ "mouse",
    grepl("^ENSMMUG", ensembl_id, ignore.case = TRUE) ~ "macaque",
    grepl("^ENSDARG", ensembl_id, ignore.case = TRUE) ~ "zebrafish",
    TRUE ~ NA_character_
  )
}
