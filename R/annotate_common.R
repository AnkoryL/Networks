#' Annotate genes with GO terms from a parsed ClueGO table
#'
#' @name annotate_with_go_terms
#'
#' @param gene_df Data frame with input genes.
#' @param go_df Data frame with GO annotations
#' @param output_folder_path Directory where the annotated table will be saved.
#' @param species_prefix Short species code used in output naming.
#'
#' @return A data frame with GO annotations merged into the input gene list.
#' @export
annotate_with_go_terms <- function(gene_df, go_df, output_folder_path, species_prefix) {
  annotated <- dplyr::left_join(
    gene_df,
    go_df,
    by = c("Ensembl_gene_id" = "ENSEMBL")
  )

  annotated_genes <- dplyr::filter(annotated, !is.na(GO_ID))
  not_mapped      <- dplyr::filter(annotated, is.na(GO_ID))

  message("Annotated genes: ", length(unique(annotated_genes$Ensembl_gene_id)))
  message("Not mapped: ", nrow(not_mapped))

  out_term_path <- file.path(
    output_folder_path,
    paste0(species_prefix, "_term_gene_table.txt")
  )

  write.table(
    annotated,
    file = out_term_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE)

  out_not_mapped <- file.path(
    output_folder_path,
    paste0(species_prefix, "_not_mapped_table.txt")
  )

  write.table(
    not_mapped[,1:3],
    file = out_not_mapped,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  return(list(
    annotated_genes = annotated_genes,
    not_mapped = nrow(not_mapped)
  ))
}


#' Annotate genes with KEGG common pathways
#'
#' @param gene_df Data frame with input genes. Must contain `Gene_symbol`.
#' @param kegg_df Data frame with KEGG annotations, containing `genes`
#'   and `pathway_id`.
#'
#' @return A list with elements:
#'   \item{annotated_genes}{Genes successfully mapped to KEGG pathways}
#'   \item{not_mapped}{Genes without KEGG annotation}
#' @export
#'
#'
annotate_with_kegg_common_pathways <- function(gene_df, species_prefix, kegg_df,output_folder_path) {
  annotated <- dplyr::left_join(
    gene_df,
    kegg_df,
    by = c("Gene_symbol" = "genes")
  )

  annotated_genes <- dplyr::filter(annotated, !is.na(pathway_id))
  not_mapped <- dplyr::filter(annotated, is.na(pathway_id))

  message("Annotated genes: ", length(unique(annotated_genes$Ensembl_gene_id)))
  message("Not mapped: ", nrow(not_mapped))

  out_not_mapped <- file.path(
    output_folder_path,
    paste0(species_prefix, "_not_mapped_table.txt")
  )

  write.table(
    not_mapped[,1:3],
    file = out_not_mapped,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )


  return(list(
    annotated_genes = annotated_genes,
    not_mapped = nrow(not_mapped)
  ))

  }
