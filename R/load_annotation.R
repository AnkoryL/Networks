#' Load internal annotation data
#'
#' @param source_type "GO" or "KEGG"
#' @param species_prefix e.g. "human", "mouse", "macaca", "zebrafish"
#' @param connection_type For KEGG: "common_KEGG_pathway" or "relation_in_KEGG_pathway"
#' @param go_type For GO: "BP", "CC", "MF", or "IP"
#'
#' @return Annotation data.frame
#' @keywords internal
load_annotation <- function(source_type, species_prefix, connection_type = NULL, go_type = NULL) {
  key <- switch(
    source_type,
    "GO"   = paste0(go_type, "_", species_prefix),
    "KEGG" = paste0(species_prefix, "_", connection_type)
  )

  annotation_path <- switch(
    key,
    # GO
    "BP_human" = system.file("extdata/go_annotations/go_bp_human.rda", package = "Networks"),
    "CC_human" = system.file("extdata/go_annotations/go_cc_human.rda", package = "Networks"),
    "MF_human" = system.file("extdata/go_annotations/go_mf_human.rda", package = "Networks"),
    "IP_human" = system.file("extdata/go_annotations/go_ip_human.rda", package = "Networks"),
    "BP_mouse" = system.file("extdata/go_annotations/go_bp_mouse.rda", package = "Networks"),
    "CC_mouse" = system.file("extdata/go_annotations/go_cc_mouse.rda", package = "Networks"),
    "MF_mouse" = system.file("extdata/go_annotations/go_mf_mouse.rda", package = "Networks"),
    "IP_mouse" = system.file("extdata/go_annotations/go_ip_mouse.rda", package = "Networks"),
    "BP_macaque" = system.file("extdata/go_annotations/go_bp_macaque.rda", package = "Networks"),
    "CC_macaque" = system.file("extdata/go_annotations/go_cc_macaque.rda", package = "Networks"),
    "MF_macaque" = system.file("extdata/go_annotations/go_mf_macaque.rda", package = "Networks"),
    "IP_macaque" = system.file("extdata/go_annotations/go_ip_macaque.rda", package = "Networks"),
    "BP_zebrafish" = system.file("extdata/go_annotations/go_bp_zebrafish.rda", package = "Networks"),
    "CC_zebrafish" = system.file("extdata/go_annotations/go_cc_zebrafish.rda", package = "Networks"),
    "MF_zebrafish" = system.file("extdata/go_annotations/go_mf_zebrafish.rda", package = "Networks"),
    "IP_zebrafish" = system.file("extdata/go_annotations/go_ip_zebrafish.rda", package = "Networks"),

    # KEGG
    "human_common_kegg_pathway" = system.file("extdata/kegg_annotations/hsa_kegg_common_KEGG_path_annotation.rda", package = "Networks"),
    "human_relation_in_kegg_pathway" = system.file("extdata/kegg_annotations/hsa_kegg_relation_in_pathway_annotation.rda", package = "Networks"),
    "mouse_common_kegg_pathway" = system.file("extdata/kegg_annotations/mouse_kegg_common_KEGG_path_annotation.rda", package = "Networks"),
    "mouse_relation_in_kegg_pathway" = system.file("extdata/kegg_annotations/mouse_kegg_relation_in_pathway_annotation.rda", package = "Networks"),
    "macaca_common_kegg_pathway" = system.file("extdata/kegg_annotations/mcc_kegg_common_KEGG_path_annotation.rda", package = "Networks"),
    "macaca_relation_in_kegg_pathway" = system.file("extdata/kegg_annotations/mcc_kegg_relation_in_pathway_annotation.rda", package = "Networks"),
    "zebrafish_common_kegg_pathway" = system.file("extdata/kegg_annotations/zebrafish_kegg_common_KEGG_path_annotation.rda", package = "Networks"),
    "zebrafish_relation_in_kegg_pathway" = system.file("extdata/kegg_annotations/zebrafish_kegg_relation_in_pathway_annotation.rda", package = "Networks"),


    stop(sprintf("No internal annotation found for %s / %s / %s / %s",
                 source_type, species_prefix, connection_type, go_type))
  )


  env <- new.env()
  load(annotation_path, envir = env)

  return(env[[ls(env)[1]]])

}
