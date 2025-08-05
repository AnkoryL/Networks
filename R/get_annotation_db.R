#' @importFrom AnnotationDbi select keys
#' @import org.Hs.eg.db org.Mm.eg.db org.Mmu.eg.db org.Dr.eg.db
NULL

get_annotation_db <- function(species_prefix) {
  switch(
    species_prefix,
    "human"     = org.Hs.eg.db,
    "mouse"     = org.Mm.eg.db,
    "macaque"   = org.Mmu.eg.db,
    "zebrafish" = org.Dr.eg.db,
    stop("Unsupported species prefix.")
  )
}
