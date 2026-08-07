#' User interface for Networks application
#'
#' @noRd

app_ui <- function(request) {

  fluidPage(

    titlePanel("Networks"),

    tabsetPanel(

      tabPanel(
        "Common GO Term",
        go_ui("go")
      ),

      tabPanel(
        "Common KEGG pathway",
        kegg_common_ui("kegg_common")
      ),

      tabPanel(
        "Relation in KEGG pathway",
        kegg_relation_ui("kegg_relation")
      ),

      tabPanel(
        "Ortholog Conversion",
        ortholog_ui("ortholog")
      )

    )

  )

}
