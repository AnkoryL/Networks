#' Server logic for Networks application
#'
#' @noRd

app_server <- function(input, output, session) {


  go_server("go")

  kegg_common_server("kegg_common")

  kegg_relation_server("kegg_relation")

  ortholog_server("ortholog")



  session$onSessionEnded(function() {

    log_message("Session ended")

    stopApp()

  })


}
