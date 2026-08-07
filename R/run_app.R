#' Run Networks application
#'
#' @return Opens Networks Shiny application
#' @export

run_app <- function() {

  shiny::runApp(shiny::shinyApp(
    ui = app_ui,
    server = app_server
  ),
  launch.browser = TRUE
  )

}
