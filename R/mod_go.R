#' User interface for Common GO Term analysis module
#'
#' Creates the Shiny UI elements for GO term network analysis,
#' including gene input, analysis parameters, result table,
#' network visualization and download button.
#'
#' @noRd
#'
go_ui <- function(id) {

  ns <- NS(id)

  tabPanel(

    "Common GO Term",

    sidebarLayout(
      sidebarPanel(
        fileInput(ns("genes"), "Genes file"),

        selectInput(
          ns("species"),
          "Species",
          choices = c("human", "mouse", "macaque", "zebrafish")
        ),

        selectInput(
          ns("type"),
          "GO type",
          choices = c("BP", "CC", "MF")
        ),

        numericInput(
          ns("level_min"),
          "GO level min",
          value = 5,
          min = 1,
          max = 30
        ),

        numericInput(
          ns("level_max"),
          "GO level max",
          value = 20,
          min = 1,
          max = 30
        ),

        numericInput(
          ns("threshold"),
          "Threshold",
          value = 5,
          min = 1
        ),

        selectInput(
          ns("ortholog"),
          "Ortholog conversion",
          choices = c("none", "human", "mouse", "macaque", "zebrafish")
        ),

        actionButton(
          ns("run"),
          "Run",
          class = "btn-primary",
          style = "width:100%;"
        ),

        br(),
        br(),

        downloadButton(
          ns("download"),
          "Download results",
          class = "btn-success",
          style = "width:100%;"
        )
      ),

      mainPanel(

        h4("Interaction Table:"),

        shinycssloaders::withSpinner(
          DT::DTOutput(ns("table")),
          type = 6,
          color = "#3c8dbc"
        ),

        hr(),

        h4("Network Graph:"),

        div(
          style = "width:100%; text-align:center;",

          shinycssloaders::withSpinner(
            imageOutput(ns("network"), height = "auto"),
            type = 6,
            color = "#3c8dbc"
          )
        )
      )
    )
  )
}

#' Server logic for Common GO Term analysis module
#'
#' Runs GO term network analysis and displays results.
#'
#' @noRd

go_server <- function(id) {

  moduleServer(id, function(input, output, session) {

    rv <- reactiveValues(
      go_dir = NULL,
      trigger = 0
    )

    observeEvent(input$run, {

      showNotification("GO analysis started")
      log_message("GO analysis started")


      rv$go_dir <- file.path(
        tempdir(),
        paste0("go_", as.integer(Sys.time()))
      )

      dir.create(
        rv$go_dir,
        recursive = TRUE
      )


      tryCatch({

        common_go_term(
          genes_list_path = input$genes$datapath,
          species_prefix = input$species,
          go_type = input$type,
          output_folder_path = rv$go_dir,
          level_from = input$level_min,
          level_to = input$level_max,
          threshold = input$threshold,
          use_ortholog = input$ortholog
        )


        log_message(
          paste0(
            "GO analysis parameters:",
            " species = ", input$species,
            ", GO type = ", input$type,
            ", threshold = ", input$threshold,
            ", levels = ",
            input$level_min,
            "-",
            input$level_max,
            ", ortholog = ",
            input$ortholog
          )
        )


        rv$trigger <- rv$trigger + 1

        showNotification(
          "Analysis completed successfully!"
        )

        log_message(
          "GO analysis completed successfully"
        )


      }, error = function(e) {

        showNotification(
          paste("Error in pipeline:", e$message),
          type = "error",
          duration = NULL
        )

        log_message(
          paste(
            "GO analysis failed:",
            conditionMessage(e)
          )
        )

      })

    })



    output$table <- DT::renderDataTable({
      req(rv$go_dir)

      files <- list.files(
        rv$go_dir,
        full.names = TRUE
      )

      table_file <- files[
        grepl(
          "gene_interaction_output_table_for_cytoscape",
          files
        )
      ]

      req(length(table_file) == 1)

        df <- read.table(
        table_file,
        sep = "\t",
        header = TRUE
      )


      DT::datatable(
        df,
        options = list(
          pageLength = 5,
          lengthMenu = c(5, 10, 25, 50),
          scrollY = "250px",
          scrollCollapse = TRUE
        ),
        rownames = FALSE
      )

    })



    output$network <- renderImage({

      req(rv$go_dir)
      rv$trigger

      img_file <- file.path(
        rv$go_dir,
        "gene_network_plot.png"
      )

      req(file.exists(img_file))


      list(
        src = img_file,
        contentType = "image/png",
        width = "100%",
        height = "100%",
        alt = "network"
      )


    }, deleteFile = FALSE)



    output$download <- downloadHandler(

      filename = function() {
        paste0(
          "go_results_",
          Sys.Date(),
          ".zip"
        )
      },


      content = function(file) {

        zip::zipr(
          zipfile = file,
          files = list.files(
            rv$go_dir,
            full.names = TRUE
          )
        )

      }

    )
    session$onSessionEnded(function() {

      unlink(isolate(rv$go_dir), recursive = TRUE, force = TRUE)

    })
  })



}
