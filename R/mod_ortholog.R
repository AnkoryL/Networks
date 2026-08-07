#' User interface for ortholog conversion module
#'
#' @noRd

ortholog_ui <- function(id) {

  ns <- NS(id)

  tagList(

    sidebarLayout(

      sidebarPanel(

        fileInput(ns("genes"), "Genes file"),

        selectInput(
          ns("from_species"),
          "From",
          choices = c("human", "mouse", "macaque", "zebrafish")
        ),

        selectInput(
          ns("to_species"),
          "To",
          choices = c("human", "mouse", "macaque", "zebrafish")
        ),

        actionButton(
          ns("run"),
          "Run",
          class = "btn-primary",
          style = "width: 100%;"
        ),

        br(), br(),

        downloadButton(
          ns("download"),
          "Download results",
          class = "btn-success",
          style = "width: 100%;"
        )

      ),


      mainPanel(

        h4("Ortholog Table:"),

        shinycssloaders::withSpinner(
          DT::DTOutput(ns("table")),
          type = 6,
          color = "#3c8dbc"
        )

      )

    )

  )

}




#' Server logic for ortholog conversion module
#'
#' @noRd

ortholog_server <- function(id) {

  moduleServer(id, function(input, output, session) {

    rv <- reactiveValues(
      ortholog_dir = NULL
    )


    observeEvent(input$run, {

      showNotification("Ortholog conversion started")
      log_message("Ortholog conversion started")

      rv$ortholog_dir <- file.path(
        tempdir(),
        paste0("ortholog_", as.integer(Sys.time()))
      )

      dir.create(rv$ortholog_dir, recursive = TRUE)



      tryCatch({
        convert_genes_to_orthologs(
          genes_list_path = input$genes$datapath,
          output_folder_path = rv$ortholog_dir,
          from_species = input$from_species,
          to_species = input$to_species
        )


        log_message(
          paste0(
            "Ortholog conversion parameters:",
            " species from = ", input$from_species,
            ", species to = ", input$to_species
          )
        )

        showNotification("Conversion completed successfully!")
        log_message("Ortholog conversion completed successfully")

      }, error = function(e) {


        showNotification(paste("Error in pipeline:", e$message),
          type = "error",
          duration = NULL
        )

        log_message(paste("Ortholog conversion failed:", conditionMessage(e))
        )


      })


    })




    output$table <- DT::renderDataTable({


      req(rv$ortholog_dir)

      files <- list.files(
        rv$ortholog_dir,
        full.names = TRUE
      )


      table_file <- files[
        grepl(
          "orthologs_list_filtered_by_input_genes.txt",
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


    output$download <- downloadHandler(


      filename = function() {
        paste0("ortholog_results_", Sys.Date(), ".zip")
      },

      content = function(file) {


        zip::zipr(
          zipfile = file,
          files = list.files(
            rv$ortholog_dir,
            full.names = TRUE
          )
        )


      }

    )

    session$onSessionEnded(function() {

      unlink(isolate(rv$ortholog_dir), recursive = TRUE, force = TRUE)


    })

  })

}
