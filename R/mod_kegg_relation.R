#' User interface for KEGG relation analysis module
#'
#' @noRd

kegg_relation_ui <- function(id) {

  ns <- NS(id)

  tagList(

    sidebarLayout(

      sidebarPanel(

        fileInput(ns("genes"), "Genes file"),

        selectInput(
          ns("species"),
          "Species",
          choices = c("human", "mouse", "macaque", "zebrafish")
        ),

        selectInput(
          ns("subtype_filter"),
          "Filter",
          choices = c(NULL, "compound")
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

        h4("Interaction Table:"),

        shinycssloaders::withSpinner(
          DT::DTOutput(ns("table")),
          type = 6,
          color = "#3c8dbc"
        ),

        hr(),

        h4("Network Graph:"),

        div(
          style = "width: 100%; text-align: center;",

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



#' Server logic for KEGG relation analysis module
#'
#' @noRd

kegg_relation_server <- function(id) {

  moduleServer(id, function(input, output, session) {


    rv <- reactiveValues(
      kegg_relation_dir = NULL,
      trigger = 0
    )


    observeEvent(input$run, {

      showNotification("KEGG relation analysis started")
      log_message("KEGG relation analysis started")


      rv$kegg_relation_dir <- file.path(
        tempdir(),
        paste0("kegg_relation_", as.integer(Sys.time()))
      )

      dir.create(rv$kegg_relation_dir, recursive = TRUE)


      tryCatch({

        relation_in_kegg_pathway(
          genes_list_path = input$genes$datapath,
          species_prefix = input$species,
          output_folder_path = rv$kegg_relation_dir,
          subtype_filter = input$subtype_filter,
          use_ortholog = input$ortholog
        )


        log_message(
          paste0(
            "KEGG relation analysis parameters:",
            " species = ", input$species,
            ", filter = ", input$subtype_filter,
            ", ortholog = ", input$ortholog
          )
        )


        rv$trigger <- rv$trigger + 1

        showNotification("Analysis completed successfully!")
        log_message("KEGG relation analysis completed successfully")


      }, error = function(e) {

        showNotification(
          paste("Error in pipeline:", e$message),
          type = "error",
          duration = NULL
        )

        log_message(
          paste("KEGG relation analysis failed:", conditionMessage(e))
        )

      })

    })



    output$table <- DT::renderDataTable({

      req(rv$kegg_relation_dir)


      files <- list.files(rv$kegg_relation_dir, full.names = TRUE)

      table_file <- files[
        grepl("gene_interaction_output_table_for_cytoscape", files)
      ]

      req(length(table_file) == 1)


      df <- read.table(table_file, sep = "\t", header = TRUE)


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

      req(rv$kegg_relation_dir)

      rv$trigger


      img_file <- file.path(
        rv$kegg_relation_dir,
        "gene_relation_network_plot.png"
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
        paste0("kegg_relation_results_", Sys.Date(), ".zip")
      },

      content = function(file) {

        zip::zipr(
          zipfile = file,
          files = list.files(rv$kegg_relation_dir, full.names = TRUE)
        )

      }
    )
    session$onSessionEnded(function() {

      unlink(isolate(rv$kegg_relation_dir), recursive = TRUE, force = TRUE)

    })

  })


}
