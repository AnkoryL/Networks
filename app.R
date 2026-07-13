suppressWarnings(library(shiny))
suppressWarnings((library(DT)))

devtools::load_all()


ui <- fluidPage(

  titlePanel("Networks"),

  tabsetPanel(

#     tabPanel(
#       "Common GO Term",
#
#       sidebarLayout(
# sidebarPanel(fileInput(
#             "go_genes",
#             "Genes file"),
#
#           selectInput(
#             "go_species",
#             "Species",
#             choices = c(
#               "human",
#               "mouse",
#               "macaque",
#               "zebrafish"
#             )
#           ),
#
#           selectInput(
#             "go_type",
#             "GO type",
#             choices = c(
#               "BP",
#               "CC",
#               "MF"
#             )
#           ),
#
#           numericInput(
#             "go_level_min",
#             "GO level min",
#             value = 5,
#             min = 1,
# 			max = 30
#           ),
#
# 		  numericInput(
#             "go_level_max",
#             "GO level max",
#             value = 20,
#             min = 1,
# 			max = 30
#           ),
#
#           numericInput(
#             "go_threshold",
#             "Threshold",
#             value = 5,
#             min = 1
#           ),
#
#           selectInput(
#             "go_ortholog",
#             "Ortholog conversion",
#             choices = c(
#               "none",
#               "human",
# 			        "mouse",
#               "macaque",
#               "zebrafish"
#             )
#           ),
#
#
#           actionButton(
#             "run_go",
#             "Run"
#           )
#         ),
#
#         mainPanel(
#
#           DTOutput("go_table"),
#           br(),
#           imageOutput("go_network"),
#           br(),
#           downloadButton("download_go", "Download results")
#         )
#       )
#     ),

    tabPanel(
      "Common GO Term",
      sidebarLayout(
        sidebarPanel(
          fileInput("go_genes", "Genes file"),
          selectInput("go_species", "Species", choices = c("human", "mouse", "macaque", "zebrafish")),
          selectInput("go_type", "GO type", choices = c("BP", "CC", "MF")),
          numericInput("go_level_min", "GO level min", value = 5, min = 1, max = 30),
          numericInput("go_level_max", "GO level max", value = 20, min = 1, max = 30),
          numericInput("go_threshold", "Threshold", value = 5, min = 1),
          selectInput("go_ortholog", "Ortholog conversion", choices = c("none", "human", "mouse", "macaque", "zebrafish")),

          actionButton("run_go", "Run", class = "btn-primary", style = "width: 100%;"),
          br(), br(),

          downloadButton("download_go", "Download results", class = "btn-success", style = "width: 100%;")
        ),

        mainPanel(
          h4("Interaction Table:"),
          shinycssloaders::withSpinner(
          DTOutput("go_table"),
          type = 6,
          color = "#3c8dbc"
          ),

          hr(),


          h4("Network Graph:"),
          div(
            style = "max-width: 100%; height: auto; max-height: 600px; overflow: hidden; text-align: center;",
            shinycssloaders::withSpinner(
              imageOutput("go_network", height = "auto"),
              type = 6,
              color = "#3c8dbc"
            )
          )
        )
      )
    ),

    tabPanel(
      "Common KEGG pathway",

      sidebarLayout(
        sidebarPanel(
          fileInput("kegg_common_genes", "Genes file"),
          selectInput("kegg_common_species","Species", choices = c("human","mouse","macaque", "zebrafish")),
          numericInput("kegg_common_threshold", "Threshold", value = 5, min = 1),
		  selectInput("kegg_common_ortholog","Ortholog conversion", choices = c("none","human","mouse", "macaque","zebrafish")),
          actionButton("run_kegg_common", "Run")),

		  mainPanel(
		  DTOutput("kegg_common_table"),
		  br(),
		  imageOutput("kegg_common_network"),
		  br(),
		  downloadButton("download_kegg_common", "Download results")
		)
      )
    ),

    tabPanel(
      "Relation in KEGG pathway",

      sidebarLayout(

        sidebarPanel(

          fileInput(
            "kegg_relation_genes",
            "Genes file"
          ),

          selectInput(
            "kegg_relation_species",
            "Species",
            choices = c(
              "human",
              "mouse",
              "macaque",
              "zebrafish"
            )
          ),

          numericInput(
            "kegg_relation_threshold",
            "Threshold",
            value = 5,
            min = 1
          ),

		  selectInput(
            "kegg_relation_ortholog",
            "Ortholog conversion",
            choices = c(
              "none",
              "human",
			  "mouse",
              "macaque",
              "zebrafish"
            )
          ),

          actionButton(
            "run_kegg_relation",
            "Run"
          )
        ),

        mainPanel(

          DTOutput("kegg_relation_table"),

          imageOutput("kegg_relation_network")
        )
      )
    ),

    tabPanel(
      "Ortholog Conversion",

      sidebarLayout(

        sidebarPanel(

          fileInput(
            "ortholog_genes",
            "Genes file"
          ),

          selectInput(
            "from_species",
            "From",
            choices = c(
              "human",
              "mouse",
              "macaque",
              "zebrafish"
            )
          ),

          selectInput(
            "to_species",
            "To",
            choices = c(
              "human",
              "mouse",
              "macaque",
              "zebrafish"
            )
          ),

          actionButton(
            "run_ortholog",
            "Run"
          )
        ),

        mainPanel(

          DTOutput("ortholog_table")
        )
      )
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    out_dir = NULL
  )
  observeEvent(input$run_go, {

    showNotification(
      "GO analysis started"
    )

    rv$out_dir <- file.path(tempdir(), paste0("go_", as.integer(Sys.time())))
    dir.create(rv$out_dir, recursive = TRUE)


      tryCatch({
        common_go_term(
          genes_list_path = input$go_genes$datapath,
          species_prefix = input$go_species,
          go_type = input$go_type,
          output_folder_path = rv$out_dir,
          level_from = input$go_level_min,
          level_to = input$go_level_max,
          threshold = input$go_threshold,
          use_ortholog = input$go_ortholog
        )

        rv$trigger <- rv$trigger + 1
        showNotification("Analysis completed successfully!", type = "default")

      }, error = function(e) {
        showNotification(paste("Error in pipeline:", e$message), type = "error", duration = NULL)
      })

    showNotification("Done")
    })


  output$go_table <- DT::renderDataTable({

    req(rv$out_dir)

    files <- list.files(rv$out_dir, full.names = TRUE)
    table_file <- files[grepl("gene_interaction_output_table_for_cytoscape", files)]

    req(length(table_file) == 1)

    read.table(table_file, sep = "\t", header = TRUE)
  })

  # output$go_network <- renderImage({
  #
  #   req(rv$out_dir)
  #
  #   img_file <- file.path(rv$out_dir, "gene_network_plot.png")
  #   req(file.exists(img_file))
  #
  #   list(
  #     src = img_file,
  #     contentType = "image/png",
  #     alt = "network"
  #   )
  # }, deleteFile = FALSE)

  output$go_network <- renderImage({
    req(rv$out_dir)
    rv$trigger

    img_file <- file.path(rv$out_dir, "gene_network_plot.png")
    req(file.exists(img_file))

    list(
      src = img_file,
      contentType = "image/png",
      width = "100%",
      height = "100%",
      alt = "network"
    )
  }, deleteFile = FALSE)

  output$download_go <- downloadHandler(

    filename = function() {
      paste0("go_results_", Sys.Date(), ".zip")
    },

    content = function(file) {

      zip::zipr(
        zipfile = file,
        files = list.files(rv$out_dir, full.names = TRUE)
      )
    }
  )

  session$onSessionEnded(function() {
    unlink(rv$out_dir, recursive = TRUE, force = TRUE)
  })

  onStop(function() {
    unlink("work", recursive = TRUE, force = TRUE)
  })


	  observeEvent(input$run_kegg_common, {

		showNotification(
		  "KEGG analysis started"
		)

	  })

	  observeEvent(input$run_kegg_relation, {

		showNotification(
		  "KEGG analysis started"
		)

	  })

	  observeEvent(input$run_ortholog, {

		showNotification(
		  "Ortholog conversion started"
		)

	  })
	onStop(function() {
	  unlink("work", recursive = TRUE, force = TRUE)
	})

}

runApp(
  shinyApp(ui, server),
  launch.browser = TRUE
)

