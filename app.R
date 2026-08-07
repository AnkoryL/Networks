# .libPaths(file.path(getwd(), "library"))
app_dir <- normalizePath(getwd())

portable_lib <- file.path(
  dirname(app_dir),
  "library"
)

.libPaths(c(portable_lib, .libPaths()))
# app_dir <- normalizePath(getwd())
# .libPaths(
#   c(file.path(app_dir, "library"),
#     file.path(app_dir, "portable-r-4.6.1-win-x64", "library")))

suppressPackageStartupMessages({
    library(shiny)
    library(DT)
    library(Networks)
})

log_dir <- file.path(getwd(), "logs")

if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
}

log_file <- file.path(
    log_dir,
    paste0("Networks_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".log")
)

log_message <- function(message) {
    cat(
        sprintf(
            "[%s] %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            message
        ),
        file = log_file,
        append = TRUE
    )
}

log_message("Application started")

ui <- fluidPage(

  titlePanel("Networks"),

  tabsetPanel(

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
            style = "width: 100%; text-align: center;",
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

          actionButton("run_kegg_common", "Run", class = "btn-primary", style = "width: 100%;"),
          br(), br(),

          downloadButton("download_kegg_common", "Download results", class = "btn-success", style = "width: 100%;")
        ),

        mainPanel(
          h4("Interaction Table:"),
          shinycssloaders::withSpinner(
            DTOutput("kegg_common_table"),
            type = 6,
            color = "#3c8dbc"
          ),

          hr(),

          h4("Network Graph:"),
          div(
            style = "width: 100%; text-align: center;",
            shinycssloaders::withSpinner(
              imageOutput("kegg_common_network", height = "auto"),
              type = 6,
              color = "#3c8dbc"
            )
          )
        )
      )
    ),


    tabPanel(
      "Relation in KEGG pathway",

      sidebarLayout(
        sidebarPanel(
          fileInput("kegg_relation_genes","Genes file"),
          selectInput("kegg_relation_species","Species",choices = c("human","mouse","macaque","zebrafish")),
          selectInput("kegg_relation_subtype_filter", "Filter", choices = c(NULL,"compound")),
          selectInput("kegg_relation_ortholog","Ortholog conversion",choices = c("none","human","mouse","macaque","zebrafish")),

          actionButton("run_kegg_relation", "Run", class = "btn-primary", style = "width: 100%;"),
          br(), br(),

          downloadButton("download_kegg_relation", "Download results", class = "btn-success", style = "width: 100%;")
        ),

        mainPanel(
          h4("Interaction Table:"),
          shinycssloaders::withSpinner(
            DTOutput("kegg_relation_table"),
            type = 6,
            color = "#3c8dbc"
          ),

          hr(),

          h4("Network Graph:"),
          div(
            style = "width: 100%; text-align: center;",
            shinycssloaders::withSpinner(
              imageOutput("kegg_relation_network", height = "auto"),
              type = 6,
              color = "#3c8dbc"
            )
          )
        )
      )
    ),


    tabPanel(
      "Ortholog Conversion",

      sidebarLayout(
        sidebarPanel(
          fileInput("ortholog_genes","Genes file"),
          selectInput("from_species","From",choices = c("human","mouse","macaque","zebrafish")),
          selectInput("to_species","To",choices = c("human","mouse","macaque","zebrafish")),

          actionButton("run_ortholog", "Run", class = "btn-primary", style = "width: 100%;"),
          br(), br(),

          downloadButton("download_ortholog", "Download results", class = "btn-success", style = "width: 100%;")
        ),

        mainPanel(
          h4("Ortholog Table:"),
          shinycssloaders::withSpinner(
            DTOutput("ortholog_table"),
            type = 6,
            color = "#3c8dbc"
          )
        )
      )
    )
  )
)



server <- function(input, output, session) {

  rv <- reactiveValues(
    go_dir = NULL,
    kegg_common_dir = NULL,
    kegg_relation_dir = NULL,
    ortholog_dir = NULL,

    go_trigger = 0,
    kegg_common_trigger = 0,
    kegg_relation_trigger = 0,
    ortholog_trigger = 0
  )

  observeEvent(input$run_go, {

    showNotification( "GO analysis started")
	log_message("GO analysis started")

    rv$go_dir <- file.path(tempdir(), paste0("go_", as.integer(Sys.time())))
    dir.create(rv$go_dir, recursive = TRUE)


      tryCatch({
        common_go_term(
          genes_list_path = input$go_genes$datapath,
          species_prefix = input$go_species,
          go_type = input$go_type,
          output_folder_path = rv$go_dir,
          level_from = input$go_level_min,
          level_to = input$go_level_max,
          threshold = input$go_threshold,
          use_ortholog = input$go_ortholog
        )
	log_message(
    paste0(
        "GO analysis parameters:",
        " species = ", input$go_species,
        ", GO type = ", input$go_type,
        ", threshold = ", input$threshold,
		", levels = ", input$go_level_min, "-", input$go_level_max,
		", ortholog = ", input$go_ortholog
    )
)
        rv$go_trigger <- rv$go_trigger + 1
        showNotification("Analysis completed successfully!", type = "default")
		log_message("GO analysis completed successfully")
      }, error = function(e) {
        showNotification(paste("Error in pipeline:", e$message), type = "error", duration = NULL)
		log_message(paste("GO analysis failed:", conditionMessage(e)))
      })

    showNotification("Done")
    })


  output$go_table <- DT::renderDataTable({

    req(rv$go_dir)

    files <- list.files(rv$go_dir, full.names = TRUE)
    table_file <- files[grepl("gene_interaction_output_table_for_cytoscape", files)]

    req(length(table_file) == 1)


    df <- read.table(table_file, sep = "\t", header = TRUE)

    DT::datatable(df,options = list(pageLength = 5,lengthMenu = c(5, 10, 25, 50), scrollY = "250px", scrollCollapse = TRUE),rownames = FALSE)

  })

  output$go_network <- renderImage({
    req(rv$go_dir)
    rv$go_trigger

    img_file <- file.path(rv$go_dir, "gene_network_plot.png")
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
        files = list.files(rv$go_dir, full.names = TRUE)
      )
    }
  )



	  observeEvent(input$run_kegg_common, {

		showNotification("KEGG common analysis started")
		log_message("KEGG common analysis started")

	    rv$kegg_common_dir <- file.path(tempdir(), paste0("kegg_common_", as.integer(Sys.time())))
	    dir.create(rv$kegg_common_dir, recursive = TRUE)


	    tryCatch({
	      common_kegg_pathway(
	        genes_list_path = input$kegg_common_genes$datapath,
	        species_prefix = input$kegg_common_species,
	        output_folder_path = rv$kegg_common_dir,
	        threshold = input$kegg_common_threshold,
	        use_ortholog = input$kegg_common_ortholog
	      )

		log_message(
		paste0(
			"common kegg analysis parameters:",
			" species = ", input$kegg_common_species,
			", threshold = ", input$kegg_common_threshold,
			", ortholog = ", input$kegg_common_ortholog
			)
		)

	      rv$kegg_common_trigger <- rv$kegg_common_trigger + 1
	      showNotification("Analysis completed successfully!", type = "default")
		  log_message("kegg analysis completed successfully")

	    }, error = function(e) {
	      showNotification(paste("Error in pipeline:", e$message), type = "error", duration = NULL)
		  log_message(paste("kegg common analysis failed:", conditionMessage(e)))
	    })

	    showNotification("Done")
	  })


	  output$kegg_common_table <- DT::renderDataTable({

	    req(rv$kegg_common_dir)

	    files <- list.files(rv$kegg_common_dir, full.names = TRUE)
	    table_file <- files[grepl("gene_interaction_output_table_for_cytoscape", files)]

	    req(length(table_file) == 1)


	    df <- read.table(table_file, sep = "\t", header = TRUE)

	    DT::datatable(df,options = list(pageLength = 5,lengthMenu = c(5, 10, 25, 50), scrollY = "250px", scrollCollapse = TRUE),rownames = FALSE)

	  })

	  output$kegg_common_network <- renderImage({
	    req(rv$kegg_common_dir)
	    rv$kegg_common_trigger

	    img_file <- file.path(rv$kegg_common_dir, "gene_network_plot.png")
	    req(file.exists(img_file))

	    list(
	      src = img_file,
	      contentType = "image/png",
	      width = "100%",
	      height = "100%",
	      alt = "network"
	    )
	  }, deleteFile = FALSE)

	  output$download_kegg_common <- downloadHandler(

	    filename = function() {
	      paste0("kegg_common_results_", Sys.Date(), ".zip")
	    },

	    content = function(file) {

	      zip::zipr(
	        zipfile = file,
	        files = list.files(rv$kegg_common_dir, full.names = TRUE)
	      )
	    }
	  )



	  observeEvent(input$run_kegg_relation, {

		showNotification("KEGG relation analysis started")
		log_message("KEGG relation analysis started")

	    rv$kegg_relation_dir <- file.path(tempdir(), paste0("kegg_relation_", as.integer(Sys.time())))
	    dir.create(rv$kegg_relation_dir, recursive = TRUE)


	    tryCatch({
	      relation_in_kegg_pathway(
	        genes_list_path = input$kegg_relation_genes$datapath,
	        species_prefix = input$kegg_relation_species,
	        output_folder_path = rv$kegg_relation_dir,
	        subtype_filter = input$kegg_relation_subtype_filter,
	        use_ortholog = input$kegg_relation_ortholog
	      )

		log_message(
		paste0(
			"kegg kegg analysis parameters:",
			" species = ", input$kegg_relation_species,
			", ortholog = ", input$kegg_relation_ortholog
			)
		)

	      rv$kegg_relation_trigger <- rv$kegg_relation_trigger + 1
	      showNotification("Analysis completed successfully!", type = "default")
		  log_message("kegg relation analysis completed successfully")

	    }, error = function(e) {
	      showNotification(paste("Error in pipeline:", e$message), type = "error", duration = NULL)
		  log_message(paste("kegg relation analysis failed:", conditionMessage(e)))

	    })

	    showNotification("Done")
	  })


	  output$kegg_relation_table <- DT::renderDataTable({

	    req(rv$kegg_relation_dir)

	    files <- list.files(rv$kegg_relation_dir, full.names = TRUE)
	    table_file <- files[grepl("gene_interaction_output_table_for_cytoscape", files)]

	    req(length(table_file) == 1)


	    df <- read.table(table_file, sep = "\t", header = TRUE)

	    DT::datatable(df,options = list(pageLength = 5,lengthMenu = c(5, 10, 25, 50), scrollY = "250px", scrollCollapse = TRUE),rownames = FALSE)

	  })

	  output$kegg_relation_network <- renderImage({
	    req(rv$kegg_relation_dir)
	    rv$kegg_relation_trigger

	    img_file <- file.path(rv$kegg_relation_dir, "gene_relation_network_plot.png")
	    req(file.exists(img_file))

	    list(
	      src = img_file,
	      contentType = "image/png",
	      width = "100%",
	      height = "100%",
	      alt = "network"
	    )
	  }, deleteFile = FALSE)

	  output$download_kegg_relation <- downloadHandler(

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



	  observeEvent(input$run_ortholog, {

		showNotification("Ortholog conversion started")
		log_message("Ortholog conversion started")
	    rv$ortholog_dir <- file.path(tempdir(), paste0("ortholog_", as.integer(Sys.time())))
	    dir.create(rv$ortholog_dir, recursive = TRUE)


	    tryCatch({
	      convert_genes_to_orthologs(
	        genes_list_path = input$ortholog_genes$datapath,
	        output_folder_path = rv$ortholog_dir,
	        from_species = input$from_species,
	        to_species = input$to_species
	      )

		log_message(
		paste0(
			"ortholog conversion parameters:",
			" species from = ", input$from_species,
			", species to = ", input$to_species
			)
		)

	      rv$ortholog_trigger <- rv$ortholog_trigger + 1
	      showNotification("Conversion completed successfully!", type = "default")
		  log_message("Conversion completed successfully")

	    }, error = function(e) {
	      showNotification(paste("Error in pipeline:", e$message), type = "error", duration = NULL)
		  log_message(paste("Ortholog conversion failed:", conditionMessage(e)))
	    })

	    showNotification("Done")
	  })


	  output$ortholog_table <- DT::renderDataTable({

	    req(rv$ortholog_dir)

	    files <- list.files(rv$ortholog_dir, full.names = TRUE)
	    table_file <- files[grepl("orthologs_list_filtered_by_input_genes.txt", files)]

	    req(length(table_file) == 1)


	    df <- read.table(table_file, sep = "\t", header = TRUE)

	    DT::datatable(df,options = list(pageLength = 5,lengthMenu = c(5, 10, 25, 50), scrollY = "250px", scrollCollapse = TRUE),rownames = FALSE)

	  })


	  output$download_ortholog <- downloadHandler(

	    filename = function() {
	      paste0("ortholog_results_", Sys.Date(), ".zip")
	    },

	    content = function(file) {

	      zip::zipr(
	        zipfile = file,
	        files = list.files(rv$ortholog_dir, full.names = TRUE)
	      )
	    }
	  )


session$onSessionEnded(function() {

    unlink(isolate(rv$go_dir), recursive = TRUE, force = TRUE)
    unlink(isolate(rv$kegg_common_dir), recursive = TRUE, force = TRUE)
    unlink(isolate(rv$kegg_relation_dir), recursive = TRUE, force = TRUE)
    unlink(isolate(rv$ortholog_dir), recursive = TRUE, force = TRUE)

	log_message("Application session ended")

    stopApp()
})


}



runApp(
  shinyApp(ui, server),
  launch.browser = TRUE
)

