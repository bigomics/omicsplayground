#' Dataset List Sub-module
#'
#' Scans PGX.DIR for .pgx files, renders selectable list

copilot_panel_datasets_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 8px;",
    shiny::actionButton(ns("refresh"), "Refresh", class = "btn-sm btn-outline-secondary mb-2"),
    DT::dataTableOutput(ns("dataset_table"), height = "auto")
  )
}

copilot_panel_datasets_server <- function(id, pgx_dir) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Currently selected dataset path
    selected_dataset <- shiny::reactiveVal(NULL)

    ## Scan pgx_dir for .pgx files
    dataset_list <- shiny::reactiveVal(character(0))

    scan_datasets <- function() {
      if (is.null(pgx_dir) || !dir.exists(pgx_dir)) return(character(0))
      files <- list.files(pgx_dir, pattern = "\\.pgx$", full.names = FALSE)
      files
    }

    ## Initial scan + refresh button
    shiny::observe({
      dataset_list(scan_datasets())
    })

    shiny::observeEvent(input$refresh, {
      dataset_list(scan_datasets())
    })

    ## Render dataset table
    output$dataset_table <- DT::renderDataTable({
      files <- dataset_list()
      if (length(files) == 0) {
        df <- data.frame(Dataset = "No datasets found")
      } else {
        df <- data.frame(Dataset = tools::file_path_sans_ext(files))
      }
      DT::datatable(df,
        selection = "single",
        rownames = FALSE,
        options = list(
          pageLength = 20,
          dom = "ft",
          scrollY = "calc(100vh - 280px)",
          scrollCollapse = TRUE
        )
      )
    })

    ## On row click: update selected dataset
    shiny::observeEvent(input$dataset_table_rows_selected, {
      idx <- input$dataset_table_rows_selected
      files <- dataset_list()
      if (!is.null(idx) && idx <= length(files)) {
        selected_dataset(file.path(pgx_dir, files[idx]))
      }
    })

    ## Return selected dataset path
    selected_dataset
  })
}
