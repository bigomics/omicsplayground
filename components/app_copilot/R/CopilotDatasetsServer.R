# CopilotDatasetsServer.R — Dataset picker panel for the Copilot board.
#
# Pure-filesystem module: scans `pgx_dir` for .pgx files and emits the
# selected full path. The orchestrator owns PGX loading and Agent
# rebinding — this panel does not touch playbase or omicsagentovi.

#' Copilot datasets panel UI.
#'
#' @param id Shiny module id.
#' @export
CopilotDatasetsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 8px;",
    shiny::actionButton(
      ns("refresh"), "Refresh",
      class = "btn-sm btn-outline-secondary mb-2"
    ),
    DT::dataTableOutput(ns("dataset_table"), height = "auto")
  )
}

#' Copilot datasets panel server.
#'
#' @param id Shiny module id.
#' @param pgx_dir character(1) path to a directory containing .pgx files.
#'   Plain scalar, not a reactive.
#' @return `list(selected = reactive(<full path or NULL>))`.
#' @export
CopilotDatasetsServer <- function(id, pgx_dir) {
  shiny::moduleServer(id, function(input, output, session) {

    .selected  <- shiny::reactiveVal(NULL)
    .file_list <- shiny::reactiveVal(character(0))

    .scan <- function() {
      if (is.null(pgx_dir) || !dir.exists(pgx_dir)) return(character(0))
      list.files(pgx_dir, pattern = "\\.pgx$", full.names = FALSE)
    }

    shiny::observe({ .file_list(.scan()) })

    shiny::observeEvent(input$refresh, {
      .selected(NULL)
      .file_list(.scan())
    })

    shiny::observeEvent(input$dataset_table_rows_selected, {
      idx   <- input$dataset_table_rows_selected
      files <- .file_list()
      if (!is.null(idx) && idx <= length(files)) {
        .selected(file.path(pgx_dir, files[[idx]]))
      }
    })

    output$dataset_table <- DT::renderDataTable({
      files <- .file_list()
      df <- if (length(files) == 0L) {
        data.frame(Dataset = "No datasets found")
      } else {
        data.frame(Dataset = tools::file_path_sans_ext(files))
      }
      DT::datatable(
        df,
        selection = "single",
        rownames  = FALSE,
        options   = list(
          pageLength     = 20,
          dom            = "ft",
          scrollY        = "calc(100vh - 280px)",
          scrollCollapse = TRUE
        )
      )
    })

    list(selected = shiny::reactive(.selected()))
  })
}
