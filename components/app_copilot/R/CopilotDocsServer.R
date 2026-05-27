# CopilotDocsServer.R — Document upload/list panel for the Copilot board.
#
# Pure-filesystem module: stores uploads in `docs_dir`. Does not call
# omicsagentovi — the agent picks up the directory lazily via its bindings.

#' Copilot docs panel UI.
#' @param id Shiny module id.
#' @export
CopilotDocsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 12px 0px;",
    shiny::fileInput(
      ns("upload"),
      #label = "Upload documents",
      label = NULL,
      accept   = c(".txt", ".pdf", ".md"),
      buttonLabel = "Upload...",
      multiple = TRUE,
      width    = "100%"
    ),
    DT::dataTableOutput(ns("docs_table"), height = "auto"),
    shiny::uiOutput(ns("delete_ui"))
  )
}

#' Copilot docs panel server.
#'
#' @param id Shiny module id.
#' @param docs_dir character(1) path to the docs directory. Must exist
#'   before the module is constructed.
#' @return `list(selected, doc_count, on_delete)`.
#' @export
CopilotDocsServer <- function(id, docs_dir) {
  shiny::moduleServer(id, function(input, output, session) {

    .refresh   <- shiny::reactiveVal(0L)
    .selected  <- shiny::reactiveVal(NULL)
    .on_delete <- shiny::reactiveVal(NULL)

    .doc_files <- shiny::reactive({
      .refresh()
      files <- list.files(docs_dir, full.names = TRUE)
      if (!length(files))
        return(data.frame(Name = character(0), Size = character(0),
                          stringsAsFactors = FALSE))
      info <- file.info(files)
      data.frame(
        Name = basename(files),
        Size = paste0(round(info$size / 1024, 1), " KB"),
        stringsAsFactors = FALSE
      )
    })

    shiny::observeEvent(input$upload, {
      shiny::req(input$upload)
      for (i in seq_len(nrow(input$upload))) {
        file.copy(
          input$upload$datapath[[i]],
          file.path(docs_dir, input$upload$name[[i]]),
          overwrite = TRUE
        )
      }
      .selected(NULL)
      .refresh(.refresh() + 1L)
    })

    shiny::observeEvent(input$docs_table_rows_selected, {
      df  <- .doc_files()
      idx <- input$docs_table_rows_selected
      if (!is.null(idx) && idx <= nrow(df)) {
        .selected(df$Name[[idx]])
      } else {
        .selected(NULL)
      }
    })

    shiny::observeEvent(input$delete_doc, {
      df  <- .doc_files()
      idx <- input$docs_table_rows_selected
      shiny::req(!is.null(idx), idx <= nrow(df))
      fpath <- file.path(docs_dir, df$Name[[idx]])
      if (file.exists(fpath)) file.remove(fpath)
      .on_delete(fpath)
      .selected(NULL)
      .refresh(.refresh() + 1L)
    })

    output$docs_table <- DT::renderDataTable({
      df <- .doc_files()
      if (nrow(df) == 0L) {
        return(NULL)
        df <- data.frame(Name = "No documents uploaded", Size = "")
      }
      DT::datatable(
        df,
        selection = "single",
        rownames  = FALSE,
        options   = list(
          pageLength     = 10,
          dom            = "t",
          scrollY        = "200px",
          scrollCollapse = TRUE
        )
      )
    })

    output$delete_ui <- shiny::renderUI({
      if (is.null(input$docs_table_rows_selected) ||
          nrow(.doc_files()) == 0L) return(NULL)
      shiny::actionButton(
        session$ns("delete_doc"),
        "Delete selected",
        class = "btn-sm btn-outline-danger mt-2"
      )
    })

    list(
      selected  = shiny::reactive(.selected()),
      doc_count = shiny::reactive(nrow(.doc_files())),
      on_delete = .on_delete
    )
  })
}
