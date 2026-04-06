#' Doc Upload/List Sub-module
#'
#' Upload .txt/.pdf files to the repo-root docs_dir (set in global.R as
#' DOCS.DIR and threaded from server.R → CopilotBoardServer).

copilot_panel_docs_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 8px;",
    shiny::fileInput(ns("upload"), "Upload documents",
      accept = c(".txt", ".pdf"),
      multiple = TRUE,
      width = "100%"
    ),
    DT::dataTableOutput(ns("docs_table"), height = "auto"),
    shiny::uiOutput(ns("delete_ui"))
  )
}

copilot_panel_docs_server <- function(id, docs_dir) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Reactive trigger for refreshing file list
    refresh_trigger <- shiny::reactiveVal(0)

    ## List files in docs_dir
    doc_files <- shiny::reactive({
      refresh_trigger()  # depend on trigger
      files <- list.files(docs_dir, full.names = TRUE)
      if (length(files) == 0) return(data.frame())
      info <- file.info(files)
      data.frame(
        Name = basename(files),
        Size = paste0(round(info$size / 1024, 1), " KB"),
        stringsAsFactors = FALSE
      )
    })

    ## Handle file upload
    shiny::observeEvent(input$upload, {
      req(input$upload)
      for (i in seq_len(nrow(input$upload))) {
        file.copy(
          input$upload$datapath[i],
          file.path(docs_dir, input$upload$name[i]),
          overwrite = TRUE
        )
      }
      refresh_trigger(refresh_trigger() + 1)
    })

    ## Render docs table
    output$docs_table <- DT::renderDataTable({
      df <- doc_files()
      if (nrow(df) == 0) {
        df <- data.frame(Name = "No documents uploaded", Size = "")
      }
      DT::datatable(df,
        selection = "single",
        rownames = FALSE,
        options = list(
          pageLength = 10,
          dom = "t",
          scrollY = "200px",
          scrollCollapse = TRUE
        )
      )
    })

    ## Delete button — shown when a row is selected
    output$delete_ui <- shiny::renderUI({
      idx <- input$docs_table_rows_selected
      if (is.null(idx) || nrow(doc_files()) == 0) return(NULL)
      shiny::actionButton(ns("delete_doc"), "Delete selected",
        class = "btn-sm btn-outline-danger mt-2"
      )
    })

    ## Handle delete
    shiny::observeEvent(input$delete_doc, {
      idx <- input$docs_table_rows_selected
      req(idx)
      df <- doc_files()
      if (idx <= nrow(df)) {
        file_path <- file.path(docs_dir, df$Name[idx])
        if (file.exists(file_path)) {
          file.remove(file_path)
        }
        refresh_trigger(refresh_trigger() + 1)
      }
    })

    ## Return reactive list of doc paths
    shiny::reactive({
      df <- doc_files()
      if (nrow(df) == 0) return(character(0))
      file.path(docs_dir, df$Name)
    })
  })
}
