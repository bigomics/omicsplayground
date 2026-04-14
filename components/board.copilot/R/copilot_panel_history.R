#' Session History Sub-module
#'
#' Lists persisted Copilot sessions (via omicsagentovi::ovi_sessions +
#' ovi_session_meta), supports "+ New chat", row-click restore, and
#' delete-selected. Returns three reactives to the parent server:
#'
#'   - restore_request  reactiveVal(session_id) — set on row click
#'   - new_chat_request reactiveVal(int)        — incremented on New chat
#'   - delete_request   reactiveVal(session_id) — set on delete confirm

copilot_panel_history_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 8px;",
    shiny::actionButton(
      ns("new_chat"), "+ New chat",
      class = "btn-sm btn-primary mb-2",
      width = "100%"
    ),
    DT::dataTableOutput(ns("history_table"), height = "auto"),
    shiny::uiOutput(ns("delete_ui"))
  )
}

copilot_format_relative_time <- function(iso) {
  if (is.null(iso) || !nzchar(iso)) return("")
  t <- suppressWarnings(as.POSIXct(iso, tz = "UTC"))
  if (is.na(t)) return(iso)
  secs <- as.numeric(difftime(Sys.time(), t, units = "secs"))
  if (secs < 60)      return("just now")
  if (secs < 3600)    return(paste0(floor(secs / 60), "m ago"))
  if (secs < 86400)   return(paste0(floor(secs / 3600), "h ago"))
  if (secs < 604800)  return(paste0(floor(secs / 86400), "d ago"))
  format(t, "%Y-%m-%d")
}

copilot_panel_history_server <- function(id, chat_store, refresh_trigger = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    session_dir <- chat_store@state$session_dir

    restore_request  <- shiny::reactiveVal(NULL)
    new_chat_request <- shiny::reactiveVal(0)
    delete_request   <- shiny::reactiveVal(NULL)

    sessions_df <- shiny::reactive({
      if (!is.null(refresh_trigger)) refresh_trigger()
      ids <- tryCatch(
        omicsagentovi::ovi_sessions(session_dir = session_dir),
        error = function(e) character(0)
      )
      if (length(ids) == 0) {
        return(data.frame(
          id = character(0), Title = character(0),
          Datasets = character(0), Turns = integer(0),
          Updated = character(0), updated_at = character(0),
          stringsAsFactors = FALSE
        ))
      }
      metas <- lapply(ids, function(i) {
        tryCatch(
          omicsagentovi::ovi_session_meta(i, session_dir = session_dir),
          error = function(e) NULL
        )
      })
      keep <- !vapply(metas, is.null, logical(1))
      ids <- ids[keep]; metas <- metas[keep]
      if (length(ids) == 0) {
        return(data.frame(
          id = character(0), Title = character(0),
          Datasets = character(0), Turns = integer(0),
          Updated = character(0), updated_at = character(0),
          stringsAsFactors = FALSE
        ))
      }
      df <- data.frame(
        id       = ids,
        Title    = vapply(metas, function(m) {
          t <- m$title %||% ""; if (!nzchar(t)) "(untitled)" else t
        }, character(1)),
        Datasets = vapply(metas, function(m) {
          ds <- m$datasets %||% character(0)
          if (length(ds) == 0) "" else paste(ds, collapse = ", ")
        }, character(1)),
        Turns    = vapply(metas, function(m) as.integer(m$turn_count %||% 0L), integer(1)),
        Updated  = vapply(metas, function(m) copilot_format_relative_time(m$updated_at %||% ""), character(1)),
        updated_at = vapply(metas, function(m) m$updated_at %||% "", character(1)),
        stringsAsFactors = FALSE
      )
      df[order(df$updated_at, decreasing = TRUE), , drop = FALSE]
    })

    output$history_table <- DT::renderDataTable({
      df <- sessions_df()
      if (nrow(df) == 0) {
        empty <- data.frame(Title = "No conversations yet",
                            Datasets = "", Turns = "", Updated = "")
        return(DT::datatable(empty, selection = "none", rownames = FALSE,
                             options = list(dom = "t", paging = FALSE)))
      }
      DT::datatable(
        df[, c("Title", "Datasets", "Updated"), drop = FALSE],
        selection = "single",
        rownames = FALSE,
        options = list(
          pageLength = 20,
          dom = "ft",
          scrollY = "calc(100vh - 320px)",
          scrollCollapse = TRUE
        )
      )
    })

    ## Row click → restore request
    shiny::observeEvent(input$history_table_rows_selected, {
      idx <- input$history_table_rows_selected
      df <- sessions_df()
      if (!is.null(idx) && idx <= nrow(df)) {
        restore_request(df$id[idx])
      }
    })

    ## Delete-selected UI
    output$delete_ui <- shiny::renderUI({
      idx <- input$history_table_rows_selected
      if (is.null(idx) || nrow(sessions_df()) == 0) return(NULL)
      shiny::actionButton(ns("delete_session"), "Delete selected",
        class = "btn-sm btn-outline-danger mt-2"
      )
    })

    shiny::observeEvent(input$delete_session, {
      idx <- input$history_table_rows_selected
      df <- sessions_df()
      if (!is.null(idx) && idx <= nrow(df)) {
        delete_request(df$id[idx])
      }
    })

    shiny::observeEvent(input$new_chat, {
      new_chat_request(new_chat_request() + 1L)
    })

    list(
      restore_request  = restore_request,
      new_chat_request = new_chat_request,
      delete_request   = delete_request
    )
  })
}
