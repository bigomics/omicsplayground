# CopilotHistoryServer.R — Saved-session history panel for the Copilot board.
#
# Lists sessions from `omicsagentovi::ovi_sessions(session_dir)`. Emits
# restore/delete requests as edge-triggered reactiveVals. The orchestrator
# is responsible for actually invoking restore/delete and for bumping
# `history_invalidation_tick` so the panel refreshes.

# ---- Local null-coalescing operator ----------------------------------------
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

#' Render a relative-time label for an ISO-8601 UTC timestamp.
#'
#' Kept in this file for now; future migration to a shared helper is a
#' separate step.
#' @param iso character(1) ISO timestamp, or NULL/"".
#' @return character(1).
copilot_format_relative_time <- function(iso) {
  if (is.null(iso) || !nzchar(iso)) return("")
  t <- suppressWarnings(
    as.POSIXct(sub("Z$", "", iso), format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  )
  if (is.na(t)) return(iso)
  secs <- as.numeric(difftime(Sys.time(), t, units = "secs"))
  if (secs < 60)     return("just now")
  if (secs < 3600)   return(paste0(floor(secs / 60),    "m ago"))
  if (secs < 86400)  return(paste0(floor(secs / 3600),  "h ago"))
  if (secs < 604800) return(paste0(floor(secs / 86400), "d ago"))
  format(t, "%Y-%m-%d")
}

#' Copilot history panel UI.
#' @param id Shiny module id.
#' @export
CopilotHistoryUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    style = "padding: 8px;",
    DT::dataTableOutput(ns("history_table"), height = "auto"),
    shiny::uiOutput(ns("delete_ui"))
  )
}

#' Copilot history panel server.
#'
#' @param id Shiny module id.
#' @param session_dir character(1) path to the SessionStore directory.
#' @param history_invalidation_tick Optional reactive (or reactiveVal). When
#'   provided, the panel takes a dependency on it so the session list
#'   refreshes after save or delete. Renamed from `refresh_trigger`.
#' @return `list(selected, on_restore, on_delete)`.
#'   - `selected` reactive() returning the currently highlighted session_id.
#'   - `on_restore` / `on_delete` reactiveVals are edge-triggered;
#'      the caller is expected to read them, act, then reset to NULL.
#' @export
CopilotHistoryServer <- function(id, session_dir,
                                 history_invalidation_tick = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    .selected   <- shiny::reactiveVal(NULL)
    .on_restore <- shiny::reactiveVal(NULL)
    .on_delete  <- shiny::reactiveVal(NULL)

    .sessions <- shiny::reactive({
      if (!is.null(history_invalidation_tick)) history_invalidation_tick()
      df <- tryCatch(
        omicsagentovi::ovi_sessions(session_dir = session_dir),
        error = function(e) data.frame()
      )
      if (!is.data.frame(df) || nrow(df) == 0L) {
        return(data.frame(
          session_id   = character(0),
          dataset_name = character(0),
          n_turns      = integer(0),
          updated_at   = numeric(0),
          stringsAsFactors = FALSE
        ))
      }
      df
    })

    .display_df <- shiny::reactive({
      df <- .sessions()
      if (nrow(df) == 0L) {
        return(data.frame(Title = character(0),
                          Turns = integer(0),
                          Updated = character(0),
                          stringsAsFactors = FALSE))
      }
      title <- ifelse(is.na(df$dataset_name) | !nzchar(df$dataset_name),
                      "(no dataset)", df$dataset_name)
      updated <- vapply(df$updated_at, function(ts) {
        copilot_format_relative_time(
          format(as.POSIXct(ts, origin = "1970-01-01", tz = "UTC"),
                 "%Y-%m-%dT%H:%M:%SZ")
        )
      }, character(1))
      data.frame(
        Title   = title,
        Turns   = df$n_turns,
        Updated = updated,
        stringsAsFactors = FALSE
      )
    })

    output$history_table <- DT::renderDataTable({
      df <- .display_df()
      if (nrow(df) == 0L) {
        empty <- data.frame(Title = "No conversations yet",
                            Turns = "",
                            Updated = "")
        return(DT::datatable(empty, selection = "none", rownames = FALSE,
                             options = list(dom = "t", paging = FALSE)))
      }
      DT::datatable(
        df,
        selection = "single",
        rownames  = FALSE,
        options   = list(
          pageLength     = 20,
          dom            = "ft",
          scrollY        = "calc(100vh - 320px)",
          scrollCollapse = TRUE
        )
      )
    })

    output$delete_ui <- shiny::renderUI({
      if (is.null(input$history_table_rows_selected) ||
          nrow(.sessions()) == 0L) return(NULL)
      shiny::actionButton(
        session$ns("delete_session"),
        "Delete selected",
        class = "btn-sm btn-outline-danger mt-2"
      )
    })

    shiny::observeEvent(input$history_table_rows_selected, {
      df  <- .sessions()
      idx <- input$history_table_rows_selected
      if (!is.null(idx) && idx <= nrow(df)) {
        sid <- df$session_id[[idx]]
        .selected(sid)
        .on_restore(sid)
      }
    })

    shiny::observeEvent(input$delete_session, {
      df  <- .sessions()
      idx <- input$history_table_rows_selected
      if (!is.null(idx) && idx <= nrow(df)) {
        .on_delete(df$session_id[[idx]])
      }
    })

    list(
      selected   = shiny::reactive(.selected()),
      on_restore = .on_restore,
      on_delete  = .on_delete
    )
  })
}
