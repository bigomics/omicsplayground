# CopilotHistoryServer.R — Saved-session history panel for the Copilot board.
#
# Lists sessions from `omicsagentovi::ovi_sessions(session_dir)`. Emits
# restore/delete requests as edge-triggered reactiveVals. The orchestrator
# is responsible for actually invoking restore/delete and for bumping
# `history_invalidation_tick` so the panel refreshes.

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
#' @param db_name character(1) filename of the SessionStore database within
#'   `session_dir`. Defaults to `"sessions.sqlite"`; the orchestrator passes
#'   the per-user name so briefings read the same db the chat is saved to.
#' @param history_invalidation_tick Optional reactive (or reactiveVal). When
#'   provided, the panel takes a dependency on it so the session list
#'   refreshes after save or delete. Renamed from `refresh_trigger`.
#' @return `list(selected, on_restore, on_delete)`.
#'   - `selected` reactive() returning the currently highlighted session_id.
#'   - `on_restore` / `on_delete` reactiveVals are edge-triggered;
#'      the caller is expected to read them, act, then reset to NULL.
#' @export
CopilotHistoryServer <- function(id, session_dir,
                                 db_name = "sessions.sqlite",
                                 history_invalidation_tick = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    .selected   <- shiny::reactiveVal(NULL)
    .on_restore <- shiny::reactiveVal(NULL)
    .on_delete  <- shiny::reactiveVal(NULL)

    .sessions <- shiny::reactive({
      if (!is.null(history_invalidation_tick)) history_invalidation_tick()
      df <- tryCatch(
        omicsagentovi::ovi_sessions(session_dir = session_dir, db_name = db_name),
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

    # First-user-message lookup. ovi_sessions doesn't expose transcript
    # contents, so we read transcript_records directly via DBI. WAL mode lets
    # us share the SessionStore file safely (read-only open per refresh).
    .briefings <- shiny::reactive({
      df <- .sessions()
      if (nrow(df) == 0L) return(character(0))
      db_path <- file.path(session_dir, db_name)
      if (!file.exists(db_path)) return(rep("", nrow(df)))
      tryCatch({
        con <- DBI::dbConnect(RSQLite::SQLite(), db_path, flags = RSQLite::SQLITE_RO)
        on.exit(DBI::dbDisconnect(con), add = TRUE)
        ids <- df$session_id
        placeholders <- paste(rep("?", length(ids)), collapse = ",")
        # Prefer content_text_visible (user-typed text without injected
        # context preamble); fall back to content_text for legacy rows
        # written before omicsagentovi 0.5.0 added the visible split.
        rows <- DBI::dbGetQuery(
          con,
          sprintf(
            "SELECT t.session_id,
                    COALESCE(t.content_text_visible, t.content_text) AS content_text
               FROM transcript_records t
              INNER JOIN (
                SELECT session_id, MIN(idx) AS min_idx
                  FROM transcript_records
                 WHERE role = 'user'
                   AND content_text IS NOT NULL
                   AND length(content_text) > 0
                   AND session_id IN (%s)
                 GROUP BY session_id
              ) m ON t.session_id = m.session_id AND t.idx = m.min_idx",
            placeholders
          ),
          params = as.list(ids)
        )
        first_msg <- setNames(rows$content_text, rows$session_id)
        out <- vapply(ids, function(sid) {
          msg <- first_msg[[sid]] %||% ""
          if (!nzchar(msg)) return("")
          msg <- gsub("\\s+", " ", trimws(msg))
          if (nchar(msg) <= 80L) msg else paste0(substr(msg, 1L, 80L), "…")
        }, character(1))
        out
      }, error = function(e) {
        log_info("copilot.history.briefing_query_failed", msg = conditionMessage(e))
        rep("", nrow(df))
      })
    })

    .display_df <- shiny::reactive({
      df <- .sessions()
      if (nrow(df) == 0L) {
        return(data.frame(Dataset = character(0),
                          Briefing = character(0),
                          Updated = character(0),
                          stringsAsFactors = FALSE))
      }
      dataset <- ifelse(is.na(df$dataset_name) | !nzchar(df$dataset_name),
                        "(no dataset)", df$dataset_name)
      briefing <- .briefings()
      if (length(briefing) != nrow(df)) briefing <- rep("", nrow(df))
      updated <- vapply(df$updated_at, function(ts) {
        copilot_format_relative_time(
          format(as.POSIXct(ts, origin = "1970-01-01", tz = "UTC"),
                 "%Y-%m-%dT%H:%M:%SZ")
        )
      }, character(1))
      data.frame(
        Dataset  = dataset,
        Briefing = briefing,
        Updated  = updated,
        stringsAsFactors = FALSE
      )
    })

    output$history_table <- DT::renderDataTable({
      df <- .display_df()
      if (nrow(df) == 0L) {
        empty <- data.frame(Dataset  = "No conversations yet",
                            Briefing = "",
                            Updated  = "")
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
