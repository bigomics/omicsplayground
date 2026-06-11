# CopilotReportsServer.R — Report context/tool controls for the Copilot board.

.copilot_report_pgx <- function(pgx) {
  value <- if (is.function(pgx)) pgx() else pgx
  if (inherits(value, "reactivevalues")) {
    # Snapshot here, without isolate(), so Shiny tracks additions/updates to
    # the global pgx reactiveValues. The generic normalizer isolates its own
    # fallback snapshot, which is right for non-reactive callers but would
    # hide this panel from later pgx$ai updates.
    value <- tryCatch(
      shiny::reactiveValuesToList(value, all.names = TRUE),
      error = function(e) NULL
    )
  }
  copilot_normalize_pgx(value, source = "reports_server")
}

copilot_report_label <- function(pgx, slot) {
  slot <- tryCatch(as.character(slot)[[1L]], error = function(e) "")
  if (!nzchar(slot)) return("")
  static <- c(
    combined = "Summary",
    de = "Differential Expression",
    pathways = "Enrichment",
    wgcna = "WGCNA",
    wgcna_mox = "moxWGCNA",
    mofa = "MOFA"
  )
  if (startsWith(slot, "drugs_")) {
    return(ai_report_drug_label(pgx, slot))
  }
  label <- static[[slot]]
  if (is.null(label)) slot else label
}

.copilot_report_default_selection <- function(slots) {
  if (!length(slots)) return(character(0))
  if ("combined" %in% slots) "combined" else slots
}

.copilot_report_state_token <- function(pgx) {
  base <- ai_report_dataset_token(pgx)
  slots <- ai_report_slots(pgx)
  if (!length(slots)) return(base)

  report_bits <- vapply(slots, function(slot) {
    entry <- ai_report_get(pgx, slot)
    report <- if (is.null(entry)) "" else entry$report
    paste(
      slot,
      nchar(report, type = "chars"),
      substr(report, 1L, 80L),
      substr(report, max(1L, nchar(report, type = "chars") - 79L),
             nchar(report, type = "chars")),
      sep = ":"
    )
  }, character(1))
  paste(c(base, report_bits), collapse = "|")
}

#' Copilot reports/context panel UI.
#' @param id Shiny module id.
#' @export
CopilotReportsUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::div(
    class = "copilot-context-panel",
    style = "padding: 12px 2px;",
    shiny::uiOutput(ns("reports_ui")),
    shiny::tags$hr(style = "margin: 10px 0;"),
    shiny::checkboxInput(
      ns("tools_enabled"),
      label = "Connect analysis tools",
      value = TRUE
    )
  )
}

#' Copilot reports/context panel server.
#'
#' Lists precomputed `pgx$ai` report slots and tracks which selected reports
#' have already been injected into a user turn.
#'
#' @param id Shiny module id.
#' @param pgx Reactive or value containing the current PGX.
#' @return list(selected_reports, mark_consumed, tools_enabled).
#' @export
CopilotReportsServer <- function(id, pgx) {
  shiny::moduleServer(id, function(input, output, session) {
    consumed <- shiny::reactiveVal(character(0))
    dataset_token <- shiny::reactiveVal("")
    refresh_tick <- shiny::reactiveVal(0L)

    refresh <- function() {
      refresh_tick(shiny::isolate(refresh_tick()) + 1L)
      invisible(NULL)
    }

    pgx_value <- shiny::reactive({
      refresh_tick()
      .copilot_report_pgx(pgx)
    })

    slots <- shiny::reactive({
      value <- pgx_value()
      if (is.null(value)) return(character(0))
      ai_report_slots(value)
    })

    shiny::observeEvent(pgx_value(), {
      value <- pgx_value()
      token <- .copilot_report_state_token(value)
      if (!identical(token, shiny::isolate(dataset_token()))) {
        dataset_token(token)
        consumed(character(0))
      }
    }, ignoreInit = FALSE)

    .input_id <- function(i) paste0("report_", i)

    selected_reports <- shiny::reactive({
      available <- slots()
      if (!length(available)) return(character(0))
      used <- intersect(consumed(), available)
      keep <- vapply(seq_along(available), function(i) {
        slot <- available[[i]]
        if (slot %in% used) return(FALSE)
        val <- input[[.input_id(i)]]
        if (is.null(val)) {
          return(slot %in% .copilot_report_default_selection(available))
        }
        isTRUE(val)
      }, logical(1))
      available[keep]
    })

    mark_consumed <- function(slots) {
      slots <- tryCatch(as.character(slots), error = function(e) character(0))
      slots <- slots[!is.na(slots) & nzchar(slots)]
      if (!length(slots)) return(invisible(NULL))
      consumed(unique(c(consumed(), slots)))
      invisible(NULL)
    }

    tools_enabled <- shiny::reactive({
      is.null(input$tools_enabled) || isTRUE(input$tools_enabled)
    })

    output$reports_ui <- shiny::renderUI({
      value <- pgx_value()
      available <- slots()
      if (is.null(value) || !length(available)) {
        return(shiny::tags$div(
          class = "text-muted small",
          "No precomputed AI reports found."
        ))
      }

      used <- intersect(consumed(), available)
      shiny::tags$div(
        shiny::tags$div(
          class = "small text-muted",
          "Precomputed reports"
        ),
        shiny::tags$ul(
          style = paste(
            "list-style: none;",
            "padding-left: 0;",
            "margin: 6px 0 0 0;"
          ),
          lapply(seq_along(available), function(i) {
            slot <- available[[i]]
            id <- .input_id(i)
            already_used <- slot %in% used
            current <- input[[id]]
            checked <- if (already_used) {
              TRUE
            } else if (is.null(current)) {
              slot %in% .copilot_report_default_selection(available)
            } else {
              isTRUE(current)
            }
            shiny::tags$li(
              style = if (already_used) {
                "opacity: 0.55;"
              } else {
                NULL
              },
              shiny::tags$label(
                class = "checkbox-inline",
                style = "display: block; font-weight: normal;",
                shiny::tags$input(
                  id = session$ns(id),
                  type = "checkbox",
                  checked = if (checked) "checked" else NULL,
                  disabled = if (already_used) "disabled" else NULL
                ),
                shiny::tags$span(
                  style = "margin-left: 6px;",
                  copilot_report_label(value, slot)
                )
              )
            )
          })
        )
      )
    })

    list(
      selected_reports = selected_reports,
      mark_consumed = mark_consumed,
      tools_enabled = tools_enabled,
      refresh = refresh
    )
  })
}
