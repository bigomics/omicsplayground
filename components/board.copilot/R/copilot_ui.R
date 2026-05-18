# copilot_ui.R — Copilot Board UI Shell
#
# Three-column bslib layout. Left: datasets/history/docs tabset.
# Centre: CopilotChatUI module (Phase 4). Right: evidence placeholder.
#
# Phase 6: this file is renamed to CopilotBoardUI.R.

#' Copilot Board UI Shell
#'
#' Three-column bslib layout: datasets/history/docs | chat | evidence.
#' Left column is internal (no sidebar settings exposed).
#'
#' @param id Module namespace id.
#' @return A `bslib::layout_columns` tag.
#' @export
CopilotBoardUI <- function(id) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = c(3, 5, 4),

    # ---- Left column: datasets / history / docs ----
    bslib::card(
      full_screen = FALSE,
      bslib::card_header(
        shiny::actionButton(
          ns("new_chat"),
          label = "New chat",
          icon  = shiny::icon("plus"),
          class = "btn-sm btn-outline-secondary w-100"
        )
      ),
      bslib::navset_underline(
        bslib::nav_panel(
          "Datasets",
          shiny::uiOutput(ns("datasets"))   # TODO(phase 6): CopilotDatasetsUI(ns("datasets"))
        ),
        bslib::nav_panel(
          "History",
          shiny::uiOutput(ns("history"))    # TODO(phase 6): CopilotHistoryUI(ns("history"))
        ),
        bslib::nav_panel(
          "Docs",
          shiny::uiOutput(ns("docs"))       # TODO(phase 6): CopilotDocsUI(ns("docs"))
        )
      )
    ),

    # ---- Centre column: chat (Phase 4 — owned by CopilotChatUI module) ----
    CopilotChatUI(ns("chat")),

    # ---- Right column: evidence ----
    shiny::uiOutput(ns("evidence"))  # TODO(phase 5): CopilotEvidenceUI(ns("evidence"))
  )
}

#' Copilot Board Sidebar Inputs (empty — left column is internal)
#'
#' @param id Module namespace id.
#' @return An empty `shiny::div()`.
#' @export
CopilotBoardInputs <- function(id) {
  shiny::div()
}
