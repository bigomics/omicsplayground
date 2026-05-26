# CopilotBoardUI.R — Copilot Board UI Shell
#
# Three-column bslib layout. Left: datasets/history/docs tabset.
# Centre: CopilotChatUI module. Right: CopilotEvidenceUI module.

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
    style = "height: calc(100vh - 80px);",

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
        bslib::nav_panel("Datasets", CopilotDatasetsUI(ns("datasets"))),
        bslib::nav_panel("History",  CopilotHistoryUI(ns("history"))),
        bslib::nav_panel("Docs",     CopilotDocsUI(ns("docs")))
      )
    ),

    # ---- Centre column: chat ----
    CopilotChatUI(ns("chat")),

    # ---- Right column: evidence ----
    CopilotEvidenceUI(ns("evidence"))
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
