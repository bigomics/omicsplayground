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

  ui <- bslib::layout_columns(
    col_widths = c(3, 5, 4),
    style = "height: calc(100vh - 80px);",

    # ---- Left column: datasets / history / docs ----
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(2,1),
      bslib::card(
        height = "calc(100% - 40px)",
        shiny::actionButton(
          ns("new_chat"),
          label = "New chat",
          icon  = shiny::icon("plus"),
          class = "btn-sm btn-outline-secondary w-100",
          style = "margin-bottom: 8px;"
        ),
        bslib::navset_underline(
          #bslib::nav_panel("Datasets", CopilotDatasetsUI(ns("datasets"))),
          bslib::nav_panel("Dataset",
            shiny::uiOutput(ns("dataset_info"))
          ),
          bslib::nav_panel("Settings",
            CopilotChatSettings(ns("chat"))
          ),
          bslib::nav_panel("History",
            CopilotHistoryUI(ns("history"))
          )
        )
      ),
      bslib::card(
        bslib::navset_underline(
          bslib::nav_panel("Sources", CopilotDocsUI(ns("docs")))
        )
      )
    ),

    # ---- Centre column: chat ----
    CopilotChatUI(ns("chat")),

    # ---- Right column: evidence ----
    CopilotEvidenceUI(ns("evidence"))
  )


  board <- OmicsBoardUI(
    id = ns("board"),
    #title = "AI Copilot",
    title = "Obi-Two Copilot",    
    info = FALSE,
    ui
  )
  return(board)
}

#' Copilot Board Sidebar Inputs (empty — left column is internal)
#'
#' @param id Module namespace id.
#' @return An empty `shiny::div()`.
#' @export
CopilotBoardInputs <- function(id) {
  shiny::div()
}
