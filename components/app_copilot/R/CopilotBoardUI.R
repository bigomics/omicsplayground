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

    # ---- Left column: new-chat button, on top of docs / reports ----
    # "New chat" lives here (own row) rather than inside the
    # Dataset/Settings/History card on the right, which already has little
    # vertical room to spare. Plain flexbox instead of layout_columns
    # row_heights, so the button keeps its intrinsic size instead of being
    # stretched to fill a grid row.
    shiny::div(
      style = "height: 100%; display: flex; flex-direction: column; gap: 8px;",
      shiny::actionButton(
        ns("new_chat"),
        label = "New chat",
        icon  = shiny::icon("plus"),
        class = "btn-sm btn-outline-secondary w-100",
        style = "flex: 0 0 auto;"
      ),
      shiny::div(
        style = "flex: 1 1 auto; min-height: 0;",
        # Header mirrors the "Plot output" card header exactly (same
        # `as.card_item(fillRow(height = "33px", class = "plotmodule-header"))`
        # shape, NOT card_header) rather than a navset_underline nav_panel —
        # a single-tab nav_panel rendered as a clickable-looking tab next to
        # "Documents" / "Reports", which read as a visual button and
        # confused reviewers. as.card_item is required for the header to
        # span the card edge-to-edge — a plain child gets wrapped in a
        # padded card_body, which insets the header and its border-bottom.
        bslib::card(
          height = "100%",
          bslib::as.card_item(shiny::fillRow(
            flex = 1,
            class = "plotmodule-header",
            height = "33px",
            shiny::div(
              class = "plotmodule-title",
              style = "white-space: nowrap; overflow: hidden; text-overflow: clip;",
              "Workspace"
            )
          )),
          bslib::navset_underline(
            bslib::nav_panel("Documents", CopilotDocsUI(ns("docs"))),
            bslib::nav_panel("Reports", CopilotReportsUI(ns("reports"))),
            bslib::nav_panel("Settings", CopilotChatSettings(ns("chat"))),
            bslib::nav_panel("History", CopilotHistoryUI(ns("history")))
          )
        )
      )
    ),

    # ---- Centre column: chat ----
    CopilotChatUI(ns("chat")),

    # ---- Right column: dataset info, on top of evidence ----
    # "Dataset" is the only remaining item here now that Settings/History
    # moved left, so it gets the same plain "Plot output"-style header as
    # Context instead of a navset_underline — a lone tab would reintroduce
    # the clickable-looking-button problem we just fixed on the left.
    bslib::layout_columns(
      col_widths = 12,
      row_heights = c(4, 7),
      bslib::card(
        height = "100%",
        bslib::as.card_item(shiny::fillRow(
          flex = 1,
          class = "plotmodule-header",
          height = "33px",
          shiny::div(
            class = "plotmodule-title",
            style = "white-space: nowrap; overflow: hidden; text-overflow: clip;",
            "Dataset"
          )
        )),
        shiny::uiOutput(ns("dataset_info"))
      ),
      CopilotEvidenceUI(ns("evidence"))
    )
  )


  board <- OmicsBoardUI(
    id = ns("board"),
    #title = "AI Copilot",
    title = "Obi",
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
