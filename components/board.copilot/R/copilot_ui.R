#' Copilot Board UI
#'
#' Three-column layout: datasets/history/docs | chat | evidence

CopilotBoardInputs <- function(id) {
  ns <- shiny::NS(id)
  ## Copilot has no sidebar settings — left column is internal
  shiny::div()
}

CopilotBoardUI <- function(id) {
  ns <- shiny::NS(id)

  ## Left column: tabbed navigation for datasets, history, docs
  left_panel <- bslib::card(
    full_screen = FALSE,
    class = "border-0",
    style = "height: calc(100vh - 80px); overflow-y: auto;",
    bslib::navset_underline(
      id = ns("left_tabs"),
      bslib::nav_panel(
        "Datasets",
        copilot_panel_datasets_ui(ns("datasets"))
      ),
      bslib::nav_panel(
        "History",
        copilot_panel_history_ui(ns("history"))
      ),
      bslib::nav_panel(
        "Docs",
        copilot_panel_docs_ui(ns("docs"))
      )
    )
  )

  ## Center column: chat interface
  center_panel <- bslib::card(
    class = "border-0",
    style = "height: calc(100vh - 80px);",
    bslib::card_header(
      class = "d-flex justify-content-between align-items-center",
      shiny::span("Copilot Chat", class = "fw-bold"),
      shiny::div(
        class = "d-flex gap-2 align-items-center",
        shiny::selectInput(
          inputId = ns("tier"),
          label   = NULL,
          choices = NULL,
          width   = "180px"
        )
      )
    ),
    shiny::div(
      class = "p-2",
      style = "flex: 1;",
      shinychat::chat_ui(ns("chat"), width = "100%", height = "calc(100vh - 240px)", fill = TRUE)
    ),
    bslib::card_footer(
      class = "d-flex flex-wrap gap-1 p-2",
      shiny::actionButton(ns("ask_describe"), "Describe experiment", class = "btn-sm btn-outline-primary"),
      shiny::actionButton(ns("ask_findings"), "Summarize findings", class = "btn-sm btn-outline-primary"),
      shiny::actionButton(ns("ask_pathways"), "Top pathways", class = "btn-sm btn-outline-primary"),
      shiny::actionButton(ns("ask_biomarkers"), "Top biomarkers", class = "btn-sm btn-outline-primary"),
      shiny::actionButton(ns("ask_plot"), "Show PCA plot", class = "btn-sm btn-outline-primary")
    )
  )

  ## Right column: evidence panel
  right_panel <- shiny::div(
    style = "height: calc(100vh - 80px); overflow-y: auto;",
    copilot_panel_evidence_ui(ns("evidence"))
  )

  ## Three-column layout
  bslib::layout_columns(
    col_widths = c(3, 5, 4),
    style = "height: calc(100vh - 80px); padding: 8px;",
    fill = TRUE,
    left_panel,
    center_panel,
    right_panel
  )
}
