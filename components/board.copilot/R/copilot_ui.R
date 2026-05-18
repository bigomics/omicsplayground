# copilot_ui.R — Copilot Board UI Shell
#
# Three-column bslib layout. Left: datasets/history/docs tabset.
# Centre: chat region. Right: evidence placeholder.
#
# Phase 1: module slots are empty uiOutput() placeholders. Later phases slot
# in CopilotDatasetsUI, CopilotHistoryUI, CopilotDocsUI, CopilotEvidenceUI,
# and extract the chat region into CopilotChatUI (Phase 4).
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

    # ---- Centre column: chat ----
    bslib::card(
      bslib::card_header(
        shiny::div(
          class = "d-flex align-items-center gap-2",
          shiny::span("Copilot Chat", class = "me-auto fw-semibold"),
          shiny::selectInput(
            ns("tier"),
            label   = NULL,
            choices = .COPILOT_TIER_IDS,
            width   = "auto"
          )
        )
      ),
      # TODO(phase 4): replace with CopilotChatUI(ns("chat")) once chat module exists
      shinychat::chat_ui(
        ns("chat"),
        height       = "100%",
        fill         = TRUE,
        placeholder  = "Ask a question about your data…"
      ),
      bslib::card_footer(
        shiny::div(
          class = "d-flex flex-wrap gap-2",
          shiny::actionButton(ns("ask_describe"),   "Describe data",     class = "btn-sm btn-outline-primary"),
          shiny::actionButton(ns("ask_findings"),   "Key findings",      class = "btn-sm btn-outline-primary"),
          shiny::actionButton(ns("ask_pathways"),   "Top pathways",      class = "btn-sm btn-outline-primary"),
          shiny::actionButton(ns("ask_biomarkers"), "Biomarker genes",   class = "btn-sm btn-outline-primary"),
          shiny::actionButton(ns("ask_plot"),       "Show a plot",       class = "btn-sm btn-outline-primary")
        )
      ),
      shinyjs::hidden(
        shiny::actionButton(
          ns("stop_btn"),
          label = "Stop",
          icon  = shiny::icon("stop"),
          class = "btn-danger btn-sm position-absolute"
        )
      )
    ),

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
