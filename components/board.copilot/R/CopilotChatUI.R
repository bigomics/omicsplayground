# CopilotChatUI.R — Chat Module UI
#
# Centre column of the copilot board: chat header + tier selector, chat body
# (shinychat::chat_ui), footer with example buttons and hidden stop button.
#
# See .active_plans/refactor_copilot/chat/specs.md §CopilotChatUI.

#' Copilot Chat UI
#'
#' bslib::card composing the chat region. The tier selectInput is empty by
#' default; the chat server populates it via `updateSelectInput` driven by
#' the `tier_choices` reactive.
#'
#' @param id Module namespace id (same id passed to CopilotChatServer).
#' @return A `bslib::card`.
#' @export
CopilotChatUI <- function(id) {
  ns <- shiny::NS(id)

  bslib::card(
    bslib::card_header(
      shiny::div(
        class = "d-flex align-items-center gap-2",
        shiny::span("Copilot Chat", class = "me-auto fw-semibold"),
        shiny::selectInput(
          ns("tier"),
          label   = NULL,
          choices = NULL,
          width   = "auto"
        )
      )
    ),
    shinychat::chat_ui(
      ns("chat"),
      height      = "100%",
      fill        = TRUE,
      placeholder = "Ask a question about your data…"
    ),
    bslib::card_footer(
      shiny::div(
        class = "d-flex flex-wrap gap-2 align-items-center",
        shiny::actionButton(ns("ask_describe"),   "Describe data",   class = "btn-sm btn-outline-primary"),
        shiny::actionButton(ns("ask_findings"),   "Key findings",    class = "btn-sm btn-outline-primary"),
        shiny::actionButton(ns("ask_pathways"),   "Top pathways",    class = "btn-sm btn-outline-primary"),
        shiny::actionButton(ns("ask_biomarkers"), "Biomarker genes", class = "btn-sm btn-outline-primary"),
        shiny::actionButton(ns("ask_plot"),       "Show a plot",     class = "btn-sm btn-outline-primary"),
        shinyjs::hidden(
          shiny::actionButton(
            ns("stop_btn"),
            label = "Stop",
            icon  = shiny::icon("stop"),
            class = "btn-sm btn-danger ms-auto"
          )
        )
      )
    )
  )
}
