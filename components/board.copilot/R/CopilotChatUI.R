# CopilotChatUI.R — Chat Module UI
#
# Centre column of the copilot board: chat header + tier selector, chat body
# (shinychat::chat_ui). Height-constrained so the chat region scrolls
# internally rather than pushing the outer page taller.
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
    style = "height: calc(100vh - 80px);",
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
    shiny::div(
      style = "flex: 1; min-height: 0;",
      shinychat::chat_ui(
        ns("chat"),
        height      = "calc(100vh - 240px)",
        fill        = TRUE,
        placeholder = "Ask a question about your data…"
      )
    )
  )
}
