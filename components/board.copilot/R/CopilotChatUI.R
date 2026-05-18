# CopilotChatUI.R — Chat Module UI
#
# Centre column of the copilot board: chat header + tier selector, chat body
# (shinychat::chat_ui). Height-constrained so the chat region scrolls
# internally rather than pushing the outer page taller.
#
# See .active_plans/refactor_copilot/chat/specs.md §CopilotChatUI.

# CSS lives inline so it travels with the module. The stop button is rendered
# at the same absolute position shinychat uses for its send button
# (.shiny-chat-btn-send — bottom: 6px; right: 8px inside shiny-chat-input),
# so when run_status == "streaming" the chat server hides the send button via
# CSS class and shows the stop button in the same spot, giving the
# send→stop morph UX without forking shinychat.
.COPILOT_CHAT_CSS <- "
.copilot-chat-card { display: flex; flex-direction: column; }
.copilot-chat-body { flex: 1 1 auto; min-height: 0; display: flex; }
.copilot-chat-body > shiny-chat-container { width: 100%; }
.copilot-chat-stop-wrap {
  position: absolute;
  bottom: 6px;
  right: 8px;
  z-index: 1000;
  pointer-events: auto;
}
.copilot-chat-stop-wrap .btn {
  background-color: var(--bs-danger, #dc3545);
  color: white;
  border: none;
  width: 28px;
  height: 28px;
  padding: 0;
  border-radius: 50%;
  display: inline-flex;
  align-items: center;
  justify-content: center;
  line-height: 1;
}
.copilot-chat-stop-wrap .btn:hover {
  background-color: var(--bs-danger, #b02a37);
}
.copilot-chat-tier-wrap {
  position: absolute;
  bottom: 6px;
  left: 8px;
  z-index: 1000;
  pointer-events: auto;
}
.copilot-chat-tier-wrap .selectize-control,
.copilot-chat-tier-wrap .form-group { margin: 0; }
.copilot-chat-tier-wrap .selectize-input {
  min-width: 140px;
  font-size: 0.85em;
  padding-block: 2px;
}
/* When streaming, hide shinychat's send button so the stop button reads
   as a morph rather than a sibling. */
.copilot-streaming shiny-chat-input .shiny-chat-btn-send { display: none; }
"

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
    class = "copilot-chat-card",
    style = "height: calc(100vh - 80px);",
    shiny::tags$head(shiny::tags$style(shiny::HTML(.COPILOT_CHAT_CSS))),
    shiny::div(
      class = "copilot-chat-body position-relative",
      shinychat::chat_ui(
        ns("chat"),
        height      = "100%",
        fill        = TRUE,
        placeholder = "Ask a question about your data…"
      ),
      # Tier selector — bottom-left, mirroring the stop button at bottom-right.
      shiny::div(
        id    = ns("tier_wrap"),
        class = "copilot-chat-tier-wrap",
        shiny::selectInput(
          ns("tier"),
          label   = NULL,
          choices = NULL,
          width   = "auto"
        )
      ),
      # Overlaid stop button — same coordinates as shinychat's send button.
      # Hidden by default; the chat server toggles visibility on run_status.
      shinyjs::hidden(
        shiny::div(
          id    = ns("stop_btn_wrap"),
          class = "copilot-chat-stop-wrap",
          shiny::actionButton(
            ns("stop_btn"),
            label = NULL,
            icon  = shiny::icon("stop"),
            title = "Stop"
          )
        )
      )
    )
  )
}
