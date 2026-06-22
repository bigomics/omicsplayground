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
  left: 16px;
  z-index: 1000;
  pointer-events: auto;
}
/* Grow shinychat's textarea so the tier trigger overlay has room to sit
   inside the input pill without colliding with the placeholder text. */
shiny-chat-input textarea { min-height: 42px !important; }
.copilot-chat-tier-wrap .copilot-tier-trigger {
  color: var(--bs-secondary-color, #6c757d);
  text-decoration: none;
  font-size: 0.85em;
  padding: 2px 6px;
  border-radius: 4px;
}
.copilot-chat-tier-wrap .copilot-tier-trigger:hover {
  background-color: var(--bs-tertiary-bg, rgba(0,0,0,0.05));
  text-decoration: none;
}
/* When streaming, hide shinychat's send button so the stop button reads
   as a morph rather than a sibling. */
.copilot-streaming shiny-chat-input .shiny-chat-btn-send { display: none; }
/* Starter suggestions — pushed inline as an assistant chat message holding
   <ul class='copilot-starter-list'><li class='suggestion submit'>…</li>…</ul>.
   shinychat handles auto-submit on click; this CSS just stacks one per row
   and styles the pills. */
.copilot-starter-list {
  list-style: none;
  padding: 0;
  margin: 4px 0 0 0;
  display: flex;
  flex-direction: column;
  align-items: flex-start;
  gap: 6px;
}
.copilot-starter-list li.suggestion {
  background-color: var(--bs-tertiary-bg, rgba(0,0,0,0.04));
  color: var(--bs-secondary-color, #6c757d);
  border: 1px solid var(--bs-border-color, rgba(0,0,0,0.1));
  border-radius: 16px;
  font-size: 0.85em;
  padding: 2px 12px;
  line-height: 1.4;
  cursor: pointer;
  display: inline-block;
}
.copilot-starter-list li.suggestion:hover {
  background-color: var(--bs-secondary-bg, rgba(0,0,0,0.08));
  color: var(--bs-body-color, #212529);
}
"


#' Copilot Chat UI
#'
#' bslib::card composing the chat region. The tier popover trigger is a
#' borderless actionLink; the server populates choices via `updateRadioButtons`
#' and renders the active label via `tier_label` uiOutput.
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
    shiny::tags$div(
      class = "alert alert-warning alert-dismissible fade show",
      role = "alert",
      shiny::HTML(paste0(
        "<b>AI-generated content.</b> This experimental assistant may produce ",
        "incomplete or inaccurate interpretations. Verify important findings ",
        "independently; prompts and selected dataset context may be processed ",
        "by the configured LLM provider."
      )),
      shiny::tags$button(
        type = "button",
        class = "btn-close",
        `data-bs-dismiss` = "alert",
        `aria-label` = "Close"
      )
    ),
    shiny::div(
      class = "copilot-chat-body position-relative",
      shinychat::chat_ui(
        ns("chat"),
        width       = "100%",
        height      = "100%",
        fill        = TRUE,
        placeholder = "Ask a question about your data…"
      ),
      # Tier selector — bottom-left, mirroring the stop button at bottom-right.
      # Borderless actionLink trigger opens a popover with radioButtons so the
      # selection lives in the chat module's own namespace (fixing the dead
      # observer that previously read input$tier from the board namespace).
      ## shiny::div(
      ##   id    = ns("tier_wrap"),
      ##   class = "copilot-chat-tier-wrap",
      ##   bslib::popover(
      ##     shiny::actionLink(
      ##       ns("tier_btn"),
      ##       label = shiny::uiOutput(ns("tier_label"), inline = TRUE),
      ##       class = "copilot-tier-trigger"
      ##     ),
      ##     shiny::radioButtons(
      ##       ns("tier_choice"),
      ##       label        = NULL,
      ##       choiceNames  = unname(vapply(COPILOT_TIERS, copilot_tier_label, character(1))),
      ##       choiceValues = unname(COPILOT_TIERS),
      ##       selected     = COPILOT_TIERS[[1]]
      ##     ),
      ##     placement = "top",
      ##     id        = ns("tier_pop")
      ##   )
      ## ),
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

#' Copilot chat settings panel (tier + style controls)
#'
#' Renders the user-facing knobs for the chat module. Hosts two radio groups
#' plus an optional free-text style box. Input ids — `tier_choice`,
#' `style_choice`, `custom_text` — are namespaced under the chat module so
#' `CopilotChatServer()` observes them directly.
#'
#' The custom-style textbox is hidden via `shinyjs::hidden` and toggled on
#' when `style_choice == "custom"`. A small caption under the style radio
#' reminds the user that the choice applies to the next message.
#' @export
CopilotChatSettings <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::br(),
    shiny::radioButtons(
      ns("tier_choice"), "Model tier:",
      choiceNames  = unname(vapply(COPILOT_TIERS, copilot_tier_label, character(1))),
      choiceValues = unname(COPILOT_TIERS),
      selected     = COPILOT_TIERS[[1]],
      inline       = TRUE
    ),
    shiny::br(),
    shiny::radioButtons(
      ns("style_choice"), "Answer style:",
      choiceNames  = unname(COPILOT_STYLE_LABELS),
      choiceValues = unname(COPILOT_STYLES),
      selected     = COPILOT_STYLES[[1]],
      inline       = TRUE
    ),
    shinyjs::hidden(
      shiny::div(
        id = ns("custom_text_wrap"),
        shiny::textAreaInput(
          ns("custom_text"),
          label       = "Custom style instructions:",
          value       = "",
          rows        = 3,
          width       = "100%",
          placeholder = "e.g. Be punchy. Use one-line answers."
        )
      )
    ),
    shiny::tags$small(
      class = "text-muted",
      "Style changes apply to the next message."
    )
  )
}
