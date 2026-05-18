# CopilotChatServer.R — Chat Module Server
#
# Owns all shinychat I/O for the copilot board. The orchestrator pushes events
# into the `chat_event` reactiveVal; this module observes them and dispatches
# to shinychat::chat_append / chat_append_message / chat_clear. The module
# never constructs events itself — that is the orchestrator's job.
#
# See .active_plans/refactor_copilot/chat/specs.md for the full contract.

# ---- Null-coalescing operator ----
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Example questions ----
.COPILOT_EXAMPLE_QUESTIONS <- list(
  describe   = "Describe my experiment and the main comparisons being made",
  findings   = "Summarize the main findings of this experiment",
  pathways   = "List the top enriched pathways and explain how they relate to the experiment",
  biomarkers = "Show the top candidate biomarkers for this experiment",
  plot       = "Show me a PCA plot of the samples"
)

#' Copilot Chat Module Server
#'
#' Owns shinychat I/O. Observes `chat_event` and dispatches `post`, `clear`,
#' `stream`, and `replay` events to shinychat. Forwards user input and
#' example-button clicks to `on_user_message`.
#'
#' @param id Module namespace id.
#' @param on_user_message function(text) — called on every user submission
#'   (text input or example button). Should normally dispatch a
#'   `run_request_ask(text)` to the run controller.
#' @param chat_event reactive() or reactiveVal() returning list|NULL — event
#'   bus written by the orchestrator. Shapes: post / clear / stream / replay.
#' @param run_status Optional reactive(character) — drives stop_btn visibility
#'   and example-button disabling when value is `"streaming"`.
#' @param tier_choices Optional reactive() returning named character vector —
#'   pushed into the tier selectInput on change.
#'
#' @return list(user_input, stream_done, last_error, on_tool_request, push_event=NULL)
#' @export
CopilotChatServer <- function(
  id,
  on_user_message,
  chat_event,
  run_status   = NULL,
  tier_choices = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Private state ----
    user_input_rv  <- shiny::reactiveVal("")
    stream_done_rv <- shiny::reactiveVal(0L)
    last_error_rv  <- shiny::reactiveVal(NULL)

    # ---- run_status: disable example buttons + toggle stop_btn ----
    if (!is.null(run_status)) {
      shiny::observe({
        streaming <- identical(run_status(), "streaming")
        if (streaming) {
          shinyjs::disable("ask_describe");   shinyjs::disable("ask_findings")
          shinyjs::disable("ask_pathways");   shinyjs::disable("ask_biomarkers")
          shinyjs::disable("ask_plot")
          shinyjs::show("stop_btn")
        } else {
          shinyjs::enable("ask_describe");    shinyjs::enable("ask_findings")
          shinyjs::enable("ask_pathways");    shinyjs::enable("ask_biomarkers")
          shinyjs::enable("ask_plot")
          shinyjs::hide("stop_btn")
        }
      })
    }

    # ---- tier choices update ----
    if (!is.null(tier_choices)) {
      shiny::observe({
        ch <- tier_choices()
        shiny::req(ch)
        shiny::updateSelectInput(session, "tier", choices = ch)
      })
    }

    # ---- Internal post helper ----
    .post <- function(role, text) {
      shinychat::chat_append_message(
        "chat",
        list(role = role, content = text),
        chunk = FALSE
      )
    }

    # ---- chat_event observer ----
    shiny::observeEvent(chat_event(), ignoreNULL = TRUE, {
      event <- chat_event()
      if (is.null(event) || is.null(event$type)) return()
      switch(event$type,
        post = {
          .post(event$role %||% "assistant", event$text %||% "")
        },
        clear = {
          shinychat::chat_clear("chat")
        },
        stream = {
          result <- shinychat::chat_append("chat", event$async_gen)
          promises::then(result,
            onFulfilled = function(v) {
              stream_done_rv(shiny::isolate(stream_done_rv()) + 1L)
            },
            onRejected = function(e) {
              msg <- conditionMessage(e)
              last_error_rv(msg)
              .post("assistant", copilot_msg("error_prefix", msg = msg))
              stream_done_rv(shiny::isolate(stream_done_rv()) + 1L)
            }
          )
        },
        replay = {
          records <- event$records %||% list()
          for (rec in records) {
            tryCatch(
              shinychat::chat_append_message(
                "chat",
                list(role = rec@role, content = rec@content_text),
                chunk = FALSE
              ),
              error = function(e) {
                log_info("copilot.chat.replay_record_failed",
                          msg = conditionMessage(e))
              }
            )
          }
        }
      )
    })

    # ---- User text input ----
    shiny::observeEvent(input$chat_user_input, {
      shiny::req(input$chat_user_input)
      text <- trimws(input$chat_user_input)
      shiny::req(nzchar(text))
      user_input_rv(text)
      on_user_message(text)
    })

    # ---- Example buttons ----
    .ask <- function(text) {
      user_input_rv(text)
      on_user_message(text)
    }
    shiny::observeEvent(input$ask_describe,   .ask(.COPILOT_EXAMPLE_QUESTIONS$describe))
    shiny::observeEvent(input$ask_findings,   .ask(.COPILOT_EXAMPLE_QUESTIONS$findings))
    shiny::observeEvent(input$ask_pathways,   .ask(.COPILOT_EXAMPLE_QUESTIONS$pathways))
    shiny::observeEvent(input$ask_biomarkers, .ask(.COPILOT_EXAMPLE_QUESTIONS$biomarkers))
    shiny::observeEvent(input$ask_plot,       .ask(.COPILOT_EXAMPLE_QUESTIONS$plot))

    # ---- on_tool_request: inline collapsible tool marker ----
    # Passed by run_controller as on_tool_request argument to agent_prompt_stream.
    # Fires once per tool start, during streaming.
    .on_tool_request <- function(req) {
      args_json <- tryCatch(
        jsonlite::toJSON(req$arguments, auto_unbox = TRUE, pretty = TRUE),
        error = function(e) "{}"
      )
      marker <- paste0(
        '<details><summary>\U0001F527 ', req$name %||% "tool", '</summary>',
        '<pre>', args_json, '</pre>',
        '</details>'
      )
      tryCatch(
        shinychat::chat_append_message(
          "chat",
          list(role = "assistant", content = marker),
          chunk = TRUE
        ),
        error = function(e) {
          log_info("copilot.chat.tool_marker_failed", msg = conditionMessage(e))
        }
      )
    }

    # ---- Public API ----
    list(
      user_input      = shiny::reactive(user_input_rv()),
      stream_done     = shiny::reactive(stream_done_rv()),
      last_error      = shiny::reactive(last_error_rv()),
      on_tool_request = .on_tool_request,
      push_event      = NULL
    )
  })
}
