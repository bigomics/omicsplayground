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

#' Copilot Chat Module Server
#'
#' Owns shinychat I/O. Observes `chat_event` and dispatches `post`, `clear`,
#' `stream`, and `replay` events to shinychat. Forwards user input to
#' `on_user_message`.
#'
#' @param id Module namespace id.
#' @param on_user_message function(text) — called on every user submission.
#'   Should normally dispatch a `run_request_ask(text)` to the run controller.
#' @param chat_event reactive() or reactiveVal() returning list|NULL — event
#'   bus written by the orchestrator. Shapes: post / clear / reset / stream / replay.
#'   `reset` combines clear + post in a single flush-safe write:
#'   `list(type = "reset", role = "assistant", text = <greeting>)`.
#' @param run_status Optional reactive(character) — drives the send→stop
#'   morph: when `run_status() == "streaming"`, shinychat's send button is
#'   hidden via CSS and an overlaid stop button is shown.
#' @param on_abort Optional function(reason) — invoked when the stop button is
#'   clicked. Should normally dispatch a `run_request_abort()` to the run
#'   controller. If NULL the stop button still renders but does nothing.
#' @param tier_choices Optional reactive() returning named character vector —
#'   pushed into the tier selectInput on change.
#'
#' @return list(user_input, stream_done, last_error, on_tool_request, push_event=NULL)
#'   chat_event shapes handled: post / clear / reset / stream / replay.
#' @export
CopilotChatServer <- function(
  id,
  on_user_message,
  chat_event,
  run_status   = NULL,
  on_abort     = NULL,
  tier_choices = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Private state ----
    user_input_rv  <- shiny::reactiveVal("")
    stream_done_rv <- shiny::reactiveVal(0L)
    last_error_rv  <- shiny::reactiveVal(NULL)

    # ---- tier choices update ----
    if (!is.null(tier_choices)) {
      shiny::observe({
        ch <- tier_choices()
        shiny::req(ch)
        shiny::updateSelectInput(session, "tier", choices = ch)
      })
    }

    # ---- run_status: morph send button into stop button while streaming ----
    # When run_status() == "streaming" we add the .copilot-streaming class on
    # body (CSS hides shinychat's send button) and show the overlaid stop
    # button at the same coordinates. On any other state, reverse.
    if (!is.null(run_status)) {
      shiny::observe({
        streaming <- identical(run_status(), "streaming")
        if (streaming) {
          shinyjs::addClass(selector = "body", class = "copilot-streaming")
          shinyjs::show("stop_btn_wrap")
        } else {
          shinyjs::removeClass(selector = "body", class = "copilot-streaming")
          shinyjs::hide("stop_btn_wrap")
        }
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
        reset = {
          shinychat::chat_clear("chat")
          text <- event$text %||% ""
          if (nzchar(text)) {
            shinychat::chat_append_message(
              "chat",
              list(role = event$role %||% "assistant", content = text),
              chunk = FALSE
            )
          }
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

    # ---- Stop button (overlay on shinychat's send button while streaming) ----
    shiny::observeEvent(input$stop_btn, {
      log_info("copilot.chat.stop_clicked",
               value = as.integer(input$stop_btn %||% 0L),
               has_on_abort = is.function(on_abort))
      if (is.function(on_abort)) on_abort("User stopped the run")
    }, ignoreInit = TRUE)

    # ---- on_tool_request: inline collapsible tool marker ----
    # Passed by run_controller as on_tool_request argument to agent_prompt_stream.
    # Fires once per tool start, during streaming.
    .on_tool_request <- function(req) {
      log_trace("copilot.chat.tool_request",
        name     = req$name %||% "unknown",
        n_args   = length(req$arguments),
        arg_keys = paste(names(req$arguments) %||% character(0), collapse = ","))

      args_text <- tryCatch({
        args <- req$arguments
        if (length(args) == 0L) {
          "(no arguments)"
        } else {
          lines <- mapply(function(key, val) {
            formatted <- if (is.list(val) || (is.vector(val) && length(val) > 1L)) {
              if (length(val) > 5L) {
                short <- paste(format(val[seq_len(5L)]), collapse = ", ")
                paste0(short, ", …")
              } else {
                paste(format(val), collapse = ", ")
              }
            } else {
              format(val)
            }
            paste0(key, " = ", formatted)
          }, names(args), args, SIMPLIFY = TRUE, USE.NAMES = FALSE)
          paste(lines, collapse = "\n")
        }
      }, error = function(e) {
        tryCatch(
          jsonlite::toJSON(req$arguments, auto_unbox = TRUE, pretty = TRUE),
          error = function(e2) "(error formatting arguments)"
        )
      })

      marker <- paste0(
        '<details><summary>\U0001F527 ', req$name %||% "tool", '</summary>',
        '<pre>', args_text, '</pre>',
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
