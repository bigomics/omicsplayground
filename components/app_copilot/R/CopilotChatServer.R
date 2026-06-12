# CopilotChatServer.R — Chat Module Server
#
# Owns all shinychat I/O for the copilot board. The orchestrator pushes events
# into the `chat_event` reactiveVal; this module observes them and dispatches
# to shinychat::chat_append / chat_append_message / chat_clear. The module
# never constructs events itself — that is the orchestrator's job.
#
# See .active_plans/refactor_copilot/chat/specs.md for the full contract.

#' Copilot Chat Module Server
#'
#' Owns shinychat I/O. Observes `chat_event` and dispatches `post`, `clear`,
#' `reset`, and `replay` events to shinychat. Forwards user input to
#' `on_user_message`.
#'
#' @param id Module namespace id.
#' @param on_user_message function(text) -- called on every user submission.
#'   Should normally dispatch a `run_request_ask(text)` to the run controller.
#' @param chat_event reactive() or reactiveVal() returning list|NULL -- event
#'   bus written by the orchestrator. Shapes: post / clear / reset / stream / replay.
#'   `reset` combines clear + post in a single flush-safe write:
#'   `list(type = "reset", role = "assistant", text = <greeting>)`.
#' @param run_status Optional reactive(character) -- drives the send/stop
#'   morph: when `run_status() == "streaming"`, shinychat's send button is
#'   hidden via CSS and an overlaid stop button is shown.
#' @param on_abort Optional function(reason) -- invoked when the stop button is
#'   clicked. Should normally dispatch a `run_request_abort()` to the run
#'   controller. If NULL the stop button still renders but does nothing.
#' @param tier_choices Optional reactive() returning named character vector --
#'   pushed into the tier radioButtons on change.
#' @param current_tier Optional reactive()/reactiveVal() returning the active
#'   tier id -- used to keep the radio selection and the trigger label in sync
#'   with the orchestrator's state.
#' @param starters Optional reactive() returning a character vector of starter
#'   questions to render as a one-shot button strip above the chat input.
#'   The strip is dismissed via `shinyjs::hide` on the first user message
#'   (typed or button) and re-shown on a `reset` chat event. When NULL or
#'   empty, no strip is rendered.
#'
#' @return list(user_input, stream_done, last_error, on_tool_request,
#'   push_event=NULL, tier_clicked) where `tier_clicked` is a reactiveVal
#'   that updates whenever the user picks a new tier in the popover.
#'   chat_event shapes handled: post / clear / reset / stream / replay.
#' @export
CopilotChatServer <- function(
  id,
  on_user_message,
  chat_event,
  run_status    = NULL,
  on_abort      = NULL,
  tier_choices  = NULL,
  current_tier  = NULL,
  style_choices = NULL,
  current_style = NULL,
  starters      = NULL
) {
  shiny::moduleServer(id, function(input, output, session) {

    # ---- Private state ----
    user_input_rv  <- shiny::reactiveVal("")
    stream_done_rv <- shiny::reactiveVal(0L)
    last_error_rv  <- shiny::reactiveVal(NULL)
    tier_clicked   <- shiny::reactiveVal(NULL)
    style_clicked  <- shiny::reactiveVal(NULL)
    custom_text_rv <- shiny::reactiveVal("")

    # ---- tier choices update ----
    if (!is.null(tier_choices)) {
      shiny::observe({
        ch <- tier_choices()
        shiny::req(ch)
        shiny::updateRadioButtons(
          session, "tier_choice",
          choiceNames  = names(ch),
          choiceValues = unname(ch),
          selected     = shiny::isolate(
            if (!is.null(current_tier)) current_tier() else character(0)
          )
        )
      })
    }

    # ---- tier label output ----
    output$tier_label <- shiny::renderUI({
      shiny::req(current_tier, tier_choices)
      cur <- current_tier()
      if (is.null(cur) || !nzchar(cur)) return(NULL)
      ch  <- tier_choices()
      lbl <- names(ch)[match(cur, unname(ch))]
      if (is.na(lbl)) lbl <- cur
      shiny::tagList(shiny::span(lbl), shiny::span(" ▾", class = "copilot-tier-chevron"))
    })

    # ---- tier_choice observer ----
    shiny::observeEvent(input$tier_choice, {
      shiny::req(input$tier_choice)
      tier_clicked(input$tier_choice)
    }, ignoreInit = TRUE)

    # ---- style choices update ----
    if (!is.null(style_choices)) {
      shiny::observe({
        ch <- style_choices()
        shiny::req(ch)
        shiny::updateRadioButtons(
          session, "style_choice",
          choiceNames  = unname(ch),
          choiceValues = names(ch),
          selected     = shiny::isolate(
            if (!is.null(current_style)) current_style() else character(0)
          )
        )
      })
    }

    # ---- style_choice observer ----
    shiny::observeEvent(input$style_choice, {
      shiny::req(input$style_choice)
      style_clicked(input$style_choice)
      # Toggle the custom-style textbox visibility.
      if (identical(input$style_choice, "custom")) {
        shinyjs::show("custom_text_wrap")
      } else {
        shinyjs::hide("custom_text_wrap")
      }
    }, ignoreInit = FALSE)

    # ---- custom_text observer (debounced via reactive) ----
    shiny::observeEvent(input$custom_text, {
      custom_text_rv(input$custom_text %||% "")
    }, ignoreInit = TRUE, ignoreNULL = FALSE)

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

    # ---- Starter suggestions ----
    # Pushed inline as an assistant chat message after every greeting:
    #   <ul class='copilot-starter-list'>
    #     <li class='suggestion submit'>...</li>
    #   </ul>
    # shinychat handles the click → auto-submit into input$chat_user_input,
    # which flows through the existing user-input observer.
    .starter_bubble <- function() {
      if (is.null(starters)) return("")
      qs <- tryCatch(starters(), error = function(e) NULL)
      if (is.null(qs) || length(qs) == 0L) return("")
      qs <- as.character(qs)
      escaped <- htmltools::htmlEscape(qs)
      items <- paste0(
        "<li class='suggestion submit'>", escaped, "</li>",
        collapse = ""
      )
      paste0("<ul class='copilot-starter-list-DISABLED'>", items, "</ul>")
    }

    .post_starters <- function() {
      html <- shiny::isolate(.starter_bubble())
      if (!nzchar(html)) return(invisible())
      shinychat::chat_append_message(
        "chat",
        list(role = "assistant", content = html),
        chunk = FALSE
      )
      log_trace("copilot.starter.posted", n = length(shiny::isolate(starters())))
    }

    # ---- Internal post helper ----
    .post <- function(role, text) {
      shinychat::chat_append_message(
        "chat",
        list(role = role, content = text),
        chunk = FALSE
      )
    }

    # First-greeting one-shot: the board emits the initial greeting as a
    # `post` event from its onFlushed hook. We can't distinguish "greeting"
    # from other assistant posts, so we tag the first ever event as the
    # greeting and append starters once.
    first_event_seen <- shiny::reactiveVal(FALSE)

    # ---- chat_event observer ----
    shiny::observeEvent(chat_event(), ignoreNULL = TRUE, {
      event <- chat_event()
      if (is.null(event) || is.null(event$type)) return()
      is_first <- !shiny::isolate(first_event_seen())
      if (is_first) first_event_seen(TRUE)
      switch(event$type,
        post = {
          .post(event$role %||% "assistant", event$text %||% "")
          if (is_first) .post_starters()
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
          # Re-seed starter pills for the new chat.
          .post_starters()
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
          # Atomic clear+replay: the restore controller sets clear_first
          # so the "Restoring…" placeholder is dropped in the same observer
          # firing as the records append. A separate `clear` event would
          # coalesce with this one at the reactiveVal layer and be lost.
          if (isTRUE(event$clear_first)) {
            shinychat::chat_clear("chat")
          }
          # Prefer content_text_visible (no injected preamble) so the user
          # sees only what they typed on replay; fall back to content_text
          # for records saved before omicsagentovi 0.5.0 added the split.
          .visible <- function(r) {
            v <- tryCatch(r@content_text_visible, error = function(e) NULL)
            if (!is.null(v) && length(v) && !is.na(v) && nzchar(v)) v
            else r@content_text
          }

          mode <- tryCatch(copilot_replay_mode(), error = function(e) "single")
          if (identical(mode, "batch") &&
              !is.null(session) &&
              is.function(session$sendCustomMessage) &&
              length(records) > 0L) {
            payload <- lapply(records, function(r) {
              list(role = as.character(r@role),
                   content = as.character(.visible(r) %||% ""))
            })
            tryCatch(
              session$sendCustomMessage("copilot-chat-batch-append", list(
                id         = session$ns("chat"),
                messages   = payload,
                batch_size = 32L
              )),
              error = function(e) {
                log_info("copilot.chat.replay_batch_failed",
                          msg = conditionMessage(e))
              }
            )
          } else {
            for (rec in records) {
              tryCatch(
                shinychat::chat_append_message(
                  "chat",
                  list(role = rec@role, content = .visible(rec)),
                  chunk = FALSE
                ),
                error = function(e) {
                  log_info("copilot.chat.replay_record_failed",
                            msg = conditionMessage(e))
                }
              )
            }
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
      copilot_debug_timing("chat.tool_request.enter",
        name = req$name %||% "unknown",
        n_args = length(req$arguments))
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
      copilot_debug_timing("chat.tool_request.after_format",
        name = req$name %||% "unknown",
        args_nchar = nchar(args_text, type = "chars"))

      marker <- paste0(
        '<details><summary>\U0001F527 ', req$name %||% "tool", '</summary>',
        '<pre>', args_text, '</pre>',
        '</details>'
      )
      copilot_debug_timing("chat.tool_request.before_append",
        name = req$name %||% "unknown",
        marker_nchar = nchar(marker, type = "chars"))
      tryCatch(
        shinychat::chat_append_message(
          "chat",
          list(role = "assistant", content = marker),
          chunk = TRUE
        ),
        error = function(e) {
          log_info("copilot.chat.tool_marker_failed", msg = conditionMessage(e))
          copilot_debug_timing("chat.tool_request.append_error",
            name = req$name %||% "unknown",
            msg = conditionMessage(e))
        }
      )
      copilot_debug_timing("chat.tool_request.exit",
        name = req$name %||% "unknown")
    }

    # ---- Public API ----
    list(
      user_input      = shiny::reactive(user_input_rv()),
      stream_done     = shiny::reactive(stream_done_rv()),
      last_error      = shiny::reactive(last_error_rv()),
      on_tool_request = .on_tool_request,
      push_event      = NULL,
      tier_clicked    = tier_clicked,
      style_clicked   = style_clicked,
      custom_text     = shiny::reactive(custom_text_rv())
    )
  })
}
