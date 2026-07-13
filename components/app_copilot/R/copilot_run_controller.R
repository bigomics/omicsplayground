# copilot_run_controller.R — Run Lifecycle Controller
#
# Unified dispatch for ask / new_chat / tier / abort. Owns `run_status`
# transitions. Calls `agent_prompt_stream` with `on_done` and
# `on_tool_request` callbacks. Constructs fresh Agents on new_chat / tier
# reset and on first dataset load.
#
# See .active_plans/refactor_copilot/run_controller/specs.md for the full
# contract.

# ---- AI provider / credential resolution (shared with restore controller) ----
#
# The AI Features settings tab (AppSettingsBoard) writes the user's provider
# choice and a session-only credential closure into `session$userData` via
# setUserOption(). `session$userData` is shared across every module session,
# so the copilot reads here exactly what the settings module wrote — the same
# cross-module channel get_ai_model() already relies on.
#
# Returns the provider id (default "bigomics") and the nullary credential
# closure `function() key` (or NULL for BigOmics / no key). BigOmics keeps the
# env-var path (credentials NULL), so existing deployments are unaffected.

#' @noRd
.copilot_ai_provider <- function(session) {
  list(
    provider    = getUserOption(session, "ai_provider") %||% "bigomics",
    credentials = get_ai_credentials(session)
  )
}

# Build the provider/credentials/model args to splice into omicsagentovi::Agent()
# for a fresh build at `tier`.
#
#   - BigOmics (or NULL provider): today's exact behaviour — the tier path only,
#     no credentials, no model. The managed backend resolves the tier to
#     gpt-5.4-nano / gpt-5.4-mini via omicsagentovi's own tier table.
#   - BYOK provider with a per-tier menu model selected: a provider-prefixed
#     model id (the catalog stores bare ids for BYOK providers) plus the user's
#     credential closure. Agent infers the provider from the model prefix.
#   - BYOK provider with an empty/unset menu: fall back to the provider-aware
#     tier table (Agent resolves `tier` for `provider`), still with the key.
#
#' @noRd
.copilot_agent_build_args <- function(session, tier) {
  sel <- .copilot_ai_provider(session)
  if (identical(sel$provider, "bigomics")) {
    return(list(tier = tier))
  }
  menu_key   <- if (identical(tier, "copilot-deep")) "llm_copilot_deep"
                else "llm_copilot_balanced"
  menu_model <- getUserOption(session, menu_key)
  if (!is.null(menu_model) && nzchar(menu_model)) {
    list(model = paste0(sel$provider, ":", menu_model),
         credentials = sel$credentials)
  } else {
    list(tier = tier, provider = sel$provider, credentials = sel$credentials)
  }
}

# ---- Request constructors (pure helpers — no reactives) ----

#' @export
run_request_ask <- function(text, show_user_msg = TRUE) {
  list(kind = "ask", text = text, show_user_msg = isTRUE(show_user_msg))
}

#' @export
run_request_new_chat <- function() {
  list(kind = "new_chat")
}

#' @export
run_request_tier <- function(new_tier) {
  list(kind = "tier", new_tier = new_tier)
}

#' @export
run_request_abort <- function(reason = "User stopped the run") {
  list(kind = "abort", reason = reason)
}

# --------------------------------------------------------------------------
# Factory
# --------------------------------------------------------------------------

#' Copilot Run Controller
#'
#' Returns a single `dispatch(request)` entry point. All run-related UI
#' events (chat input, example buttons, new chat, tier change, stop) flow
#' through this dispatch.
#'
#' @param agent reactiveVal(Agent|NULL) — read + write.
#' @param run_status reactiveVal(character) — owned by this controller after
#'   construction. Transitions: "idle" -> "streaming" -> result$status.
#' @param tier reactiveVal(character) — current tier.
#' @param save_ctrl list returned by copilot_save_controller().
#' @param restore_ctrl list returned by copilot_restore_controller(), or NULL.
#'   Currently unused — restore is dispatched directly from the orchestrator;
#'   kept here for API symmetry.
#' @param evidence list with $clear() or NULL. NULL until Phase 5.
#' @param store SessionStore — for tier-change pre-save.
#' @param chat_event reactiveVal — chat event bus. The run controller pushes
#'   `post`, `clear`, and `stream` events here.
#' @param chat_on_tool_request function(req) — owned by chat module; passed
#'   to agent_prompt_stream as on_tool_request.
#' @param pgx reactiveValues — global PGX.
#' @param pgx_dir character — directory containing .pgx files.
#' @param docs_dir character — directory for uploaded docs.
#' @param pgx_loaded_event reactiveVal — for bindings construction.
#' @param maxturns integer or Inf — refuse asks beyond this count.
#' @param session Shiny session.
#'
#' @return list(dispatch, apply_dataset)
#' @export
copilot_run_controller <- function(
  agent,
  run_status,
  tier,
  save_ctrl,
  restore_ctrl     = NULL,
  evidence         = NULL,
  store,
  chat_event,
  chat_on_tool_request,
  pgx,
  pgx_dir,
  docs_dir,
  pgx_loaded_event,
  maxturns         = Inf,
  session,
  style            = NULL,
  custom           = NULL,
  report_context   = NULL,
  doc_context      = NULL,
  tools_enabled    = NULL
) {

  # ---- Follow-up suggestion generator (side LLM call, runtime-optional) ----
  # Returns NULL if `omicsai` is not installed — downstream code guards on
  # this and skips the feature silently.
  followup_gen <- tryCatch(make_followup_generator(), error = function(e) NULL)

  # ---- System-prompt resolver ----
  # Single source for the system prompt used at every Agent() construction
  # site (.do_reset and apply_dataset). Reads the live style/custom reactives
  # when wired, falls back to the option-driven default otherwise. A build
  # failure (bad fragment, etc.) is logged and falls back so we never block
  # Agent construction on a prompt assembly issue.
  .resolve_system_prompt <- function() {
    if (is.null(style)) return(copilot_system_prompt())
    tryCatch(
      omicsagentovi::ovi_build_system_prompt(
        style  = shiny::isolate(style()),
        custom = if (!is.null(custom)) shiny::isolate(custom()) else NULL
      ),
      error = function(e) {
        log_info("copilot.run.prompt_build_failed", msg = conditionMessage(e))
        copilot_system_prompt()
      }
    )
  }

  # ---- Bindings factory (closes over session + host integration handles) ----
  .build_bindings <- function() {
    build_run_bindings(
      session          = session,
      evidence_api     = if (!is.null(evidence))
                           list(append_artifact = evidence$append_artifact)
                         else NULL,
      docs_dir         = docs_dir,
      data_dir         = pgx_dir,
      pgx_loaded_event = pgx_loaded_event
    )
  }

  .tools_are_enabled <- function() {
    if (is.null(tools_enabled)) return(TRUE)
    value <- tryCatch(shiny::isolate(tools_enabled()), error = function(e) TRUE)
    isTRUE(value)
  }

  .sync_runtime_controls <- function(agent_value) {
    if (is.null(agent_value)) return(agent_value)
    tryCatch({
      agent_value@runtime$tools_enabled <- .tools_are_enabled()
      agent_value
    }, error = function(e) {
      log_info("copilot.run.runtime_sync_failed", msg = conditionMessage(e))
      agent_value
    })
  }

  .selected_report_slots <- function() {
    if (is.null(report_context) ||
        !is.function(report_context[["selected_reports"]])) {
      return(character(0))
    }
    slots <- tryCatch(
      shiny::isolate(report_context$selected_reports()),
      error = function(e) character(0)
    )
    slots <- tryCatch(as.character(slots), error = function(e) character(0))
    slots[!is.na(slots) & nzchar(slots)]
  }

  .mark_reports_consumed <- function(slots) {
    if (is.null(report_context) ||
        !is.function(report_context[["mark_consumed"]])) {
      return(invisible(NULL))
    }
    tryCatch(
      report_context$mark_consumed(slots),
      error = function(e) {
        log_info("copilot.run.report_consume_failed", msg = conditionMessage(e))
      }
    )
    invisible(NULL)
  }

  .selected_doc_names <- function() {
    if (is.null(doc_context) ||
        !is.function(doc_context[["selected_docs"]])) {
      return(character(0))
    }
    docs <- tryCatch(
      shiny::isolate(doc_context$selected_docs()),
      error = function(e) character(0)
    )
    docs <- tryCatch(as.character(docs), error = function(e) character(0))
    docs[!is.na(docs) & nzchar(docs)]
  }

  .mark_docs_consumed <- function(docs) {
    if (is.null(doc_context) ||
        !is.function(doc_context[["mark_consumed"]])) {
      return(invisible(NULL))
    }
    tryCatch(
      doc_context$mark_consumed(docs),
      error = function(e) {
        log_info("copilot.run.doc_consume_failed", msg = conditionMessage(e))
      }
    )
    invisible(NULL)
  }

  .reset_context_panels <- function() {
    if (!is.null(report_context) &&
        is.function(report_context[["reset_consumed"]])) {
      tryCatch(report_context$reset_consumed(), error = function(e) {
        log_info("copilot.run.report_reset_failed", msg = conditionMessage(e))
      })
    }
    if (!is.null(doc_context) &&
        is.function(doc_context[["reset_consumed"]])) {
      tryCatch(doc_context$reset_consumed(), error = function(e) {
        log_info("copilot.run.doc_reset_failed", msg = conditionMessage(e))
      })
    }
    invisible(NULL)
  }

  # ---- Session id generator ----
  # AgentSession's default session_id is NA_character_, which makes the
  # sessions/turns/transcript_records save fail with NOT NULL constraint
  # errors. The package leaves session-id minting to the host, so we
  # generate one at every fresh-agent construction site (first dataset
  # load, new_chat reset, tier change).
  .new_session_id <- function() {
    paste0(
      "sess-",
      format(Sys.time(), "%Y%m%d-%H%M%S"),
      "-",
      paste(sample(c(0:9, letters), 8, replace = TRUE), collapse = "")
    )
  }

  # ---- Current pgx extractor ----
  # The global `pgx` is a reactiveValues whose `$name` field signals "a
  # On reset (new_chat / tier change) we rebuild an Agent and need the
  # currently-loaded PGX in the same shape the existing agent stored —
  # a plain list with class "pgx", not the reactiveValues handle (tools
  # run outside the reactive domain, and reactiveValues breaks
  # downstream consumers like .ovi_pgx_name and omicspgxmcp plot tools).
  # The current agent already holds the materialised PGX in its context,
  # set there by an earlier `apply_dataset` / `agent_set_pgx`, so prefer
  # that path. As a fallback (no current agent yet), materialise from the
  # reactiveValues directly.
  .current_pgx <- function() {
    current <- shiny::isolate(agent())
    if (!is.null(current)) {
      ctx_pgx <- tryCatch(current@context@pgx, error = function(e) NULL)
      if (!is.null(ctx_pgx)) {
        # Defensive normalisation: if an upstream code path (restore,
        # legacy plumbing) leaked a reactiveValues into the agent's
        # context, this is the last line of defence before .pgx_check
        # blows up on the next omicspgx call.
        return(copilot_normalize_pgx(ctx_pgx,
                                     source = "run_controller/.current_pgx/agent"))
      }
    }
    # Fallback: snapshot the global pgx reactiveValues through the
    # centralised normaliser.
    copilot_normalize_pgx(pgx, source = "run_controller/.current_pgx/global")
  }

  # Construct a fresh Agent at `the_tier` for `pgx_val`, splicing in the
  # user's provider/credentials/model via .copilot_agent_build_args. Returns
  # NULL (and notifies) on failure. Shared by the dispatch reset path and the
  # first-load apply_dataset path.
  .build_new_agent <- function(the_tier, pgx_val, bindings) {
    tryCatch(
      do.call(omicsagentovi::Agent, c(
        .copilot_agent_build_args(session, the_tier),
        list(
          system_prompt = .resolve_system_prompt(),
          context       = omicsagentovi::RunContext(pgx = pgx_val),
          session       = omicsagentovi::AgentSession(session_id = .new_session_id()),
          bindings      = bindings
        )
      )),
      error = function(e) {
        log_info("copilot.run.agent_construct_failed", msg = conditionMessage(e))
        shiny::showNotification(
          copilot_msg("agent_failed", msg = conditionMessage(e)),
          type = "error", session = session
        )
        NULL
      }
    )
  }

  # ------------------------------------------------------------------------
  # .run_ask
  # ------------------------------------------------------------------------
  .run_ask <- function(request) {
    current <- shiny::isolate(agent())
    if (is.null(current)) {
      chat_event(list(type = "post", role = "assistant",
                      text = copilot_msg("no_dataset")))
      return(invisible(NULL))
    }

    # Max-turns guard
    if (is.finite(maxturns)) {
      records <- tryCatch(
        omicsagentovi::session_transcript(current@session, view = "user"),
        error = function(e) list()
      )
      n_user <- sum(vapply(records,
                           function(r) identical(r@role, "user"),
                           logical(1)))
      if (n_user >= maxturns) {
        chat_event(list(type = "post", role = "assistant",
          text = copilot_msg("max_turns", n = maxturns)))
        return(invisible(NULL))
      }
    }

    current <- .sync_runtime_controls(current)
    selected_reports <- .selected_report_slots()
    if (length(selected_reports)) {
      staged <- .copilot_stage_ai_report_context(current, selected_reports)
      current <- .sync_runtime_controls(staged$agent)
      if (isTRUE(staged$staged) && length(staged$slots)) {
        .mark_reports_consumed(staged$slots)
      }
    }
    selected_docs <- .selected_doc_names()
    if (length(selected_docs)) {
      doc_staged <- .copilot_stage_user_docs_context(
        current, docs_dir, selected_docs
      )
      current <- .sync_runtime_controls(doc_staged$agent)
      if (isTRUE(doc_staged$staged) && length(doc_staged$docs)) {
        .mark_docs_consumed(doc_staged$docs)
      }
    }
    agent(current)

    # Optionally post user message into chat
    if (isTRUE(request$show_user_msg)) {
      chat_event(list(type = "post", role = "user", text = request$text))
    }

    gen <- tryCatch(
      omicsagentovi::agent_prompt_stream(
        current,
        request$text,
        on_tool_request = chat_on_tool_request,
        on_done = function(result) {
          if (!is.null(result$agent)) agent(result$agent)
          run_status(result$status %||% "completed")
          updated <- shiny::isolate(agent())
          if (!is.null(updated) &&
              omicsagentovi::session_is_dirty(updated@session)) {
            save_ctrl$on_run_settled()
          }

          # ---- Follow-up suggestion bubble ---------------------------------
          # Gates: completed status, non-empty text, generator available.
          # On any failure, silently skip — never block the chat.
          status_ok <- identical(result$status %||% "completed", "completed") &&
                       !is.null(result$text) && nzchar(result$text)
          if (!status_ok) {
            log_trace("copilot.followup.skipped", reason = "non_completed_or_empty")
            return(invisible())
          }
          if (is.null(followup_gen)) {
            log_trace("copilot.followup.skipped", reason = "no_generator")
            return(invisible())
          }
          # Snapshot what reports will be attached on the user's NEXT click.
          # `.selected_report_slots()` reads the report panel's
          # ticked-minus-consumed set; by on_done time the slots used in
          # THIS turn have already been moved into `consumed` (see
          # .mark_reports_consumed call earlier in .run_ask), so this returns
          # the post-turn attached set — exactly what the helper should
          # ground report-related suggestions on.
          attached_slots <- tryCatch(.selected_report_slots(),
                                     error = function(e) character(0))
          followup_payload <- tryCatch(
            build_followup_payload(
              agent                 = result$agent,
              last_text             = result$text,
              attached_report_slots = attached_slots
            ),
            error = function(e) {
              log_info("copilot.followup.failed",
                       phase = "payload_build",
                       msg = conditionMessage(e))
              list(last_text = result$text)
            }
          )
          p <- tryCatch(followup_gen$generate(followup_payload),
                        error = function(e) NULL)
          if (is.null(p)) {
            log_info("copilot.followup.failed", phase = "generate_call")
            return(invisible())
          }
          promises::then(p,
            onFulfilled = function(qs) {
              if (length(qs) == 0L) {
                log_trace("copilot.followup.skipped", reason = "empty_result")
                return(invisible())
              }
              chat_event(list(
                type = "post",
                role = "assistant",
                text = format_followup_bubble(qs)
              ))
            },
            onRejected = function(e) {
              log_info("copilot.followup.failed",
                       phase = "promise_reject",
                       msg = conditionMessage(e))
            }
          )
        }
      ),
      error = function(e) {
        log_info("copilot.run.stream_invoke_failed", msg = conditionMessage(e))
        run_status("failed")
        chat_event(list(type = "post", role = "assistant",
          text = copilot_friendly_error(e)))
        NULL
      }
    )
    if (is.null(gen)) return(invisible(NULL))

    run_status("streaming")
    chat_event(list(type = "stream", async_gen = gen))
    invisible(NULL)
  }

  # ------------------------------------------------------------------------
  # .run_abort
  # ------------------------------------------------------------------------
  .run_abort <- function(request) {
    current <- shiny::isolate(agent())
    if (is.null(current)) {
      log_info("copilot.run.abort_skipped", reason = "no_agent")
      return(invisible(NULL))
    }
    rs_status <- tryCatch(
      omicsagentovi::run_state_status(current@run_state),
      error = function(e) "unknown"
    )
    log_info("copilot.run.abort_invoked",
             reason     = request$reason %||% "User stopped the run",
             run_status = shiny::isolate(run_status()),
             rs_status  = rs_status)
    ok <- tryCatch(
      omicsagentovi::agent_request_abort(current, request$reason %||% "User stopped the run"),
      error = function(e) {
        log_info("copilot.run.abort_failed", msg = conditionMessage(e))
        FALSE
      }
    )
    log_info("copilot.run.abort_result", ok = isTRUE(ok))
    # run_status will be flipped to "aborted" by the stream's on_done callback.
    invisible(NULL)
  }

  # ------------------------------------------------------------------------
  # .do_reset — shared body for new_chat and tier change
  # ------------------------------------------------------------------------
  .do_reset <- function(new_tier = NULL) {
    current <- shiny::isolate(agent())
    # Pre-save synchronously if dirty (Q4 — full reset, mustn't drop work)
    if (!is.null(current) &&
        omicsagentovi::session_is_dirty(current@session)) {
      saved <- tryCatch(
        omicsagentovi::session_save(
          store, current,
          style  = if (!is.null(style))  shiny::isolate(style())  else NULL,
          custom = if (!is.null(custom)) shiny::isolate(custom()) else NULL
        ),
        error = function(e) {
          log_info("copilot.run.presave_failed", msg = conditionMessage(e))
          current
        }
      )
      agent(saved)
    }

    the_tier <- new_tier %||% shiny::isolate(tier())
    tier(the_tier)

    bindings <- .build_bindings()
    pgx_val  <- .current_pgx()

    # Carry the prior agent's dataset locator into the new agent so the new
    # session knows the dataset name + path on disk even when the LLM has not
    # called manage_pgx yet. Without this, session_save persists
    # dataset_name = NA and history shows "(no dataset)".
    prior_locator <- if (!is.null(current)) {
      tryCatch(current@context@dataset_locator,
               error = function(e) list(name = NA_character_,
                                        path = NA_character_,
                                        data_dir = NA_character_))
    } else {
      list(name = NA_character_, path = NA_character_, data_dir = NA_character_)
    }

    new_agent <- NULL
    if (!is.null(pgx_val)) {
      new_agent <- .build_new_agent(the_tier, pgx_val, bindings)
      # Install the dataset locator + current_dataset memory block so the
      # agent is aware of the dataset on the very first turn.
      if (!is.null(new_agent)) {
        new_agent <- tryCatch(
          omicsagentovi::agent_set_pgx(
            new_agent,
            pgx          = pgx_val,
            dataset_name = prior_locator$name,
            dataset_path = prior_locator$path,
            data_dir     = prior_locator$data_dir
          ),
          error = function(e) {
            log_info("copilot.run.set_pgx_after_reset_failed", msg = conditionMessage(e))
            new_agent
          }
        )
        # Stage context blocks (current_dataset, future providers) so the
        # LLM sees host-prepared context on the first user turn.
        new_agent <- .sync_runtime_controls(new_agent)
        new_agent <- .copilot_stage_context_blocks(new_agent)
      }
    }
    agent(new_agent)

    if (!is.null(evidence)) {
      tryCatch(evidence$clear(), error = function(e) NULL)
    }

    greeting <- if (is.null(new_agent)) {
      copilot_msg("greeting")
    } else {
      copilot_msg("greeting_active")
    }
    chat_event(list(type = "reset", role = "assistant", text = greeting))

    # New chat / tier change starts a fresh conversation — context tickboxes
    # should re-arm so previously-staged reports/docs can be sent again.
    .reset_context_panels()

    run_status("idle")
    invisible(NULL)
  }

  # ------------------------------------------------------------------------
  # .run_new_chat / .run_tier
  # ------------------------------------------------------------------------
  .run_new_chat <- function(request) .do_reset(new_tier = NULL)

  .run_tier <- function(request) {
    new_tier <- request$new_tier
    if (is.null(new_tier) || identical(new_tier, shiny::isolate(tier()))) {
      return(invisible(NULL))
    }
    .do_reset(new_tier = new_tier)
  }

  # ------------------------------------------------------------------------
  # dispatch
  # ------------------------------------------------------------------------
  dispatch <- function(request) {
    if (!is.list(request) || is.null(request$kind)) {
      log_info("copilot.run.bad_request")
      return(invisible(NULL))
    }

    # Guard: while streaming, only abort is permitted.
    if (identical(shiny::isolate(run_status()), "streaming") &&
        !identical(request$kind, "abort")) {
      log_info("copilot.run.dispatch_busy", kind = request$kind)
      return(invisible(NULL))
    }

    switch(request$kind,
      ask      = .run_ask(request),
      abort    = .run_abort(request),
      new_chat = .run_new_chat(request),
      tier     = .run_tier(request),
      {
        log_info("copilot.run.unknown_kind", kind = request$kind)
        invisible(NULL)
      }
    )
  }

  # ------------------------------------------------------------------------
  # apply_dataset — shared dataset-change controller
  # ------------------------------------------------------------------------
  apply_dataset <- function(pgx_val, name, path, data_dir) {
    # Funnel every caller through the normaliser so a stray reactiveValues
    # can never reach RunContext / agent_set_pgx.
    pgx_val <- copilot_normalize_pgx(pgx_val,
                                     source = "run_controller/apply_dataset")
    if (is.null(pgx_val)) return(invisible(NULL))

    current <- shiny::isolate(agent())
    if (is.null(current)) {
      # First dataset load: construct a fresh Agent at the current tier,
      # then install the dataset locator + current_dataset memory block so
      # the agent is aware of the dataset on the very first turn and the
      # next save records the real dataset name (not "(no dataset)").
      bindings <- .build_bindings()
      the_tier <- shiny::isolate(tier())
      new_agent <- .build_new_agent(the_tier, pgx_val, bindings)
      if (!is.null(new_agent)) {
        new_agent <- tryCatch(
          omicsagentovi::agent_set_pgx(
            new_agent,
            pgx          = pgx_val,
            dataset_name = name,
            dataset_path = path,
            data_dir     = data_dir
          ),
          error = function(e) {
            log_info("copilot.run.set_pgx_at_first_load_failed",
                     msg = conditionMessage(e))
            new_agent
          }
        )
        new_agent <- .sync_runtime_controls(new_agent)
        new_agent <- .copilot_stage_context_blocks(new_agent)
      }
      agent(new_agent)
    } else {
      # Switch dataset on existing Agent. Chat preserved (Q3).
      updated <- tryCatch(
        omicsagentovi::agent_set_pgx(
          current,
          pgx          = pgx_val,
          dataset_name = name,
          dataset_path = path,
          data_dir     = data_dir
        ),
        error = function(e) {
          log_info("copilot.run.set_pgx_failed", msg = conditionMessage(e))
          shiny::showNotification(
            copilot_msg("switch_failed", msg = conditionMessage(e)),
            type = "error", session = session
          )
          NULL
        }
      )
      if (!is.null(updated)) {
        # Re-stage context blocks so the next turn sees the new dataset.
        updated <- .sync_runtime_controls(updated)
        updated <- .copilot_stage_context_blocks(updated)
        agent(updated)
      }
    }
    invisible(NULL)
  }

  # ---- Public surface ----
  list(
    dispatch        = dispatch,
    apply_dataset   = apply_dataset,
    # Settle the run lifecycle after a mid-stream promise rejection.
    # Called by CopilotChatServer's onRejected handler (wired via
    # CopilotBoardServer) so run_status leaves "streaming" and dispatch
    # unblocks for the next user turn.
    on_stream_error = function() run_status("failed")
  )
}
