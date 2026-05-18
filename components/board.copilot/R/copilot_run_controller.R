# copilot_run_controller.R — Run Lifecycle Controller
#
# Unified dispatch for ask / new_chat / tier / abort. Owns `run_status`
# transitions. Calls `agent_prompt_stream` with `on_done` and
# `on_tool_request` callbacks. Constructs fresh Agents on new_chat / tier
# reset and on first dataset load.
#
# See .active_plans/refactor_copilot/run_controller/specs.md for the full
# contract.

# ---- Null-coalescing operator ----
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Local minimal logger (TODO phase 6: replace with copilot_logger.R log_info) ----
.log_info <- function(event, ...) {
  args <- list(...)
  kv <- if (length(args)) paste(names(args), unlist(lapply(args, as.character)),
                                sep = "=", collapse = " ") else ""
  message("[CopilotRun] ", event, if (nzchar(kv)) paste0(" ", kv) else "")
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
  session
) {

  # ---- Bindings factory (closes over session + host integration handles) ----
  .build_bindings <- function() {
    build_run_bindings(
      session          = session,
      evidence_api     = if (!is.null(evidence)) evidence$api else NULL,
      docs_dir         = docs_dir,
      data_dir         = pgx_dir,
      pgx_loaded_event = pgx_loaded_event
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
      if (!is.null(ctx_pgx)) return(ctx_pgx)
    }
    has_data <- tryCatch(
      !is.null(shiny::isolate(pgx$name)) && !is.null(shiny::isolate(pgx$X)),
      error = function(e) FALSE
    )
    if (!isTRUE(has_data)) return(NULL)
    snapshot <- shiny::reactiveValuesToList(pgx)
    class(snapshot) <- unique(c("pgx", class(snapshot)))
    snapshot
  }

  # ------------------------------------------------------------------------
  # .run_ask
  # ------------------------------------------------------------------------
  .run_ask <- function(request) {
    current <- shiny::isolate(agent())
    if (is.null(current)) {
      chat_event(list(type = "post", role = "assistant",
                      text = "Load a dataset first, then I can answer questions."))
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
          text = paste0("This session has reached the maximum of ", maxturns,
                        " user turns. Start a new chat to continue.")))
        return(invisible(NULL))
      }
    }

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
        }
      ),
      error = function(e) {
        .log_info("copilot.run.stream_invoke_failed", msg = conditionMessage(e))
        run_status("failed")
        chat_event(list(type = "post", role = "assistant",
          text = paste("Copilot error:", conditionMessage(e))))
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
    if (is.null(current)) return(invisible(NULL))
    tryCatch(
      omicsagentovi::agent_request_abort(current, request$reason %||% "User stopped the run"),
      error = function(e) {
        .log_info("copilot.run.abort_failed", msg = conditionMessage(e))
      }
    )
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
        omicsagentovi::session_save(store, current),
        error = function(e) {
          .log_info("copilot.run.presave_failed", msg = conditionMessage(e))
          current
        }
      )
      agent(saved)
    }

    the_tier <- new_tier %||% shiny::isolate(tier())
    tier(the_tier)

    bindings <- .build_bindings()
    pgx_val  <- .current_pgx()

    new_agent <- NULL
    if (!is.null(pgx_val)) {
      new_agent <- tryCatch(
        omicsagentovi::Agent(
          tier     = the_tier,
          context  = omicsagentovi::RunContext(pgx = pgx_val),
          bindings = bindings
        ),
        error = function(e) {
          .log_info("copilot.run.agent_construct_failed", msg = conditionMessage(e))
          shiny::showNotification(
            paste("Copilot: failed to create agent —", conditionMessage(e)),
            type = "error", session = session
          )
          NULL
        }
      )
    }
    agent(new_agent)

    if (!is.null(evidence)) {
      tryCatch(evidence$clear(), error = function(e) NULL)
    }

    chat_event(list(type = "clear"))
    greeting <- if (is.null(new_agent)) {
      # TODO(phase 6): copilot_msg("greeting.no_dataset")
      "Hi — load a dataset and ask me anything about your experiment."
    } else {
      # TODO(phase 6): copilot_msg("greeting")
      "Hi — what would you like to explore?"
    }
    chat_event(list(type = "post", role = "assistant", text = greeting))

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
      .log_info("copilot.run.bad_request")
      return(invisible(NULL))
    }

    # Guard: while streaming, only abort is permitted.
    if (identical(shiny::isolate(run_status()), "streaming") &&
        !identical(request$kind, "abort")) {
      .log_info("copilot.run.dispatch_busy", kind = request$kind)
      return(invisible(NULL))
    }

    switch(request$kind,
      ask      = .run_ask(request),
      abort    = .run_abort(request),
      new_chat = .run_new_chat(request),
      tier     = .run_tier(request),
      {
        .log_info("copilot.run.unknown_kind", kind = request$kind)
        invisible(NULL)
      }
    )
  }

  # ------------------------------------------------------------------------
  # apply_dataset — shared dataset-change controller
  # ------------------------------------------------------------------------
  apply_dataset <- function(pgx_val, name, path, data_dir) {
    if (is.null(pgx_val)) return(invisible(NULL))

    current <- shiny::isolate(agent())
    if (is.null(current)) {
      # First dataset load: construct a fresh Agent at the current tier.
      bindings <- .build_bindings()
      new_agent <- tryCatch(
        omicsagentovi::Agent(
          tier     = shiny::isolate(tier()),
          context  = omicsagentovi::RunContext(pgx = pgx_val),
          bindings = bindings
        ),
        error = function(e) {
          .log_info("copilot.run.agent_construct_failed", msg = conditionMessage(e))
          shiny::showNotification(
            paste("Copilot: failed to create agent —", conditionMessage(e)),
            type = "error", session = session
          )
          NULL
        }
      )
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
          .log_info("copilot.run.set_pgx_failed", msg = conditionMessage(e))
          shiny::showNotification(
            paste("Copilot: failed to switch dataset —", conditionMessage(e)),
            type = "error", session = session
          )
          NULL
        }
      )
      if (!is.null(updated)) agent(updated)
    }
    invisible(NULL)
  }

  # ---- Public surface ----
  list(
    dispatch      = dispatch,
    apply_dataset = apply_dataset
  )
}
