# copilot_restore_controller.R — Restore Lifecycle Controller
#
# Owns the async restore lifecycle: session loading, PGX injection, and
# transcript replay. Called once at board init via CopilotBoardServer.
#
# Returns list(start = function(session_id), status = reactive).
#
# Two deviations from the original spec:
#   1. `busy` arg replaced by `run_status` — `run_status()` is the truth;
#      guard checks `identical(run_status(), "streaming")`.
#   2. `bindings` arg replaced by `bindings_factory` — called inside start()
#      so bindings closures capture the current reactive environment.
#
# Phase 4 refactor: shinychat I/O is no longer called directly. The
# controller now writes events into the chat module's `chat_event`
# reactiveVal; CopilotChatServer dispatches them to shinychat.
#
# See .active_plans/refactor_copilot/restore_controller/specs.md for the
# full contract.

# --------------------------------------------------------------------------
# Internal helpers (file-local, prefixed `.`)
# --------------------------------------------------------------------------

#' Replay the visible transcript by pushing a single `replay` event.
#'
#' Fetches `session_transcript(agent@session, view = "user")`, filters out
#' records with empty `content_text`, and pushes one `replay` event into
#' `chat_event`. The chat module's replay handler does the per-record
#' append.
#'
#' @param agent      An `omicsagentovi::Agent`.
#' @param chat_event reactiveVal — chat module event bus writer.
#' @return Integer count of records included in the replay payload.
.replay_transcript <- function(agent, chat_event) {
  records <- omicsagentovi::session_transcript(agent@session, view = "user")
  # Drop empty-content records (keep filter semantics from Phase 3).
  records <- Filter(function(r) nzchar(r@content_text), records)
  # Always emit the event (even with zero records) and instruct the chat
  # module to clear first — this drops the "Restoring…" placeholder line
  # before the historical messages land. clear_first is consumed by
  # CopilotChatServer's replay handler in the same observer firing, so
  # there is no chance of the reactiveVal write coalescing away the clear.
  chat_event(list(type = "replay", clear_first = TRUE, records = records))
  length(records)
}

# --------------------------------------------------------------------------
# Factory
# --------------------------------------------------------------------------

#' Copilot Restore Controller
#'
#' Owns the async restore lifecycle. Constructs an `ExtendedTask` wrapping
#' `omicsagentovi::ovi_restore()`, wires its result observer, and handles
#' PGX injection and transcript replay on the Shiny thread.
#'
#' @param store             A `SessionStore` — constructed by the orchestrator.
#' @param agent             `reactiveVal(Agent|NULL)` — written only on
#'   successful restore; never touched on failure (previous agent preserved).
#' @param run_status        `reactiveVal(character)` — read-only; checked on
#'   `start()`. Guards against restoring while the agent is streaming.
#' @param restore_inflight  `reactiveVal(NULL|character)` — owned by this
#'   controller. `NULL` = idle; `"<session_id>"` = restoring.
#' @param bindings_factory  `function() -> RunBindings`. Called inside
#'   `start()` to build a fresh `RunBindings` per restore attempt.
#' @param local_pgx         `reactive()` returning current PGX or `NULL`.
#'   Evaluated on the Shiny thread after the task resolves.
#' @param data_dir          Character scalar — passed to `agent_set_pgx`.
#' @param evidence          List with `$clear()` or `NULL`. Guard for NULL.
#'   Phase 5 wires a real evidence handle here.
#' @param chat_event        `reactiveVal` — chat module event bus writer.
#'   Receives `clear`, `post`, and `replay` events.
#' @param session           Shiny session object — required for
#'   `ExtendedTask` scoping and `showNotification`.
#'
#' @return `list(start = function(session_id), status = reactive)`.
#'   `status()` values: `"idle"` / `"restoring"` / `"attaching_dataset"` /
#'   `"replaying"` / `"failed"`.
#' @export
copilot_restore_controller <- function(
  store,
  agent,
  run_status,
  restore_inflight,
  bindings_factory,
  local_pgx,
  data_dir,
  evidence,
  chat_event,
  session
) {

  # ---- Internal state ----
  restore_status <- shiny::reactiveVal("idle")

  # ---- ExtendedTask ----
  # `provider` + `credentials` are captured on the Shiny thread by start()
  # and passed in as task args, never read from the reactive inside the
  # worker (the future-worker gotcha). `credentials` is a plain value-capturing
  # closure (or NULL), so it serialises into the future cleanly.
  restore_task <- shiny::ExtendedTask$new(
    function(store, session_id, bindings, provider, credentials) {
      promises::future_promise({
        omicsagentovi::ovi_restore(
          session_id  = session_id,
          session_dir = store@session_dir,
          bindings    = bindings,
          restore_pgx = "never",
          db_name     = basename(store@db_path),
          provider    = provider,
          credentials = credentials
        )
      }, seed = TRUE)
    })

  # ---- Failure handler ----
  .handle_restore_failure <- function(err, session_id, previous_agent) {
    msg <- if (!is.null(err)) conditionMessage(err) else "Restore failed (unknown error)."

    # Classify by message pattern (package uses plain stop(), no typed conditions)
    is_fatal <- !is.null(err) &&
      grepl("model.*not resolvable", msg, ignore.case = TRUE)

    if (is_fatal) {
      shiny::showNotification(msg, type = "error", session = session)
    } else {
      shiny::showNotification(msg, type = "warning", session = session)
    }

    chat_event(list(
      type = "post", role = "assistant",
      text = copilot_msg("restore_failed")
    ))

    restore_inflight(NULL)
    restore_status("failed")

    # Do NOT touch agent() — previous_agent stays active by virtue of not
    # being overwritten. (previous_agent arg retained for future use / tests.)
    invisible(NULL)
  }

  # ---- Shared completion handler (sync and async paths) ----
  # Performs PGX inject, context-block re-stage, agent commit, and replay.
  # `restore_mode` is logged for diagnostic separation of fast path vs future.
  .complete_restore <- function(restored, session_id, restore_mode = "async") {
    if (is.null(restored)) return(invisible(NULL))

    raw_pgx <- tryCatch(shiny::isolate(local_pgx()), error = function(e) NULL)
    # Read $name BEFORE normalisation: on a reactiveValues the field
    # lookup still works, and after normalisation we still want the same
    # name; doing it here keeps NA-handling symmetric with the old behaviour.
    pgx_name <- tryCatch(
      shiny::isolate(raw_pgx$name),
      error = function(e) NA_character_
    )

    # Funnel the global reactiveValues through the centralised normaliser
    # so agent@context@pgx is ALWAYS a classed list (or NULL). Bypassing
    # this (as the pre-fix code did) leaks a `reactivevalues` object into
    # the agent context and trips omicspgx::.pgx_check on the next call.
    current_pgx <- copilot_normalize_pgx(raw_pgx,
                                         source = "restore_controller/.complete_restore")

    if (!is.null(current_pgx)) {
      restore_status("attaching_dataset")
      started <- Sys.time()
      restored <- tryCatch(
        omicsagentovi::agent_set_pgx(
          restored,
          current_pgx,
          dataset_name = pgx_name %||% NA_character_,
          data_dir     = data_dir
        ),
        error = function(e) {
          log_info("copilot.restore_pgx_inject_failed",
            msg = conditionMessage(e))
          restored
        }
      )
      elapsed_ms <- round(
        as.numeric(difftime(Sys.time(), started, units = "secs")) * 1000, 1)
      log_info("copilot.restore_pgx_injected",
        dataset    = pgx_name %||% "unknown",
        elapsed_ms = elapsed_ms)
    }

    restored <- .copilot_stage_context_blocks(restored)

    agent(restored)
    restore_status("replaying")

    n_replayed <- .replay_transcript(restored, chat_event)

    restore_inflight(NULL)
    restore_status("idle")

    log_info("copilot.restore_complete",
      session_id   = session_id %||% "unknown",
      n_replayed   = n_replayed,
      has_pgx      = !is.null(current_pgx),
      restore_mode = restore_mode)

    invisible(NULL)
  }

  # ---- Result observer (async path) ----
  shiny::observeEvent(restore_task$result(), {
    session_id <- shiny::isolate(restore_inflight())
    restored <- tryCatch(
      restore_task$result(),
      error = function(e) {
        .handle_restore_failure(e, session_id, NULL)
        NULL
      }
    )
    if (is.null(restored)) return()
    .complete_restore(restored, session_id, restore_mode = "async")
  }, ignoreNULL = TRUE)

  # ---- start() — sole public entry point ----
  start <- function(session_id) {
    # Guard: refuse while agent is actively streaming
    if (identical(run_status(), "streaming")) {
      shiny::showNotification(
        copilot_msg("restore_busy"),
        type = "warning",
        session = session
      )
      return(invisible(NULL))
    }

    # Guard: another restore already in flight
    if (!is.null(restore_inflight())) {
      return(invisible(NULL))
    }

    previous_agent <- shiny::isolate(agent())

    # Capture the user's AI provider + credential closure on the Shiny thread
    # (never inside the future). BigOmics -> provider "bigomics" + NULL
    # credentials, which ovi_restore treats exactly as today (env-var path,
    # openai tier fallback). A BYOK provider re-resolves the saved tier for
    # that provider and applies the user key.
    ai_sel   <- .copilot_ai_provider(session)
    provider <- ai_sel$provider
    cred     <- ai_sel$credentials

    restore_inflight(session_id)
    restore_status("restoring")

    # Clear evidence panel if wired (Phase 5)
    if (!is.null(evidence)) evidence$clear()

    # Clear chat AND post placeholder atomically. Two back-to-back writes
    # to a reactiveVal in the same flush cycle coalesce — only the second
    # value reaches the observer, which dropped the "clear" and left the
    # previous chat in place. The `reset` event type bundles clear+post
    # in a single observer firing (see CopilotChatServer chat_event
    # handler). Restoring without this is exactly the "appends to current
    # conversation" bug.
    chat_event(list(type = "reset", role = "assistant",
                    text = copilot_msg("restore_started")))

    has_local_pgx <- !is.null(tryCatch(local_pgx(), error = function(e) NULL))
    log_info("copilot.restore_start",
      session_id    = session_id,
      restore_mode  = if (has_local_pgx) "sync" else "async")

    bindings <- tryCatch(bindings_factory(), error = function(e) {
      .handle_restore_failure(e, session_id, previous_agent)
      NULL
    })
    if (is.null(bindings)) return(invisible(NULL))

    if (has_local_pgx) {
      # Sync fast path: a PGX is already attached so there is no future-side
      # heavy work to wait on. Restore inside onFlushed so the "Restoring..."
      # placeholder is painted first, then run synchronously on the main
      # thread — avoids the future/promise queue + observer-fire roundtrip.
      session$onFlushed(function() {
        started <- Sys.time()
        restored <- tryCatch(
          omicsagentovi::ovi_restore(
            session_id  = session_id,
            session_dir = store@session_dir,
            bindings    = bindings,
            restore_pgx = "never",
            db_name     = basename(store@db_path),
            provider    = provider,
            credentials = cred
          ),
          error = function(e) {
            log_info("copilot.restore_sync_failed", msg = conditionMessage(e))
            .handle_restore_failure(e, session_id, previous_agent)
            NULL
          }
        )
        if (is.null(restored)) return()
        restore_ms <- round(
          as.numeric(difftime(Sys.time(), started, units = "secs")) * 1000, 1)
        log_info("copilot.restore_worker_done",
          session_id = session_id, restore_ms = restore_ms,
          restore_mode = "sync")
        .complete_restore(restored, session_id, restore_mode = "sync")
      }, once = TRUE)
      return(invisible(NULL))
    }

    tryCatch(
      {
        restore_task$invoke(store, session_id, bindings, provider, cred)
      },
      error = function(e) {
        log_info("copilot.restore_invoke_failed",
          msg = conditionMessage(e))
        .handle_restore_failure(e, session_id, previous_agent)
      }
    )

    invisible(NULL)
  }

  # ---- Public surface ----
  list(
    start  = start,
    status = shiny::reactive(restore_status())
  )
}
