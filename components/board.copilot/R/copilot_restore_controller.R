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
# See .active_plans/refactor_copilot/restore_controller/specs.md for the
# full contract and failure mode table.

# ---- Null-coalescing operator ----
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Local minimal logger (TODO phase 6: replace with copilot_logger.R log_info) ----
.log_info_restore <- function(event, ...) {
  args <- list(...)
  kv <- if (length(args)) paste(names(args), unlist(lapply(args, as.character)),
                                sep = "=", collapse = " ") else ""
  message("[CopilotRestore] ", event, if (nzchar(kv)) paste0(" ", kv) else "")
}

# --------------------------------------------------------------------------
# Internal helpers (file-local, prefixed `.`)
# --------------------------------------------------------------------------

#' Replay the visible transcript into shinychat
#'
#' Iterates `session_transcript(agent@session, view = "user")` and calls
#' `shinychat::chat_append_message` for each record with non-empty
#' `content_text`. Returns count of records actually appended.
#'
#' @param agent   An `omicsagentovi::Agent`.
#' @param chat_ns Character scalar — namespaced shinychat id.
#' @return Integer count of records appended.
.replay_transcript <- function(agent, chat_ns) {
  records <- omicsagentovi::session_transcript(agent@session, view = "user")
  n <- 0L
  for (rec in records) {
    if (!nzchar(rec@content_text)) next
    tryCatch(
      {
        shinychat::chat_append_message(
          chat_ns,
          list(role = rec@role, content = rec@content_text),
          chunk = FALSE
        )
        n <- n + 1L
      },
      error = function(e) {
        .log_info_restore("copilot.replay_record_failed",
          idx = rec@idx,
          msg = conditionMessage(e))
      }
    )
  }
  n
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
#' @param chat_ns           Character — namespaced shinychat id
#'   (e.g. `session$ns("chat")`).
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
  chat_ns,
  session
) {

  # ---- Internal state ----
  restore_status <- shiny::reactiveVal("idle")

  # ---- ExtendedTask ----
  restore_task <- shiny::ExtendedTask$new(function(store, session_id, bindings) {
    promises::future_promise({
      omicsagentovi::ovi_restore(
        session_id  = session_id,
        session_dir = store@session_dir,
        bindings    = bindings,
        restore_pgx = "never"
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

    # Placeholder chat message — TODO phase 6: replace via copilot_msg()
    shinychat::chat_append(
      chat_ns,
      "Restore failed. You can continue with the current session or pick another from history."
    )

    restore_inflight(NULL)
    restore_status("failed")

    # Do NOT touch agent() — previous_agent stays active by virtue of not
    # being overwritten. (previous_agent arg retained for future use / tests.)
    invisible(NULL)
  }

  # ---- Result observer ----
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

    # Inject current PGX on the Shiny thread (no future-boundary serialisation)
    current_pgx <- tryCatch(shiny::isolate(local_pgx()), error = function(e) NULL)

    if (!is.null(current_pgx)) {
      restore_status("attaching_dataset")
      started <- Sys.time()
      restored <- tryCatch(
        omicsagentovi::agent_set_pgx(
          restored,
          current_pgx,
          dataset_name = current_pgx$name %||% NA_character_,
          data_dir     = data_dir
        ),
        error = function(e) {
          .log_info_restore("copilot.restore_pgx_inject_failed",
            msg = conditionMessage(e))
          restored  # continue without PGX rather than failing the whole restore
        }
      )
      elapsed_ms <- round(
        as.numeric(difftime(Sys.time(), started, units = "secs")) * 1000, 1)
      .log_info_restore("copilot.restore_pgx_injected",
        dataset    = current_pgx$name %||% "unknown",
        elapsed_ms = elapsed_ms)
    }

    agent(restored)
    restore_status("replaying")

    n_replayed <- .replay_transcript(restored, chat_ns)

    restore_inflight(NULL)
    restore_status("idle")

    .log_info_restore("copilot.restore_complete",
      session_id = session_id %||% "unknown",
      n_replayed = n_replayed,
      has_pgx    = !is.null(current_pgx))

  }, ignoreNULL = TRUE)

  # ---- start() — sole public entry point ----
  start <- function(session_id) {
    # Guard: refuse while agent is actively streaming
    if (identical(run_status(), "streaming")) {
      shiny::showNotification(
        "Restore is unavailable while the agent is responding.",
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

    restore_inflight(session_id)
    restore_status("restoring")

    # Clear evidence panel if wired (Phase 5)
    if (!is.null(evidence)) evidence$clear()

    # Clear chat and show placeholder — TODO phase 6: replace via copilot_msg()
    shinychat::chat_clear(chat_ns)
    shinychat::chat_append(chat_ns, "Restoring previous session…")

    .log_info_restore("copilot.restore_start", session_id = session_id)

    tryCatch(
      {
        # Build fresh bindings per restore attempt
        bindings <- bindings_factory()
        restore_task$invoke(store, session_id, bindings)
      },
      error = function(e) {
        .log_info_restore("copilot.restore_invoke_failed",
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
