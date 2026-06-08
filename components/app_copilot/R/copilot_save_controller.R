# copilot_save_controller.R — Persistence Lifecycle Controller
#
# Owns the session_save path. Called once at board init. Exposes:
#   on_run_settled()  — call after each run settles; saves if dirty
#   status            — reactive "idle" | "saving" | "failed"
#
# Does NOT own: agent_rv, busy state, session dirty/saved generations.
# Those belong to the orchestrator and the package respectively.
#
# Concurrency model
# -----------------
# The save runs synchronously on the main thread, but is deferred to the next
# `session$onFlushed(once = TRUE)` so the user has already seen the agent's
# response by the time the SQLite write begins.
#
# The earlier design wrapped `session_save()` in `promises::future_promise()`
# under an `ExtendedTask`. That crashed with `bad_weak_ptr` after the first
# turn because two of the captured objects hold C++ resources that do not
# survive `future` serialization into a worker process:
#   1. SessionStore's DBI/RSQLite connection — a raw SQLite pointer.
#   2. The Agent's @bindings closures, which capture the live Shiny session
#      handle (httpuv/later weak refs) via session$wrapFunction().
# Either alone is enough to break a future-based save.
#
# session_save() is fast and IO-bound (a single BEGIN IMMEDIATE write), so
# running it on the main thread after the flush cycle is the right tradeoff.
# See .active_plans/refactor_copilot/save_controller/specs.md for the original
# contract; this implementation preserves the public surface and concurrency
# semantics (skip-if-running, session-id mismatch guard) without `future`.

# --------------------------------------------------------------------------
# Internal prune helper
# --------------------------------------------------------------------------

# Delete oldest sessions when the store exceeds max_history entries.
# ovi_sessions() returns rows ordered by updated_at DESC so the tail is oldest.
.prune_sessions <- function(store, max_history) {
  rows <- omicsagentovi::ovi_sessions(session_dir = store@session_dir)
  if (nrow(rows) <= max_history) return(invisible(NULL))
  n_delete <- nrow(rows) - max_history
  to_delete <- tail(rows$session_id, n_delete)
  for (sid in to_delete) {
    tryCatch(
      omicsagentovi::ovi_session_delete(sid, session_dir = store@session_dir),
      error = function(e) {
        log_info("copilot.prune_failed", sid = sid, msg = conditionMessage(e))
      }
    )
  }
  invisible(NULL)
}

# --------------------------------------------------------------------------
# Factory
# --------------------------------------------------------------------------

#' Copilot Save Controller
#'
#' Owns the persistence lifecycle. Defers `omicsagentovi::session_save()` to
#' the next flushed tick so it runs after the current reactive flush cycle —
#' the user sees the streamed response before the SQLite write begins.
#'
#' Tier-change pre-save is synchronous and is done inline in the run
#' controller; this factory does NOT provide that path.
#'
#' @param store             A `SessionStore` — constructed by the orchestrator.
#' @param agent             `reactiveVal(Agent|NULL)` — read + write.
#' @param history_invalidation_tick `reactiveVal(integer)` — incremented after
#'   each successful save.
#' @param session           Shiny session object — used for `onFlushed`.
#' @param max_history       Maximum saved sessions to keep after pruning.
#'   Defaults to `copilot_max_history()` (reads `copilot.max_history` option,
#'   falls back to `.COPILOT_MAX_HISTORY`).
#'
#' @return `list(on_run_settled = function(), status = reactive)`
#' @export
copilot_save_controller <- function(
  store,
  agent,
  history_invalidation_tick,
  session,
  max_history = copilot_max_history(),
  style       = NULL,
  custom      = NULL
) {

  # ---- Internal state ----
  save_status <- shiny::reactiveVal("idle")  # "idle" | "saving" | "failed"

  # ---- Save request ----
  .request_save <- function(reason) {
    current <- shiny::isolate(agent())
    if (is.null(current)) return(invisible(NULL))
    if (!omicsagentovi::session_is_dirty(current@session)) return(invisible(NULL))

    if (identical(shiny::isolate(save_status()), "saving")) {
      # A save is already queued (or just running). Skip — the next settled
      # run will find the same dirty state and save again.
      log_info("copilot.save_skipped_busy",
        reason = reason,
        dirty_gen = current@session@dirty_generation)
      return(invisible(NULL))
    }

    save_status("saving")
    session$onFlushed(function() .do_save(reason), once = TRUE)
    invisible(NULL)
  }

  # ---- Actual save body — runs on the main thread after flush ----
  .do_save <- function(reason) {
    current <- shiny::isolate(agent())
    if (is.null(current)) {
      save_status("idle")
      return(invisible(NULL))
    }

    updated_agent <- tryCatch(
      omicsagentovi::session_save(
        store, current,
        style  = if (!is.null(style))  shiny::isolate(style())  else NULL,
        custom = if (!is.null(custom)) shiny::isolate(custom()) else NULL
      ),
      error = function(e) {
        # Surface to both the structured log (telemetry) and the R console
        # (server-side dev/ops visibility). The end user is intentionally
        # not notified — session_save failures are infra-level and the chat
        # itself remains usable.
        log_info("copilot.save_failed",
          reason = reason, msg = conditionMessage(e))
        message(sprintf("[copilot.save_failed] reason=%s msg=%s",
                        reason, conditionMessage(e)))
        save_status("failed")
        NULL
      }
    )
    if (is.null(updated_agent)) return(invisible(NULL))

    # Only write back if the active session hasn't been replaced while we
    # were preparing the save (e.g. tier change between scheduling and
    # onFlushed firing — narrow window but worth guarding).
    latest <- shiny::isolate(agent())
    if (!is.null(latest) &&
        identical(latest@session@session_id, updated_agent@session@session_id)) {
      agent(updated_agent)
    }

    .prune_sessions(store, max_history)
    history_invalidation_tick(shiny::isolate(history_invalidation_tick()) + 1L)
    save_status("idle")

    log_info("copilot.save_complete",
      session_id = updated_agent@session@session_id,
      saved_gen  = updated_agent@session@saved_generation)
    invisible(NULL)
  }

  # ---- Public surface ----
  list(
    on_run_settled = function() .request_save("run_settled"),
    status         = shiny::reactive(save_status())
  )
}
