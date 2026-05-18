# copilot_save_controller.R — Persistence Lifecycle Controller
#
# Owns the async session_save path. Called once at board init. Exposes:
#   on_run_settled()  — call after each run settles; fires async save if dirty
#   status            — reactive "idle" | "saving" | "failed"
#
# Does NOT own: agent_rv, busy state, session dirty/saved generations.
# Those belong to the orchestrator and the package respectively.
#
# See .active_plans/refactor_copilot/save_controller/specs.md for the full
# contract.

# --------------------------------------------------------------------------
# Internal prune helper
# --------------------------------------------------------------------------

# Delete oldest sessions when the store exceeds max_history entries.
# ovi_sessions() returns rows ordered by updated_at DESC so the tail is oldest.
.prune_sessions <- function(store, max_history) {
  rows <- omicsagentovi::ovi_sessions(session_dir = store@session_dir)
  if (nrow(rows) <= max_history) return(invisible(NULL))
  n_delete <- nrow(rows) - max_history
  # tail() of a DESC-ordered frame gives the oldest rows
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
#' Owns the async persistence lifecycle. Constructs an `ExtendedTask` wrapping
#' `omicsagentovi::session_save()` and wires its result observer. The controller
#' is the only board code that calls `session_save()` on the normal async path.
#'
#' Tier-change pre-save is synchronous and is done inline in the run controller;
#' this factory does NOT provide that path.
#'
#' @param store             A `SessionStore` — constructed by the orchestrator.
#' @param agent             `reactiveVal(Agent|NULL)` — read + write.
#' @param history_invalidation_tick `reactiveVal(integer)` — incremented after
#'   each successful save.
#' @param session           Shiny session object — required for `ExtendedTask`
#'   scoping.
#' @param max_history       Maximum saved sessions to keep after pruning.
#'   Defaults to `getOption("copilot.max_history", 20L)`.
#'
#' @return `list(on_run_settled = function(), status = reactive)`
#' @export
copilot_save_controller <- function(
  store,
  agent,
  history_invalidation_tick,
  session,
  max_history = getOption("copilot.max_history", 20L)
) {

  # ---- Internal state ----
  save_status <- shiny::reactiveVal("idle")  # "idle" | "saving" | "failed"

  # ---- ExtendedTask ----
  save_task <- shiny::ExtendedTask$new(function(store, agent) {
    promises::future_promise({
      omicsagentovi::session_save(store, agent)
    }, seed = TRUE)
  })

  # ---- Async save request helper ----
  .request_async_save <- function(reason) {
    current <- shiny::isolate(agent())
    if (is.null(current)) return(invisible(NULL))
    if (!omicsagentovi::session_is_dirty(current@session)) return(invisible(NULL))

    if (identical(save_task$status(), "running")) {
      # Skip — the next settled run will find the same dirty state and save again.
      log_info("copilot.save_skipped_busy",
        reason = reason,
        dirty_gen = current@session@dirty_generation)
      return(invisible(NULL))
    }

    save_status("saving")
    tryCatch(
      save_task$invoke(store, current),
      error = function(e) {
        save_status("failed")
        log_info("copilot.save_invoke_failed", msg = conditionMessage(e))
      }
    )
    invisible(NULL)
  }

  # ---- Result observer ----
  shiny::observeEvent(save_task$result(), {
    updated_agent <- tryCatch(save_task$result(), error = function(e) NULL)

    if (is.null(updated_agent)) {
      save_status("failed")
      log_info("copilot.save_failed")
      return()
    }

    # Only write back if the active session hasn't been replaced while we
    # were saving (e.g. tier change or new chat mid-flight).
    current <- shiny::isolate(agent())
    if (!is.null(current) &&
        identical(current@session@session_id, updated_agent@session@session_id)) {
      agent(updated_agent)
    }

    .prune_sessions(store, max_history)
    history_invalidation_tick(shiny::isolate(history_invalidation_tick()) + 1L)
    save_status("idle")

    log_info("copilot.save_complete",
      session_id = updated_agent@session@session_id,
      saved_gen  = updated_agent@session@saved_generation)

  }, ignoreNULL = TRUE)

  # ---- Public surface ----
  list(
    on_run_settled = function() .request_async_save("run_settled"),
    status         = shiny::reactive(save_status())
  )
}
