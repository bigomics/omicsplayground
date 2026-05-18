# copilot_bindings.R — RunBindings factory for the Copilot board
#
# This is a pure adapter: no reactives, no shiny::observe, no UI.
# Called once per Agent construction (and once per tier change when a fresh
# Agent is needed). Returns a typed RunBindings S7 object the package can consume.
#
# All callback signatures verified against actual call sites in omicsagentovi:
#   plot_callback:     tool-show-plot.R:79     — function(pgx, plot_type, args, artifact)
#   notification_sink: tool-manage-pgx.R:202   — function(level, payload)
#   progress_callback: NEVER CALLED            — stub only

# ---- Tier IDs (from contract_audit.md §1 / copilot-models.R:12-25) ----
.COPILOT_TIER_IDS <- c("copilot-default", "copilot-fast", "copilot-deep")

#' Build RunBindings for the Copilot board
#'
#' Creates a typed `omicsagentovi::RunBindings` S7 object that wires host
#' integration callbacks and filesystem locations into the agent.
#'
#' @param session Shiny session object, or NULL (used for reactive-domain
#'   wrapping of plot_callback). NULL -> callbacks that need a session are NULL.
#' @param evidence_api List with at minimum `$append_artifact(pgx, plot_type, args, artifact)`.
#'   Returned by `CopilotEvidenceServer`. NULL -> plot_callback set to NULL.
#' @param docs_dir character(0 or 1). Passed through to RunBindings; NULL/"" -> character(0).
#' @param data_dir character(0 or 1). Same normalisation as docs_dir.
#' @param pgx_loaded_event reactiveVal handle (or bare function, or NULL).
#'   The notification_sink writes the "pgx_loaded" payload here when the
#'   agent's manage_pgx tool loads a dataset. A Shiny-thread observer then
#'   calls agent_set_pgx(). NULL -> no-op for pgx_loaded events.
#'
#' @return An `omicsagentovi::RunBindings` S7 object.
#'
#' @export
build_run_bindings <- function(
  session,
  evidence_api,
  docs_dir         = character(0),
  data_dir         = character(0),
  pgx_loaded_event = NULL
) {
  # ---- Factory-time validation ----
  if (!is.null(evidence_api) && !is.function(evidence_api[["append_artifact"]])) {
    stop("evidence_api must have $append_artifact as a function", call. = FALSE)
  }
  if (!is.null(pgx_loaded_event) &&
      !is.function(pgx_loaded_event) &&
      !inherits(pgx_loaded_event, "reactiveVal")) {
    stop(
      "pgx_loaded_event must be a reactiveVal, a function, or NULL",
      call. = FALSE
    )
  }

  # ---- Normalise directory slots ----
  docs_dir <- .normalise_dir(docs_dir)
  data_dir <- .normalise_dir(data_dir)

  # ---- plot_callback ----
  # Phase 1: forward raw arguments to evidence_api$append_artifact.
  # TODO(phase 5): render plot from recipe here using copilot_build_plot(),
  #   copilot_detect_plot_kind(), copilot_parse_features() before forwarding.
  plot_cb <- if (is.null(evidence_api)) {
    NULL
  } else {
    impl <- function(pgx, plot_type, args, artifact) {
      evidence_api$append_artifact(
        pgx       = pgx,
        plot_type = plot_type,
        args      = args,
        artifact  = artifact
      )
    }
    # Wrap for reactive-domain safety (callback fires from tool-execution thread)
    if (!is.null(session) && is.function(session[["wrapFunction"]])) {
      session$wrapFunction(impl)
    } else if (!is.null(session)) {
      function(...) shiny::withReactiveDomain(session, impl(...))
    } else {
      impl
    }
  }

  # ---- notification_sink ----
  # Trampoline pattern: write payload to pgx_loaded_event reactiveVal so the
  # Shiny thread observer can call agent_set_pgx(). Do NOT call agent_set_pgx
  # here — this closure executes on the tool-execution thread.
  #
  # Actual call site (tool-manage-pgx.R:202):
  #   notification_sink("pgx_loaded", list(pgx, dataset_name, name_arg, data_dir))
  notification_cb <- if (is.null(session)) {
    NULL
  } else {
    function(level, payload) {
      if (identical(level, "pgx_loaded")) {
        if (!is.null(pgx_loaded_event)) {
          pgx_loaded_event(payload)  # write to reactiveVal — safe from any thread
        }
      }
      # All other levels: no-op (forward-compatible placeholder).
      # TODO(phase 6): route to copilot_logger when shared/ arrives.
      invisible(NULL)
    }
  }

  # ---- progress_callback ----
  # Slot is validated in RunBindings but never called by the package (commit
  # 51ad59d). Keep as a no-op stub rather than NULL so the slot is always a
  # function and the contract can be wired later without touching RunBindings.
  progress_cb <- function(...) invisible(NULL)

  # ---- Construct and return ----
  omicsagentovi::RunBindings(
    plot_callback     = plot_cb,
    progress_callback = progress_cb,
    notification_sink = notification_cb,
    docs_dir          = docs_dir,
    data_dir          = data_dir
  )
}

# ---- Internal helpers ----

# Normalise a directory slot: NULL, NA, or "" -> character(0); else pass through.
.normalise_dir <- function(path) {
  if (is.null(path) || (length(path) == 1L && (is.na(path) || !nzchar(path)))) {
    return(character(0))
  }
  path
}
