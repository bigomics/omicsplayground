# copilot_bindings.R — RunBindings factory for the Copilot board
#
# This is a pure adapter: no reactives, no shiny::observe, no UI.
# Called once per Agent construction (and once per tier change when a fresh
# Agent is needed). Returns a typed RunBindings S7 object the package can consume.
#
# All callback signatures verified against actual call sites in omicsagentovi:
#   plot_callback:     tool-show-plot.R        — function(pgx, plot_type, args, artifact, plot_result = NULL)
#   notification_sink: tool-manage-pgx.R:202   — function(level, payload)
#   progress_callback: NEVER CALLED            — stub only

#' Build RunBindings for the Copilot board
#'
#' Creates a typed `omicsagentovi::RunBindings` S7 object that wires host
#' integration callbacks and filesystem locations into the agent.
#'
#' @param session Shiny session object, or NULL (used for reactive-domain
#'   wrapping of plot_callback). NULL -> callbacks that need a session are NULL.
#' @param evidence_api List with at minimum `$append_artifact(record)` (single-arg,
#'   the fully-rendered record list). Returned by `CopilotEvidenceServer`. NULL ->
#'   plot_callback set to NULL.
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
  # Receives (pgx, plot_type, args, artifact, plot_result) from the package.
  # Builds the plot from the recipe using the board-side renderer, then hands a
  # fully-rendered record to evidence_api$append_artifact(record) (single-arg).
  # Errors are UI failures, not run failures (overview §5): caught locally and
  # surfaced via showNotification without propagating to the package.
  plot_cb <- if (is.null(evidence_api)) {
    NULL
  } else {
    impl <- function(pgx, plot_type, args, artifact, plot_result = NULL) {
      tryCatch({
        has_plot_result <- is.list(plot_result) && !is.null(plot_result$plot)
        if (has_plot_result) {
          plot_obj <- plot_result$plot
          if (!is.null(plot_result$plot_type)) {
            plot_type <- plot_result$plot_type
          }
        } else {
          plot_obj <- copilot_build_plot(pgx, plot_type, args)
        }
        kind     <- copilot_detect_plot_kind(plot_obj)
        if (is.null(kind)) {
          stop("unrecognised plot class: ", paste(class(plot_obj), collapse = " "))
        }

        prerendered_path <- NULL
        if (identical(kind, "ggplot")) {
          prerendered_path <- tryCatch(
            copilot_prerender_ggplot(plot_obj),
            error = function(e) {
              message("[CopilotBindings] prerender_ggplot failed: ", conditionMessage(e))
              NULL
            }
          )
        } else if (identical(kind, "plotly")) {
          plot_obj <- copilot_prerender_plotly(plot_obj)
        }

        label <- trimws(paste(
          plot_type,
          args$contrast   %||%
          args$collection %||%
          "",
          sep = " "
        ))

        record <- list(
          kind             = kind,
          plot             = plot_obj,
          prerendered_path = prerendered_path,
          plot_type        = plot_type,
          args             = args,
          artifact         = artifact,
          label            = label,
          timestamp        = Sys.time()
        )
        evidence_api$append_artifact(record)
      }, error = function(e) {
        # Artifact failures are UI failures, not run failures (overview §5).
        # Surface as a warning notification; do not propagate so the stream
        # on_done handler can still complete normally.
        # Guard: showNotification needs a live Shiny session — fall back to
        # message() in headless/test contexts where session is NULL or lacks
        # sendNotification.
        has_session <- !is.null(session) &&
          is.function(tryCatch(session$sendNotification, error = function(e2) NULL))
        if (has_session) {
          shiny::showNotification(
            copilot_msg("plot_failed", msg = conditionMessage(e)),
            type    = "warning",
            session = session
          )
        } else {
          log_info("copilot.bindings.plot_callback_failed",
                   msg = conditionMessage(e))
        }
      })
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
