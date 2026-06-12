# copilot_logger.R — Structured key=value logging for the Copilot board.
#
# Depends on:
#   - `info()` from components/app/R/utils/utils.R
#   - `.OPT_COPILOT_TRACE` from copilot_options.R
#
# Format: `[CopilotBoard] <event> key1=val1 key2=val2 ...`. Keys ending in
# `_ms` are rounded to 1 decimal; keys ending in `_s` to 3 decimals. NULL,
# NA, and zero-length values render as `"NULL"`, `"NA"`, `"(empty)"`.

#' Emit a structured INFO log line.
#'
#' @param event character(1) dotted event name (e.g. `copilot.run.start`).
#' @param ... Named key=value pairs.
#' @export
log_info <- function(event, ...) {
  .copilot_log_emit(event, ..., type = "INFO")
  invisible(NULL)
}

#' Emit a structured DEBUG log line, gated on `getOption(.OPT_COPILOT_TRACE)`.
#' @inheritParams log_info
#' @export
log_trace <- function(event, ...) {
  if (!isTRUE(getOption(.OPT_COPILOT_TRACE, FALSE))) return(invisible(NULL))
  .copilot_log_emit(event, ..., type = "DBUG")
  invisible(NULL)
}

#' Capture a Shiny-side timing probe for live debugging.
#'
#' Writes a compact event into `.GlobalEnv$copilot_timing_trace` and mirrors the
#' latest snapshot to `.GlobalEnv$trace` for quick MCPR inspection. If a
#' host-provided `do_global_export()` helper exists, it is called as well.
#'
#' @param event character(1) dotted event name.
#' @param ... Named values to attach to the probe.
#' @export
copilot_debug_timing <- function(event, ...) {
  fields <- list(...)
  snapshot <- list(
    event = as.character(event)[[1L]],
    wall_time = Sys.time(),
    proc_elapsed = unname(proc.time()[["elapsed"]]),
    fields = fields
  )

  previous <- get0("copilot_timing_trace", envir = .GlobalEnv,
                   ifnotfound = list())
  assign("copilot_timing_trace", c(previous, list(snapshot)),
         envir = .GlobalEnv)

  exported <- NULL
  if (exists("do_global_export", mode = "function", inherits = TRUE)) {
    exported <- tryCatch(
      do_global_export(),
      error = function(e) {
        list(error = conditionMessage(e), snapshot = snapshot)
      }
    )
  }

  # DEBUG: timings
  trace <- if (is.null(exported)) snapshot else exported
  assign("trace", trace, envir = .GlobalEnv)

  log_trace("copilot.debug_timing", event = event)
  invisible(snapshot)
}

# ---- Internal emit ---------------------------------------------------------

.copilot_log_emit <- function(event, ..., type = "INFO") {
  pairs <- list(...)
  kv <- if (length(pairs)) {
    paste(mapply(.fmt_kv, names(pairs), pairs,
                 SIMPLIFY = TRUE, USE.NAMES = FALSE),
          collapse = " ")
  } else {
    ""
  }
  msg <- trimws(paste("[CopilotBoard]", event, kv))
  if (exists("info", mode = "function")) {
    info(msg, type = type)
  } else {
    # Fallback for unit-test contexts where components/app/ is not loaded.
    message("[", type, "] ", msg)
  }
  invisible(NULL)
}

# NA-safe scalar formatter for log values.
.fmt_val <- function(x) {
  if (is.null(x))            return("NULL")
  if (length(x) == 0L)       return("(empty)")
  if (length(x) == 1L && is.na(x)) return("NA")
  x
}

# Format a single key=value pair, with duration rounding by key suffix.
.fmt_kv <- function(key, val) {
  if (is.numeric(val) && length(val) == 1L && !is.na(val)) {
    if (endsWith(key, "_ms")) val <- round(val, 1)
    else if (endsWith(key, "_s"))  val <- round(val, 3)
  }
  paste0(key, "=", .fmt_val(val))
}
