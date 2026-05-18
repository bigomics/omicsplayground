# copilot_options.R — Shared option constants and helpers for the Copilot board.
#
# Pure R, no imports. Loaded before copilot_logger.R (which consumes
# `.OPT_COPILOT_TRACE`) and before copilot_server.R (which consumes
# `COPILOT_TIERS` / `copilot_tier_label()`).

# ---- Authoritative tier enumeration (from contract_audit.md §1) -------------
#' @export
COPILOT_TIERS <- c("copilot-default", "copilot-fast", "copilot-deep")

# ---- Option name constants --------------------------------------------------
#' @export
.OPT_COPILOT_TRACE       <- "copilot.trace"
#' @export
.OPT_COPILOT_REPLAY_MODE <- "copilot.restore_replay_mode"

# ---- Numeric limits ---------------------------------------------------------
#' @export
.COPILOT_MAX_HISTORY <- 100L
#' @export
.COPILOT_MAX_TURNS   <- Inf

# ---- Functions --------------------------------------------------------------

#' Is the Copilot board enabled?
#' @export
is_copilot_enabled <- function() {
  isTRUE(getOption("ENABLE_CHIRP", FALSE))
}

#' Human-readable label for a tier id.
#' Unknown tiers pass through as-is.
#' @export
copilot_tier_label <- function(tier) {
  switch(tier,
    `copilot-default` = "Balanced (Default)",
    `copilot-fast`    = "Fast",
    `copilot-deep`    = "Deep Think",
    tier
  )
}

#' Is the `copilot.trace` option enabled?
#' @export
copilot_trace_enabled <- function() {
  isTRUE(getOption(.OPT_COPILOT_TRACE, FALSE))
}

#' Max history (cap on sessions kept). Reads `copilot.max_history` with default.
#' @export
copilot_max_history <- function(default = .COPILOT_MAX_HISTORY) {
  as.integer(getOption("copilot.max_history", default))
}

#' Restore replay mode — "single" (default) or "batch".
#' Unknown values clamp to "single".
#' @export
copilot_replay_mode <- function() {
  mode <- getOption(.OPT_COPILOT_REPLAY_MODE, "single")
  if (!identical(mode, "batch")) "single" else "batch"
}
