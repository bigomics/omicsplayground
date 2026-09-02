# copilot_options.R — Shared option constants and helpers for the Copilot board.
#
# Pure R, no imports. Loaded before copilot_logger.R (which consumes
# `.OPT_COPILOT_TRACE`) and before copilot_server.R (which consumes
# `COPILOT_TIERS` / `copilot_tier_label()`).

# All constants in this block are internal module config — referenced by
# other files inside components/app_copilot/ via source-time visibility.
# They are intentionally NOT @export'd so that downstream callers cannot
# bind to them as if they were API; rename freely.

# ---- Authoritative tier enumeration (from contract_audit.md §1) -------------
COPILOT_TIERS <- c("copilot-default", "copilot-deep")

# ---- Style preset enumeration ---------------------------------------------
# Mirrors `omicsagentovi::ovi_prompt_styles()`. Single source of truth for
# the style-radio choices in the Copilot chat settings.
COPILOT_STYLES <- c("concise", "technical", "custom")

# Human-readable labels for the style radio (mirrors `ovi_prompt_styles()`).
COPILOT_STYLE_LABELS <- c(
  concise   = "Concise",
  technical = "Technical",
  custom    = "Custom"
)

# ---- Starter questions (chat-column suggestion strip, mimics edgy) ---------
# Hardcoded starter questions rendered as a one-shot button strip in the chat
# column. Strip is dismissed (`shinyjs::hide`) on first user message, then
# re-shown on `chat_event` type `reset` (new chat).
COPILOT_STARTERS <- c(
  "Summarize main findings",
  "What pathways are involved?",
  "Show top biomarkers",
  "Find references",
  "Show a volcano plot"
)

# ---- Option name constants --------------------------------------------------
.OPT_COPILOT_TRACE          <- "copilot.trace"
.OPT_COPILOT_REPLAY_MODE    <- "copilot.restore_replay_mode"
.OPT_COPILOT_SHOW_THINKING  <- "copilot.show_thinking"
.OPT_COPILOT_SHOW_TOOLS     <- "copilot.show_tools"

# ---- Reasoning-trace markers ------------------------------------------------
# Wrapped around each reasoning run by `agent_prompt_stream()` when the
# "Show reasoning trace" switch is on, so the trace reads as an aside rather
# than as the opening of the answer. The closing "\n\n" is what puts the
# answer in its own paragraph.
#
# Styling is inline rather than a class in `.COPILOT_CHAT_CSS`: shinychat
# renders message bodies inside its own web component, so a module stylesheet
# is not guaranteed to reach them. The colour falls back to Bootstrap's muted
# grey when the CSS variable is unavailable.
.COPILOT_THINKING_PREFIX <- paste0(
  "<span style=\"color: var(--bs-secondary-color, #6c757d); ",
  "font-style: italic;\">["
)
.COPILOT_THINKING_SUFFIX <- "]</span>\n\n"

# ---- Numeric limits ---------------------------------------------------------
.COPILOT_MAX_HISTORY <- 100L

# ---- Required in-house packages --------------------------------------------
# Missing any of these hides the Copilot nav entry and skips the server
# wiring. Install from GitHub: `remotes::install_github("bigomics/<pkg>")`.
COPILOT_REQUIRED_PKGS <- c(
  "omicsai", "omicsagentovi", "omicspgx", "omicspgxmcp", "omicsplots"
)

#' Check that every package in `COPILOT_REQUIRED_PKGS` is installed.
#'
#' Logs a single multi-line message listing the missing packages and the
#' install command. Returns `TRUE` only when every package is available.
#' @export
copilot_packages_ok <- function() {
  missing <- Filter(
    function(p) !requireNamespace(p, quietly = TRUE),
    COPILOT_REQUIRED_PKGS
  )
  if (length(missing) > 0L) {
    lines <- c(
      "[copilot] disabled - missing in-house packages:",
      sprintf("  remotes::install_github(\"bigomics/%s\")", missing)
    )
    message(paste(lines, collapse = "\n"))
    return(FALSE)
  }
  TRUE
}

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
    `copilot-default` = "Balanced",
    `copilot-deep`    = "Deep",
    `copilot-fast`    = "Balanced",
    tier
  )
}

#' Human-readable label for a style id.
#' Unknown styles pass through as-is.
#' @export
copilot_style_label <- function(style) {
  lbl <- COPILOT_STYLE_LABELS[[style]]
  if (is.null(lbl)) style else lbl
}

#' Is the `copilot.trace` option enabled?
#' @export
copilot_trace_enabled <- function() {
  isTRUE(getOption(.OPT_COPILOT_TRACE, FALSE))
}

#' Default for the "Show reasoning" switch.
#'
#' The switch itself is per-session and user-visible (chat settings panel);
#' this only seeds its initial value, so a deployment can turn the reasoning
#' trace on by default with
#' `options(copilot.show_thinking = TRUE)` in `global.R`.
#'
#' Reasoning is always requested from the model (both Copilot tiers set
#' `reasoning_effort = "medium"`); this governs display only.
#'
#' @return Logical scalar.
#' @export
copilot_show_thinking_default <- function() {
  isTRUE(getOption(.OPT_COPILOT_SHOW_THINKING, FALSE))
}

#' Default for the "Show tool calls" switch.
#'
#' Defaults to TRUE -- the markers have always been shown, and they are the
#' only signal that a long pause is a tool running rather than a stall.
#'
#' Display only: the switch decides whether the `on_tool_request` callback is
#' handed to `agent_prompt_stream()`, and that callback is a pure
#' notification. Tools still run either way -- what actually disconnects them
#' is the Reports panel's `tools_enabled` flag, which is a different control.
#'
#' Turn the markers off by default with
#' `options(copilot.show_tools = FALSE)` in `global.R`.
#'
#' @return Logical scalar.
#' @export
copilot_show_tools_default <- function() {
  isTRUE(getOption(.OPT_COPILOT_SHOW_TOOLS, FALSE))
}

#' System prompt for Agent construction
#'
#' Returns the string passed to `system_prompt` when constructing a fresh
#' `omicsagentovi::Agent`. Reads the `copilot.system_prompt` R option first;
#' falls back to the default fragment-assembled prompt via
#' `omicsagentovi::ovi_build_system_prompt("concise")`.
#'
#' To use a deployment-specific override, set the option before the app
#' starts (e.g. in `global.R`):
#' `options(copilot.system_prompt = omicsagentovi::ovi_build_system_prompt("technical"))`
#'
#' @return Character scalar system prompt.
#' @export
copilot_system_prompt <- function() {
  getOption("copilot.system_prompt", omicsagentovi::ovi_build_system_prompt("concise"))
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
