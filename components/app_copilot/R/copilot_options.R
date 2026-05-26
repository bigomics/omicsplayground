# copilot_options.R — Shared option constants and helpers for the Copilot board.
#
# Pure R, no imports. Loaded before copilot_logger.R (which consumes
# `.OPT_COPILOT_TRACE`) and before copilot_server.R (which consumes
# `COPILOT_TIERS` / `copilot_tier_label()`).

# ---- Authoritative tier enumeration (from contract_audit.md §1) -------------
#' @export
COPILOT_TIERS <- c("copilot-default", "copilot-fast", "copilot-deep")

# ---- Starter questions (chat-column suggestion strip, mimics edgy) ---------
#' Hardcoded starter questions rendered as a one-shot button strip in the chat
#' column. Strip is dismissed (`shinyjs::hide`) on first user message, then
#' re-shown on `chat_event` type `reset` (new chat).
#' @export
COPILOT_STARTERS <- c(
  "Summarize main findings",
  "What pathways are involved?",
  "Show top biomarkers",
  "Find references",
  "Get differential expression top genes"
)

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

# ---- Required in-house packages --------------------------------------------
#' Required in-house packages for the Copilot to boot.
#' Missing any one of these hides the Copilot nav entry and skips the server
#' wiring. Install from GitHub: `remotes::install_github("bigomics/<pkg>")`.
#' @export
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

#' System prompt for Agent construction
#'
#' Returns the string passed to `system_prompt` when constructing a fresh
#' `omicsagentovi::Agent`. Reads the `copilot.system_prompt` R option first;
#' falls back to the default bundled prompt via
#' `omicsagentovi::ovi_system_prompt()`.
#'
#' To use a board-specific or deployment-specific prompt, set the option
#' before the app starts (e.g. in `global.R`):
#' `options(copilot.system_prompt = omicsagentovi::ovi_system_prompt("agent_system_report"))`
#'
#' @return Character scalar system prompt, or `NA_character_` if no bundled
#'   file is found and no option is set.
#' @export
copilot_system_prompt <- function() {
  getOption("copilot.system_prompt", omicsagentovi::ovi_system_prompt())
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
