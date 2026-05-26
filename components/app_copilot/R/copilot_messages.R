# copilot_messages.R — User-facing message templates for the Copilot board.
#
# Pure R, no imports. The `{{key}}` placeholders are substituted via
# `substitute_template()`. Unknown keys are left in place.

#' Substitute `{{key}}` placeholders in a template using a named params list.
#'
#' @param template character(1).
#' @param params Named list of replacements; values coerced via `as.character`.
#' @return character(1).
#' @export
substitute_template <- function(template, params) {
  if (length(params) == 0L) return(template)
  Reduce(
    function(s, key) gsub(paste0("{{", key, "}}"),
                          as.character(params[[key]]),
                          s, fixed = TRUE),
    names(params),
    init = template
  )
}

#' Canonical message templates for the Copilot board.
#'
#' All user-facing strings live here so product wording is changed in one place
#' and so call-sites read as `copilot_msg("key", ...)` rather than inline prose.
#' @export
COPILOT_MSG <- list(
  greeting         = "Hi — load a dataset and ask me anything about your experiment.",
  greeting_active  = "Hi — what would you like to explore?",
  restore_started  = "Restoring previous session…",
  restore_failed   = "Restore failed. You can continue with the current session or pick another from history.",
  no_dataset       = "Load a dataset first, then I can answer questions.",
  max_turns        = "This session has reached the maximum of {{n}} user turns. Start a new chat to continue.",
  plot_failed      = "Plot rendering failed: {{msg}}",
  agent_failed     = "Copilot: failed to create agent — {{msg}}",
  switch_failed    = "Copilot: failed to switch dataset — {{msg}}",
  error_prefix     = "Copilot error: {{msg}}",
  restore_busy     = "Restore is unavailable while the agent is responding.",
  model_gone       = "ovi_restore: cannot restore session — model is not resolvable. Pick a different model.",
  no_session       = "ovi_restore: no saved session with that id."
)

#' Render a `COPILOT_MSG` template with optional named replacements.
#'
#' @param key character(1) name in `COPILOT_MSG`.
#' @param ... Named arguments passed as the params list.
#' @return character(1). Returns the key itself if it is not registered.
#' @export
copilot_msg <- function(key, ...) {
  template <- COPILOT_MSG[[key]]
  if (is.null(template)) return(key)
  substitute_template(template, list(...))
}
