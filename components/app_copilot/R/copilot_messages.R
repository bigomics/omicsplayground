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
  no_session       = "ovi_restore: no saved session with that id.",
  # Provider-specific error categories (shown to users instead of raw API text)
  error_quota      = "You've exceeded your current AI provider quota. Please check your plan and billing, or contact us.",
  error_auth       = "Your API key was rejected. Please re-check the key in Settings ▸ AI Features.",
  error_rate       = "The AI provider is rate-limiting requests. Please wait a moment and try again.",
  error_generic    = "Copilot couldn't complete the request. Please try again, or contact us."
)

#' Map a provider error to a friendly user-facing message.
#'
#' Matches known provider failure patterns (quota, auth, rate-limit) and
#' returns the corresponding `COPILOT_MSG` template. Falls back to
#' `error_generic` for unrecognised errors. Raw error text is NOT included
#' in the returned string — log it separately via `log_info()`.
#'
#' @param e_or_msg A condition object or character(1) message string.
#' @return character(1) user-facing message from `COPILOT_MSG`.
#' @export
copilot_friendly_error <- function(e_or_msg) {
  msg <- if (is.character(e_or_msg)) e_or_msg else conditionMessage(e_or_msg)
  msg_lower <- tolower(msg)
  # Quota / billing exhausted
  if (grepl("insufficient_quota|exceeded.*quota|billing|credit", msg_lower) ||
      (grepl("429", msg) && grepl("quota|billing|credit", msg_lower))) {
    return(copilot_msg("error_quota"))
  }
  # Auth / invalid key. Note: no bare `api.*key` — provider rate-limit
  # messages often link to ".../api-keys", which must not be read as an
  # auth failure (rate-limit is categorised below).
  if (grepl("401|invalid.?api.?key|authentication|unauthorized", msg_lower)) {
    return(copilot_msg("error_auth"))
  }
  # Rate limiting (429 not already matched by quota above)
  if (grepl("429|rate.?limit|too many request", msg_lower)) {
    return(copilot_msg("error_rate"))
  }
  copilot_msg("error_generic")
}

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
