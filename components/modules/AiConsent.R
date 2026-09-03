##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## AI data-sharing consent — per-user, persisted to disk.
##
## Every other AI setting (provider, credentials, the four model menus) lives
## only in `session$userData` and resets on logout. That is fine for a
## preference; it is not fine for a consent. A consent that silently reverts to
## its default on the next login is not a consent, so this one field gets real
## per-user persistence, modelled on the plot colour theme
## (components/ui/ui-ColorDefaults.R): a small JSON file in the user's dir,
## written on change and read back on login.
##
## Default is FALSE — opt-in, not opt-out. A user who never opens Settings has
## not consented, and their calls go out under omicsai's private routing
## defaults (see omicsai::omicsai_privacy_extra).
##
## Scope note: the switch governs *conversation* data, i.e. Copilot chats. AI
## reports, infographics and card summaries are not conversations and always
## run under the private default; they never read this flag. AI usage telemetry
## (components/modules/AiTelemetry.R) records token counts and cost only, never
## prompt or completion text, and is operational data that runs regardless.

AI_CONSENT_FILE <- "ai_consent.json"

## Bump when the consent wording changes materially enough that a previously
## recorded opt-in should no longer be considered informed. `load_ai_consent()`
## treats a stored record with an older version as "no consent", which
## re-prompts the user by returning the switch to its default position.
AI_CONSENT_VERSION <- 1L

#' Persist the AI data-sharing consent for a user
#'
#' @param user_dir The user's personal directory (`auth$user_dir`). A NULL or
#'   non-existent directory is a no-op — anonymous/no-auth deployments have no
#'   place to write, and their consent stays session-scoped.
#' @param consent Logical scalar; coerced with `isTRUE()`.
#' @return invisible path written, or invisible NULL.
save_ai_consent <- function(user_dir, consent) {
  if (is.null(user_dir) || !nzchar(user_dir) || !dir.exists(user_dir)) {
    return(invisible(NULL))
  }
  path <- file.path(user_dir, AI_CONSENT_FILE)
  record <- list(
    share_data   = isTRUE(consent),
    version      = AI_CONSENT_VERSION,
    consented_at = as.numeric(Sys.time())
  )
  tryCatch(
    {
      jsonlite::write_json(record, path, auto_unbox = TRUE, pretty = TRUE)
      invisible(path)
    },
    error = function(e) {
      warning("[save_ai_consent] could not write ", path, ": ",
              conditionMessage(e), call. = FALSE)
      invisible(NULL)
    }
  )
}

#' Read the persisted AI data-sharing consent for a user
#'
#' Fails closed: a missing file, an unreadable file, or a record written under
#' an older `AI_CONSENT_VERSION` all read as FALSE.
#'
#' @param user_dir The user's personal directory (`auth$user_dir`).
#' @return Logical scalar.
load_ai_consent <- function(user_dir) {
  if (is.null(user_dir) || !nzchar(user_dir)) {
    return(FALSE)
  }
  path <- file.path(user_dir, AI_CONSENT_FILE)
  if (!file.exists(path)) {
    return(FALSE)
  }
  record <- tryCatch(
    jsonlite::read_json(path, simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(record) || !isTRUE(record$share_data)) {
    return(FALSE)
  }
  ## Stale consent wording — treat as not given so the user is asked again.
  isTRUE(as.integer(record$version %||% 0L) >= AI_CONSENT_VERSION)
}

#' Read the current session's AI data-sharing consent
#'
#' The session copy is authoritative at runtime; the JSON file is what seeds it
#' at login (see the `auth$logged` observer in appsettings_server.R).
#'
#' @param parent_session Parent Shiny session.
#' @return Logical scalar; FALSE when unset.
get_ai_data_consent <- function(parent_session) {
  isTRUE(getUserOption(parent_session, "ai_share_data"))
}
