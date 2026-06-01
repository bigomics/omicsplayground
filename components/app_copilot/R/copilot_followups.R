# copilot_followups.R — Side LLM call that produces 2 follow-up question
# suggestions after every completed assistant turn.
#
# Architecture: small side call via ../omicsai (NOT the main omicsagentovi
# Agent) so the follow-up prompt does not pollute the agent's transcript,
# turn budget, or tool surface. Cheap, fast, runtime-optional.
#
# Failure modes are silent — if omicsai is missing, the generator factory
# returns NULL; if the LLM call rejects, the promise resolves to character(0).
# Caller (copilot_run_controller$on_done) treats empty results as "skip".

# ---- Prompt template (markdown-free, n is interpolated) --------------------
.COPILOT_FOLLOWUP_TEMPLATE <- paste0(
  "Given the assistant response below, suggest {n} short follow-up questions ",
  "a biologist exploring this data might ask next. ",
  "Return a clean numbered list, no markdown, no preamble, no explanations. ",
  "Each question must be a single line, end with a question mark, and be ",
  "concrete (not 'Tell me more').\n\n",
  "Assistant response:\n{last_text}\n"
)

#' Factory: side LLM follow-up generator
#'
#' Returns a list with a single `generate(last_text)` method that yields a
#' promise resolving to a `character()` vector of 0..n follow-up questions.
#'
#' @param model Provider-qualified model id. Default
#'   `"groq:openai/gpt-oss-120b"` matches `../omicsai/R/model_registry.R`.
#' @param n Integer, number of follow-up questions to request (default 2L).
#' @param config Optional pre-built `omicsai_config`. If supplied, `model` is
#'   ignored.
#'
#' @return A list `list(generate = function(last_text) <promise>)`, or `NULL`
#'   if the `omicsai` package is not installed. Callers must guard with
#'   `is.null(...)`.
#' @export
make_followup_generator <- function(
  model  = getOption("copilot.followup_model", "groq:openai/gpt-oss-120b"),
  n      = 2L,
  config = NULL
) {
  if (!requireNamespace("omicsai", quietly = TRUE)) return(NULL)
  if (!requireNamespace("promises", quietly = TRUE)) return(NULL)

  cfg <- if (!is.null(config)) config else {
    tryCatch(omicsai::omicsai_config(model = model),
             error = function(e) NULL)
  }
  if (is.null(cfg)) return(NULL)

  n_ <- as.integer(n)

  generate <- function(last_text) {
    if (is.null(last_text) || !nzchar(last_text)) {
      return(promises::promise_resolve(character(0)))
    }
    promises::future_promise({
      out <- tryCatch(
        omicsai::omicsai_gen_text(
          template = .COPILOT_FOLLOWUP_TEMPLATE,
          params   = list(n = n_, last_text = last_text),
          config   = cfg
        ),
        error = function(e) NULL
      )
      raw <- tryCatch(as.character(out$text %||% out), error = function(e) "")
      if (length(raw) == 0L || !nzchar(raw[[1]])) return(character(0))
      parse_followup_list(raw[[1]], n = n_)
    })
  }

  list(generate = generate)
}

#' Parse a raw LLM response into 0..n follow-up question strings.
#'
#' Accepts numbered (`1. `, `2) `), dashed (`- `), or bulleted (`* `, `• `)
#' lists. Strips inline markdown (`**bold**`, `*italic*`, `_under_`, backticks).
#' Returns at most `n` matches; may return fewer.
#'
#' @param raw_text Single character string from the LLM.
#' @param n Integer cap (default 2L).
#' @return `character()` of length 0..n.
#' @export
parse_followup_list <- function(raw_text, n = 2L) {
  if (is.null(raw_text) || length(raw_text) == 0L) return(character(0))
  raw <- as.character(raw_text)[[1]]
  if (!nzchar(raw)) return(character(0))

  lines <- strsplit(raw, "\r?\n", perl = TRUE)[[1]]
  pat   <- "^\\s*(?:[0-9]+[.)]|[-*•])\\s+(.+?)\\s*$"

  matches <- regmatches(lines, regexec(pat, lines, perl = TRUE))
  out <- vapply(matches, function(m) {
    if (length(m) >= 2L) m[[2]] else NA_character_
  }, character(1))
  out <- out[!is.na(out) & nzchar(out)]

  # Strip inline markdown: **bold**, *italic*, _under_, `code`.
  out <- gsub("\\*\\*(.*?)\\*\\*", "\\1", out, perl = TRUE)
  out <- gsub("\\*(.+?)\\*",       "\\1", out, perl = TRUE)
  out <- gsub("_(.+?)_",           "\\1", out, perl = TRUE)
  out <- gsub("`(.+?)`",           "\\1", out, perl = TRUE)
  out <- trimws(out)
  out <- out[nzchar(out)]

  n_ <- as.integer(n)
  if (length(out) > n_) out <- out[seq_len(n_)]
  out
}

#' Render follow-up questions as a shinychat suggestion bubble.
#'
#' Produces `<ul><li class='suggestion submit'>...</li></ul>`. shinychat
#' auto-submits the `<li>` text into `input$chat_user_input` on click.
#'
#' HTML-escapes every question via `htmltools::htmlEscape` — LLM output is
#' untrusted and may contain `<`, `>`, `&`, quotes, or partial tags.
#'
#' @param questions `character()` from `parse_followup_list`.
#' @return Character scalar; empty string if `length(questions) == 0L`.
#' @export
format_followup_bubble <- function(questions) {
  if (length(questions) == 0L) return("")
  escaped <- htmltools::htmlEscape(questions)
  items <- paste0(
    "<li class='suggestion submit'>", escaped, "</li>",
    collapse = ""
  )
  paste0("<ul>", items, "</ul>")
}
