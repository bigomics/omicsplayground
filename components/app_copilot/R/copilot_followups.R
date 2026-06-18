# copilot_followups.R — Side LLM call producing follow-up question
# suggestions after every completed assistant turn.
#
# Architecture: small side call via ../omicsai (NOT the main omicsagentovi
# Agent) so the follow-up prompt does not pollute the agent's transcript,
# turn budget, or tool surface. Cheap, fast, runtime-optional.
#
# Input shape: structured payload built by build_followup_payload()
# (in copilot_followups_inputs.R). The helper sees a strict subset of the
# main agent's knowledge:
#   - last_text          (assistant's final reply)
#   - dataset_context    (same string the main agent gets via current_dataset)
#   - recent_tool_outputs (≤2 trimmed ContentToolResult bodies)
#   - tool_catalog       (live ovi_tools() metadata)
#   - available_reports  (pgx$ai entries with non-empty $report)
#
# Failure modes are silent — if omicsai is missing, the generator factory
# returns NULL; if the LLM call rejects, the promise resolves to character(0).
# Caller (copilot_run_controller$on_done) treats empty results as "skip".

# ---- Prompt template ------------------------------------------------------

.COPILOT_FOLLOWUP_TEMPLATE <- paste(
  "You are suggesting {n} short follow-up questions for a biologist exploring",
  "an omics dataset in a chat interface. The questions are clickable and will",
  "be submitted verbatim as the user's next message to a main analysis agent.",
  "",
  "Write each question in the user's own voice, addressing the agent directly",
  "in the second person — as a request or query the user would type. Start",
  "with verbs like 'Can you', 'Show me', 'What', 'Which', 'How', 'List'.",
  "NEVER use survey/third-person framings like 'Would you like to see...',",
  "'Are you interested in...', 'Do you want to explore...' — those phrasings",
  "treat the user as the answerer, but here the user IS the asker.",
  "",
  "# Dataset",
  "{dataset_context}",
  "",
  "# Available tools (the agent can call these)",
  "{tool_catalog}",
  "",
  "# Attached precomputed reports (TICKED by the user — agent will see",
  "# these on the next user message; empty means the agent will NOT have",
  "# any precomputed report body for the next turn)",
  "{available_reports}",
  "",
  "# Recent tool outputs (most recent last; may be empty)",
  "{recent_tool_outputs}",
  "",
  "# Last assistant reply",
  "{last_text}",
  "",
  "# Output rules",
  "1. Exactly {n} questions, numbered '1.' and '2.', one per line, ending with '?'.",
  "2. Phrase each as a direct request from the user to the agent. Examples of",
  "   the RIGHT voice: 'Can you show me a volcano plot for act48h_vs_notact?',",
  "   'Which L1000 drugs reverse the activation signature?', 'List the top",
  "   genes in module ME1.' Examples of the WRONG voice: 'Would you like to",
  "   see a volcano plot...?', 'Are you interested in exploring...?'.",
  "3. Natural-language only. NEVER output JSON, tool names like `query_de`, or",
  "   function arguments. If a recent tool output suggests",
  "   `omicspgxmcp.show_omics_plot(...)`, translate to plain English",
  "   (e.g. 'Can you show a volcano plot for <contrast>?').",
  "4. PLOT-SUGGESTION HARD RULE — read carefully:",
  "   * IF the `Recent tool outputs` block above is NOT empty (i.e. it does",
  "     NOT say '(none)'): you MAY suggest a plot, but ONLY using",
  "     contrast / feature / module / pathway identifiers that appear",
  "     verbatim in those tool outputs.",
  "   * IF the `Recent tool outputs` block IS empty (shows '(none)'):",
  "     DO NOT SUGGEST ANY PLOTTING FOLLOW-UP. Forbidden phrases include",
  "     'volcano plot', 'heatmap', 'barplot', 'UMAP', 'PCA plot',",
  "     'box plot', 'scatter plot', 'show me a plot of', 'visualize',",
  "     'show_omics_plot', or any request whose answer is a chart.",
  "     Reason: the agent has no verified identifiers (contrast names,",
  "     gene names, module ids) and will either fail or fabricate arguments.",
  "     Instead, suggest DISCOVERY questions whose first hop is itself a",
  "     tool call that returns concrete identifiers, e.g.:",
  "       - 'Can you list the available contrasts?'",
  "       - 'Which WGCNA modules are largest?'",
  "       - 'What are the top differentially expressed genes?'",
  "       - 'Which pathways are most enriched?'",
  "       - 'What does the <report-name> report cover?'",
  "5. REPORTS HARD RULE:",
  "   * Only suggest follow-ups about a precomputed report if that report",
  "     appears in the `Attached precomputed reports` block above. Reports",
  "     not listed there are NOT attached and the agent has no body text",
  "     for them — never suggest 'summarize the WGCNA report' if WGCNA is",
  "     not in the attached block.",
  "   * IF BOTH the `Recent tool outputs` block AND the",
  "     `Attached precomputed reports` block are empty (both show '(none",
  "     available)' / '(none)'): your suggestions MUST be pure",
  "     tool-discovery questions — single-step queries the agent can",
  "     answer by calling one of the listed tools (e.g. 'Can you list the",
  "     available contrasts?', 'Which pathways are enriched?'). No plots,",
  "     no report references, no fabricated entity names.",
  "6. Only reference entities (gene names, contrast names, modules, pathways,",
  "   report names) that appear in the blocks above. Never invent them.",
  "7. Do not ask about wet-lab metadata, sample collection, or sequencing",
  "   protocol unless that exact information is in the Dataset block.",
  "8. No markdown, no preamble, no explanations.",
  "",
  "Questions:",
  sep = "\n"
)

# Render the prompt by substituting payload blocks. Missing blocks are
# substituted as the literal string "(none available)" so the LLM gets a
# consistent surface and can condition cleanly on the empty case.
render_followup_prompt <- function(payload, n = 2L) {
  fmt_outputs <- function(items) {
    if (!length(items)) return("(none)")
    parts <- vapply(items, function(x) {
      nm_raw <- x$name; if (is.null(nm_raw) || !nzchar(nm_raw)) nm_raw <- "tool"
      tx_raw <- x$text; if (is.null(tx_raw)) tx_raw <- ""
      sprintf("## %s\n%s",
              as.character(nm_raw)[[1L]],
              as.character(tx_raw)[[1L]])
    }, character(1L))
    paste(parts, collapse = "\n\n")
  }
  none_if_empty <- function(s) if (is.null(s) || !nzchar(s)) "(none available)" else s

  out <- .COPILOT_FOLLOWUP_TEMPLATE
  out <- gsub("{n}",                  as.character(as.integer(n)),                 out, fixed = TRUE)
  out <- gsub("{dataset_context}",    none_if_empty(payload$dataset_context),      out, fixed = TRUE)
  out <- gsub("{tool_catalog}",       none_if_empty(payload$tool_catalog),         out, fixed = TRUE)
  out <- gsub("{available_reports}",  none_if_empty(payload$available_reports),    out, fixed = TRUE)
  out <- gsub("{recent_tool_outputs}", fmt_outputs(payload$recent_tool_outputs),   out, fixed = TRUE)
  out <- gsub("{last_text}",          none_if_empty(payload$last_text),            out, fixed = TRUE)
  out
}

#' Factory: side LLM follow-up generator
#'
#' Returns a list with a single `generate(payload)` method that yields a
#' promise resolving to a `character()` vector of 0..n follow-up questions.
#'
#' @param model Provider-qualified model id. Default
#'   `"groq:openai/gpt-oss-120b"` matches `../omicsai/R/model_registry.R`.
#' @param n Integer, number of follow-up questions to request (default 2L).
#' @param config Optional pre-built `omicsai_config`. If supplied, `model` is
#'   ignored.
#'
#' @return A list `list(generate = function(payload) <promise>)`, or `NULL`
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

  generate <- function(payload) {
    # Accept legacy plain-string input for backward compat (treats it as
    # last_text-only); caller in copilot_run_controller now passes a list.
    if (is.character(payload)) {
      payload <- list(last_text = payload)
    }
    last_text_raw <- payload$last_text
    if (is.null(last_text_raw)) last_text_raw <- ""
    last_text <- as.character(last_text_raw)[[1L]]
    if (is.null(last_text) || is.na(last_text) || !nzchar(last_text)) {
      return(promises::promise_resolve(character(0)))
    }
    prompt_text <- render_followup_prompt(payload, n = n_)
    promises::future_promise({
      out <- tryCatch(
        omicsai::omicsai_gen_text(
          template = prompt_text,
          params   = list(),
          config   = cfg
        ),
        error = function(e) NULL
      )
      raw_src <- if (is.list(out) && !is.null(out$text)) out$text else out
      raw <- tryCatch(as.character(raw_src), error = function(e) "")
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
