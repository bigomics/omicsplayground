# copilot_context_blocks.R — Host-side context-block injection registry.
#
# The chat board stages named text blocks into the agent's pending registry
# via omicsagentovi::agent_inject_text_block() before user turns where the
# LLM needs fresh host-side context (current dataset, AI reports, etc.).
# The package wraps each staged block as <name>text</name> in the next
# user-turn entry point, then clears the registry (consume-once).
#
# Add a provider to .COPILOT_CONTEXT_PROVIDERS to ship a new block — the
# stage helper picks it up automatically. Providers receive the live Agent
# and return either a non-empty character scalar or NULL (skipped).
#
# Stage points (callers of .copilot_stage_context_blocks):
#   - copilot_run_controller.R apply_dataset() first-load
#   - copilot_run_controller.R apply_dataset() switch
#   - copilot_run_controller.R .do_reset() (new chat / tier change)
#   - copilot_restore_controller.R result observer (after restore)
#
# Token cost: blocks are consume-once and only re-stage on dataset change
# or session reset. The composed preamble persists in chat$get_turns() for
# the remainder of the conversation, so each provider should be brief —
# tens of lines, not kilobytes — or rely on prompt-caching at the provider.
# ---- Block providers -----------------------------------------------------

# Each entry: name (string) -> function(agent) -> character(1) | NULL.
.COPILOT_CONTEXT_PROVIDERS <- list(

  current_dataset = function(agent) {
    pgx <- tryCatch(agent@context@pgx, error = function(e) NULL)
    if (is.null(pgx)) return(NULL)

    name <- tryCatch(agent@context@dataset_name, error = function(e) "")
    if (is.null(name) || is.na(name)) name <- ""

    n_samples <- tryCatch(nrow(pgx$samples), error = function(e) NA_integer_)
    n_genes   <- tryCatch(nrow(pgx$X),       error = function(e) NA_integer_)
    contrasts <- tryCatch(names(pgx$contrasts), error = function(e) character(0))
    organism  <- tryCatch(pgx$organism,       error = function(e) NA_character_)
    desc      <- tryCatch(pgx$description,    error = function(e) NA_character_)

    .fmt_field <- function(label, value) {
      if (is.null(value)) return(NULL)
      if (length(value) == 0L) return(NULL)
      if (all(is.na(value))) return(NULL)
      paste0(label, ": ", paste(value, collapse = ", "))
    }
    lines <- c(
      .fmt_field("name",      if (nzchar(name)) name else NULL),
      .fmt_field("organism",  if (!is.na(organism) && nzchar(organism)) organism else NULL),
      .fmt_field("samples",   if (!is.na(n_samples)) n_samples else NULL),
      .fmt_field("genes",     if (!is.na(n_genes))   n_genes   else NULL),
      .fmt_field("contrasts", if (length(contrasts)) contrasts else NULL),
      .fmt_field("description", if (!is.na(desc) && nzchar(desc)) desc else NULL)
    )
    if (!length(lines)) return(NULL)
    paste(lines, collapse = "\n")
  }

  # ---- Future providers ----
  # ai_report  = function(agent) { ... }   # WGCNA / Report board summaries
  # docs       = function(agent) { ... }   # uploaded-doc inventory
  # user_pref  = function(agent) { ... }   # persisted preferences
)

# ---- Stager --------------------------------------------------------------

# Walks .COPILOT_CONTEXT_PROVIDERS, calls each, and stages non-empty blocks
# onto `agent` via omicsagentovi::agent_inject_text_block. Returns the
# updated Agent. NULL agent passes through.
#
# Errors inside a single provider are caught and logged; one bad block does
# not block the others. The package wraps each staged block as
# <name>text</name> when the next user-turn entry point fires.
.copilot_stage_context_blocks <- function(agent) {
  if (is.null(agent)) return(NULL)
  staged <- character(0)
  for (nm in names(.COPILOT_CONTEXT_PROVIDERS)) {
    val <- tryCatch(
      .COPILOT_CONTEXT_PROVIDERS[[nm]](agent),
      error = function(e) {
        log_info("copilot.context.provider_failed",
                 name = nm, msg = conditionMessage(e))
        NULL
      }
    )
    if (is.null(val) || !length(val) || !nzchar(as.character(val)[[1L]])) next
    agent <- tryCatch(
      omicsagentovi::agent_inject_text_block(
        agent, name = nm, text = as.character(val)[[1L]]
      ),
      error = function(e) {
        log_info("copilot.context.inject_failed",
                 name = nm, msg = conditionMessage(e))
        agent
      }
    )
    staged <- c(staged, nm)
  }
  if (length(staged)) {
    log_trace("copilot.context.staged", names = paste(staged, collapse = ","))
  }
  agent
}
