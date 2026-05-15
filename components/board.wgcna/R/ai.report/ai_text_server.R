## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# AI Text Server for WGCNA
# =============================================================================
# Shiny module that handles text generation for both summary and report modes.
# All pure data extraction lives in ai_report_data.R.
#
# Prompt assembly uses omicsai structured prompt classes (summary_prompt,
# report_prompt) with frag() for deferred template loading.

#' WGCNA AI Text Server
#'
#' Handles text generation for both summary and report modes.
#' In summary mode, generates a single-module summary via omicsai_gen_text().
#' In report mode, generates a full integrated report via omicsai_gen_text()
#' with structured prompt assembly.
#'
#' @param id Module namespace ID
#' @param wgcna Reactive returning WGCNA results object
#' @param pgx PGX object (non-reactive)
#' @param controls List of control reactives (from ai_report_controls_server)
#' @param parent_session Parent Shiny session (for getUserOption)
#'
#' @return List with reactives: text (current mode), report_text (report-only)
wgcna_ai_text_server <- function(id, wgcna, pgx, controls, parent_session,
                                  progress_reactive = NULL) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: prompt template paths ----

    summary_data_path   <- file.path(BOARD_PROMPTS_DIR, "wgcna_summary_data.md")
    interpretation_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_interpretation.md")

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)
    deep_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- Data tables cache (passed to diagram for supporting context) ----
    report_data_tables <- shiny::reactiveVal(NULL)

    # ---- Deep Report state (cost-cap warning surfaced via response prefix) ----
    deep_cost_warning <- shiny::reactiveVal(NULL)

    ## v03 Deep Report cost cap (USD, post-hoc, soft warning).
    ## Pricing per 1M tokens for the copilot-deep tier (gpt-5.4-mini class).
    ## Mirrors tmp/.../v03_iterate.R:69-78.
    DEEP_COST_CAP_USD <- 0.50
    .deep_price <- list(in_fresh = 0.25, in_cached = 0.025, out_ = 2.00)
    .deep_usd_cost <- function(u) {
      na0 <- function(x) if (is.null(x) || is.na(x)) 0L else as.integer(x)
      (na0(u$input_tokens_fresh)  / 1e6) * .deep_price$in_fresh +
      (na0(u$input_tokens_cached) / 1e6) * .deep_price$in_cached +
      (na0(u$output_tokens)       / 1e6) * .deep_price$out_
    }

    # ---- SUMMARY MODE: per-module summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
    {
      mode <- controls$mode()
      if (mode != "summary" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      module <- controls$selected_module()
      shiny::req(w, module)
      model <- get_ai_model(parent_session)

      style <- "short_summary"
      data_params <- wgcna_build_summary_params(w, module, pgx)
      organism <- pgx$organism %||% NULL

      p <- omicsai::summary_prompt(
        role    = omicsai::frag("system_base"),
        task    = omicsai::frag(paste0("text/", style)),
        species = omicsai::omicsai_species_prompt(organism),
        context = omicsai::frag(interpretation_path),
        data    = omicsai::frag(summary_data_path, data_params)
      )
      bp <- omicsai::build_prompt(p)

      ## Cache the prompt for instant toggle
      summary_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system)
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
        }
      )
      result$text
    }, ignoreNULL = FALSE)

    # ---- REPORT MODE: single-shot integrated report ----

    ai_report <- shiny::eventReactive(controls$trigger(), {
      mode <- controls$mode()
      if (mode != "report" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      shiny::req(w)
      model <- get_ai_model(parent_session)

      # Use coordinator progress if available, otherwise create local
      p <- if (!is.null(progress_reactive)) progress_reactive() else NULL
      local_progress <- is.null(p)
      if (local_progress) {
        p <- shiny::Progress$new()
      }

      ## Step 1: Build structured data tables
      message(sprintf("[INFO][%s] --- [AI-REPORT] extracting module data...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Extracting module data...", value = 0.05)
      tables <- wgcna_build_report_tables(w, pgx)
      report_data_tables(tables$text)

      ## Step 2: Build structured prompt
      ## v03: tables$text carries the dual-trait data block; the module-tier
      ## label is already rendered per row by wgcna_build_report_tables, so
      ## the prior wgcna_rank_modules() injection is no longer needed.
      message(sprintf("[INFO][%s] --- [AI-REPORT] assembling prompt...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Assembling prompt...", value = 0.10)
      board_rules_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_report_rules.md")
      organism <- pgx$organism %||% NULL

      ## Pre-render methods/context (reused in Step 6 as deterministic appendix)
      methods_text <- wgcna_build_methods(w, pgx)

      rp <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path),
        board_rules = omicsai::frag(board_rules_path),
        data        = tables$text
      )
      bp <- omicsai::build_prompt(rp)

      ## Cache the full prompt for instant toggle
      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      ## Step 5: Generate report via single LLM call
      ## Use omicsai_gen_text (NOT gen_report — that would double-load text/report.md)
      message(sprintf("[INFO][%s] --- [AI-REPORT] generating AI report (model: %s)...", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), model))
      if (!is.null(p)) p$set(message = "Generating AI report...", value = 0.20)

      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system, max_tokens = AI_BUDGETS$report)
      cache <- omicsai::omicsai_cache_init("mem")
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg, cache = cache),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
        }
      )

      ## Safety net: prepend H1 if missing (same as gen_report does)
      if (!grepl("^# ", result$text)) {
        result$text <- paste0("# Analysis Report\n\n", result$text)
      }

      ## Step 6: Append deterministic methods section
      full_report <- paste(result$text, methods_text, sep = "\n\n")

      message(sprintf("[INFO][%s] --- [AI-REPORT] report text complete", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Report text complete!", value = 0.65)
      # Don't close progress — coordinator closes after diagram completes
      if (local_progress && !is.null(p)) p$close()
      full_report
    }, ignoreNULL = FALSE)

    # ---- DEEP REPORT MODE: agentic report with literature retrieval ----
    ##
    ## Shares wgcna_build_report_tables + wgcna_build_methods with the
    ## Report branch. Differences:
    ##   - `board_rules` slot swapped from wgcna_report_rules.md to
    ##     wgcna_agent_skills.md (Phase A/B agent skills + section rules
    ##     for the Deep variant). `context` stays at wgcna_interpretation.md.
    ##   - Dispatch via omicsagentovi::Agent + agent_prompt() instead of
    ##     omicsai::omicsai_gen_text. The 6-arg constructor signature
    ##     mirrors tmp/.../v03_iterate.R:318-327 (the verified campaign
    ##     call); tools, max_turns, system_prompt, chat_args are explicit.
    ##   - Post-hoc cost cap ($0.50 default): if exceeded, the agent
    ##     output is still rendered but prefixed with a visible warning
    ##     (per Phase 5 brief — soften the plan's hard halt to a warn).
    ##   - Tool-call trace is intentionally NOT surfaced to the UI: keeping
    ##     the agent's pipeline opaque to the end user; the agent's text
    ##     output already carries Vancouver citations and a Bibliography.
    ##
    ## NOTE: synchronous agent_prompt() is used (not _async). The
    ## campaign verified sync, and the existing Report branch is also
    ## sync (omicsai_gen_text blocks). Wrapping this in an ExtendedTask
    ## for true non-blocking dispatch is out of scope for Phase 5 — see
    ## the worker report.
    ai_deep_report <- shiny::eventReactive(controls$deep_trigger(), {
      if (controls$deep_trigger() < 1) return(NULL)

      w <- wgcna()
      shiny::req(w)

      p <- shiny::Progress$new(session)
      on.exit(try(p$close(), silent = TRUE), add = TRUE)

      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] extracting module data...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      p$set(message = "Extracting module data...", value = 0.05)
      tables <- wgcna_build_report_tables(w, pgx)
      report_data_tables(tables$text)

      p$set(message = "Assembling Deep Report prompt...", value = 0.10)
      agent_skills_path <- file.path(BOARD_PROMPTS_DIR, "wgcna_agent_skills.md")
      organism <- pgx$organism %||% NULL

      methods_text <- wgcna_build_methods(w, pgx)

      ## Slot swap: same `context` as Report, board_rules points at the
      ## Deep agent-skills file.
      rp <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path),
        board_rules = omicsai::frag(agent_skills_path),
        data        = tables$text
      )
      bp <- omicsai::build_prompt(rp)

      deep_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      ## Construct agent — verified 6-arg signature from the campaign
      ## (tmp/prompt_optim_aligment/candidate_agentic/wgcna/v03_iterate.R:318-327).
      p$set(message = "Initializing Deep Report agent (literature retrieval)...", value = 0.20)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] constructing agent...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      agent <- omicsagentovi::Agent(
        tier          = "copilot-deep",
        tools         = c("search_context", "read_context", "query_wgcna",
                          "manage_pgx", "query_gene_info"),
        context       = omicsagentovi::ovi_context(
                          pgx        = pgx,
                          session_id = paste0(session$token, "_wgcna_deep")
                        ),
        max_turns     = 50L,
        system_prompt = bp$system,
        chat_args     = list(api_args = list(reasoning = list(summary = "auto")))
      )

      deep_cost_warning(NULL)

      p$set(message = "Running Deep Report agent (this can take 1–3 minutes)...", value = 0.30)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] running agent_prompt (max 50 turns)...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      response <- tryCatch(
        omicsagentovi::agent_prompt(agent, bp$board),
        error = function(e) {
          shiny::validate(shiny::need(FALSE,
            .aicards_friendly_error(conditionMessage(e))))
        }
      )

      ## Post-hoc cost cap check (soft — render output, prepend warning).
      ## In-house package contract: agent_usage_summary must be present;
      ## a missing/erroring call is a real bug, not graceful-degradation.
      usage <- omicsagentovi::agent_usage_summary(agent)
      cost <- .deep_usd_cost(usage)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] usage: in_fresh=%d in_cached=%d out=%d est_cost=$%.4f",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      usage$input_tokens_fresh %||% 0L,
                      usage$input_tokens_cached %||% 0L,
                      usage$output_tokens %||% 0L,
                      cost))
      if (isTRUE(cost > DEEP_COST_CAP_USD)) {
        warn_msg <- sprintf(
          "Deep Report cost (~$%.2f) exceeded the $%.2f cap — partial output below.",
          cost, DEEP_COST_CAP_USD)
        deep_cost_warning(warn_msg)
        response <- paste0("> **Warning:** ", warn_msg, "\n\n", response)
      }

      ## Safety net: prepend H1 if missing
      if (!grepl("^# ", response)) {
        response <- paste0("# Analysis Report\n\n", response)
      }

      ## Append deterministic methods section (same as Report branch)
      full_report <- paste(response, methods_text, sep = "\n\n")

      p$set(message = "Deep Report complete!", value = 1.0)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] complete (%d turns)",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      length(turns)))
      full_report
    }, ignoreNULL = FALSE)

    # ---- Unified text content (instant toggle between result and prompt) ----
    ##
    ## Deep Report has its own trigger, independent of the mode radio. We
    ## remember the last branch fired so toggling Show Prompt while
    ## viewing a Deep Report shows the deep prompt (not whichever mode
    ## the radio happens to be on).
    last_deep <- shiny::reactiveVal(FALSE)
    shiny::observeEvent(controls$deep_trigger(), {
      if (controls$deep_trigger() >= 1) last_deep(TRUE)
    }, ignoreInit = TRUE)
    shiny::observeEvent(controls$trigger(), {
      if (controls$trigger() >= 1) last_deep(FALSE)
    }, ignoreInit = TRUE)

    ai_text <- shiny::reactive({
      mode <- controls$mode()
      show <- isTRUE(controls$show_prompt())

      if (isTRUE(last_deep())) {
        if (show) deep_prompt_cache() else ai_deep_report()
      } else if (mode == "summary") {
        if (show) summary_prompt_cache() else ai_summary()
      } else {
        if (show) report_prompt_cache() else ai_report()
      }
    })

    list(
      text = ai_text,
      report_text = ai_report,
      deep_text = ai_deep_report,
      deep_cost_warning = deep_cost_warning,
      last_deep = shiny::reactive(last_deep()),
      report_data_tables = report_data_tables
    )
  })
}
