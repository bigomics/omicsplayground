##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# MOFA AI Text Server
# =============================================================================
# Shiny module that handles text generation for Summary, Report, and Deep
# Report modes. All pure data extraction lives in mofa_ai_report_data.R.
#
# Prompt assembly uses omicsai structured prompt classes (summary_prompt,
# report_prompt) with frag() for deferred template loading. Deep Report
# dispatches via omicsagentovi::Agent with the MOFA tool surface.

#' MOFA AI Text Server
#'
#' @param id Module namespace ID.
#' @param mofa_reactive Reactive returning the MOFA results object.
#' @param pgx PGX object (non-reactive).
#' @param controls List of control reactives (from
#'   `multiomics_ai_report_controls_server`).
#' @param parent_session Parent Shiny session (for `get_ai_model`).
#'
#' @return List with reactives: `text` (current mode output, possibly the
#'   prompt cache when Show Prompt is on), `report_text`, `deep_text`,
#'   `deep_cost_warning`, `last_deep`.
mofa_ai_text_server <- function(id, mofa_reactive, pgx, controls, parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Prompt template paths ----
    ## Summary mode uses the outer wrapper (mofa_summary.md) which
    ## embeds the inner FACTOR_TEMPLATE via {{factor_detail}}. Report
    ## mode renders the full report data block directly.
    summary_data_path   <- file.path(MOFA_PROMPTS_DIR, "mofa_summary.md")
    interpretation_path <- file.path(MOFA_PROMPTS_DIR, "mofa_interpretation.md")
    report_data_path    <- file.path(MOFA_PROMPTS_DIR, "mofa_report_data.md")
    report_rules_path   <- file.path(MOFA_PROMPTS_DIR, "mofa_report_rules.md")
    agent_skills_path   <- file.path(MOFA_PROMPTS_DIR, "mofa_agent_skills.md")

    # ---- Prompt caches (for instant Show-Prompt toggle) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache  <- shiny::reactiveVal(NULL)
    deep_prompt_cache    <- shiny::reactiveVal(NULL)

    # ---- Data tables cache (passed to diagram/image cards) ----
    report_data_tables   <- shiny::reactiveVal(NULL)

    # ---- Deep Report cost-cap state (soft, post-hoc) ----
    deep_cost_warning    <- shiny::reactiveVal(NULL)
    DEEP_COST_CAP_USD    <- 0.50
    .deep_price <- list(in_fresh = 0.25, in_cached = 0.025, out_ = 2.00)
    .deep_usd_cost <- function(u) {
      na0 <- function(x) if (is.null(x) || is.na(x)) 0L else as.integer(x)
      (na0(u$input_tokens_fresh)  / 1e6) * .deep_price$in_fresh +
      (na0(u$input_tokens_cached) / 1e6) * .deep_price$in_cached +
      (na0(u$output_tokens)       / 1e6) * .deep_price$out_
    }

    # ---- SUMMARY MODE: per-factor summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- get_ai_model(parent_session)
        m <- mofa_reactive()
        factor_name <- controls$selected_module()
        shiny::req(m, factor_name)

        params <- mofa_build_summary_params(m, factor_name, pgx, ntop = 12L)
        shiny::req(params)

        organism <- pgx$organism %||% NULL
        style <- controls$summary_style() %||% "short_summary"

        p <- omicsai::summary_prompt(
          role    = omicsai::frag("system_base"),
          task    = omicsai::frag(paste0("text/", style)),
          species = omicsai::omicsai_species_prompt(organism),
          context = omicsai::frag(interpretation_path,
                                  list(experiment = params$experiment)),
          data    = omicsai::frag(summary_data_path, params)
        )
        bp <- omicsai::build_prompt(p)

        summary_prompt_cache(paste0(
          "# SYSTEM PROMPT\n\n", bp$system,
          "\n\n---\n\n",
          "# USER MESSAGE\n\n", bp$board
        ))

        cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system)
        result <- tryCatch(
          omicsai::omicsai_gen_text(bp$board, config = cfg),
          error = function(e) {
            shiny::validate(shiny::need(FALSE,
              .aicards_friendly_error(conditionMessage(e))))
          }
        )
        result$text
      },
      ignoreNULL = FALSE
    )

    # ---- REPORT MODE: single-shot integrated report ----

    ai_report <- shiny::eventReactive(controls$trigger(), {
      if (controls$mode() != "report" || controls$trigger() < 1) return(NULL)

      model <- get_ai_model(parent_session)
      m <- mofa_reactive()
      shiny::req(m)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      progress$set(message = "Extracting factor data...", value = 0.1)
      tables <- mofa_build_report_tables(m, pgx, n_factors = 8L, ntop = 10L)
      report_data_tables(tables$text)

      organism <- pgx$organism %||% NULL
      experiment_label <- multiomics_ai_experiment_label(pgx, m$experiment)

      methods_context <- mofa_build_methods(m, pgx)

      p <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path,
                                    list(experiment = experiment_label)),
        board_rules = omicsai::frag(report_rules_path),
        data        = tables$text
      )
      bp <- omicsai::build_prompt(p)

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      progress$set(message = "Generating report...", value = 0.3)
      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system,
                                     max_tokens = AI_BUDGETS$report)
      cache <- omicsai::omicsai_cache_init("mem")
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg, cache = cache),
        error = function(e) {
          shiny::validate(shiny::need(FALSE,
            .aicards_friendly_error(conditionMessage(e))))
        }
      )

      if (!grepl("^# ", result$text)) {
        result$text <- paste0("# Analysis Report\n\n", result$text)
      }

      full_report <- paste(result$text, methods_context, sep = "\n\n")
      progress$set(message = "Done!", value = 1)
      full_report
    }, ignoreNULL = FALSE)

    # ---- DEEP REPORT MODE: agentic dispatch with literature retrieval ----

    ai_deep_report <- shiny::eventReactive(controls$deep_trigger(), {
      if (controls$deep_trigger() < 1) return(NULL)

      m <- mofa_reactive()
      shiny::req(m)

      p <- shiny::Progress$new(session)
      on.exit(try(p$close(), silent = TRUE), add = TRUE)

      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] extracting MOFA data...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      p$set(message = "Extracting MOFA data...", value = 0.05)
      tables <- mofa_build_report_tables(m, pgx, n_factors = 8L, ntop = 10L)
      report_data_tables(tables$text)

      p$set(message = "Assembling Deep Report prompt...", value = 0.10)
      organism <- pgx$organism %||% NULL
      experiment_label <- multiomics_ai_experiment_label(pgx, m$experiment)
      methods_text <- mofa_build_methods(m, pgx)

      ## Slot swap: same `context` as Report, board_rules points at the
      ## Deep agent-skills file.
      rp <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path,
                                    list(experiment = experiment_label)),
        board_rules = omicsai::frag(agent_skills_path),
        data        = tables$text
      )
      bp <- omicsai::build_prompt(rp)

      deep_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      p$set(message = "Initializing Deep Report agent (literature retrieval)...",
            value = 0.20)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] constructing agent...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      agent <- omicsagentovi::Agent(
        tier          = "copilot-deep",
        tools         = c("query_mofa", "query_mofa_factor",
                          "query_gene_info",
                          "search_context", "read_context"),
        context       = omicsagentovi::RunContext(pgx = pgx),
        session       = omicsagentovi::AgentSession(
                          session_id = paste0(session$token, "_mofa_deep")
                        ),
        bindings      = omicsagentovi::RunBindings(data_dir = PGX.DIR),
        max_turns     = 50L,
        system_prompt = bp$system
      )
      agent <- omicsagentovi::agent_set_pgx(
        agent,
        pgx          = pgx,
        dataset_name = as.character(pgx$name %||% ""),
        dataset_path = NULL,
        data_dir     = PGX.DIR
      )

      deep_cost_warning(NULL)

      p$set(message = "Running Deep Report agent (this can take 1-3 minutes)...",
            value = 0.30)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] running agent_prompt (max 50 turns)...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      agent_result <- tryCatch(
        omicsagentovi::agent_prompt(agent, bp$board),
        error = function(e) {
          shiny::validate(shiny::need(FALSE,
            .aicards_friendly_error(conditionMessage(e))))
        }
      )
      response <- agent_result$text %||% ""

      ## Post-hoc cost cap check (soft — render output, prepend warning).
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

      if (!grepl("^# ", response)) {
        response <- paste0("# Analysis Report\n\n", response)
      }

      full_report <- paste(response, methods_text, sep = "\n\n")
      p$set(message = "Deep Report complete!", value = 1.0)
      full_report
    }, ignoreNULL = FALSE)

    # ---- Unified text content (instant toggle between result and prompt) ----
    ##
    ## Deep Report has its own trigger, independent of the mode radio.
    ## Remember the last branch fired so toggling Show Prompt while viewing
    ## a Deep Report shows the deep prompt (not whichever mode the radio
    ## happens to be on).
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
      text               = ai_text,
      report_text        = ai_report,
      deep_text          = ai_deep_report,
      deep_cost_warning  = deep_cost_warning,
      last_deep          = shiny::reactive(last_deep()),
      report_data_tables = report_data_tables
    )
  })
}
