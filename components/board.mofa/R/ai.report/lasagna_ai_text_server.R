##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# LASAGNA AI Text Server
# =============================================================================
# Shiny module that handles text generation for Summary, Report, and Deep
# Report modes. All pure data extraction lives in lasagna_ai_report_data.R.
#
# Deep Report dispatches via omicsagentovi::Agent with the LASAGNA tool
# surface (query_de, query_pathways, query_correlation, query_gene_info,
# search_context, read_context).

#' LASAGNA AI Text Server
#'
#' @return List with reactives: `text` (current mode output), `report_text`,
#'   `deep_text`, `deep_cost_warning`, `last_deep`.
lasagna_ai_text_server <- function(id,
                                   graph_data_reactive,
                                   contrast_reactive,
                                   pgx,
                                   controls,
                                   parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Prompt template paths ----
    summary_data_path   <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_summary.md")
    interpretation_path <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_interpretation.md")
    report_rules_path   <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_report_rules.md")
    agent_skills_path   <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_agent_skills.md")

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

    # ---- SUMMARY MODE: lead-module summary for the active contrast ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module(),
           contrast_reactive()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- get_ai_model(parent_session)
        contrast <- controls$selected_module() %||% contrast_reactive()
        res <- graph_data_reactive()
        shiny::req(contrast, res)

        params <- lasagna_build_summary_params(res, contrast, pgx, ntop = 12L)
        shiny::req(params)

        organism <- pgx$organism %||% NULL
        style <- controls$summary_style() %||% "short_summary"

        p <- omicsai::summary_prompt(
          role    = omicsai::frag("system_base"),
          task    = omicsai::frag(paste0("text/", style)),
          species = omicsai::omicsai_species_prompt(organism),
          context = omicsai::frag(interpretation_path),
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
      contrast <- contrast_reactive()
      res <- graph_data_reactive()
      shiny::req(contrast, res)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      progress$set(message = "Extracting contrast data...", value = 0.1)
      tables <- lasagna_build_report_tables(res, contrast, pgx,
                                            n_modules = 8L, ntop = 10L)
      report_data_tables(tables$text)

      organism <- pgx$organism %||% NULL
      methods_context <- lasagna_build_methods(pgx, contrast)

      p <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path),
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
      cfg <- omicsai::omicsai_config(model = model,
                                     system_prompt = bp$system,
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

      contrast <- contrast_reactive()
      res <- graph_data_reactive()
      shiny::req(contrast, res)

      p <- shiny::Progress$new(session)
      on.exit(try(p$close(), silent = TRUE), add = TRUE)

      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] extracting LASAGNA data...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      p$set(message = "Extracting LASAGNA data...", value = 0.05)
      tables <- lasagna_build_report_tables(res, contrast, pgx,
                                            n_modules = 8L, ntop = 10L)
      report_data_tables(tables$text)

      p$set(message = "Assembling Deep Report prompt...", value = 0.10)
      organism <- pgx$organism %||% NULL
      methods_text <- lasagna_build_methods(pgx, contrast)

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

      p$set(message = "Initializing Deep Report agent (literature retrieval)...",
            value = 0.20)
      message(sprintf("[INFO][%s] --- [AI-DEEP-REPORT] constructing agent...",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

      agent <- omicsagentovi::Agent(
        tier          = "copilot-deep",
        tools         = c("query_de", "query_pathways",
                          "query_correlation", "query_gene_info",
                          "search_context", "read_context"),
        context       = omicsagentovi::RunContext(pgx = pgx),
        session       = omicsagentovi::AgentSession(
                          session_id = paste0(session$token, "_lasagna_deep")
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

    # ---- Unified text content (instant Show-Prompt toggle) ----
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
