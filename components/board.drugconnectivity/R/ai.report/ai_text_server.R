##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# AI Text Server for DrugConnectivity
# =============================================================================
# Shiny module that handles text generation for both summary and report modes.
# All pure data extraction lives in ai_report_data.R.
#
# Prompt assembly uses omicsai structured prompt classes (summary_prompt,
# report_prompt) with frag() for deferred template loading.

#' Drug Connectivity AI Text Server
#'
#' Handles text generation for both summary and report modes.
#' In summary mode, generates a single-contrast summary via omicsai_gen_text().
#' In report mode, generates a full integrated report via omicsai_gen_text()
#' with structured prompt assembly. Fixes a previous double-system-prompt bug
#' by using report_prompt() + omicsai_gen_text() instead of omicsai_gen_report().
#'
#' @param id Module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param controls List of control reactives (from ai_report_controls_server)
#' @param method_reactive Reactive returning the currently selected analysis method
#' @param contrast_reactive Reactive returning the currently active contrast
#' @param annotated_only_reactive Reactive returning whether to show annotated drugs only
#' @param parent_session Parent Shiny session (for getUserOption)
#'
#' @return List with reactives: text (current mode), report_text (report-only)
drugconnectivity_ai_text_server <- function(id,
                                            pgx,
                                            controls,
                                            method_reactive,
                                            contrast_reactive,
                                            annotated_only_reactive,
                                            parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: prompt template paths ----

    summary_data_path <- file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_summary_data.md")
    report_data_path  <- file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_report_data.md")
    context_path      <- file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_interpretation.md")

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache  <- shiny::reactiveVal(NULL)

    # ---- Data tables cache (passed to diagram for supporting context) ----
    report_data_tables <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-contrast summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module(), method_reactive(), contrast_reactive()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model          <- get_ai_model(parent_session)
        method         <- method_reactive()
        fallback_contrast <- contrast_reactive()
        contrast       <- controls$selected_module() %||% fallback_contrast
        only_annotated <- isTRUE(annotated_only_reactive())
        shiny::req(method, contrast)

        at <- dc_analysis_type_info(method)

        params <- drugconnectivity_build_summary_params(
          pgx = pgx,
          contrast = contrast,
          method = method,
          only_annotated = only_annotated,
          ntop = 12L
        )
        shiny::req(params)

        style    <- controls$summary_style() %||% "short_summary"
        organism <- pgx$organism %||% NULL

        p <- omicsai::summary_prompt(
          role    = omicsai::frag("system_base"),
          task    = omicsai::frag(paste0("text/", style)),
          species = omicsai::omicsai_species_prompt(organism),
          context = omicsai::frag(context_path, list(
            experiment = params$experiment %||% "omics experiment"
          )),
          data    = omicsai::frag(summary_data_path, params)
        )
        bp <- omicsai::build_prompt(p)

        ## Cache the prompt for instant toggle
        summary_prompt_cache(paste0(
          "# SYSTEM PROMPT\n\n", bp$system,
          "\n\n---\n\n",
          "# USER MESSAGE\n\n", bp$board
        ))

        cfg    <- omicsai::omicsai_config(model = model, system_prompt = bp$system)
        result <- tryCatch(
          omicsai::omicsai_gen_text(bp$board, config = cfg),
          error = function(e) {
            shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
          }
        )
        result$text
      },
      ignoreNULL = FALSE
    )

    # ---- REPORT MODE: single-shot integrated report ----

    ai_report <- shiny::eventReactive(controls$trigger(), {
      if (controls$mode() != "report" || controls$trigger() < 1) return(NULL)

      model          <- get_ai_model(parent_session)
      method         <- method_reactive()
      only_annotated <- isTRUE(annotated_only_reactive())
      shiny::req(method)

      at <- dc_analysis_type_info(method)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## Step 1: Build structured data tables
      progress$set(message = "Extracting drug data...", value = 0.1)
      tables <- drugconnectivity_build_report_tables(
        pgx = pgx,
        method = method,
        only_annotated = only_annotated,
        max_contexts = 8L,
        ntop = 6L
      )
      report_data_tables(tables$text)

      ## Step 2: Assemble data content (domain data only, no instructions)
      data_content <- tables$text

      ## Step 3: Build structured prompt
      board_rules_path <- file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_report_rules.md")
      organism <- pgx$organism %||% NULL

      ## Pre-render methods/context (reused in Step 5 as deterministic appendix)
      methods_context <- drugconnectivity_build_methods(pgx = pgx, method = method)

      p <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_interpretation.md"),
          list(experiment = pgx$name %||% "omics experiment")),
        board_rules = omicsai::frag(board_rules_path),
        data        = data_content
      )
      bp <- omicsai::build_prompt(p)

      ## Cache the full prompt for instant toggle
      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      ## Step 4: Generate report via single LLM call
      ## Use omicsai_gen_text (NOT gen_report — that would double-load text/report.md)
      progress$set(message = "Generating report...", value = 0.3)

      cfg    <- omicsai::omicsai_config(model = model, system_prompt = bp$system, max_tokens = AI_BUDGETS$report)
      cache  <- omicsai::omicsai_cache_init("mem")
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg, cache = cache),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
        }
      )

      ## Safety net: prepend H1 if missing
      if (!grepl("^# ", result$text)) {
        result$text <- paste0("# Analysis Report\n\n", result$text)
      }

      ## Step 5: Append deterministic methods section
      full_report <- paste(result$text, methods_context, sep = "\n\n")

      progress$set(message = "Done!", value = 1)
      full_report
    }, ignoreNULL = FALSE)

    # ---- Unified text content (instant toggle between result and prompt) ----

    ai_text <- shiny::reactive({
      mode <- controls$mode()
      show <- isTRUE(controls$show_prompt())

      if (mode == "summary") {
        if (show) summary_prompt_cache() else ai_summary()
      } else {
        if (show) report_prompt_cache() else ai_report()
      }
    })

    list(text = ai_text, report_text = ai_report, report_data_tables = report_data_tables)
  })
}
