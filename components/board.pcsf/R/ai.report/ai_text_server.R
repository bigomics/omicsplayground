##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# =============================================================================
# AI Text Server for PCSF
# =============================================================================
# Shiny module that handles text generation for both summary and report modes.
# All pure data extraction lives in ai_report_data.R.
#
# Prompt assembly uses omicsai structured prompt classes (summary_prompt,
# report_prompt) with frag() for deferred template loading.

#' PCSF AI Text Server
#'
#' Handles text generation for both summary and report modes.
#' In summary mode, generates a single-contrast summary via omicsai_gen_text().
#' In report mode, generates a full integrated report via omicsai_gen_text()
#' with structured prompt assembly.
#'
#' @param id Module namespace ID
#' @param pgx PGX object (non-reactive)
#' @param controls List of control reactives (from ai_report_controls_server)
#' @param contrast_reactive Reactive returning the currently active contrast
#' @param current_graph_reactive Reactive returning the current PCSF graph
#' @param current_table_reactive Reactive returning the current centrality table
#' @param pcsf_params_reactive Reactive returning PCSF algorithm parameters
#' @param parent_session Parent Shiny session (for getUserOption)
#'
#' @return List with reactives: text (current mode), report_text (report-only)
pcsf_ai_text_server <- function(id,
                                pgx,
                                controls,
                                contrast_reactive,
                                current_graph_reactive,
                                current_table_reactive,
                                pcsf_params_reactive,
                                parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: prompt template paths ----

    summary_data_path   <- file.path(PCSF_PROMPTS_DIR, "pcsf_summary_data.md")
    interpretation_path <- file.path(PCSF_PROMPTS_DIR, "pcsf_interpretation.md")
    report_data_path    <- file.path(PCSF_PROMPTS_DIR, "pcsf_report_data.md")

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-contrast summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module(), contrast_reactive()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- get_ai_model(parent_session)
        fallback_contrast <- contrast_reactive()
        contrast <- controls$selected_module() %||% fallback_contrast
        shiny::req(contrast)

        active_contrast <- fallback_contrast %||% ""
        params <- if (!is.null(active_contrast) && nzchar(active_contrast) && identical(contrast, active_contrast)) {
          pcsf_build_summary_params(
            pgx = pgx,
            contrast = contrast,
            pcsf_graph = current_graph_reactive(),
            centrality_table = current_table_reactive(),
            ntop = 12L,
            pcsf_params = pcsf_params_reactive()
          )
        } else {
          pcsf_build_summary_params(
            pgx = pgx,
            contrast = contrast,
            pcsf_graph = NULL,
            centrality_table = NULL,
            ntop = 12L,
            pcsf_params = pcsf_params_reactive()
          )
        }
        shiny::req(params)

        style    <- controls$summary_style() %||% "short_summary"
        organism <- pgx$organism %||% NULL

        p <- omicsai::summary_prompt(
          role    = omicsai::frag("system_base"),
          task    = omicsai::frag(paste0("text/", style)),
          species = omicsai::omicsai_species_prompt(organism),
          context = omicsai::frag(interpretation_path),
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

      model <- get_ai_model(parent_session)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## Step 1: Build structured data tables
      progress$set(message = "Extracting network data...", value = 0.1)
      tables <- pcsf_build_report_tables(
        pgx = pgx,
        contrasts = pcsf_get_contrasts(pgx),
        max_contexts = 6L,
        ntop = 10L,
        ntop_network = 750L,
        pcsf_params = pcsf_params_reactive()
      )

      ## Step 2: Assemble data content (domain data only, no instructions)
      data_content <- tables$text

      ## Step 3: Build structured prompt
      board_rules_path <- file.path(PCSF_PROMPTS_DIR, "pcsf_report_rules.md")
      organism <- pgx$organism %||% NULL

      ## Pre-render methods/context (reused in Step 5 as deterministic appendix)
      methods_context <- pcsf_build_methods(pgx = pgx)

      p <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path),
        board_rules = omicsai::frag(board_rules_path),
        data        = omicsai::frag(report_data_path, list(report_tables = data_content))
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

      cfg    <- omicsai::omicsai_config(model = model, system_prompt = bp$system, max_tokens = 8192L)
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

    list(text = ai_text, report_text = ai_report)
  })
}
