## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

# =============================================================================
# AI Text Server for WGCNA
# =============================================================================
# Shiny module that handles text generation for both summary and report modes.
# All pure data extraction lives in ai_report_data.R.

#' WGCNA AI Text Server
#'
#' Handles text generation for both summary and report modes.
#' In summary mode, generates a single-module summary via omicsai_gen_text().
#' In report mode, generates multi-module summaries then integrates via omicsai_create_report().
#'
#' @param id Module namespace ID
#' @param wgcna Reactive returning WGCNA results object
#' @param pgx PGX object (non-reactive)
#' @param controls List of control reactives (from ai_report_controls_server)
#' @param parent_session Parent Shiny session (for getUserOption)
#'
#' @return List with reactives: text (current mode), report_text (report-only)
wgcna_ai_text_server <- function(id, wgcna, pgx, controls, parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: templates, context, model ----

    ai_summary_template <- omicsai::omicsai_load_template(
      file.path(OPG, "components/board.wgcna/prompts/wgcna_summary_template.md")
    )

    ai_context_template <- omicsai::omicsai_load_template(
      file.path(OPG, "components/board.wgcna/prompts/WGCNA_methods.md")
    )

    ai_context <- shiny::reactive({
      w <- wgcna()
      shiny::req(w)
      omicsai::omicsai_substitute_template(
        ai_context_template,
        list(experiment = w$experiment %||% "omics experiment")
      )
    })

    ai_model <- shiny::reactive({
      m <- getUserOption(parent_session, "llm_model")
      shiny::req(m, m != "")
      m
    })

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-module summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
    {
      mode <- controls$mode()
      if (mode != "summary" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      module <- controls$selected_module()
      shiny::req(w, module)
      model <- ai_model()
      shiny::req(model)

      style <- controls$summary_style() %||% "short"
      params <- wgcna_build_summary_params(w, module, pgx)
      params$methods_context <- ai_context()
      params$style_instructions <- omicsai::omicsai_instructions(paste0("format_", style))

      ## Cache the prompt for instant toggle
      prompt <- omicsai::omicsai_substitute_template(ai_summary_template, params)
      summary_prompt_cache(paste0("# Summary Prompt\n\n", prompt))

      result <- omicsai::omicsai_gen_text(
        template = ai_summary_template,
        params = params,
        config = omicsai::omicsai_config(model = model %||% Sys.getenv("OMICS_AI_MODEL", "ollama:llama3.2"))
      )
      result$text
    }, ignoreNULL = FALSE)

    # ---- REPORT MODE: single-shot integrated report ----

    ai_report <- shiny::eventReactive(controls$trigger(), {
      mode <- controls$mode()
      if (mode != "report" || controls$trigger() < 1) return(NULL)

      w <- wgcna()
      shiny::req(w)
      model <- ai_model()
      shiny::req(model)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## Step 1: Build structured data tables
      progress$set(message = "Extracting module data...", value = 0.1)
      tables <- wgcna_build_report_tables(w, pgx)

      ## Step 2: Classify module signal strength
      progress$set(message = "Classifying modules...", value = 0.2)
      ranking <- wgcna_rank_modules(w)

      ## Step 3: Load board rules
      board_rules <- paste(readLines(
        file.path(BOARD_PROMPTS_DIR, "wgcna_report_rules.md"),
        warn = FALSE
      ), collapse = "\n")

      ## Step 4: Species context (graceful: empty string if unknown)
      organism <- pgx$organism %||% NULL
      species_text <- omicsai::omicsai_species_prompt(organism)

      ## Step 5: Assemble user message (board owns all domain logic)
      user_message <- omicsai::collapse_lines(
        board_rules,
        species_text,
        "---",
        "## Module Signal Classification",
        omicsai::omicsai_format_ranking(ranking),
        "---",
        "## Input Data",
        tables$text,
        sep = "\n\n"
      )

      ## Cache the full prompt for instant toggle
      sys_prompt <- tryCatch({
        fp <- omicsai::omicsai_prompt_path("report_format.md")
        txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
        omicsai::omicsai_substitute_template(txt, list(max_words = "1500"))
      }, error = function(e) "(report_format.md not found)")

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", sys_prompt,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", user_message
      ))

      ## Step 6: Generate report via single LLM call
      progress$set(message = "Generating report...", value = 0.3)

      result <- omicsai::omicsai_gen_report(
        template = user_message,
        model = model
      )

      ## Step 7: Append deterministic methods section
      methods <- wgcna_build_methods(w, pgx)
      full_report <- paste(result$text, methods, sep = "\n\n")

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

    list(
      text = ai_text,
      report_text = ai_report
    )
  })
}
