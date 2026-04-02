##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_ai_text_server <- function(id, mwgcna, pgx, controls, parent_session) {
  moduleServer(id, function(input, output, session) {
    summary_data_path <- file.path(MULTIWGCNA_PROMPTS_DIR, "multiwgcna_summary_data.md")
    interpretation_path <- file.path(MULTIWGCNA_PROMPTS_DIR, "multiwgcna_interpretation.md")

    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)
    report_data_tables <- shiny::reactiveVal(NULL)

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        w <- mwgcna()
        module_ref <- controls$selected_module()
        shiny::req(w, module_ref)

        params <- multiwgcna_build_summary_params(w, module_ref, pgx)
        shiny::req(params)
        model <- get_ai_model(parent_session)
        style <- controls$summary_style() %||% "short_summary"
        organism <- pgx$organism %||% NULL

        p <- omicsai::summary_prompt(
          role = omicsai::frag("system_base"),
          task = omicsai::frag(paste0("text/", style)),
          species = omicsai::omicsai_species_prompt(organism),
          context = omicsai::frag(interpretation_path),
          data = omicsai::frag(summary_data_path, params)
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
            shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
          }
        )
        result$text
      },
      ignoreNULL = FALSE
    )

    ai_report <- shiny::eventReactive(controls$trigger(), {
      if (controls$mode() != "report" || controls$trigger() < 1) return(NULL)

      w <- mwgcna()
      shiny::req(w)
      model <- get_ai_model(parent_session)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      progress$set(message = "Extracting module data...", value = 0.1)
      tables <- multiwgcna_build_report_tables(w, pgx)
      report_data_tables(tables$text)

      progress$set(message = "Classifying modules...", value = 0.15)
      ranking <- multiwgcna_rank_modules(w, pgx)

      data_content <- omicsai::collapse_lines(
        "## Module Signal Classification",
        omicsai::omicsai_format_ranking(ranking),
        "---",
        "## Input Data",
        tables$text,
        sep = "\n\n"
      )

      board_rules_path <- file.path(MULTIWGCNA_PROMPTS_DIR, "multiwgcna_report_rules.md")
      organism <- pgx$organism %||% NULL
      methods_text <- multiwgcna_build_methods(w, pgx)

      rp <- omicsai::report_prompt(
        role = omicsai::frag("system_base"),
        task = omicsai::frag("text/report"),
        species = omicsai::omicsai_species_prompt(organism),
        context = omicsai::frag(interpretation_path),
        board_rules = omicsai::frag(board_rules_path),
        data = data_content
      )
      bp <- omicsai::build_prompt(rp)

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      progress$set(message = "Generating AI report...", value = 0.25)
      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system, max_tokens = 8192L)
      cache <- omicsai::omicsai_cache_init("mem")
      result <- tryCatch(
        omicsai::omicsai_gen_text(bp$board, config = cfg, cache = cache),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, .aicards_friendly_error(conditionMessage(e))))
        }
      )

      if (!grepl("^# ", result$text)) {
        result$text <- paste0("# Analysis Report\n\n", result$text)
      }

      progress$set(message = "Done!", value = 1)
      paste(result$text, methods_text, sep = "\n\n")
    }, ignoreNULL = FALSE)

    ai_text <- shiny::reactive({
      show <- isTRUE(controls$show_prompt())
      if (controls$mode() == "summary") {
        if (show) summary_prompt_cache() else ai_summary()
      } else {
        if (show) report_prompt_cache() else ai_report()
      }
    })

    list(
      text = ai_text,
      report_text = ai_report,
      report_data_tables = report_data_tables
    )
  })
}
