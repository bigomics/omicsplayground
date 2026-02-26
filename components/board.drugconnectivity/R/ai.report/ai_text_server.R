##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Drug Connectivity AI Text Server
drugconnectivity_ai_text_server <- function(id,
                                            pgx,
                                            controls,
                                            method_reactive,
                                            contrast_reactive,
                                            annotated_only_reactive,
                                            parent_session) {
  moduleServer(id, function(input, output, session) {
    summary_template <- omicsai::omicsai_load_template(
      file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_ai_report_summary_template.md")
    )

    report_rules <- omicsai::omicsai_load_template(
      file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "drugconnectivity_report_rules.md")
    )

    methods_context <- omicsai::omicsai_load_template(
      file.path(DRUGCONNECTIVITY_PROMPTS_DIR, "DrugConnectivity_methods.md")
    )

    ai_model <- shiny::reactive({
      m <- getUserOption(parent_session, "llm_model")
      shiny::req(m, m != "")
      m
    })

    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module(), method_reactive(), contrast_reactive()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)
        model <- ai_model()
        method <- method_reactive()
        fallback_contrast <- contrast_reactive()
        contrast <- controls$selected_module() %||% fallback_contrast
        only_annotated <- isTRUE(annotated_only_reactive())
        shiny::req(model, method, contrast)
        at <- dc_analysis_type_info(method)

        params <- drugconnectivity_build_summary_params(
          pgx = pgx,
          contrast = contrast,
          method = method,
          only_annotated = only_annotated,
          ntop = 12L
        )
        shiny::req(params)

        params$methods_context <- omicsai::omicsai_substitute_template(
          methods_context,
          list(
            experiment = params$experiment,
            analysis_type = at$analysis_type,
            analysis_type_description = at$analysis_type_description
          )
        )
        params$style_instructions <- omicsai::omicsai_instructions(
          paste0("format_", controls$summary_style() %||% "short")
        )

        prompt <- omicsai::omicsai_substitute_template(summary_template, params)
        summary_prompt_cache(paste0("# Summary Prompt\n\n", prompt))

        out <- omicsai::omicsai_gen_text(
          template = summary_template,
          params = params,
          config = omicsai::omicsai_config(model = model)
        )
        out$text
      },
      ignoreNULL = FALSE
    )

    ai_report <- shiny::eventReactive(controls$trigger(), {
      if (controls$mode() != "report" || controls$trigger() < 1) return(NULL)

      model <- ai_model()
      method <- method_reactive()
      only_annotated <- isTRUE(annotated_only_reactive())
      shiny::req(model, method)
      at <- dc_analysis_type_info(method)

      tables <- drugconnectivity_build_report_tables(
        pgx = pgx,
        method = method,
        only_annotated = only_annotated,
        max_contexts = 8L,
        ntop = 10L
      )

      species_text <- omicsai::omicsai_species_prompt(pgx$organism %||% NULL)

      sys_prompt <- tryCatch({
        fp <- omicsai::omicsai_prompt_path("report_format.md")
        txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
        omicsai::omicsai_substitute_template(txt, list(max_words = "1500"))
      }, error = function(e) "(report_format.md not found)")

      analysis_type_block <- omicsai::omicsai_substitute_template(
        "ANALYSIS TYPE: {{analysis_type}}\n{{analysis_type_description}}",
        list(
          analysis_type = at$analysis_type,
          analysis_type_description = at$analysis_type_description
        )
      )

      user_message <- omicsai::collapse_lines(
        omicsai::omicsai_substitute_template(
          methods_context,
          list(
            experiment = pgx$name %||% pgx$description %||% "omics experiment",
            analysis_type = at$analysis_type,
            analysis_type_description = at$analysis_type_description
          )
        ),
        analysis_type_block,
        report_rules,
        "---",
        species_text,
        "---",
        tables$text,
        sep = "\n\n"
      )

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", sys_prompt,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", user_message
      ))

      report_input <- omicsai::collapse_lines(
        sys_prompt,
        "---",
        user_message,
        sep = "\n\n"
      )

      result <- omicsai::omicsai_gen_report(
        template = report_input,
        model = model,
        max_words = 1500L
      )

      methods <- drugconnectivity_build_methods(pgx = pgx, method = method)
      paste(result$text, methods, sep = "\n\n")
    }, ignoreNULL = FALSE)

    ai_text <- shiny::reactive({
      show <- isTRUE(controls$show_prompt())
      if (controls$mode() == "summary") {
        if (show) summary_prompt_cache() else ai_summary()
      } else {
        if (show) report_prompt_cache() else ai_report()
      }
    })

    list(text = ai_text, report_text = ai_report)
  })
}
