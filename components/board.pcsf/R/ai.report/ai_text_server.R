##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' PCSF AI Text Server
pcsf_ai_text_server <- function(id,
                                pgx,
                                controls,
                                contrast_reactive,
                                current_graph_reactive,
                                current_table_reactive,
                                pcsf_params_reactive,
                                parent_session) {
  moduleServer(id, function(input, output, session) {
    summary_template <- omicsai::omicsai_load_template(
      file.path(PCSF_PROMPTS_DIR, "pcsf_ai_report_summary_template.md")
    )

    report_rules <- omicsai::omicsai_load_template(
      file.path(PCSF_PROMPTS_DIR, "pcsf_report_rules.md")
    )

    methods_context <- omicsai::omicsai_load_template(
      file.path(PCSF_PROMPTS_DIR, "PCSF_methods.md")
    )

    ai_model <- shiny::reactive({
      m <- getUserOption(parent_session, "llm_model")
      shiny::req(m, m != "")
      m
    })

    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module(), contrast_reactive()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- ai_model()
        fallback_contrast <- contrast_reactive()
        contrast <- controls$selected_module() %||% fallback_contrast
        shiny::req(model, contrast)

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

        params$methods_context <- omicsai::omicsai_substitute_template(
          methods_context,
          list(experiment = params$experiment)
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
      shiny::req(model)

      tables <- pcsf_build_report_tables(
        pgx = pgx,
        contrasts = pcsf_get_contrasts(pgx),
        max_contexts = 6L,
        ntop = 10L,
        ntop_network = 750L,
        pcsf_params = pcsf_params_reactive()
      )

      user_message <- omicsai::collapse_lines(
        omicsai::omicsai_substitute_template(
          methods_context,
          list(experiment = pgx$name %||% pgx$description %||% "omics experiment")
        ),
        report_rules,
        omicsai::omicsai_species_prompt(pgx$organism %||% NULL),
        "---",
        tables$text,
        sep = "\n\n"
      )

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

      result <- omicsai::omicsai_gen_report(
        template = user_message,
        model = model
      )

      methods <- pcsf_build_methods(pgx = pgx)
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
