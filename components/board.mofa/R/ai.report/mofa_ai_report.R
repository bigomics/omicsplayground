##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

mofa_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)

  show_prompt_input <- shiny::checkboxInput(
    shiny::NS(ns("controls"), "show_prompt"),
    "Show prompt",
    FALSE
  )

  multiomics_ai_report_layout_ui(
    ns("layout"),
    text_title = "AI MOFA Report",
    diagram_title = "Factor Mechanism Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input
  )
}

mofa_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  multiomics_ai_report_controls_ui(ns("controls"), module_label = "Factor:")
}

mofa_diagram_style <- function() {
  list(
    node_styles = list(
      factor = list(bg = "#E0F2FE", border = "#1E3A8A", shape = "ellipse"),
      trait = list(bg = "#DCFCE7", border = "#166534", shape = "box"),
      pathway = list(bg = "#FCE7F3", border = "#9D174D", shape = "diamond"),
      feature = list(bg = "#FEF3C7", border = "#92400E", shape = "box")
    ),
    edge_styles = list(
      associates = list(color = "#2563EB", dashes = FALSE, arrows = "arrow", arrow_on = TRUE, width = 2),
      enriches = list(color = "#BE185D", dashes = FALSE, arrows = "arrow", arrow_on = TRUE, width = 2),
      loads = list(color = "#92400E", dashes = TRUE, arrows = "arrow", arrow_on = TRUE, width = 1.5)
    ),
    show_genes_in_label = FALSE
  )
}

mofa_build_diagram_prompt <- function(report_text, organism, board_root) {
  style <- mofa_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")

  layers <- list()
  layers[[1]] <- omicsai::omicsai_instructions("diagram_network")
  layers[[2]] <- tryCatch({
    tpl <- omicsai::omicsai_load_template("prompts/mofa/diagram_mofa_rules.md", root = board_root)
    omicsai::omicsai_substitute_template(
      tpl,
      list(
        node_names = fmt_names(names(style$node_styles)),
        link_names = fmt_names(names(style$edge_styles))
      )
    )
  }, error = function(e) "")
  layers[[3]] <- tryCatch(omicsai::omicsai_species_prompt(organism), error = function(e) "")
  layers[[4]] <- paste("## AI Report\n\n", report_text)
  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n---\n\n")
}

mofa_build_image_prompt <- function(report_text, organism, diagram_edgelist = NULL) {
  layers <- list()
  layers[[1]] <- omicsai::omicsai_image_species_visual(organism)
  layers[[2]] <- paste(
    "<report>",
    omicsai::omicsai_strip_report_noise(report_text),
    "</report>",
    sep = "\n"
  )

  if (!is.null(diagram_edgelist)) {
    edges <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(edges)) layers[[3]] <- paste("<diagram>", edges, "</diagram>", sep = "\n")
  }

  layers <- layers[nzchar(layers)]
  paste(layers, collapse = "\n\n")
}

mofa_ai_text_server <- function(id, mofa_reactive, pgx, controls, parent_session) {
  moduleServer(id, function(input, output, session) {
    summary_template <- omicsai::omicsai_load_template(
      file.path(MOFA_PROMPTS_DIR, "mofa_ai_report_summary_template.md")
    )

    report_rules <- omicsai::omicsai_load_template(
      file.path(MOFA_PROMPTS_DIR, "mofa_report_rules.md")
    )

    methods_context <- omicsai::omicsai_load_template(
      file.path(MOFA_PROMPTS_DIR, "MOFA_methods.md")
    )

    ai_model <- shiny::reactive({
      m <- getUserOption(parent_session, "llm_model")
      shiny::req(m, m != "")
      m
    })

    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- ai_model()
        m <- mofa_reactive()
        factor_name <- controls$selected_module()
        shiny::req(model, m, factor_name)

        params <- mofa_ai_build_summary_params(m, pgx, factor_name, ntop = 12L)
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
      m <- mofa_reactive()
      shiny::req(model, m)

      tables <- mofa_ai_build_report_tables(m, pgx, max_contexts = 8L, ntop = 10L)

      sys_prompt <- tryCatch({
        fp <- omicsai::omicsai_prompt_path("report_format.md")
        txt <- paste(readLines(fp, warn = FALSE), collapse = "\n")
        omicsai::omicsai_substitute_template(txt, list(max_words = "1500"))
      }, error = function(e) "(report_format.md not found)")

      user_message <- omicsai::collapse_lines(
        omicsai::omicsai_substitute_template(
          methods_context,
          list(experiment = m$experiment %||% pgx$name %||% pgx$description %||% "omics experiment")
        ),
        report_rules,
        omicsai::omicsai_species_prompt(pgx$organism %||% NULL),
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

      methods <- mofa_ai_build_methods(m, pgx)
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

mofa_ai_report_server <- function(id, mofa_reactive, pgx, parent_session, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    cache <- omicsai::omicsai_cache_init("mem")

    factor_choices <- shiny::reactive({
      m <- mofa_reactive()
      shiny::req(m)
      mofa_ai_factor_choices(m)
    })

    controls <- multiomics_ai_report_controls_server(
      "controls",
      module_choices = factor_choices
    )

    text_result <- mofa_ai_text_server(
      "text",
      mofa_reactive = mofa_reactive,
      pgx = pgx,
      controls = controls,
      parent_session = parent_session
    )

    diagram_result <- AiDiagramCardServer(
      "layout-diagram",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        prompt <- mofa_build_diagram_prompt(
          report_text = txt,
          organism = pgx$organism %||% "human",
          board_root = file.path(OPG, "components/board.mofa")
        )
        list(content = prompt)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        llm <- getUserOption(parent_session, "llm_model")
        omicsai::omicsai_diagram_config(
          model = llm %||% "ollama:llama3.2",
          default_regulation = "association"
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (controls$mode() == "report") controls$trigger() else 0
      }),
      style = mofa_diagram_style()
    )

    AiImageCardServer(
      "layout-infographic",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        diag <- diagram_result()
        edgelist <- if (!is.null(diag)) diag$edgelist else NULL
        prompt <- mofa_build_image_prompt(
          report_text = txt,
          organism = pgx$organism %||% "human",
          diagram_edgelist = edgelist
        )
        list(content = prompt)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        img_model <- getUserOption(parent_session, "image_model")
        omicsai::omicsai_image_config(
          model = img_model %||% "gemini-3-pro-image-preview",
          style = controls$image_style() %||% "bigomics",
          n_blocks = as.integer(controls$image_blocks() %||% 1L),
          image_size = "1K"
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (controls$mode() == "report" && isTRUE(controls$include_infographic())) {
          controls$trigger()
        } else {
          0
        }
      }),
      watermark = watermark
    )

    multiomics_ai_report_layout_server(
      "layout",
      text_reactive = text_result$text,
      filename = "mofa-ai-report",
      title = "MOFA AI Report",
      watermark = watermark
    )
  })
}
