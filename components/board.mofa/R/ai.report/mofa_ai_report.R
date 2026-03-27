##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

mofa_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)
  controls_ns <- shiny::NS(ns("controls"))

  show_prompt_input <- shiny::checkboxInput(
    controls_ns("show_prompt"),
    "Show prompt",
    FALSE
  )

  ## Infographic style/layout controls: namespaced to controls so the controls
  ## server can read them, but rendered in the image card's hamburger menu.
  ## The include_infographic checkbox lives in the sidebar (ai_report_shared.R).
  infographic_options <- shiny::tagList(
    shiny::tags$hr(),
    shiny::selectInput(
      controls_ns("image_style"),
      "Infographic Style:",
      choices = NULL,
      width = "100%"
    ),
    shiny::radioButtons(
      controls_ns("image_blocks"),
      "Layout:",
      choices = c("1 Panel" = "1", "2 Panels" = "2", "3 Panels" = "3"),
      selected = "1",
      inline = TRUE
    )
  )

  multiomics_ai_report_layout_ui(
    ns("layout"),
    text_title = "AI MOFA Report",
    diagram_title = "Factor Mechanism Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input,
    infographic_options = infographic_options
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
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  rules_path <- file.path(board_root, "prompts/mofa/mofa_diagram_rules.md")

  p <- omicsai::diagram_prompt(
    role        = omicsai::frag("system_base"),
    task        = omicsai::frag("diagram/network"),
    species     = omicsai::omicsai_species_prompt(organism),
    board_rules = omicsai::frag(rules_path, list(node_names = node_names, link_names = link_names)),
    report      = paste("## AI Report\n\n", report_text)
  )
  omicsai::build_prompt(p)
}

mofa_build_image_prompt <- function(report_text, organism,
                                    diagram_edgelist = NULL,
                                    style_name = NULL,
                                    n_blocks = 1L) {
  species_img <- omicsai::omicsai_image_species_visual(organism)
  clean <- omicsai::omicsai_strip_report_noise(report_text)

  edge_text <- NULL
  if (!is.null(diagram_edgelist)) {
    dot <- omicsai::omicsai_edgelist_to_text(diagram_edgelist)
    if (nzchar(dot)) edge_text <- dot
  }

  style_frag  <- if (!is.null(style_name) && nzchar(style_name)) {
    omicsai::frag(paste0("image/styles/", style_name))
  }
  blocks_frag <- omicsai::frag(paste0("image/blocks/blocks_", as.integer(n_blocks)))

  p <- omicsai::image_prompt(
    role             = omicsai::frag("system_base"),
    task             = omicsai::frag("image/infographic", params = list(board_name = "MOFA")),
    species          = species_img,
    style            = style_frag,
    blocks           = blocks_frag,
    report           = clean,
    diagram_edgelist = edge_text
  )
  omicsai::build_prompt(p)
}

mofa_ai_text_server <- function(id, mofa_reactive, pgx, controls, parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: prompt template paths ----

    summary_data_path   <- file.path(MOFA_PROMPTS_DIR, "mofa_summary_data.md")
    interpretation_path <- file.path(MOFA_PROMPTS_DIR, "mofa_interpretation.md")

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-factor summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- get_ai_model(parent_session)
        m <- mofa_reactive()
        factor_name <- controls$selected_module()
        shiny::req(m, factor_name)

        params <- mofa_ai_build_summary_params(m, pgx, factor_name, ntop = 12L)
        shiny::req(params)

        organism <- pgx$organism %||% NULL
        style <- controls$summary_style() %||% "short_summary"

        p <- omicsai::summary_prompt(
          role    = omicsai::frag("system_base"),
          task    = omicsai::frag(paste0("text/", style)),
          species = omicsai::omicsai_species_prompt(organism),
          context = omicsai::frag(interpretation_path, list(experiment = params$experiment)),
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
      m <- mofa_reactive()
      shiny::req(m)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## Step 1: Build structured data tables
      progress$set(message = "Extracting factor data...", value = 0.1)
      tables <- mofa_ai_build_report_tables(m, pgx, max_contexts = 8L, ntop = 10L)

      ## Step 2: Build template parameters from structured data tables
      report_data_path <- file.path(MOFA_PROMPTS_DIR, "mofa_report_data.md")
      data_tables <- tables[c(
        "experiment", "organism", "n_factors",
        "factor_ranking", "factor_details"
      )]

      ## Step 3: Build structured prompt
      board_rules_path <- file.path(MOFA_PROMPTS_DIR, "mofa_report_rules.md")
      organism <- pgx$organism %||% NULL

      ## Pre-render methods/context (reused in Step 5 as deterministic appendix)
      methods_context <- mofa_ai_build_methods(m, pgx)

      interpretation_path <- file.path(MOFA_PROMPTS_DIR, "mofa_interpretation.md")
      experiment_label <- mofa_reactive()$experiment %||% pgx$name %||% pgx$description %||% "omics experiment"

      p <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path, list(experiment = experiment_label)),
        board_rules = omicsai::frag(board_rules_path),
        data        = omicsai::frag(report_data_path, data_tables)
      )
      bp <- omicsai::build_prompt(p)

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      ## Step 4: Generate report via single LLM call
      ## Use omicsai_gen_text (NOT gen_report — that would double-load text/report.md)
      progress$set(message = "Generating report...", value = 0.3)

      cfg <- omicsai::omicsai_config(model = model, system_prompt = bp$system, max_tokens = 8192L)
      cache <- omicsai::omicsai_cache_init("mem")
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
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.mofa")
        bp <- mofa_build_diagram_prompt(txt, organism, board_root)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.mofa")
        bp <- mofa_build_diagram_prompt(txt, organism, board_root)
        llm <- get_ai_model(parent_session)
        make_llm_diagram_config(llm,
          system_prompt = bp$system,
          default_regulation = "associates",
          node_styles = mofa_diagram_style()$node_styles,
          edge_styles = mofa_diagram_style()$edge_styles
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
        organism <- pgx$organism %||% "human"
        img_style <- controls$image_style() %||% "bigomics"
        img_blocks <- as.integer(controls$image_blocks() %||% 1L)
        bp <- mofa_build_image_prompt(txt, organism, edgelist,
                                      style_name = img_style, n_blocks = img_blocks)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        diag <- diagram_result()
        edgelist <- if (!is.null(diag)) diag$edgelist else NULL
        organism <- pgx$organism %||% "human"
        img_style <- controls$image_style() %||% "bigomics"
        img_blocks <- as.integer(controls$image_blocks() %||% 1L)
        bp <- mofa_build_image_prompt(txt, organism, edgelist,
                                      style_name = img_style, n_blocks = img_blocks)
        img_model <- getUserOption(parent_session, "image_model")
        shiny::validate(shiny::need(
          !is.null(img_model) && nzchar(img_model),
          "No image model configured. Please select a model in Settings."
        ))
        omicsai::omicsai_image_config(
          model = img_model,
          system_prompt = bp$system,
          style = img_style,
          n_blocks = img_blocks,
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
