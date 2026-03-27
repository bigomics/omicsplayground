##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

lasagna_ai_report_ui <- function(id) {
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
    text_title = "AI LASAGNA Report",
    diagram_title = "Layer-Interaction Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input,
    infographic_options = infographic_options
  )
}

lasagna_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  multiomics_ai_report_controls_ui(ns("controls"), module_label = "Contrast:")
}

lasagna_diagram_style <- function() {
  list(
    node_styles = list(
      contrast = list(bg = "#E0F2FE", border = "#1E3A8A", shape = "ellipse"),
      layer = list(bg = "#DCFCE7", border = "#166534", shape = "box"),
      feature = list(bg = "#FEF3C7", border = "#92400E", shape = "box")
    ),
    edge_styles = list(
      cross_layer = list(color = "#2563EB", dashes = FALSE, arrows = "arrow", arrow_on = TRUE, width = 2),
      intra_layer = list(color = "#B45309", dashes = TRUE, arrows = "arrow", arrow_on = FALSE, width = 1.5),
      association = list(color = "#0F766E", dashes = FALSE, arrows = "arrow", arrow_on = TRUE, width = 2)
    ),
    show_genes_in_label = FALSE
  )
}

lasagna_build_diagram_prompt <- function(report_text, organism, board_root) {
  style <- lasagna_diagram_style()
  fmt_names <- function(nms) paste(sprintf("- `%s`", nms), collapse = "\n")
  node_names <- fmt_names(names(style$node_styles))
  link_names <- fmt_names(names(style$edge_styles))

  rules_path <- file.path(board_root, "prompts/lasagna/lasagna_diagram_rules.md")

  p <- omicsai::diagram_prompt(
    role        = omicsai::frag("system_base"),
    task        = omicsai::frag("diagram/network"),
    species     = omicsai::omicsai_species_prompt(organism),
    board_rules = omicsai::frag(rules_path, list(node_names = node_names, link_names = link_names)),
    report      = paste("## AI Report\n\n", report_text)
  )
  omicsai::build_prompt(p)
}

lasagna_build_image_prompt <- function(report_text, organism,
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
    task             = omicsai::frag("image/infographic", params = list(board_name = "LASAGNA")),
    species          = species_img,
    style            = style_frag,
    blocks           = blocks_frag,
    report           = clean,
    diagram_edgelist = edge_text
  )
  omicsai::build_prompt(p)
}

lasagna_ai_text_server <- function(id,
                                   graph_data_reactive,
                                   contrast_reactive,
                                   pgx,
                                   controls,
                                   parent_session) {
  moduleServer(id, function(input, output, session) {

    # ---- Shared: prompt template paths ----

    summary_data_path   <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_summary_data.md")
    interpretation_path <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_interpretation.md")

    # ---- Prompt caches (for instant toggle without regeneration) ----
    summary_prompt_cache <- shiny::reactiveVal(NULL)
    report_prompt_cache <- shiny::reactiveVal(NULL)

    # ---- SUMMARY MODE: per-contrast summary ----

    ai_summary <- shiny::eventReactive(
      list(controls$trigger(), controls$selected_module(), contrast_reactive()),
      {
        if (controls$mode() != "summary" || controls$trigger() < 1) return(NULL)

        model <- get_ai_model(parent_session)
        contrast <- controls$selected_module() %||% contrast_reactive()
        res <- graph_data_reactive()
        shiny::req(contrast, res)

        params <- lasagna_ai_build_summary_params(res, contrast, pgx, ntop = 12L)
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
      contrast <- contrast_reactive()
      res <- graph_data_reactive()
      shiny::req(contrast, res)

      progress <- shiny::Progress$new()
      on.exit(progress$close())

      ## Step 1: Build structured data tables
      progress$set(message = "Extracting contrast data...", value = 0.1)
      tables <- lasagna_ai_build_report_tables(res, contrast, pgx, ntop = 12L)

      ## Step 2: Build template parameters from structured data tables
      report_data_path <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_report_data.md")
      data_tables <- tables[c(
        "experiment", "contrast", "network_snapshot",
        "layer_participation", "top_nodes_table",
        "top_edges_table", "caveats"
      )]

      ## Step 3: Build structured prompt
      board_rules_path <- file.path(LASAGNA_PROMPTS_DIR, "lasagna_report_rules.md")
      organism <- pgx$organism %||% NULL

      ## Pre-render methods/context (reused as deterministic appendix)
      methods_context <- lasagna_ai_build_methods(pgx, contrast)

      p <- omicsai::report_prompt(
        role        = omicsai::frag("system_base"),
        task        = omicsai::frag("text/report"),
        species     = omicsai::omicsai_species_prompt(organism),
        context     = omicsai::frag(interpretation_path),
        board_rules = omicsai::frag(board_rules_path),
        data        = omicsai::frag(report_data_path, data_tables)
      )
      bp <- omicsai::build_prompt(p)

      report_prompt_cache(paste0(
        "# SYSTEM PROMPT\n\n", bp$system,
        "\n\n---\n\n",
        "# USER MESSAGE\n\n", bp$board
      ))

      ## Step 3: Generate report via single LLM call
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

      ## Step 4: Append deterministic methods section
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

lasagna_ai_report_server <- function(id,
                                     graph_data_reactive,
                                     contrast_reactive,
                                     contrast_choices_reactive,
                                     pgx,
                                     parent_session,
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    cache <- omicsai::omicsai_cache_init("mem")

    controls <- multiomics_ai_report_controls_server(
      "controls",
      module_choices = contrast_choices_reactive
    )

    text_result <- lasagna_ai_text_server(
      "text",
      graph_data_reactive = graph_data_reactive,
      contrast_reactive = contrast_reactive,
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
        bp <- lasagna_build_diagram_prompt(txt, organism, board_root)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.mofa")
        bp <- lasagna_build_diagram_prompt(txt, organism, board_root)
        llm <- get_ai_model(parent_session)
        make_llm_diagram_config(llm,
          system_prompt = bp$system,
          default_regulation = "cross_layer",
          node_styles = lasagna_diagram_style()$node_styles,
          edge_styles = lasagna_diagram_style()$edge_styles
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (controls$mode() == "report") controls$trigger() else 0
      }),
      style = lasagna_diagram_style()
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
        bp <- lasagna_build_image_prompt(txt, organism, edgelist,
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
        bp <- lasagna_build_image_prompt(txt, organism, edgelist,
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
      filename = "lasagna-ai-report",
      title = "LASAGNA AI Report",
      watermark = watermark
    )
  })
}
