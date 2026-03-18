##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

drugconnectivity_ai_report_server <- function(id,
                                              pgx,
                                              method_reactive,
                                              contrast_reactive,
                                              annotated_only_reactive,
                                              parent_session,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    cache <- omicsai::omicsai_cache_init("mem")

    summary_choices <- shiny::reactive({
      method <- method_reactive()
      shiny::req(method)
      dr <- dc_get_method_data(pgx, method)
      choices <- dc_get_contrasts(dr)
      fallback <- contrast_reactive()
      if (!is.null(fallback) && nzchar(fallback) && !fallback %in% choices) {
        choices <- c(fallback, choices)
      }
      choices
    })

    controls <- drugconnectivity_ai_report_controls_server(
      "controls",
      module_choices = summary_choices
    )

    text_result <- drugconnectivity_ai_text_server(
      "text",
      pgx = pgx,
      controls = controls,
      method_reactive = method_reactive,
      contrast_reactive = contrast_reactive,
      annotated_only_reactive = annotated_only_reactive,
      parent_session = parent_session
    )

    diagram_result <- AiDiagramCardServer(
      "layout-diagram",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.drugconnectivity")
        data_tables <- text_result$report_data_tables()
        bp <- drugconnectivity_build_diagram_prompt(txt, organism, board_root, data_tables = data_tables)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.drugconnectivity")
        data_tables <- text_result$report_data_tables()
        bp <- drugconnectivity_build_diagram_prompt(txt, organism, board_root, data_tables = data_tables)
        llm <- get_ai_model(parent_session)
        make_llm_diagram_config(llm,
          system_prompt = bp$system,
          default_regulation = "association",
          node_styles = drugconnectivity_diagram_style()$node_styles,
          edge_styles = drugconnectivity_diagram_style()$edge_styles
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (controls$mode() == "report") controls$trigger() else 0
      }),
      style = drugconnectivity_diagram_style()
    )

    AiImageCardServer(
      "layout-infographic",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        diag <- diagram_result()
        edgelist <- if (!is.null(diag)) diag$edgelist else NULL
        img_style <- controls$image_style() %||% "bigomics"
        img_blocks <- as.integer(controls$image_blocks() %||% 1L)
        bp <- drugconnectivity_build_image_prompt(txt, organism, edgelist,
                                                  style_name = img_style, n_blocks = img_blocks)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        diag <- diagram_result()
        edgelist <- if (!is.null(diag)) diag$edgelist else NULL
        img_style <- controls$image_style() %||% "bigomics"
        img_blocks <- as.integer(controls$image_blocks() %||% 1L)
        bp <- drugconnectivity_build_image_prompt(txt, organism, edgelist,
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

    drugconnectivity_ai_report_layout_server(
      "layout",
      text_reactive = text_result$text,
      watermark = watermark
    )
  })
}
