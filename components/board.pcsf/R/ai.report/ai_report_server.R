##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

pcsf_ai_report_server <- function(id,
                                  pgx,
                                  contrast_reactive,
                                  current_graph_reactive,
                                  current_table_reactive,
                                  pcsf_params_reactive,
                                  parent_session,
                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    cache <- omicsai::omicsai_cache_init("mem")

    summary_choices <- shiny::reactive({
      choices <- pcsf_get_contrasts(pgx)
      fallback <- contrast_reactive()
      if (!is.null(fallback) && nzchar(fallback) && !fallback %in% choices) {
        choices <- c(fallback, choices)
      }
      choices
    })

    controls <- pcsf_ai_report_controls_server(
      "controls",
      module_choices = summary_choices
    )

    text_result <- pcsf_ai_text_server(
      "text",
      pgx = pgx,
      controls = controls,
      contrast_reactive = contrast_reactive,
      current_graph_reactive = current_graph_reactive,
      current_table_reactive = current_table_reactive,
      pcsf_params_reactive = pcsf_params_reactive,
      parent_session = parent_session
    )

    diagram_result <- AiDiagramCardServer(
      "layout-diagram",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.pcsf")
        prompt <- pcsf_build_diagram_prompt(txt, organism, board_root)
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
      style = pcsf_diagram_style()
    )

    AiImageCardServer(
      "layout-infographic",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        diag <- diagram_result()
        edgelist <- if (!is.null(diag)) diag$edgelist else NULL
        prompt <- pcsf_build_image_prompt(txt, organism, edgelist)
        list(content = prompt)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        img_model <- getUserOption(parent_session, "image_model")
        omicsai::omicsai_image_config(
          model = img_model %||% "gemini-3.1-flash-image-preview",
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

    pcsf_ai_report_layout_server(
      "layout",
      text_reactive = text_result$text,
      watermark = watermark
    )
  })
}
