##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' WGCNA AI Report Server (Coordinator)
#'
#' Thin coordinator that wires controls, text, diagram, infographic,
#' and layout servers together. This is the single entry point called
#' from wgcna_server.R.
#'
#' @param id Module namespace ID
#' @param wgcna Reactive returning WGCNA results object
#' @param pgx PGX object (non-reactive)
#' @param parent_session Parent Shiny session (for getUserOption)
#' @param watermark Logical; add watermark to outputs
#'
#' @return NULL
wgcna_ai_report_server <- function(id, wgcna, pgx, parent_session, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    cache <- omicsai::omicsai_cache_init("mem")

    # Module choices for summary mode selector
    ai_module_choices <- shiny::reactive({
      w <- wgcna()
      shiny::req(w)
      setdiff(names(sort(lengths(w$me.genes), decreasing = TRUE)), "grey")
    })

    # Controls
    controls <- ai_report_controls_server("controls", module_choices = ai_module_choices)

    # Card servers
    text_result <- wgcna_ai_text_server("text", wgcna, pgx, controls, parent_session)
    diagram_result <- AiDiagramCardServer(
      "layout-diagram",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.wgcna")
        prompt <- wgcna_build_diagram_prompt(txt, organism, board_root)
        list(content = prompt)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        llm <- getUserOption(parent_session, "llm_model")
        omicsai::omicsai_diagram_config(
          model = llm %||% "ollama:llama3.2",
          default_regulation = "positive"
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (controls$mode() == "report") controls$trigger() else 0
      }),
      style = wgcna_diagram_style()
    )

    AiImageCardServer(
      "layout-infographic",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        diag <- diagram_result()
        edgelist <- if (!is.null(diag)) diag$edgelist else NULL
        prompt <- wgcna_build_image_prompt(txt, organism, edgelist)
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

    # Layout rendering
    ai_report_layout_server(
      "layout",
      text_reactive = text_result$text,
      watermark = watermark
    )
  })
}
