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

    # ── Coordinator-owned progress bar (spans text + diagram) ──
    report_progress <- shiny::reactiveVal(NULL)

    shiny::observeEvent(controls$trigger(), {
      if (controls$mode() != "report" || controls$trigger() < 1) return()
      old <- report_progress()
      if (!is.null(old)) try(old$close(), silent = TRUE)
      p <- shiny::Progress$new(session)
      p$set(message = "Starting report generation...", value = 0)
      report_progress(p)
      message(sprintf("[INFO][%s] --- [AI-REPORT] starting report generation...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    })

    # Card servers
    text_result <- wgcna_ai_text_server(
      "text", wgcna, pgx, controls, parent_session,
      progress_reactive = report_progress
    )

    # Update progress when text completes → diagram starts
    shiny::observeEvent(text_result$report_text(), {
      txt <- text_result$report_text()
      shiny::req(txt)
      p <- report_progress()
      message(sprintf("[INFO][%s] --- [AI-REPORT] generating network diagram...", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) p$set(message = "Generating network diagram...", value = 0.7)
    })

    diagram_result <- AiDiagramCardServer(
      "layout-diagram",
      params_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.wgcna")
        data_tables <- text_result$report_data_tables()
        bp <- wgcna_build_diagram_prompt(txt, organism, board_root, data_tables = data_tables)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- text_result$report_text() %||% ""
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.wgcna")
        data_tables <- text_result$report_data_tables()
        bp <- wgcna_build_diagram_prompt(txt, organism, board_root, data_tables = data_tables)
        llm <- get_ai_model(parent_session)
        make_llm_diagram_config(llm,
          system_prompt = bp$system,
          default_regulation = "positive",
          node_styles = wgcna_diagram_style()$node_styles,
          edge_styles = wgcna_diagram_style()$edge_styles
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
        img_style <- controls$image_style() %||% "bigomics"
        img_blocks <- as.integer(controls$image_blocks() %||% 1L)
        bp <- wgcna_build_image_prompt(txt, organism, edgelist,
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
        bp <- wgcna_build_image_prompt(txt, organism, edgelist,
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

    # Close progress when diagram completes
    shiny::observeEvent(diagram_result(), {
      diag <- diagram_result()
      shiny::req(diag)
      p <- report_progress()
      message(sprintf("[INFO][%s] --- [AI-REPORT] diagram complete, closing progress", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
      if (!is.null(p)) {
        p$set(message = "Done! Infographic generating in background...", value = 1)
        p$close()
        report_progress(NULL)
      }
    })

    # Layout rendering
    ai_report_layout_server(
      "layout",
      text_reactive = text_result$text,
      watermark = watermark
    )
  })
}
