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

    # Low-signal heuristic: TRUE when no module clears q < 0.05 in its
    # gse table. Deep Report retrieval is unproductive without anchors,
    # so the Deep button is disabled in this case (see controls server).
    low_signal_reactive <- shiny::reactive({
      w <- wgcna()
      if (is.null(w) || is.null(w$gse) || length(w$gse) == 0) return(FALSE)
      any_sig <- vapply(w$gse, function(g) {
        is.data.frame(g) && "q.value" %in% colnames(g) &&
          any(!is.na(g$q.value) & g$q.value < 0.05)
      }, logical(1))
      !any(any_sig)
    })

    # Controls
    controls <- ai_report_controls_server(
      "controls",
      module_choices = ai_module_choices,
      low_signal = low_signal_reactive
    )

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

    # Image / infographic card.
    # Reads from text_result$text (unified) so Deep Report's output also
    # feeds the image when the deep button fires. Gated on the opt-in
    # checkbox; trigger is the OR of Report and Deep Report triggers.
    image_text_reactive <- shiny::reactive(text_result$text() %||% "")
    image_params_reactive <- shiny::reactive({
      txt <- image_text_reactive()
      shiny::req(nzchar(txt))
      organism <- pgx$organism %||% "human"
      diag <- diagram_result()
      edgelist <- if (!is.null(diag)) diag$edgelist else NULL
      img_style <- controls$image_style() %||% "bigomics"
      img_blocks <- as.integer(controls$image_blocks() %||% 1L)
      bp <- wgcna_build_image_prompt(txt, organism, edgelist,
                                     style_name = img_style, n_blocks = img_blocks)
      list(content = bp$board, system = bp$system, style = img_style, blocks = img_blocks)
    })
    AiImageCardServer(
      "layout-infographic",
      params_reactive = shiny::reactive({
        p <- image_params_reactive()
        list(content = p$content)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        p <- image_params_reactive()
        img_model <- getUserOption(parent_session, "image_model")
        shiny::validate(shiny::need(
          !is.null(img_model) && nzchar(img_model),
          "No image model configured. Please select a model in Settings."
        ))
        omicsai::omicsai_image_config(
          model = img_model,
          system_prompt = p$system,
          style = p$style,
          n_blocks = p$blocks,
          image_size = "1K"
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (!isTRUE(controls$include_infographic())) return(0)
        ## image_trigger ticks only on Report or Deep clicks (never
        ## on Summary or on mode switches), so the card fires once
        ## per relevant Generate press and never spuriously.
        controls$image_trigger()
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
