##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' LASAGNA AI Report Tab UI
#'
#' Entry point for the AI Report tab body.
#' @param id Module namespace ID.
#' @return Shiny UI element.
lasagna_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)
  controls_ns <- shiny::NS(ns("controls"))

  show_prompt_input <- shiny::checkboxInput(
    controls_ns("show_prompt"),
    "Show prompt",
    FALSE
  )

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
    text_title        = "AI LASAGNA Report",
    diagram_title     = "Layer-Interaction Diagram",
    infographic_title = "Graphical Abstract",
    text_options      = show_prompt_input,
    infographic_options = infographic_options
  )
}

#' LASAGNA AI Report Sidebar UI
#'
#' Entry point for the AI Report sidebar controls.
#' @param id Module namespace ID.
#' @return Shiny UI element.
lasagna_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  multiomics_ai_report_controls_ui(ns("controls"), module_label = "Contrast:")
}

#' LASAGNA AI Report Server (coordinator)
#'
#' Thin coordinator that wires controls, text, diagram, infographic, and
#' layout servers together. This is the single entry point called from
#' `board_lasagna_server.R`.
#'
#' @return NULL.
lasagna_ai_report_server <- function(id,
                                     graph_data_reactive,
                                     contrast_reactive,
                                     contrast_choices_reactive,
                                     pgx,
                                     parent_session,
                                     watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    cache <- omicsai::omicsai_cache_init("mem")

    ## Low-signal heuristic: TRUE when the active contrast yields a
    ## graph with very few inter-layer edges. Deep Report retrieval is
    ## unproductive without anchors, so the Deep button is disabled.
    low_signal_reactive <- shiny::reactive({
      res <- graph_data_reactive()
      contrast <- contrast_reactive()
      if (is.null(res) || is.null(contrast)) return(FALSE)
      ctx <- tryCatch(lasagna_ai_extract_context(res, contrast, pgx, ntop = 5L),
                      error = function(e) NULL)
      if (is.null(ctx)) return(TRUE)
      isTRUE(ctx$network$n_inter_edges < 5L)
    })

    controls <- multiomics_ai_report_controls_server(
      "controls",
      module_choices = contrast_choices_reactive,
      low_signal     = low_signal_reactive
    )

    text_result <- lasagna_ai_text_server(
      "text",
      graph_data_reactive = graph_data_reactive,
      contrast_reactive   = contrast_reactive,
      pgx                 = pgx,
      controls            = controls,
      parent_session      = parent_session
    )

    ## Unified report text: picks deep_text when Deep Report was the
    ## last fired branch, otherwise the normal report_text. Never
    ## returns the prompt cache so diagram/image always receive actual
    ## LLM output.
    active_report_text <- shiny::reactive({
      if (isTRUE(text_result$last_deep())) {
        text_result$deep_text() %||% ""
      } else {
        text_result$report_text() %||% ""
      }
    })

    diagram_result <- AiDiagramCardServer(
      "layout-diagram",
      params_reactive = shiny::reactive({
        txt <- active_report_text()
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.mofa")
        bp <- lasagna_build_diagram_prompt(txt, organism, board_root)
        list(content = bp$board)
      }),
      template_reactive = shiny::reactive("{{content}}"),
      config_reactive = shiny::reactive({
        txt <- active_report_text()
        shiny::req(nzchar(txt))
        organism <- pgx$organism %||% "human"
        board_root <- file.path(OPG, "components/board.mofa")
        bp <- lasagna_build_diagram_prompt(txt, organism, board_root)
        llm <- get_ai_model(parent_session)
        make_llm_diagram_config(llm,
          system_prompt      = bp$system,
          default_regulation = "cross_layer",
          node_styles        = lasagna_diagram_style()$node_styles,
          edge_styles        = lasagna_diagram_style()$edge_styles
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (controls$mode() == "report") {
          controls$trigger()
        } else if (isTRUE(text_result$last_deep())) {
          controls$deep_trigger()
        } else {
          0
        }
      }),
      style = lasagna_diagram_style()
    )

    image_params_reactive <- shiny::reactive({
      txt <- active_report_text()
      shiny::req(nzchar(txt))
      organism <- pgx$organism %||% "human"
      diag <- diagram_result()
      edgelist <- if (!is.null(diag)) diag$edgelist else NULL
      img_style  <- controls$image_style() %||% "bigomics"
      img_blocks <- as.integer(controls$image_blocks() %||% 1L)
      bp <- lasagna_build_image_prompt(txt, organism, edgelist,
                                       style_name = img_style,
                                       n_blocks   = img_blocks)
      list(content = bp$board, system = bp$system,
           style = img_style, blocks = img_blocks)
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
          model         = img_model,
          system_prompt = p$system,
          style         = p$style,
          n_blocks      = p$blocks,
          image_size    = "1K"
        )
      }),
      cache = cache,
      trigger_reactive = shiny::reactive({
        if (!isTRUE(controls$include_infographic())) return(0)
        controls$image_trigger()
      }),
      watermark = watermark
    )

    multiomics_ai_report_layout_server(
      "layout",
      text_reactive = text_result$text,
      filename      = "lasagna-ai-report",
      title         = "LASAGNA AI Report",
      watermark     = watermark
    )

    invisible(NULL)
  })
}
