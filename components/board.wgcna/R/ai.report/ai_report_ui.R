##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' WGCNA AI Report Tab UI
#'
#' Entry point for the AI Report tab body.
#' @param id Module namespace ID
#' @return Shiny UI element
wgcna_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)

  ## "Show prompt" checkbox: namespaced to controls so the controls server
  ## can read it via input$show_prompt, but rendered in the card hamburger menu.
  show_prompt_input <- shiny::checkboxInput(
    shiny::NS(ns("controls"), "show_prompt"),
    "Show prompt",
    FALSE
  )

  ai_report_layout_ui(
    ns("layout"),
    text_title = "AI Report",
    diagram_title = "Module Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input
  )
}

#' WGCNA AI Report Sidebar UI
#'
#' Entry point for the AI Report sidebar controls.
#' @param id Module namespace ID
#' @return Shiny UI element
wgcna_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  ai_report_controls_ui(ns("controls"))
}

#' AI Report Layout - Always 3 Cards
#'
#' Creates a fixed 3-card layout: text (left 50%) + diagram and infographic
#' (right 50% stacked). The text card is shared between summary and report
#' modes — the parent server decides what content feeds into it.
#'
#' @param id Module namespace ID
#' @param text_title Title for text panel (default "AI Report")
#' @param diagram_title Title for diagram panel (default "Board Diagram")
#' @param infographic_title Title for infographic panel (default "Graphical Abstract")
#'
#' @return Shiny UI tagList
ai_report_layout_ui <- function(id,
                                text_title = "AI Report",
                                diagram_title = "Board Diagram",
                                infographic_title = "Graphical Abstract",
                                text_options = NULL) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = 12,
    height = "calc(100vh - 180px)",
    row_heights = c("auto", 1),

    # Disclaimer banner
    bs_alert(
      HTML(paste0(
        "\u26a0\ufe0f <b>Disclaimer:</b> This page contains AI-generated content. ",
        "Please verify important information independently."
      )),
      translate = FALSE
    ),

    # Always 3 cards: text (left) | diagram + infographic (right stacked)
    bslib::layout_columns(
      col_widths = c(6, 6),
      height = "calc(100vh - 220px)",

      # LEFT: Text card (shared between summary and report modes)
      PlotModuleUI(
        ns("report_text"),
        outputFunc = htmlOutput,
        title = text_title,
        options = text_options,
        caption = "AI-generated analysis content",
        height = c("100%", TABLE_HEIGHT_MODAL),
        width = c("auto", "100%"),
        download.fmt = c("pdf")
      ),

      # RIGHT: Diagram + Infographic stacked
      bslib::layout_columns(
        col_widths = 12,
        row_heights = c(1, 1),

        AiDiagramCardUI(
          ns("diagram"),
          title = diagram_title,
          caption = "AI-generated module diagram",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),

        AiImageCardUI(
          ns("infographic"),
          title = infographic_title,
          caption = "AI-generated graphical abstract",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )
}

#' AI Report Layout - Server (mode-agnostic)
#'
#' Renders the 3 cards with whatever content the parent provides.
#' No mode awareness — the parent server handles mode-switching logic.
#'
#' @param id Module namespace ID
#' @param text_reactive Reactive returning markdown text string.
#'   Fed by parent: summary text in summary mode, report text in report mode.
#' @param diagram_result_reactive Deprecated; unused
#' @param infographic_reactive Deprecated; unused
#' @param watermark Logical; add watermark to outputs
#'
#' @return NULL (no return values needed)
ai_report_layout_server <- function(id,
                                    text_reactive,
                                    diagram_result_reactive = NULL,
                                    infographic_reactive = NULL,
                                    watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    # ============================================================
    # CARD 1: Text (left panel — shared between summary and report)
    # ============================================================
    text.RENDER <- function() {
      text <- text_reactive()
      shiny::validate(shiny::need(!is.null(text), "Click 'Generate!' to create AI content"))
      shiny::div(
        class = "gene-info",
        shiny::HTML(opg_markdown_to_html(text))
      )
    }

    text.RENDER2 <- function() {
      text <- text_reactive()
      shiny::validate(shiny::need(!is.null(text), "Content not available"))
      shiny::div(
        style = "font-size: 18px;",
        shiny::HTML(opg_markdown_to_html(text))
      )
    }

    PlotModuleServer(
      "report_text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text.RENDER,
      func2 = text.RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      download.pdf = .aicards_markdown_pdf_download(
        text_reactive = shiny::reactive({
          txt <- text_reactive()
          shiny::req(!is.null(txt))
          .aicards_text_markdown(txt)
        }),
        filename = "wgcna-ai-report",
        title = "WGCNA AI Report"
      ),
      pdf.width = 8, pdf.height = 11,
      res = c(75, 100),
      add.watermark = watermark
    )

    invisible(NULL)
  })
}
