##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Multi-omics AI Report Layout UI
#'
#' Shared by MofaBoard and LasagnaBoard. Renders the fixed 3-card layout:
#' text (left 50%) + diagram and infographic (right 50% stacked).
#'
#' @param id Module namespace ID.
#' @param text_title Title for the text panel.
#' @param diagram_title Title for the diagram panel.
#' @param infographic_title Title for the infographic panel.
#' @param text_options Extra UI rendered in the text card's hamburger menu
#'   (typically the "Show prompt" checkbox).
#' @param infographic_options Extra UI rendered in the image card's
#'   hamburger menu (style / layout pickers).
#'
#' @return Shiny tagList.
multiomics_ai_report_layout_ui <- function(id,
                                           text_title = "AI Report",
                                           diagram_title = "Board Diagram",
                                           infographic_title = "Graphical Abstract",
                                           text_options = NULL,
                                           infographic_options = NULL) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = 12,
    height = "calc(100vh - 180px)",
    row_heights = c("auto", 1),

    bs_alert(
      HTML(paste0(
        "⚠️ <b>Disclaimer:</b> This page contains AI-generated content. ",
        "Please verify important information independently."
      )),
      translate = FALSE
    ),

    bslib::layout_columns(
      col_widths = c(6, 6),
      height = "calc(100vh - 220px)",

      # LEFT: Text card (shared between summary / report / deep report)
      PlotModuleUI(
        ns("report_text"),
        outputFunc = htmlOutput,
        title = text_title,
        options = text_options,
        caption = "AI-generated analysis content",
        height = c("100%", TABLE_HEIGHT_MODAL),
        width = c("auto", "100%"),
        download.fmt = c("pdf", "docx", "md")
      ),

      # RIGHT: Diagram + Infographic stacked
      bslib::layout_columns(
        col_widths = 12,
        row_heights = c(1, 1),

        AiDiagramCardUI(
          ns("diagram"),
          title = diagram_title,
          caption = "AI-generated mechanism diagram",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        ),

        AiImageCardUI(
          ns("infographic"),
          title = infographic_title,
          caption = "AI-generated graphical abstract",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%"),
          extra_options = infographic_options
        )
      )
    )
  )
}

#' Multi-omics AI Report Layout Server
#'
#' Renders the text card with whatever content the parent provides via
#' `text_reactive`. The diagram and image cards are wired by the parent
#' (via AiDiagramCardServer / AiImageCardServer); this server owns only
#' the text card and its PDF/docx/md download handlers.
#'
#' @param id Module namespace ID.
#' @param text_reactive Reactive returning markdown text string. The parent
#'   server is responsible for choosing report-vs-deep-report-vs-prompt-cache.
#' @param filename Base filename for downloads (no extension).
#' @param title Document title used in PDF/docx headers.
#' @param watermark Logical; add watermark to outputs.
#'
#' @return invisible(NULL).
multiomics_ai_report_layout_server <- function(id,
                                               text_reactive,
                                               filename = "multiomics-ai-report",
                                               title = "Multi-omics AI Report",
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    text.RENDER <- function() {
      text <- text_reactive()
      shiny::validate(shiny::need(!is.null(text), "Click 'Generate!' to create AI content"))
      shiny::div(class = "gene-info", shiny::HTML(opg_markdown_to_html(text)))
    }

    text.RENDER2 <- function() {
      text <- text_reactive()
      shiny::validate(shiny::need(!is.null(text), "Content not available"))
      shiny::div(style = "font-size: 18px;", shiny::HTML(opg_markdown_to_html(text)))
    }

    dl_text <- shiny::reactive({
      txt <- text_reactive()
      shiny::req(!is.null(txt))
      .aicards_text_markdown(txt)
    })

    PlotModuleServer(
      "report_text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text.RENDER,
      func2 = text.RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      download.fmt = c("pdf", "docx", "md"),
      download.pdf  = .aicards_download_handler(dl_text, filename, title, "pdf"),
      download.docx = .aicards_download_handler(dl_text, filename, title, "docx"),
      download.md   = .aicards_download_handler(dl_text, filename, title, "md"),
      pdf.width = 8,
      pdf.height = 11,
      res = c(75, 100),
      add.watermark = watermark
    )

    invisible(NULL)
  })
}
