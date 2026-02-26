##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

pcsf_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)

  show_prompt_input <- shiny::checkboxInput(
    shiny::NS(ns("controls"), "show_prompt"),
    "Show prompt",
    FALSE
  )

  pcsf_ai_report_layout_ui(
    ns("layout"),
    text_title = "AI PCSF Report",
    diagram_title = "PCSF Mechanism Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input
  )
}

pcsf_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  pcsf_ai_report_controls_ui(ns("controls"))
}

pcsf_ai_report_layout_ui <- function(id,
                                     text_title = "AI Report",
                                     diagram_title = "Board Diagram",
                                     infographic_title = "Graphical Abstract",
                                     text_options = NULL) {
  ns <- shiny::NS(id)

  bslib::layout_columns(
    col_widths = 12,
    height = "calc(100vh - 180px)",
    row_heights = c("auto", 1),
    bs_alert(
      HTML(paste0(
        "\u26a0\ufe0f <b>Disclaimer:</b> This page contains AI-generated content. ",
        "Please verify important information independently."
      )),
      translate = FALSE
    ),
    bslib::layout_columns(
      col_widths = c(6, 6),
      height = "calc(100vh - 220px)",
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
          width = c("auto", "100%")
        )
      )
    )
  )
}

pcsf_ai_report_layout_server <- function(id, text_reactive, watermark = FALSE) {
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
        filename = "pcsf-ai-report",
        title = "PCSF AI Report"
      ),
      pdf.width = 8,
      pdf.height = 11,
      res = c(75, 100),
      add.watermark = watermark
    )

    invisible(NULL)
  })
}
