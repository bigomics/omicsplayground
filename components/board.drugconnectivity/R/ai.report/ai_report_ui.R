##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

drugconnectivity_ai_report_ui <- function(id) {
  ns <- shiny::NS(id)
  controls_ns <- shiny::NS(ns("controls"))

  show_prompt_input <- shiny::checkboxInput(
    controls_ns("show_prompt"),
    "Show prompt",
    FALSE
  )

  ## Infographic style/layout controls: namespaced to controls so the controls
  ## server can read them, but rendered in the image card's hamburger menu.
  ## The include_infographic checkbox lives in the sidebar (ai_report_controls.R).
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

  drugconnectivity_ai_report_layout_ui(
    ns("layout"),
    text_title = "AI Drug Connectivity Report",
    diagram_title = "Drug-MOA-Target Diagram",
    infographic_title = "Graphical Abstract",
    text_options = show_prompt_input,
    infographic_options = infographic_options
  )
}

drugconnectivity_ai_report_inputs_ui <- function(id) {
  ns <- shiny::NS(id)
  drugconnectivity_ai_report_controls_ui(ns("controls"))
}

drugconnectivity_ai_report_layout_ui <- function(id,
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
        download.fmt = c("pdf", "docx", "md")
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
          width = c("auto", "100%"),
          extra_options = infographic_options
        )
      )
    )
  )
}

drugconnectivity_ai_report_layout_server <- function(id, text_reactive, watermark = FALSE) {
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
      download.pdf  = .aicards_download_handler(dl_text, "drugconnectivity-ai-report", "Drug Connectivity AI Report", "pdf"),
      download.docx = .aicards_download_handler(dl_text, "drugconnectivity-ai-report", "Drug Connectivity AI Report", "docx"),
      download.md   = .aicards_download_handler(dl_text, "drugconnectivity-ai-report", "Drug Connectivity AI Report", "md"),
      pdf.width = 8,
      pdf.height = 11,
      res = c(75, 100),
      add.watermark = watermark
    )

    invisible(NULL)
  })
}
