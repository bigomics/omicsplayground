##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiomics_ai_report_controls_ui <- function(id, module_label = "Item:") {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::radioButtons(
      ns("mode"),
      "Mode:",
      choices = c("Summary" = "summary", "Report" = "report"),
      selected = "summary",
      inline = TRUE
    ),
    shiny::actionButton(
      ns("generate_btn"),
      "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary btn-block",
      style = "margin-bottom: 10px;"
    ),
    shiny::div(
      id = ns("summary_controls"),
      shiny::selectInput(
        ns("summary_module"),
        module_label,
        choices = NULL,
        width = "100%"
      ),
      shiny::radioButtons(
        ns("summary_style"),
        "Summary Style:",
        choices = c("Short Summary" = "short", "Long Summary" = "long"),
        selected = "short",
        inline = FALSE
      )
    ),
    shinyjs::hidden(
      shiny::div(
        id = ns("image_controls"),
        shiny::checkboxInput(
          ns("include_infographic"),
          "Include infographic",
          value = FALSE
        ),
        shiny::selectInput(
          ns("image_style"),
          "Infographic Style:",
          choices = NULL,
          width = "100%"
        ),
        shiny::radioButtons(
          ns("image_blocks"),
          "Layout:",
          choices = c("1 Panel" = "1", "2 Panels" = "2", "3 Panels" = "3"),
          selected = "1",
          inline = TRUE
        )
      )
    )
  )
}

multiomics_ai_report_controls_server <- function(id, module_choices = NULL) {
  moduleServer(id, function(input, output, session) {
    shiny::observe({
      mode <- input$mode %||% "summary"
      if (mode == "summary") {
        shinyjs::show("summary_controls")
        shinyjs::hide("image_controls")
      } else {
        shinyjs::hide("summary_controls")
        shinyjs::show("image_controls")
      }
    })

    shiny::observe({
      styles <- tryCatch(
        omicsai::omicsai_available_image_styles(),
        error = function(e) "bigomics"
      )
      shiny::updateSelectInput(
        session, "image_style",
        choices = styles,
        selected = if ("bigomics" %in% styles) "bigomics" else styles[1]
      )
    })

    shiny::observe({
      if (is.null(module_choices)) return()
      choices <- module_choices()
      shiny::req(choices)
      shiny::updateSelectInput(
        session, "summary_module",
        choices = choices,
        selected = choices[1]
      )
    })

    trigger <- shiny::reactiveVal(0)
    shiny::observeEvent(input$generate_btn, {
      trigger(trigger() + 1)
    })

    list(
      trigger = shiny::reactive(trigger()),
      mode = shiny::reactive(input$mode %||% "summary"),
      summary_style = shiny::reactive(input$summary_style),
      show_prompt = shiny::reactive(input$show_prompt),
      selected_module = shiny::reactive(input$summary_module),
      include_infographic = shiny::reactive(isTRUE(input$include_infographic)),
      image_style = shiny::reactive(input$image_style %||% "bigomics"),
      image_blocks = shiny::reactive(input$image_blocks %||% "1")
    )
  })
}

multiomics_ai_report_layout_ui <- function(id,
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
        filename = filename,
        title = title
      ),
      pdf.width = 8,
      pdf.height = 11,
      res = c(75, 100),
      add.watermark = watermark
    )

    invisible(NULL)
  })
}
