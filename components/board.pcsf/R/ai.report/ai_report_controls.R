##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' PCSF AI Report Controls UI
pcsf_ai_report_controls_ui <- function(id) {
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
        "Contrast:",
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

#' PCSF AI Report Controls Server
pcsf_ai_report_controls_server <- function(id, module_choices = NULL) {
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

    trigger <- reactiveVal(0)
    observeEvent(input$generate_btn, {
      trigger(trigger() + 1)
    })

    list(
      trigger = reactive(trigger()),
      mode = reactive(input$mode %||% "summary"),
      summary_style = reactive(input$summary_style),
      show_prompt = reactive(input$show_prompt),
      selected_module = reactive(input$summary_module),
      include_infographic = reactive(isTRUE(input$include_infographic)),
      image_style = reactive(input$image_style %||% "bigomics"),
      image_blocks = reactive(input$image_blocks %||% "1")
    )
  })
}
