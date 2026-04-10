##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' AI Report Controls UI
#'
#' Standard control panel for AI report options with mode toggle.
#' Mode and Generate button are always visible. Module selector and
#' style options are only visible in Summary mode (toggled via shinyjs).
#'
#' @param id Module namespace ID
#'
#' @return Shiny tagList
ai_report_controls_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    # Always visible: Mode toggle
    shiny::radioButtons(
      ns("mode"),
      "Mode:",
      choices = c("Report" = "report", "Summary" = "summary"),
      selected = "report",
      inline = TRUE
    ),

    # Always visible: Generate button
    shiny::actionButton(
      ns("generate_btn"),
      "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary btn-block",
      style = "margin-bottom: 10px;"
    ),

    # Summary-only controls (hidden when in Report mode)
    shinyjs::hidden(
      shiny::div(
        id = ns("summary_controls"),
        shiny::selectInput(
          ns("summary_module"),
          "Module:",
          choices = NULL,
          width = "100%"
        ),

      )
    ),

    ## NOTE: "Show Prompt" checkbox is rendered in the text card's hamburger
    ## menu (see ai_report_ui.R). Infographic style/blocks controls are rendered
    ## in the image card's hamburger menu. All are namespaced to this controls
    ## module so input$* is read here.
    ## Infographics are always computed (no toggle needed).
  )
}

#' AI Report Controls Server
#'
#' Server logic for control panel. Toggles visibility of summary-only
#' controls based on mode. Accepts module_choices reactive from parent
#' to populate the module selector.
#'
#' @param id Module namespace ID
#' @param module_choices Reactive returning character vector of module names
#'
#' @return List with reactives: trigger, mode, show_prompt, selected_module
ai_report_controls_server <- function(id, module_choices = NULL) {
  moduleServer(id, function(input, output, session) {

    # Toggle mode-specific controls based on mode selection
    shiny::observe({
      mode <- input$mode %||% "report"
      if (mode == "summary") {
        shinyjs::show("summary_controls")
      } else {
        shinyjs::hide("summary_controls")
      }
    })

    # Populate image style choices on init (deferred from UI to avoid startup crash)
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

    # Update module selector when choices change
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

    # Trigger counter (increments on button click)
    trigger <- reactiveVal(0)
    observeEvent(input$generate_btn, {
      trigger(trigger() + 1)
    })

    # Return reactive values
    list(
      trigger = reactive(trigger()),
      mode = reactive(input$mode %||% "report"),
      show_prompt = reactive(input$show_prompt),
      selected_module = reactive(input$summary_module),
      include_infographic = reactive(TRUE),
      image_style = reactive(input$image_style %||% "bigomics"),
      image_blocks = reactive(input$image_blocks %||% "1")
    )
  })
}
