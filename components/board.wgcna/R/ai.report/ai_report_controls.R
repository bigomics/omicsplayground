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
      choices = c("Summary" = "summary", "Report" = "report"),
      selected = "summary",
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
    shiny::div(
      id = ns("summary_controls"),
      shiny::selectInput(
        ns("summary_module"),
        "Module:",
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

    # Report-only controls (hidden when in Summary mode)
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
    ## NOTE: "Show Prompt" checkbox is rendered in the text card's hamburger
    ## menu (see ai_report_ui.R). It is namespaced to this controls module
    ## so input$show_prompt is still read here.
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
#' @return List with reactives: trigger, mode, summary_style, show_prompt, selected_module
ai_report_controls_server <- function(id, module_choices = NULL) {
  moduleServer(id, function(input, output, session) {

    # Toggle mode-specific controls based on mode selection
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
