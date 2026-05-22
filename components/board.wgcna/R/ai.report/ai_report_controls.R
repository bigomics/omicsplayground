##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' AI Report Controls UI
#'
#' Single-button control panel: mode picks the path (Report / Summary /
#' Deep Report) and Generate fires it. Mode-specific extras (module
#' selector, infographic toggle) are shown/hidden via shinyjs.
#'
#' @param id Module namespace ID
#'
#' @return Shiny tagList
ai_report_controls_ui <- function(id) {
  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::radioButtons(
      ns("mode"),
      "Mode:",
      choices = c(
        "Summary" = "summary",
        "Report" = "report",
        "Deep Report" = "deep_report"
      ),
      selected = "report"
    ),

    shiny::actionButton(
      ns("generate_btn"),
      "Generate!",
      icon = icon("refresh"),
      class = "btn-outline-primary btn-block",
      style = "margin-bottom: 10px;"
    ),

    # Report/Deep-Report extras (hidden in Summary mode).
    shinyjs::hidden(
      shiny::div(
        id = ns("report_extras"),
        aicards_include_infographic_input(ns("include_infographic"))
      )
    ),

    # Summary-only controls (hidden in Report/Deep modes).
    shinyjs::hidden(
      shiny::div(
        id = ns("summary_controls"),
        shiny::selectInput(
          ns("summary_module"),
          "Module:",
          choices = NULL,
          width = "100%"
        )
      )
    )

    ## NOTE: "Show Prompt" checkbox is rendered in the text card's hamburger
    ## menu (see ai_report_ui.R). Infographic style/blocks controls are
    ## rendered in the image card's hamburger menu. All are namespaced to
    ## this controls module so input$* is read here.
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
#' @param low_signal Reactive returning TRUE when the active WGCNA result
#'   has weak module-trait signal (no module clears `q < 0.05` enrichment).
#'   When TRUE the Deep Report button is disabled with an explanatory
#'   tooltip — literature retrieval would be unproductive without anchors.
#'   Defaults to `reactive(FALSE)` so callers that don't pass this keep
#'   the previous always-enabled behaviour.
#'
#' @return List with reactives: trigger, deep_trigger, mode, show_prompt,
#'   selected_module, include_infographic, image_style, image_blocks
ai_report_controls_server <- function(id, module_choices = NULL,
                                      low_signal = shiny::reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {

    # Toggle mode-specific controls based on mode selection.
    # Summary shows module dropdown; Report / Deep Report show the
    # infographic checkbox. Mode block is mutually-exclusive.
    shiny::observe({
      mode <- input$mode %||% "report"
      if (mode == "summary") {
        shinyjs::show("summary_controls")
        shinyjs::hide("report_extras")
      } else {
        shinyjs::hide("summary_controls")
        shinyjs::show("report_extras")
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

    # Single Generate button; mode at click time decides which counter
    # ticks. Trigger and deep_trigger are independent so existing
    # eventReactives (Report/Summary on trigger, Deep on deep_trigger)
    # keep firing exactly once per relevant click. image_trigger is a
    # separate counter that captures "this was a Report or Deep click"
    # so the image card never fires from mode switches alone.
    trigger <- reactiveVal(0)
    deep_trigger <- reactiveVal(0)
    image_trigger <- reactiveVal(0)
    observeEvent(input$generate_btn, {
      mode <- shiny::isolate(input$mode %||% "report")
      if (mode == "deep_report") {
        deep_trigger(deep_trigger() + 1)
      } else {
        trigger(trigger() + 1)
      }
      if (mode %in% c("report", "deep_report")) {
        image_trigger(image_trigger() + 1)
      }
    })

    # Disable Generate when mode = Deep Report on low-signal datasets —
    # literature retrieval is unproductive without enrichment anchors.
    # Toggled via shinyjs; cost-hint label stays visible.
    shiny::observe({
      is_low <- isTRUE(low_signal())
      mode <- input$mode %||% "report"
      if (is_low && mode == "deep_report") {
        shinyjs::disable("generate_btn")
      } else {
        shinyjs::enable("generate_btn")
      }
    })

    # Return reactive values
    list(
      trigger = reactive(trigger()),
      deep_trigger = reactive(deep_trigger()),
      image_trigger = reactive(image_trigger()),
      mode = reactive(input$mode %||% "report"),
      show_prompt = reactive(input$show_prompt),
      selected_module = reactive(input$summary_module),
      include_infographic = reactive(isTRUE(input$include_infographic)),
      image_style = reactive(input$image_style %||% "bigomics"),
      image_blocks = reactive(input$image_blocks %||% "1")
    )
  })
}
