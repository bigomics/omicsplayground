##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Cost / ETA indicators rendered under each generate button.
## Single source of truth — bump here when pricing or model changes.
.WGCNA_AI_REPORT_COSTS <- list(
  report      = list(cost = "~$0.01", eta = "~20s"),
  deep_report = list(cost = "~$0.05", eta = "~2min")
)

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

  .cost_hint <- function(slot) {
    cfg <- .WGCNA_AI_REPORT_COSTS[[slot]]
    shiny::tags$small(
      style = "display:block; color:#888; margin-top:-6px; margin-bottom:10px; text-align:center;",
      sprintf("%s / %s", cfg$cost, cfg$eta)
    )
  }

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
      style = "margin-bottom: 2px;"
    ),
    .cost_hint("report"),

    # Deep Report button (literature-grounded agentic report).
    # Disabled on low-signal datasets via shinyjs (toggled in server).
    shiny::actionButton(
      ns("deep_generate_btn"),
      "Generate Deep Report",
      icon = icon("book"),
      class = "btn-outline-primary btn-block",
      style = "margin-bottom: 2px;"
    ),
    .cost_hint("deep_report"),

    # Opt-in infographic generation (default off to keep cost low).
    shiny::checkboxInput(
      ns("include_infographic"),
      "Include infographic",
      value = FALSE
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

    # Deep Report trigger counter (driven by deep_generate_btn).
    # Kept separate from `trigger` so the regular Report path is never
    # fired by the deep button and vice-versa.
    deep_trigger <- reactiveVal(0)
    observeEvent(input$deep_generate_btn, {
      deep_trigger(deep_trigger() + 1)
    }, ignoreNULL = TRUE, ignoreInit = TRUE)

    # Disable Deep Report on low-signal datasets — literature retrieval
    # would be unproductive without enrichment anchors. Toggled via
    # shinyjs so the visible cost-hint label remains.
    shiny::observe({
      is_low <- isTRUE(low_signal())
      if (is_low) {
        shinyjs::disable("deep_generate_btn")
      } else {
        shinyjs::enable("deep_generate_btn")
      }
    })

    # Return reactive values
    list(
      trigger = reactive(trigger()),
      deep_trigger = reactive(deep_trigger()),
      mode = reactive(input$mode %||% "report"),
      show_prompt = reactive(input$show_prompt),
      selected_module = reactive(input$summary_module),
      include_infographic = reactive(isTRUE(input$include_infographic)),
      image_style = reactive(input$image_style %||% "bigomics"),
      image_blocks = reactive(input$image_blocks %||% "1")
    )
  })
}
