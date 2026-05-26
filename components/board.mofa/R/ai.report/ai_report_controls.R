##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Multi-omics AI Report Controls UI
#'
#' Shared by MofaBoard and LasagnaBoard. Single Generate button dispatches
#' on the active mode (Summary / Report / Deep Report). The Deep Report
#' mode runs the agentic literature-retrieval pipeline; mode-specific
#' extras (module selector for Summary, infographic toggle for Report /
#' Deep) are shown/hidden via shinyjs.
#'
#' @param id Module namespace ID.
#' @param module_label Label for the per-unit picker shown in Summary mode
#'   (e.g. "Factor:" for MOFA, "Contrast:" for LASAGNA).
#'
#' @return Shiny tagList.
multiomics_ai_report_controls_ui <- function(id, module_label = "Item:") {
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
          module_label,
          choices = NULL,
          width = "100%"
        ),
        shiny::radioButtons(
          ns("summary_style"),
          "Summary Style:",
          choices = c("Short Summary" = "short_summary", "Long Summary" = "long_summary"),
          selected = "short_summary",
          inline = FALSE
        )
      )
    )

    ## NOTE: "Show prompt" checkbox is rendered in the text card's hamburger
    ## menu (see ai_report_ui.R on each board). Infographic style/blocks
    ## controls are rendered in the image card's hamburger menu. All are
    ## namespaced to this controls module so input$* is read here.
  )
}

#' Multi-omics AI Report Controls Server
#'
#' Companion server for `multiomics_ai_report_controls_ui()`. Toggles
#' visibility of mode-specific blocks and emits two independent triggers
#' so existing eventReactives downstream keep firing exactly once per
#' relevant click:
#'   - `trigger`       — Summary or Report click
#'   - `deep_trigger`  — Deep Report click
#'   - `image_trigger` — Report or Deep Report click (image card fires only
#'                       on Report/Deep, never on a bare mode switch)
#'
#' @param id Module namespace ID.
#' @param module_choices Reactive returning character vector of per-unit
#'   names (factors for MOFA, contrasts for LASAGNA) — populates the
#'   Summary mode picker.
#' @param low_signal Reactive returning TRUE when the active dataset has
#'   no anchor for literature retrieval (no enrichment passes a meaningful
#'   threshold, etc.). When TRUE the Generate button is disabled while
#'   Deep Report mode is selected. Defaults to `reactive(FALSE)`.
#'
#' @return Named list of reactives: trigger, deep_trigger, image_trigger,
#'   mode, summary_style, show_prompt, selected_module,
#'   include_infographic, image_style, image_blocks.
multiomics_ai_report_controls_server <- function(id, module_choices = NULL,
                                                 low_signal = shiny::reactive(FALSE)) {
  moduleServer(id, function(input, output, session) {

    # Toggle mode-specific controls based on mode selection. Summary shows
    # the module dropdown; Report / Deep Report show the infographic
    # checkbox. Blocks are mutually exclusive.
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

    # Populate image style choices on init (deferred from UI to avoid
    # startup crash if omicsai is not yet loaded).
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
    # keep firing exactly once per relevant click. image_trigger captures
    # "this was a Report or Deep click" so the image card never fires
    # from mode switches alone.
    trigger       <- shiny::reactiveVal(0)
    deep_trigger  <- shiny::reactiveVal(0)
    image_trigger <- shiny::reactiveVal(0)
    shiny::observeEvent(input$generate_btn, {
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

    # Disable Generate when mode = Deep Report on low-signal datasets.
    shiny::observe({
      is_low <- isTRUE(low_signal())
      mode <- input$mode %||% "report"
      if (is_low && mode == "deep_report") {
        shinyjs::disable("generate_btn")
      } else {
        shinyjs::enable("generate_btn")
      }
    })

    list(
      trigger             = shiny::reactive(trigger()),
      deep_trigger        = shiny::reactive(deep_trigger()),
      image_trigger       = shiny::reactive(image_trigger()),
      mode                = shiny::reactive(input$mode %||% "report"),
      summary_style       = shiny::reactive(input$summary_style),
      show_prompt         = shiny::reactive(input$show_prompt),
      selected_module     = shiny::reactive(input$summary_module),
      include_infographic = shiny::reactive(isTRUE(input$include_infographic)),
      image_style         = shiny::reactive(input$image_style %||% "bigomics"),
      image_blocks        = shiny::reactive(input$image_blocks %||% "1")
    )
  })
}
