# CopilotEvidenceUI.R — Evidence module UI
#
# Card layout for the evidence panel: dataset context, plot viewer,
# plot-history carousel, results table, and empty state. Conditional panels
# are driven by output flags set in CopilotEvidenceServer (suspendWhenHidden
# = FALSE so they evaluate while hidden).

#' Copilot Evidence Module UI
#'
#' Three-slot plot display (ggplot / plotly / iheatmapr) with conditional
#' visibility, a history carousel, and a results table.
#'
#' @param id Character scalar — Shiny module namespace.
#' @return A `shiny.tag` object.
#' @export
CopilotEvidenceUI <- function(id) {
  ns <- shiny::NS(id)

  ggplot_output_id    <- ns("evidence_ggplot")
  plotly_output_id    <- ns("evidence_plotly")
  iheatmapr_output_id <- ns("evidence_iheatmapr")

  shiny::div(
    style = "padding: 8px; height: 100%; overflow-y: auto;",

    # ---- Dataset context card (visible when a dataset is loaded) ----
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_dataset"), "']"),
      bslib::card(
        class = "mb-2",
        bslib::card_header("Dataset Info", class = "py-1 px-2 fw-bold"),
        bslib::card_body(
          class = "p-2",
          shiny::uiOutput(ns("dataset_info"))
        )
      )
    ),

    # ---- Plot card ----
    bslib::card(
      class = "mb-2",
      bslib::card_header("Evidence Plot", class = "py-1 px-2 fw-bold"),
      bslib::card_body(
        class = "p-2",

        # Empty state — shown when no artifact exists
        shiny::conditionalPanel(
          condition = paste0("!output['", ns("has_plot"), "']"),
          shiny::div(
            style = paste0(
              "min-height: 220px; display: flex; align-items: center;",
              "justify-content: center; text-align: center; color: #999;"
            ),
            shiny::div(
              shiny::icon("chart-area", class = "fa-2x", style = "color: #ddd;"),
              shiny::br(), shiny::br(),
              shiny::tags$span("Ask Copilot for a plot and it will appear here.")
            )
          )
        ),

        # ggplot renderer
        shiny::conditionalPanel(
          condition = paste0("output['", ns("is_ggplot_active"), "']"),
          shiny::imageOutput(ns("evidence_ggplot"), height = "auto")
        ),

        # plotly renderer
        shiny::conditionalPanel(
          condition = paste0("output['", ns("is_plotly_active"), "']"),
          plotly::plotlyOutput(ns("evidence_plotly"), height = "350px")
        ),

        # iheatmapr renderer
        shiny::conditionalPanel(
          condition = paste0("output['", ns("is_iheatmapr_active"), "']"),
          iheatmapr::iheatmaprOutput(ns("evidence_iheatmapr"), height = "350px")
        )
      )
    ),

    # ---- Plot history carousel (visible when 2+ plots in history) ----
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_history"), "']"),
      shiny::div(
        class = "mb-2",
        shiny::div(
          class = "d-flex align-items-center mb-1",
          shiny::tags$small(class = "text-muted fw-bold", "Plot History")
        ),
        shiny::div(
          style = "overflow-x: auto; white-space: nowrap; padding: 4px 0;",
          shiny::uiOutput(ns("plot_carousel"))
        )
      )
    ),

    # ---- Results table card ----
    shiny::conditionalPanel(
      condition = paste0("output['", ns("has_table"), "']"),
      bslib::card(
        class = "mb-2",
        bslib::card_header("Results", class = "py-1 px-2"),
        bslib::card_body(
          class = "p-1",
          DT::dataTableOutput(ns("evidence_table"))
        )
      )
    ),

    # ---- Empty state (no dataset loaded yet) ----
    shiny::conditionalPanel(
      condition = paste0("!output['", ns("has_dataset"), "']"),
      shiny::div(
        style = "text-align: center; padding: 40px; color: #999;",
        shiny::icon("chart-line", class = "fa-3x", style = "color: #ddd;"),
        shiny::br(), shiny::br(),
        shiny::p("Evidence and visualisations will appear here when you interact with the Copilot.")
      )
    ),

    # ---- plot_rendered JS listener (telemetry only — no save logic depends on it) ----
    shiny::tags$script(shiny::HTML(sprintf(
      "(function() {
        var ids = [%s, %s, %s];
        var inputId = '%s';
        ids.forEach(function(id) {
          $(document).on('shiny:value', '#' + id, function() {
            Shiny.setInputValue(inputId, {
              output_id: id,
              ts: new Date().toISOString(),
              nonce: Math.random()
            }, {priority: 'event'});
          });
        });
      })();",
      jsonlite::toJSON(ggplot_output_id,    auto_unbox = TRUE),
      jsonlite::toJSON(plotly_output_id,    auto_unbox = TRUE),
      jsonlite::toJSON(iheatmapr_output_id, auto_unbox = TRUE),
      ns("plot_rendered")
    )))
  )
}
