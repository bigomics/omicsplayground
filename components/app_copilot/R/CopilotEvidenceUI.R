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
    style = "padding: 0px; height: 100%; overflow-y: auto;",

    ## # ---- Dataset context card (visible when a dataset is loaded) ----
    ## shiny::conditionalPanel(
    ##   condition = paste0("output['", ns("has_dataset"), "']"),
    ##   bslib::card(
    ##     class = "mb-2",
    ##     bslib::card_header("Dataset Info", class = "py-1 px-2 fw-bold"),
    ##     bslib::card_body(
    ##       class = "p-2",
    ##       shiny::uiOutput(ns("dataset_info"))
    ##     )
    ##   )
    ## ),

    # ---- Plot card ----
    # Header mirrors the standard PlotModuleUI look: a `plotmodule-header`
    # fillRow placed via `as.card_item` (NOT card_header), with the icon
    # buttons wearing the same `btn-circle-xs` / `download-button` /
    # `zoom-button` classes so they pick up the project's existing CSS.
    bslib::card(
      class = "mb-2",
      min_height = 500,
      bslib::as.card_item(shiny::div(
        shiny::fillRow(
          flex = c(1, NA, NA),
          class = "plotmodule-header",
          height = "33px",
          shiny::div(
            class = "plotmodule-title",
            style = "white-space: nowrap; overflow: hidden; text-overflow: clip;",
            "Plot output"
          ),
          shiny::div(
            class = "download-button", title = "download",
            shiny::conditionalPanel(
              condition = paste0("output['", ns("has_plot"), "']"),
              # DropdownMenu trigger (download icon, circular) + popover body
              # holding width / height in inches and the actual downloadButton.
              # Mirrors PlotModule's pattern at ui-PlotModule.R:216-245.
              DropdownMenu(
                shiny::div(
                  style = "width: 240px;",
                  # Format selector — rendered server-side because the available
                  # choices depend on artifact kind (ggplot=PNG only;
                  # plotly/iheatmapr=HTML+PNG).
                  shiny::uiOutput(ns("dl_format_ui")),
                  shiny::div(
                    style = "display: flex; gap: 10px; margin-bottom: 12px;",
                    shiny::div(
                      style = "flex: 1;",
                      shiny::numericInput(
                        ns("dl_width"), "Width (in)",
                        value = 8, min = 1, max = 20, step = 1, width = "100%"
                      )
                    ),
                    shiny::div(
                      style = "flex: 1;",
                      shiny::numericInput(
                        ns("dl_height"), "Height (in)",
                        value = 6, min = 1, max = 20, step = 1, width = "100%"
                      )
                    )
                  ),
                  shiny::downloadButton(
                    ns("evidence_download"),
                    label = "Download",
                    class = "btn-outline-primary btn-sm w-100"
                  )
                ),
                size = "xs",
                icon = shiny::icon("download"),
                status = "default"
              )
            )
          ),
          shiny::div(
            class = "zoom-button", title = "zoom",
            shiny::conditionalPanel(
              condition = paste0("output['", ns("has_plot"), "']"),
              # data-bs-toggle modalTrigger (no R round-trip), pointed at the
              # always-mounted .popup-modal declared below. Mirrors PlotModule's
              # zoom path so the modal inherits the same global CSS and sizing.
              modalTrigger(
                ns("evidence_zoombutton"),
                ns("evidencePopup"),
                shiny::icon("up-right-and-down-left-from-center"),
                class = "btn-circle-xs"
              )
            )
          )
        )
      )),
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

    # ---- Maximize modal ----
    # Recipe copied from PlotModule (`ui-PlotModule.R:480-512`):
    #   * `modalUI(size = "fullscreen")` -> Bootstrap `.modal-fullscreen` so
    #     the dialog spans the viewport instead of getting capped by the
    #     project-global `.modal-lg { width: 640px }` override.
    #   * `.popup-modal` / `.popup-plot` wrappers pick up the global rules
    #     in `playground.css` (margin auto, fonts, footer hidden).
    #   * Four inline CSS rules per instance — `.modal-content { width: 100vw }`
    #     is the load-bearing one that stretches the content layer to viewport
    #     width; without it the plotly inside collapses to a narrow column.
    #   * Inner output is a `uiOutput` slot so only the active kind's output
    #     element is mounted, avoiding the zero-width autosize problem that
    #     `display: none` `conditionalPanel`s cause for plotly.
    shiny::div(
      class = "popup-modal",
      modalUI(
        id = ns("evidencePopup"),
        title = "Evidence Plot",
        size = "fullscreen",
        footer = NULL,
        shiny::div(
          class = "popup-plot-body",
          shiny::div(
            class = "popup-plot",
            shiny::uiOutput(ns("evidence_modal_slot"))
          )
        )
      ),
      shiny::tags$head(shiny::tags$style(shiny::HTML(sprintf(
        "#%s .modal-dialog  { width: 100%%; }
         #%s .modal-body    { min-height: calc(80vh - 100px); padding: 30px 150px; }
         #%s .modal-content { width: 100vw; }
         #%s .modal-footer  { display: none; }",
        ns("evidencePopup"), ns("evidencePopup"),
        ns("evidencePopup"), ns("evidencePopup")
      ))))
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
