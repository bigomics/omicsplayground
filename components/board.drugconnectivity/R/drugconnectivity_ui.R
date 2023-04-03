##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

DrugConnectivityInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(),
    withTooltip(shiny::selectInput(ns("dsea_contrast"), "Contrast:", choices = NULL),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    ),
    withTooltip(shiny::selectInput(ns("dsea_method"), "Analysis type:", choices = ""),
      "Select type of drug enrichment analysis: activity or sensitivity (if available).",
      placement = "top"
    ),
    shiny::hr(),
    withTooltip(
      shiny::checkboxInput(
        ns("dseatable_filter"),
        "only annotated drugs",
        FALSE
      ),
      "Show only annotated drugs."
    )
  )
}

DrugConnectivityUI <- function(id) {
  ns <- shiny::NS(id)

  div(
    boardHeader(title = "Drug Connectivity", info_link = ns("dsea_info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Drug enrichment",
        div(
          class = "row",
          div(
            class = "col-md-9",
            div(
              class = "row",
              div(
                class = "col-md-6",
                drugconnectivity_plot_enplots_ui(
                  id = ns("dsea_enplots"),
                  height = c(350, TABLE_HEIGHT_MODAL),
                  label = "a"
                )
              ),
              div(
                class = "col-md-6",
                drugconnectivity_plot_moa_ui(
                  id = ns("dsea_moaplot"),
                  height = c(350, TABLE_HEIGHT_MODAL),                  
                  label = "c"
                )
              )
            ),
            drugconnectivity_table_dsea_ui(
              ns("dsea_table"),
              height = c(300, TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )
          ),
          div(
            class = "col-md-3",
            drugconnectivity_plot_actmap_ui(
              ns("dsea_actmap"),
              height = c(700, TABLE_HEIGHT_MODAL),
              width = c("100%", "100%"),              
              label = "d"
            )
          )
        )
      ),
      shiny::tabPanel(
        "Connectivity map (beta)",
        shiny::div(
          class = "row",
          div(
            class = "col-md-5",
            drugconnectivity_plot_cmap_enplot_ui(
              id = ns("cmap_enplot"),
              label = "a"
            ),
            drugconnectivity_table_cmap_ui(
              id = ns("cmap_table"),
              height = c(380, TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )
          ),
          div(
            class = "col-md-7",
            drugconnectivity_plot_cmap_dsea_ui(
              id = ns("cmap_dsea"),
              label = "c"
            )
          )
        ) ## end of row
      )
    )
  )
}
