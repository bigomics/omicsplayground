##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FunctionalInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(
      shiny::selectInput(ns("fa_contrast"), "Contrast:",
        choices = NULL
      ),
      "Select the contrast corresponding to the comparison of interest.",
      placement = "top"
    ),
    withTooltip(
      shiny::actionLink(ns("fa_options"), "Options",
        icon = icon("cog", lib = "glyphicon")
      ),
      "Show/hide advanced options",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.fa_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::checkboxInput(
            ns("fa_filtertable"),
            "filter signficant (tables)",
            FALSE
          ),
          "Click to filter the significant entries in the tables."
        )
      )
    )
  )
}

FunctionalUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  tabs <- shiny::tabsetPanel(
    id = ns("tabs"),
    shiny::tabPanel(
      "KEGG",
      div(
        class = "row",
        div(
          class = "col-md-6",
          functional_plot_kegg_graph_ui(
            ns("kegg_graph"),
            label = "a",
            height = c(400, TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")            
          ),
          functional_table_kegg_table_ui(
            ns("kegg_table"),
            label = "b",
            height = c(330, TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        div(
          class = "col-md-6",
          functional_plot_kegg_actmap_ui(
            ns("kegg_actmap"),
            label = "c",
            height = c(750, TABLE_HEIGHT_MODAL)            
          )
        )
      )
    ),
    shiny::tabPanel(
      "GO graph",
      div(
        class = "row",
        div(
          class = "col-md-6",
          functional_plot_go_network_ui(
            id = ns("GO_network"),
            height = c(400, TABLE_HEIGHT_MODAL),            
            label = "a"
          ),
          functional_table_go_table_ui(
            id = ns("GO_table"),
            height = c(330, TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        div(
          class = "col-md-6",
          functional_plot_go_actmap_ui(
            id = ns("GO_actmap"),
            height = c(750, TABLE_HEIGHT_MODAL),           
            label = "c"
          )
        )
      )
    )
  )

  page_ui <- div(
    boardHeader(title = "Pathway Analysis", info_link = ns("fa_info")),
    tabs
  )
  return(page_ui)
}
