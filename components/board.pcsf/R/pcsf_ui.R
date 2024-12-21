##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PcsfInputs <- function(id) {
  ns <- NS(id)

  bigdash::tabSettings(
    hr(),
    br(),
    withTooltip(
      selectInput(ns("contrast"), "Select contrast:",
        choices = NULL, multiple = FALSE
      ),
      "Select contrast.",
      placement = "right"
    ),
    hr(),
    withTooltip(
      shiny::sliderInput(ns("pcsf_beta"), "Prize strength (beta):", -4, 1, 0, 1),
      "Select prize strength. Smaller beta value corresponds to lower node prizes resulting in smaller solution size. A larger beta value corresponds to higher node prizes resulting in a larger graph (more greedy solution).",
      placement = "right"
    ),
    hr(),
    br(),
    withTooltip(
      actionLink(ns("adv_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::conditionalPanel(
      "input.adv_options % 2 == 1",
      ns = ns,
      br(),
      withTooltip(
        shiny::radioButtons(ns("pcsf_ntop"), "Initial network size:",
          choices = c("S" = 150, "M" = 500, "L" = 1500),
          selected = 500, inline = TRUE
        ),
        "Select initial network size (number of top genes) for ."
      )
    )
  )
}

pcsf_module_info <- "The PCSF network analysis uses the Prize-collection Steiner Forest algorithm to determine high-confidence subnetworks of highly correlated and highly differentially expressed genes. The STRING protein-protein interaction network is used as template. The PCSF solution may be used to identify 'driver' genes that appear as hubs in the network computed."

pcsf_graph_info <- "Prize-collection Steiner Forest solution for the top differential genes using the STRING database as backbone. 'Driver' genes appear as hubs in the network computed using a page-rank centrality measure."

PcsfUI <- function(id) {
  ns <- NS(id)
  tagList(
    boardHeader(
      title = "Prize-Collecting Steiner Forest", info_link = ns("pcsf_info")
    ),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "PCSF network",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 190px)",
          bs_alert(pcsf_module_info),
          bslib::layout_columns(
            col_widths = c(6, 6),
            height = "calc(100vh - 190px)",
            pcsf_plot_network_ui(
              ns("pcsf_network"),
              caption = paste(
                "PCSF network analysis",
                "Functional analysis of biological networks using Prize-collection Steiner Forest algorithm that determines high-confidence subnetworks."
              ),
              info.text = pcsf_graph_info,
              height = c("100%", "75vh"),
              width = c("auto", "100%")
            ),
            pcsf_table_centrality_ui(
              ns("centrality_table"),
              title = "Centrality score",
              info.text = "",
              caption = "Table showing the centrality score of genes.",
              width = c("100%", "100%"),
              height = c("100%", TABLE_HEIGHT_MODAL)
            )
          )
        )
      )
    )
  )
}
