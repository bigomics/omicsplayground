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
      radioButtons(
        ns("colorby"),
        "Color by:",
        choices = c("gene.cluster", "contrast"),
        selected = "gene.cluster",
        inline = TRUE
      ),
      "Choose how to color the nodes",
      placement = "right",
      options = list(container = "body")
    ),
    conditionalPanel(
      "input.colorby == 'contrast'",
      ns = ns,
      withTooltip(
        selectInput(ns("contrast"), NULL, choices = NULL, multiple = FALSE),
        "Select contrast.",
        placement = "right",
        options = list(container = "body")
      ),
    ),
    br(),
    withTooltip(
      actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top",
      options = list(container = "body")
    )
  )
}

pcsf_module_info <- "PCSF Network Analysis. Functional analysis of biological networks using Prize-collection Steiner Forest algorithm that determines high-confidence subnetworks."

pcsf_graph_info <- "Prize-collection Steiner Forest solution for the top differential genes using the STRING database as backbone."

PcsfUI <- function(id) {
  ns <- NS(id)

  div(
    boardHeader(
      title = "PCSF",
      info_link = ns("pcsf_info")
    ),
    tabsetPanel(
      id = ns("tabs1"),
      tabPanel(
        "PCSF network",
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 190px)",
          heights_equal = "row",
          ## bs_alert(pcsf_module_info),
          bslib::layout_column_wrap(
            width = 1,
            height = "100%",
            heights_equal = "row",            
            style = htmltools::css(grid_template_columns = "5fr 7fr"),
            pcsf_plot_heatmap_ui(
              ns("pcsf_heatmap"),
              caption = paste(
                "PCSF gene modules",
                ""),
              info.text = "",
              height = c("100%", "75vh"),
              width = c("auto", "100%")
            ),
            pcsf_plot_network_ui(
              ns("pcsf_network"),
              caption = paste(
                "PCSF network analysis",
                "Functional analysis of biological networks using Prize-collection Steiner Forest algorithm that determines high-confidence subnetworks."),
              info.text = pcsf_graph_info,
              height = c("100%", "75vh"),
              width = c("auto", "100%")
            )
          )
        )
      )
    )
  )

}
