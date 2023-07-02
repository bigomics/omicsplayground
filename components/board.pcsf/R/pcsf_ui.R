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
        "Color nodes by:",
        choices = c("gene.cluster", "contrast"),
        selected = "gene.cluster",
        inline = TRUE
      ),
      "Choose how to color the nodes",
      placement = "right"
    ),
    conditionalPanel(
      "input.colorby == 'contrast'",
      ns = ns,
      withTooltip(
        selectInput(ns("contrast"), NULL, choices = NULL, multiple = FALSE),
        "Select contrast.",
        placement = "right"
      )
    ),
    hr(),
    withTooltip(
      radioButtons(
        ns("highlightby"),
        "Highlight labels by:",
        choices = c("none", "FC", "centrality"),
        selected = "centrality",
        inline = TRUE
      ),
      "Highlight labels by scaling with selection.",
      placement = "top"
    ),
    hr(),
    withTooltip(    
      shiny::sliderInput(ns("pcsf_beta"), "Solution size:", -5, 5, 0, 0.5),
      "Select contrast.",
      placement = "right"
    ),
    br(),
    withTooltip(
      actionLink(ns("adv_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::conditionalPanel(
      "input.adv_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(shiny::checkboxInput(ns("check1"), "check1", TRUE),
          "Some check.",
          placement = "top")
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
    ## "Hello",
    ##bs_alert(pcsf_module_info),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "PCSF network",
        bslib::layout_column_wrap(
          width = 1,
          height = "calc(100vh - 190px)",
          heights_equal = "row",            
          div(bs_alert(pcsf_module_info), style="margin-bottom:15px;"),
          bslib::layout_column_wrap(
            width = 1,
            height = "calc(100vh - 190px)",                   
            style = htmltools::css(grid_template_columns = "5fr 7fr"),
            pcsf_plot_heatmap_ui(
              id = ns("pcsf_heatmap"),
              caption = "PCSF gene modules",
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
