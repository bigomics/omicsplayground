##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

IntersectionInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("comparisons"), "Contrasts:", choices = NULL, multiple = TRUE),
      "Select the contrasts that you want to compare. If you select N=2 contrast a single scatterplot will be drawn. For N>=3 a scatterplot matrix will be drawn.",
      placement = "top"
    ),
    withTooltip(shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), br(),
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      withTooltip(
        shiny::radioButtons(ns("level"), "Level:",
          choices = c("gene", "geneset"), inline = TRUE
        ),
        "Select feature level: gene or geneset",
        placement = "top"
      ),
      withTooltip(shiny::selectInput(ns("filter"), "Filter:", choices = NULL, multiple = FALSE),
        "Filter features",
        placement = "top"
      ),
      shiny::conditionalPanel(
        "input.filter == '<custom>'",
        ns = ns,
        withTooltip(
          shiny::textAreaInput(ns("customlist"), NULL,
            value = NULL,
            rows = 5, placeholder = "Paste your custom gene list"
          ),
          "Paste a custom list of genes to highlight.",
          placement = "bottom"
        )
      )
    )
  )
}

IntersectionUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Pairwise scatter",
      div(
        class = "row",
        div(
          class = "col-md-7",
          intersection_scatterplot_pairs_ui(
            id = ns("scatterplot"),
            height = c("70vh", TABLE_HEIGHT_MODAL)
          )
        ),
        div(
          class = "col-md-5",
          intersection_plot_venn_diagram_ui(
            id = ns("venndiagram")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Signature clustering",
      div(
        class = "row",
        div(
          class = "col-md-7",
          foldchange_heatmap_ui(
            id = ns("FoldchangeHeatmap"),
            height = c("70vh", TABLE_HEIGHT_MODAL)            
          )
        ),
        div(
          class = "col-md-5",
          contrast_correlation_ui(
            id = ns("ctcorrplot"),
            height = c("70vh", TABLE_HEIGHT_MODAL)            
          )
        )
      ),
    ),
  )


  ## return this div
  div(
    boardHeader(
      title = "Compare signatures",
      info_link = ns("info")
    ),
    tabs
  )
}
