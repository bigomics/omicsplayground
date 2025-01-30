##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

IntersectionInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::br(),
    withTooltip(shiny::selectInput(ns("comparisons"), "Contrasts:", choices = NULL, multiple = TRUE),
      "Select the contrasts that you want to compare. If you select N=2 contrast a single scatterplot will be drawn. For N>=3 a scatterplot matrix will be drawn.",
      placement = "top"
    ),
    shiny::br(), shiny::br(),
    withTooltip(shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(),
    # shiny::conditionalPanel(
    #   "input.options % 2 == 1",
    #   ns = ns,
      withTooltip(
        shiny::radioButtons(ns("level"), "Level:",
          choices = c("gene", "geneset"), inline = TRUE
        ),
        "Select feature level: gene or geneset",
        placement = "top"
      ),
      withTooltip(shiny::selectInput(ns("filter"), "Filter:", choices = "<all>", multiple = FALSE),
        "Filter features",
        placement = "top"
      ),
      # shiny::conditionalPanel(
      #   "input.filter == '<custom>'",
      #   ns = ns,
        withTooltip(
          shiny::textAreaInput(ns("customlist"), NULL,
            value = NULL,
            rows = 5
          ),
          "Paste a custom list of genes to highlight.",
          placement = "bottom"
        )
      #)
    #)
  )
}

IntersectionUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- "calc(100vh - 200px)" ## full height of page (minus header)
  fullH.css <- "height: calc(100vh - 130px);" ## full height of page (minus header)
  halfH.css <- "height: calc(50vh - 130px);" ## half height of page

  tabs <- shiny::tabsetPanel(
    id = ns("tabs1"),
    shiny::tabPanel(
      "Pairwise scatter",
      bslib::layout_columns(
        col_widths = c(7, 5),
        height = fullH,
        intersection_scatterplot_pairs_ui(
          id = ns("scatterplot"),
          title = "Scatterplot pairs",
          info.text = "For the selected contrasts, the Pairs panel provides pairwise scatterplots for the differential expression profiles corresponding to multiple contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. When K >= 3 contrasts are selected, the figure shows a KxK scatterplot matrix. When K <= 2, the Pairs panel provides an interactive pairwise scatterplots for the differential expression profiles of the two selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over. Users can also select a number points by selecting points with the mouse, using the box selection or the lasso selection tool. Note that the selected genes will appear in input panel on the left sidebar as selection.",
          caption = "Pairwise scatterplots for two or more differential expression profiles for multiple selected contrasts.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("100%", "100%")
        ),
        bslib::layout_columns(
          col_widths = 12,
          intersection_plot_venn_diagram_ui(
            id = ns("venndiagram"),
            title = "Venn diagram",
            info.text = "The Venn diagram visualizes the number of intersecting genes between the profiles. The list of intersecting genes with further details is also reported in an interactive table below, where users can select and remove a particular contrasts from the intersection analysis.",
            caption = "Venn diagram showing the intersections between the expression profiles of selected contrasts.",
            height = c("65%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          ),
          intersection_table_venn_diagram_ui(
            id = ns("venndiagram"),
            title = "Leading-edge table",
            info.text = "Venn diagram areas can be selected via the settings menu and are represented by corresponding letters (e.g. 'ABC' represents the intersection of contrasts A, B and C). Contrast letter identifiers are provided in the Venn Diagram.",
            caption = "Table of genes in a selected intersection.",
            height = c("35%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        )
      )
    ),
    shiny::tabPanel(
      "Signature clustering",
      bslib::layout_columns(
        col_widths = c(7, 5),
        height = fullH,
        foldchange_heatmap_ui(
          id = ns("FoldchangeHeatmap"),
          title = "Foldchange heatmap",
          info.text = "The Connectivity Heatmap shows the most similar profiles as a heatmap.
                       Contrasts that are similar will be clustered close together.",
          caption = "Signature heatmap visualizing the similarity of all available contrasts
                     with each other for the top most differentially expressed genes.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", 1600)
        ),
        contrast_correlation_ui(
          id = ns("ctcorrplot"),
          title = "Contrast correlation",
          info.text = "Contrasts that are similar will be clustered close together. The numeric value in the cells correspond to the Pearson correlation coefficient between contrast signatures. Red corresponds to positive correlation and blue to negative correlation.",
          caption = "Circle heatmap showing the correlation between all available contrasts.",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width = c("auto", "100%")
        )
      )
    )
  )

  ## return this div
  div(
    boardHeader(title = "Compare signatures", info_link = ns("info")),
    tabs
  )
}
