##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FunctionalInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("fa_contrast"), "Contrast:",
                                   choices = NULL),
                "Select the contrast corresponding to the comparison of interest.",
                placement = "top"
    ),
    withTooltip(shiny::actionLink(ns("fa_options"), "Options",
                                  icon = icon("cog", lib = "glyphicon")),
                "Show/hide advanced options",
                placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.fa_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::checkboxInput(ns("fa_filtertable"),
                               "filter signficant (tables)",
                               FALSE),
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
          functional_plot_kegg_graph_ui(ns("kegg_graph"),
                                        label = "a"),
          functional_plot_kegg_actmap_ui(
            ns("kegg_actmap"),label = "c")
        ),
        div(
          class = "col-md-6",
          functional_table_kegg_table_ui(
            ns("kegg_table"),
            height = c(270, 700),
            width = c("100%", "90%")
          )
        )
      ),
      tags$div(
        HTML(strwrap("<b>(a)</b> <b>KEGG pathway map.</b> Genes are colored according to
             their upregulation (red) or downregulation (blue) in the contrast
             profile. <b>(b)</b> <b>Enrichment table</b> reporting enrichment
             score for each pathway for the selected contrast profile.
             <b>(c)</b> <b>Activation matrix</b> visualizing the activation
             levels of pathways across contrasts."))
      )
    ),
    shiny::tabPanel(
      "GO graph",
      div(
        class = "row",
        div(
          class = "col-md-6",
          functional_plot_go_network_ui(
            ns("GO_network"),
            label = "a"),
          functional_table_go_table_ui(
            ns("GO_table"),
            height = c(270, 700),
            width = c("100%", "90%")
          )
        ),
        div(
          class = "col-md-6",
          functional_plot_go_actmap_ui(ns("GO_actmap"),
                                       label = "c")
        )
      ),
      tags$div(
        HTML(strwrap("<b>(a)</b> <b>Gene Ontology graph.</b> The graph represents
             the enrichment of the GO terms as a tree structure. <b>(b)</b><b> GO
             score table.</b> The score of a GO term is the cumulative score of
             all higher order terms. <b>(c)</b> <b>Activation matrix</b>
             visualizing the enrichment of GO terms across multiple contrast
             profiles."))
      )
    )
  )

  page_ui <- div(
    boardHeader(title = "Pathway Analysis", info_link = ns("fa_info")),
    tabs
  )
  return(page_ui)
}
