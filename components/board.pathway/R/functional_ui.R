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
    ##----------------------------- WikiPathway -------------------------------
    shiny::tabPanel(
      "WikiPathways",
      bslib::layout_column_wrap(
        width = 1/2,
        style = "height: calc(100vh - 180px);",
        bslib::layout_column_wrap(
          width = 1,
          height = '100%',
          functional_plot_wikipathway_graph_ui(
            ns("wikipathway_graph"),
            title = "WikiPathway",
            info.text = "Genes are colored according to their upregulation (red) or downregulation (blue) in the contrast profile. Each pathway is scored for the selected contrast profile and reported in the table below.",
            caption = "Visualisation of the selected WikiPathway with highlighted up- and down-regulated genes.",
            info.width = "350px",
            label = "a",
            ##height = c(435, TABLE_HEIGHT_MODAL),
            height = c("60%", TABLE_HEIGHT_MODAL),            
            width = c("100%", "100%")            
          ),
          functional_table_wikipathway_ui(
            ns("wikipathway_table"),
            title = "Enrichment table",
            info.text = "Scoring is performed by considering the total number of genes in the pathway (n), the number of genes in the pathway supported by the contrast profile (k), the ratio of k/n, and the ratio of |upregulated or downregulated genes|/k. Additionally, the table contains the list of the upregulated and downregulated genes for each pathway and a q value from the Fisher’s test for the overlap.",
            caption = "Reporting enrichment score for each pathway for the selected contrast profile.",
            label = "b",
            ##height = c(340, TABLE_HEIGHT_MODAL),
            height = c("40%", TABLE_HEIGHT_MODAL),            
            width = c("100%", "100%")
          )
        ),
        functional_plot_wikipathway_actmap_ui(
          ns("wikipathway_actmap"),
          title = "Activation matrix",
          info.text = "The activation matrix facilitates the rapid perusal and detection of the similarities of certain pathways between contrasts. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
          caption = "The matrix allow visualizing the activation levels of pathways (or pathway keywords) across multiple contrast profiles.",
          label = "c",
          ##height = c(790, TABLE_HEIGHT_MODAL)
          height = c("100%", TABLE_HEIGHT_MODAL)            
        )
      )
    ),
    ##----------------------------- REACTOME -------------------------------
    shiny::tabPanel(
      "Reactome",
      bslib::layout_column_wrap(
        width = 1/2,
        style = "height: calc(100vh - 180px);",
        bslib::layout_column_wrap(
          width = 1,
          height = '100%',
          functional_plot_reactome_graph_ui(
            ns("reactome_graph"),
            title = "Reactome pathway",
            info.text = "Reactome (http://reactome.org) is a free, open-source, curated and peer-reviewed knowledge-base of biomolecular pathways. Genes are colored according to their upregulation (red) or downregulation (blue) in the contrast profile. Each pathway is scored for the selected contrast profile and reported in the table below.",
            caption = "Visualisation of the selected Reactome pathway with highlighted up (red) and down (blue) regulated genes.",
            info.width = "350px",
            label = "a",
            height = c("60%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")            
          ),
          functional_table_reactome_ui(
            ns("reactome_table"),
            title = "Enrichment table",
            info.text = "Scoring is performed by considering the total number of genes in the pathway (n), the number of genes in the pathway supported by the contrast profile (k), the ratio of k/n, and the ratio of |upregulated or downregulated genes|/k. Additionally, the table contains the list of the upregulated and downregulated genes for each pathway and a q value from the Fisher’s test for the overlap.",
            caption = "Reporting enrichment score for each pathway for the selected contrast profile.",
            label = "b",
            height = c("40%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        functional_plot_reactome_actmap_ui(
          ns("reactome_actmap"),
          title = "Activation matrix",
          info.text = "The activation matrix facilitates the rapid perusal and detection of the similarities of certain pathways between contrasts. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
          caption = "The matrix allow visualizing the activation levels of pathways (or pathway keywords) across multiple contrast profiles.",
          label = "c",
          height = c("100%", TABLE_HEIGHT_MODAL)
        )
      )
    ),
    ##----------------------------- GO GRAPH -------------------------------
    shiny::tabPanel(
      "GO graph",
      bslib::layout_column_wrap(
        width = 1/2,
        style = "height: calc(100vh - 180px);",
        bslib::layout_column_wrap(
          width = 1,
          height = '100%',
          functional_plot_go_network_ui(
            id = ns("GO_network"),
            title = "Gene Ontology graph",
            info.text = "Gene Ontology (GO) provides a computational representation of the current knowledge about roles of genes for many organisms in terms of molecular functions, cellular components and biological processes. The structure of GO can be described in terms of a graph, where each GO term is a node, and the relationships between the terms are edges between the nodes. GO is loosely hierarchical, with 'child' terms being more specialized than their 'parent' terms. The graph is interactive. You can move the graph and zoom in using the mouse.",
            caption = "The GO graph represents the enrichment of the GO terms as a tree structure.",
            height = c(435, TABLE_HEIGHT_MODAL),
            width = c("60%", "100%"),
            label = "a"
          ),
          functional_table_go_table_ui(
            id = ns("GO_table"),
            title = "GO score table",
            info.text = "The scoring of a GO term is performed by considering the cumulative score of all terms from that term to the root node. That means that GO terms that are supported by higher level terms levels are preferentially scored.",
            caption = "The GO scores table provides statistical information on the relevance of correlated GO terms. ",
            height = c("40%", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          )
        ),
        functional_plot_go_actmap_ui(
          id = ns("GO_actmap"),
          title = "Activation matrix",
          info.text = "From this figure, you can easily detect GO terms that are consistently up/down across conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
          caption = "The GO activation matrix visualizes the activation of GO terms across contrast profiles. ",
          height = c("100%", TABLE_HEIGHT_MODAL),
          width =  c("100%", "100%"),        
          label = "c"
        )
      )
    ), ## end of GO tabpanel
    ##----------------------------- Enrichment Map  -------------------------------
    shiny::tabPanel(
      "Enrichment Map (beta)",
      bslib::layout_column_wrap(
        width = 1,
        style = "height: calc(100vh - 180px);",
        functional_plot_enrichmap_ui(
          id = ns("enrichment_map"),
          title = "Enrichment Map",
          info.text = "The Enrichment Map visualizes enrichments of pathways as an enrichment map, a network representing overlaps among enriched pathways. Nodes represent gene-sets and edges represent mutual overlap; in this way, highly redundant gene-sets are grouped together as clusters, dramatically improving the capability to navigate and interpret enrichment results. Reference: 'Enrichment Map: A Network-Based Method for Gene-Set Enrichment Visualization and Interpretation' Merico D et.al, PLoS One, 2010.",
          caption = "The Enrichment Map visualizes enrichments of pathways as an enrichment map, a network representing overlaps among enriched pathways. Nodes represent gene-sets and edges represent mutual overlap.",
          # height = c("70vh", TABLE_HEIGHT_MODAL),
          height = c("100%", TABLE_HEIGHT_MODAL),            
          width = c("100%", "100%"),
          label = "a"
        )
      )
    ) ## Enrichment map tabpanel

  ) ## end of tabset panel

  page_ui <- div(
    boardHeader(title = "Pathway Analysis", info_link = ns("fa_info")),
    tabs
  )
  return(page_ui)
}
