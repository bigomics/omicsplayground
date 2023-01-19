##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FeatureMapInputs <- function(id) {
  ns <- shiny::NS(id)  ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    ## data set parameters
    withTooltip(shiny::selectInput(ns('sigvar'),'Show phenotype:', choices=NULL, multiple=FALSE),
                "To update", placement = "top"),
    shiny::br(),
    shiny::br(),
    withTooltip( shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
                 "Toggle advanced options.", placement="top"),
    shiny::br(),br(),
    shiny::conditionalPanel(
      "input.options % 2 == 1", ns=ns,
      shiny::tagList(
        tipifyR(shiny::radioButtons(ns('umap_type'),'UMAP datatype:',
                                    choices=c('logCPM','logFC'), inline=TRUE),
                "The UMAP can be computed from the normalized log-expression (logCPM), or from the log-foldchange matrix (logFC). Clustering based on logCPM is the default, but when batch/tissue effects are present the logFC might be better."),
        tipifyR(shiny::selectInput(ns('filter_genes'),'Show genes:',
                                   choices=NULL, multiple=FALSE),
                "Filter the genes to highlight on the map."),
        tipifyR(shiny::selectInput(ns('filter_gsets'),'Show genesets:',
                                   choices=NULL, multiple=FALSE),
                "Filter the genesets to highlight on the map.")
      )
    )
  )
}

FeatureMapUI <- function(id) {
  ns <- shiny::NS(id)  ## namespace
  # NOTE: Output from `featuremap_plot_gene_map_ui` is
  # [[1]] The results from `PlotModuleUI` (what we want to draw on the UI)
  # [[2]] The `ns()` function from that module
  # We need that so we can place the Gene table on a different div()
  ns2 <- featuremap_plot_gene_map_ui(ns('gene_map'))
  # NOTE: Same as above
  ns3 <- featuremap_plot_table_geneset_map_ui(ns('gsetUMAP'))

  div(
    boardHeader(title = "Cluster features", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel("Gene",
                      div(
                        class = "row",
                        div(
                          class = "col-md-6",
                          ns2[[1]]
                        ),
                        div(
                          class = "col-md-6",
                          featuremap_plot_gene_sig_ui(ns('gene_sig'))
                        )
                      ),
                      tableWidget(ns2[[2]]('geneTable'))
      ),
      shiny::tabPanel("Geneset",
                      div(
                        class = "row",
                        div(
                          class = "col-md-6",
                          ns3[[1]]
                        ),
                        div(
                          class = "col-md-6",
                          featuremap_plot_gset_sig_ui(ns('gsetSigPlots'))
                        )
                      ),
                      tableWidget(ns3[[2]]('gsetTable'))
      )
    )
  )
}