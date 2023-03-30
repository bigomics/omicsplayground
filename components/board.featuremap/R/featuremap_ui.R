##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FeatureMapInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),
    ## data set parameters
    withTooltip(shiny::selectInput(ns("sigvar"), "Show phenotype:", choices = NULL, multiple = FALSE),
      "To update",
      placement = "top"
    ),
    shiny::br(),
    shiny::br(),
    withTooltip(shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(), br(),
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        tipifyR(
          shiny::radioButtons(ns("umap_type"), "UMAP datatype:",
            choices = c("logCPM", "logFC"), inline = TRUE
          ),
          "The UMAP can be computed from the normalized log-expression (logCPM), or from the log-foldchange matrix (logFC). Clustering based on logCPM is the default, but when batch/tissue effects are present the logFC might be better."
        ),
        tipifyR(
          shiny::selectInput(ns("filter_genes"), "Show genes:",
            choices = NULL, multiple = FALSE
          ),
          "Filter the genes to highlight on the map."
        ),
        tipifyR(
          shiny::selectInput(ns("filter_gsets"), "Show genesets:",
            choices = NULL, multiple = FALSE
          ),
          "Filter the genesets to highlight on the map."
        )
      )
    )
  )
}

FeatureMapUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  div(
    boardHeader(title = "Cluster features", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Gene",
        bslib::layout_column_wrap(
          width = 1,
          heights_equal = "row",
          bslib::layout_column_wrap(
            width = 1/2,
            featuremap_plot_gene_map_ui(
                ns("geneUMAP"),
                height = c("50vh", 800)
            ),
            featuremap_plot_gene_sig_ui(
                ns("geneSigPlots"),
                height = c("50vh", 800)
            )
          ),
          featuremap_table_gene_map_ui(
              ns("geneUMAP"),
              height = c(400, 800)
          ),
          br()
        )
      ),
      shiny::tabPanel(
        "Geneset",
        bslib::layout_column_wrap(
          width = 1,
          heights_equal = "row",
          bslib::layout_column_wrap(
            width = 1/2,
            featuremap_plot_geneset_map_ui(
                ns("gsetUMAP"),
                height = c("50vh", 800)                
            ),                     
            featuremap_plot_gset_sig_ui(
                ns("gsetSigPlots"),
                height = c("50vh", 800)
            )
          ),
          featuremap_table_geneset_map_ui(
              ns("gsetUMAP"),
              height = c(400, 800)
          ),
          br()
        )
      )
    )
  )
}
