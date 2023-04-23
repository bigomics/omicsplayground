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
      "Select the phenotype to show in the signatures plot.",
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
          shiny::selectInput(ns("ref_group"), "Reference:", choices = NULL),
          "Reference group. If no group is selected the average is used as reference."
        ),
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

  height1 <- c("calc(60vh - 100px)", "70vh")
  height2 <- c("calc(40vh - 100px)", "70vh")  
  
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
                title = "Gene UMAP",
                info.text = "UMAP clustering of genes colored by standard-deviation of log-expression(sd.X), or standard-deviation of the fold-change (sd.FC). The distance metric is covariance of the gene expression. Genes that are clustered nearby have high covariance.The colour intensity threshold can be set with the Settings icon.",
                caption = "Gene UMAP coloured by level of variance. Shades of red indicate high variance.",
                height = height1,
                width = c("auto", "100%")
            ),
            featuremap_plot_gene_sig_ui(
                ns("geneSigPlots"),
                title = "Gene signatures",
                info.text = "UMAP clustering of genes colored by relative log-expression of the phenotype group. The distance metric is covariance. Genes that are clustered nearby have high covariance.",
                caption = "Gene signature maps coloured by differential expression.",
                height = height1,
                width =  c("auto", "100%") 
            )
          ),
          featuremap_table_gene_map_ui(
              ns("geneUMAP"),
              title = "Gene table",
              info.text = "The contents of this table can be subsetted by selecting (by click&drag) on the Gene map plot.",
              caption = "",
              height = height2,
              width = c("auto", "100%")
          )
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
                title = "Geneset UMAP",
                info.text = "UMAP clustering of genesets colored by standard-deviation of log-expression(sd.X), or standard-deviation of the fold-change (sd.FC). The distance metric is covariance of the geneset expression. Genesets that are clustered nearby have high covariance. The colour intensity threshold can be set with the Settings icon.",
                caption = "Geneset UMAP coloured by level of variance. Shades of red indicate high variance.",
                height = height1,
                width = c("auto", "100%")
            ),                     
            featuremap_plot_gset_sig_ui(
                ns("gsetSigPlots"),
                title = "Geneset signatures",
                info.text = "UMAP clustering of genesets colored by relative log-expression of the phenotype group. The distance metric is covariance. Genesets that are clustered nearby have high covariance.",
                caption = "Geneset signature maps coloured by differential expression.",
                height = height1,
                width = c("auto", "100%")
            )
          ),
          featuremap_table_geneset_map_ui(
              ns("gsetUMAP"),
              title = "Geneset table",
              info.text = "The contents of this table can be subsetted by selecting an area (by click&drag) on the Geneset map plot.",
              caption = "",
              height = height2,
              width = c("auto", "100%")
          )
        )
      )
    )
  )
}
