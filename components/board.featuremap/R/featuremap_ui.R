##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FeatureMapInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::br(),
    ## data set parameters
    withTooltip(
      shiny::radioButtons(
        ns("showvar"), "Show:",
        inline = TRUE,
        choices = c("phenotype", "comparisons")
      ),
      "Show gene signatures colored by phenotype conditions (relative expression)
       or by comparisons (logFC).",
      placement = "right", options = list(container = "body")
    ),
    shiny::conditionalPanel(
      "input.showvar == 'phenotype'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("sigvar"), NULL, choices = NULL, multiple = FALSE),
        "Select the phenotype conditions to show in the signatures plot.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("ref_group"), "Reference:", choices = NULL),
        "Reference group. If no group is selected the average is used as reference.",
        placement = "right", options = list(container = "body")
      )
    ),
    shiny::conditionalPanel(
      "input.showvar == 'comparisons'",
      ns = ns,
      withTooltip(
        shiny::selectizeInput(ns("selcomp"), NULL, choices = NULL, multiple = TRUE),
        "Select the comparisons to show in the signatures plot.",
        placement = "top"
      )
    ),
    hr(),
    shiny::conditionalPanel(
      "input.tabs == 'Gene'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("filter_genes"), "Filter genes:",
          choices = NULL, multiple = TRUE
        ),
        "Filter the genes to highlight on the map.",
        placement = "right", options = list(container = "body")
      )
    ),
    shiny::conditionalPanel(
      "input.tabs == 'Gene' && input.filter_genes.includes('<custom>')",
      ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("customlist"), NULL,
          value = NULL,
          rows = 5, placeholder = "Paste your custom gene list"
        ),
        "Paste a custom list of genes to highlight.",
        placement = "bottom"
      )
    ),
    shiny::conditionalPanel(
      "input.tabs == 'Geneset'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("filter_gsets"), "Filter genesets:",
          choices = NULL, multiple = TRUE
        ),
        "Filter the genesets to highlight on the map.",
        placement = "right", options = list(container = "body")
      )
    ),
    shiny::hr(),
    withTooltip(
      shiny::checkboxInput(ns("show_fulltable"), "Show full table", FALSE),
      "Show full table. Not filtered."
    )
    ## shiny::br(),
    ## shiny::br(),
    ## withTooltip(shiny::actionLink(ns("options"), "Advanced options",
    ##   icon = icon("cog", lib = "glyphicon")),
    ##   "Toggle advanced options.", placement = "top"
    ## ),
    ## shiny::br(), br(),
    ## shiny::conditionalPanel(
    ##   "input.options % 2 == 1",
    ##   ns = ns,
    ##   tagList()
    ## )
  )
}

FeatureMapUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  height1 <- c("60%", "70vh")
  height2 <- c("40%", "70vh")
  fullH <- "calc(100vh - 180px)"

  div(
    boardHeader(title = "Cluster features", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Gene",
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          bs_alert("Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). This is in contrast to sample-level clustering which clusters samples by similarity of their expression profile. Feature-level clustering allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers. By coloring the UMAP with the foldchange, one can visually compare the global effect between different conditions."),
          bslib::layout_columns(
            col_widths = c(5, 7),
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
              width = c("auto", "100%")
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
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          bs_alert("Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). This is in contrast to sample-level clustering which clusters samples by similarity of their expression profile. Feature-level clustering allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers. By coloring the UMAP with the foldchange, one can visually compare the global effect between different conditions."),
          bslib::layout_columns(
            col_widths = c(5, 7),
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
