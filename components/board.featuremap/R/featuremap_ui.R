##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FeatureMapInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
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
      ),
    ),
    shiny::conditionalPanel(
      "input.showvar == 'comparisons'",
      ns = ns,
      withTooltip(
        shiny::selectizeInput(ns("selcomp"), NULL, choices = NULL, multiple = TRUE),
        "Select the comparisons to show in the signatures plot.",
        placement = "top"
      ),
    ),
    hr(),
    shiny::conditionalPanel(
      "input.tabs == 'Gene'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("filter_genes"), tspan("Filter genes:"),
          choices = NULL, multiple = TRUE
        ),
        "Filter the genes to highlight on the map.",
        placement = "right", options = list(container = "body")
      ),
    ),
    shiny::conditionalPanel(
      "input.tabs == 'Gene' && input.filter_genes.includes('<custom>')",
      ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("customlist"), NULL,
          value = NULL,
          rows = 5
        ),
        "Paste a custom list of genes to highlight.",
        placement = "bottom"
      ),
    ),
    shiny::conditionalPanel(
      "input.tabs == 'Geneset'",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("filter_gsets"), tspan("Filter genesets:"),
          choices = NULL, multiple = TRUE
        ),
        "Filter the genesets to highlight on the map.",
        placement = "right", options = list(container = "body")
      ),
    )
  )
}

FeatureMapUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  height1 <- c("60%", "70vh")
  height2 <- c("40%", "70vh")
  fullH <- "calc(100vh - 181px)"

  div(
    boardHeader(title = "Cluster features", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "Gene",
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          row_heights = list("auto", 1),
          bs_alert("Feature-level UMAP clustering is based on pairwise co-expression between features (genes or genesets). Correlated features are clustered closer together. Feature-level clustering allows to detect gene modules, explore gene neighbourhoods, and identify potential drivers. By coloring the UMAP with the foldchange, one can visually compare the global effect between different conditions."),
          bslib::layout_columns(
            col_widths = 12,
            row_heights = list(1, 1),
            bslib::layout_columns(
              col_widths = c(5, 7),
              featuremap_plot_gene_map_ui(
                ns("geneUMAP"),
                title = "Feature UMAP",
                info.text = "UMAP clustering of features colored by standard-deviation of log-expression or fold-change which can be set using the {color by} plot setting. The color intensity threshold can be set using the {color gamma} plot setting. Additionally it is possible to select the number of labels displayed using the {nr labels} plot setting. By selecting this plot (by drag&drop) the Feature table is subset.",
                info.methods = "Clustering of features performed with Uniform Manifold Approximation and Projection (UMAP) using the top 1000 most varying features, then reduced to 50 PCA dimensions before computing the UMAP embedding. Performed using the uwot R package [1]. The distance metric is covariance of the feature expression. Features that are clustered nearby have high covariance.",
                info.references = list(
                  list(
                    "Melville J (2024) uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.",
                    "https://doi.org/10.32614/CRAN.package.uwot"
                  )
                ),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "Feature UMAP coloured by level of variance. Shades of red indicate high variance.",
                height = height1,
                width = c("auto", "100%")
              ),
              featuremap_plot_gene_sig_ui(
                ns("geneSigPlots"),
                title = "Feature signatures",
                info.text = "UMAP clustering of features colored by relative log-expression of the phenotype group.",
                info.methods = "See Feature UMAP",
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "Feature signature maps coloured by differential expression.",
                height = height1,
                width = c("auto", "100%")
              )
            ),
            featuremap_table_gene_map_ui(
              ns("geneUMAP"),
              title = "Feature table",
              info.text = "The contents of this table can be subsetted by selecting (by click&drag) on the Feature UMAP plot.",
              caption = "",
              height = height2,
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Geneset",
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          row_heights = list("auto", 1),
          bs_alert("Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). This is in contrast to sample-level clustering which clusters samples by similarity of their expression profile. Feature-level clustering allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers. By coloring the UMAP with the foldchange, one can visually compare the global effect between different conditions."),
          bslib::layout_columns(
            col_widths = 12,
            row_heights = c(1, 1),
            bslib::layout_columns(
              col_widths = c(5, 7),
              featuremap_plot_geneset_map_ui(
                ns("gsetUMAP"),
                title = "Geneset UMAP",
                info.text = "UMAP clustering of genesets colored by standard-deviation of log-expression or fold-change which can be set using the {color by} plot setting. The color intensity threshold can be set using the {color gamma} plot setting. Additionally it is possible to select the number of labels displayed using the {nr labels} plot setting. By selecting this plot (by drag&drop) the Geneset table is subset.",
                info.methods = "Clustering of genesets performed with Uniform Manifold Approximation and Projection (UMAP) using the top 1000 most varying genesets, then reduced to 50 PCA dimensions before computing the UMAP embedding. Performed using the uwot R package [1]. The distance metric is covariance of the geneset expression. Genesets that are clustered nearby have high covariance.",
                info.references = list(
                  list(
                    "Melville J (2024) uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.",
                    "https://doi.org/10.32614/CRAN.package.uwot"
                  )
                ),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "Geneset UMAP coloured by level of variance. Shades of red indicate high variance.",
                height = height1,
                width = c("auto", "100%")
              ),
              featuremap_plot_gset_sig_ui(
                ns("gsetSigPlots"),
                title = "Geneset signatures",
                "UMAP clustering of genesets colored by relative log-expression of the phenotype group.",
                info.methods = "See Geneset UMAP",
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "Geneset signature maps coloured by differential expression.",
                height = height1,
                width = c("auto", "100%")
              )
            ),
            featuremap_table_geneset_map_ui(
              ns("gsetUMAP"),
              title = "Geneset table",
              info.text = "The contents of this table can be subsetted by selecting an area (by click&drag) on the Geneset UMAP plot.",
              caption = "",
              height = height2,
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "AI Summary",
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          AiTextCardUI(
            ns("featuremapAISummary"),
            title = "AI Feature Map Summary",
            info.text = "AI-generated summary of the feature map analysis for the selected contrast.",
            caption = "AI-generated feature map summary.",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
  )
}
