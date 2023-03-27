##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ClusteringInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    #        withTooltip( shiny::actionLink(ns("clust_info"), "Tutorial", icon = shiny::icon("youtube")),
    #                "Show more information and video tutorial about this module."),
    #        shiny::hr(), shiny::br(),
    withTooltip(shiny::selectInput(ns("hm_features"), "Features:", choices = NULL, multiple = FALSE),
      "Select a family of features.",
      placement = "top"
    ),
    withTooltip(
      shiny::radioButtons(ns("hm_clustmethod"), "Layout:",
        c("default", "tsne", "pca", "umap"),
        inline = TRUE
      ),
      "Choose the layout method for clustering to visualise.",
    ),
    shiny::conditionalPanel(
      "input.hm_features == '<custom>'",
      ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("hm_customfeatures"), NULL,
          value = NULL,
          height = "150px", width = "100%",
          rows = 5, placeholder = "Paste your custom gene list"
        ),
        "Paste a custom list of genes to be used as features.",
        placement = "bottom"
      )
    ),
    shiny::conditionalPanel(
      "input.hm_features == '<contrast>'",
      ns = ns,
      tipifyR(
        shiny::selectInput(ns("hm_contrast"), NULL, choices = NULL),
        "Select contrast to be used as signature."
      )
    ),
    withTooltip(shiny::selectInput(ns("hm_group"), "Group by:", choices = NULL),
      "Group the samples by condition.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::selectInput(ns("hm_samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter the relevant samples for the analysis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::actionLink(ns("hm_options"), "Options", icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.hm_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(shiny::selectInput(ns("hm_level"), "Level:", choices = c("gene", "geneset")),
          "Specify the level analysis: gene or geneset level.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(shiny::checkboxInput(ns("hm_filterXY"), "exclude X/Y genes", FALSE),
          "Exclude genes on X/Y chromosomes.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(
          shiny::checkboxInput(
            ns("hm_filterMitoRibo"),
            "exclude mito/ribo genes", FALSE
          ),
          "Exclude mitochondrial (MT) and ribosomal protein (RPS/RPL) genes.",
          placement = "top", options = list(container = "body")
        )
      )
    )
  )
}

ClusteringUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 850 ## full height of page

  div(
    class = "row",
    ## h4("Cluster Samples"),
    boardHeader(title = "Cluster Samples", info_link = ns("clust_info")),
    div(
      class = "col-md-7",
      shiny::tabsetPanel(
        id = ns("tabs1"),
        shiny::tabPanel(
          "Heatmap",
          clustering_plot_hm_splitmap_ui(
            id = ns("hm_splitmap"),
            label = "a",
            height = fullH - 80,
            width = "100%"
          )
        ),
        shiny::tabPanel(
          "PCA/tSNE",
          clustering_plot_clustpca_ui(
            ns("PCAplot"),
            label = "",
            height = c("70vh", "70vh"),
            parent = ns
          ),
          tags$div(
            class = "caption",
            HTML("<b>PCA/tSNE plot.</b> The plot visualizes the similarity in expression of
                          samples as a scatterplot in reduced dimension (2D or 3D).
                          Samples that are similar are clustered near to each other, while samples with different
                          expression are positioned farther away. Groups of samples with similar profiles
                          will appear as <i>clusters</i> in the plot.")
          )
        ),
        shiny::tabPanel(
          "Parallel",
          clustering_plot_table_hm_parcoord_ui(
            id = ns("hm_parcoord"),
            label = "a",
            width = c("100%", "100%"),
            height = c(0.45 * fullH, 600)
          ),
          br(),
          tags$div(
            class = "caption",
            HTML("<b>Parallel Coordinates plot.</b> <b>(a)</b>The Parallel Coordinates plot displays
                            the expression levels of selected genes across all conditions.
                            On the x-axis the experimental conditions are plotted. The y-axis shows the expression level
                            of the genes grouped by condition. The colors correspond to the gene groups as
                            defined by the hierarchical clustered heatmap. <b>(b)</b>
                            Average expression of selected genes across conditions.")
          )
        )
      ),
    ),
    div(
      class = "col-md-5",
      shiny::tabsetPanel(
        id = ns("tabs2"),
        shiny::tabPanel(
          "Annotate clusters",
          clustering_plot_clusterannot_ui(
            id = ns("plots_clustannot"),
            label = "a",
            height = c(360, 600),
            width = c("100%", "100%")
          ),
          clustering_table_clustannot_ui(
            ns("tables_clustannot"),
            height = c(330, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          ),
          tags$div(
            class = "caption",
            HTML("<b>Cluster annotation.</b> <b>(a)</b> Top ranked annotation features (by correlation) for each gene cluster as defined  in the heatmap. <b>(b)</b> Table of average correlation values of annotation features, for each gene cluster.")
          )
        ),
        shiny::tabPanel(
          "Phenotypes",
          clustering_plot_phenoplot_ui(
            id = ns("clust_phenoplot"),
            label = "",
            height = c(fullH - 80, 700)
          ),
          tags$div(
            class = "caption",
            HTML("<b>Phenotype distribution.</b> The plots show the distribution of the phenotypes
                            superposed on the t-SNE clustering. Often, we can expect the t-SNE distribution to be
                            driven by the particular phenotype that is controlled by the experimental condition
                            or unwanted batch effects.")
          )
        ),
        shiny::tabPanel(
          "Feature ranking",
          clustering_plot_featurerank_ui(
            id = ns("clust_featureRank"),
            label = "",
            height = c(fullH - 80, 650),
            width = c("auto", "100%")
          ),
          tags$div(
            class = "caption",
            HTML("<b>Feature-set ranking.</b> Ranked discriminant score for top feature sets.
                            The plot ranks the discriminative power of feature sets (or gene sets) as the
                            cumulative discriminant score for all phenotype variables.")
          )
        )
      )
    )
  )
}
