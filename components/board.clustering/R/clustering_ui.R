##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ClusteringInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
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

  fullH <- 800 ## full height of page
  rowH  <- 350

  fullH <- "80vh" ## full height of full page
  rowH  <- "40vh"
  
  div(
    class = "row",
    ## h4("Cluster Samples"),
    boardHeader(title = "Cluster Samples", info_link = ns("board_info")),
    div(
      class = "col-md-7",
      shiny::tabsetPanel(
        id = ns("tabs1"),
        shiny::tabPanel(
          "Heatmap",
          clustering_plot_splitmap_ui(
            id = ns("splitmap"),
            label = "a",
            title = "Clustered Heatmap",
            caption = "Heatmap showing gene expression sorted by 2-way hierarchical clustering.",
            info.text = "In the heatmap, red corresponds to overexpression, blue to underexpression of the gene. Gene clusters are also functionally annotated in the 'Annotate clusters' panel on the right. Hierarchical clustering can be performed on gene level or gene set level expression in which users have to specify it under the {Level} dropdown list. Under the plot settings, users can split the samples by a phenotype class (e.g., tissue, cell type, or gender) using the {split by} setting. In addition, users can specify the top N = (50, 150, 500) features to be used in the heatmap. The ordering of top features is selected under {top mode}. The criteria to select the top features are: SD - features with the highest standard deviation across all the samples,specific - features that are overexpressed in each phenotype class compared to the rest, or by PCA - by principal components. Users can also choose between 'relative' or 'absolute' expression scale. Under the {cexCol} and {cexRow} settings, it is also possible to adjust the cex for the column and row labels.",
            height = c(fullH,"80vh"),
            width = "100%"
          )
        ),
        shiny::tabPanel(
          "PCA/tSNE",
          clustering_plot_clustpca_ui(
            ns("PCAplot"),
            title = "PCA/tSNE plot",
            info.text = "The PCA/tSNE panel visualizes unsupervised clustering obtained by the principal components analysis ( PCA), t-distributed stochastic embedding ( tSNE) or the Uniform Manifold Approximation and Projection (UMAP) algorithms. This plot shows the relationship (or similarity) between the samples for visual analytics, where similarity is visualized as proximity of the points. Samples that are ‘similar’ will be placed close to each other. Users can select from three different clustering approaches (default=t-SNE).",
            caption = "Clustering plot of the dataset samples.",
            label = "",
            height = c(fullH, "70vh"),
            width = c("auto", "100%"),
            parent = ns
          )
        ),
        shiny::tabPanel(
          "Parallel",
          shinyjqui::jqui_sortable(
              bslib::layout_column_wrap(
                 width = 1,                 
                 clustering_plot_parcoord_ui(
                     id = ns("parcoord"),
                     title = "Parallel coordinates",
                     info.text = "The Parallel Coordinates panel displays the expression levels of selected genes across all conditions in the analysis. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustered heatmap. The plot is interactive.",
                     caption = "The interactive Parallel Coordinates plot displays the expression levels of selected genes across all conditions.",
                     label = "a",
                     width = c("100%", "100%"),
                     height = c(rowH, TABLE_HEIGHT_MODAL)
                 ),
                 clustering_table_parcoord_ui(
                     id = ns("parcoord"),
                     title = "Selected genes",
                     info.text = "In this table, users can check mean expression values of features across the conditions for the selected genes.",
                     caption = "Table showing the expression in each sample of the  genes displayed in the Parallel Coordinates.",
                     label = "a",
                     width = c("100%", "100%"),
                     height = c(rowH, TABLE_HEIGHT_MODAL)
                 )
              ) ## layout   
          ) ## sortable
        )
    )),
    div(
      class = "col-md-5",
      shiny::tabsetPanel(
        id = ns("tabs2"),
        shiny::tabPanel(
          "Annotate clusters",
          clustering_plot_clusterannot_ui(
            id = ns("plots_clustannot"),
            title = "Functional annotation of clusters",
            info.text =  "For each cluster, functional annotation terms are ranked by correlating gene sets from more than 42 published reference databases, including well-known databases such as GO, KEGG and Gene Ontology. In the plot settings, users can specify the level and reference set to be used under the Reference level and Reference set settings, respectively.",
            caption = "Top ranked annotation features (by correlation) for each gene cluster as defined in the heatmap.",
            label = "a",
            height = c("45vh", TABLE_HEIGHT_MODAL),
            width = c("100%", "100%")
          ),
          clustering_table_clustannot_ui(
            ns("tables_clustannot"),
            title = "Annotation scores",
            info.text = "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings.",
            caption = "Average correlation values of annotation terms, for each gene cluster.",
            height = c(280, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        shiny::tabPanel(
          "Phenotypes",
          clustering_plot_phenoplot_ui(
            id = ns("clust_phenoplot"),
            title = "Phenotype distribution",
            info.text = "This figure visualizes the distribution of the available phenotype data. The plots show the distribution of the phenotypes superposed on the t-SNE clustering. You can choose to put the group labels in the figure or as separate legend in the plot settings",
            caption = "t-SNE clustering plot of phenotype distribution for the current samples.",
            label = "",
            height = c(fullH, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        ),
        shiny::tabPanel(
          "Feature ranking",
          clustering_plot_featurerank_ui(
            id = ns("clust_featureRank"),
            title = "Feature-set ranking",
            info.text = "Ranked discriminant score for top feature sets. The plot ranks the discriminitive power of the feature set (genes) as a cumulative discriminant score for all phenotype variables. In this way, we can find which feature set (or gene family/set) can explain the variance in the data the best. Correlation-based discriminative power is calculated as the average '(1-cor)' between the groups. Thus, a feature set is highly discriminative if the between-group correlation is low. P-value based scoring is computed as the average negative log p-value from the ANOVA. The 'meta' method combines the score of the former methods in a multiplicative manner.",
            caption = "Ranked discriminant score for top feature sets.",
            label = "",
            height = c(fullH, TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    )
  )
}
