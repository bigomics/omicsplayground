##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ClusteringInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace

  settings_items1 <- tagList(
    withTooltip(shiny::selectInput(ns("selected_phenotypes"), "Show phenotypes:", choices = NULL, multiple = TRUE),
      "Select phenotypes to show in heatmap and phenotype distribution plots.",
      placement = "top"
    ),
    hr(),
    withTooltip(
      shiny::radioButtons(
        ns("hm_splitby"), "Split samples by:", inline = TRUE,
        #
        choices = c("none", "phenotype", "gene")
      ),
      "Split the samples by phenotype or expression level of a gene.",
      placement = "right", options = list(container = "body")
    ),
    shiny::conditionalPanel(
      "input.hm_splitby != 'none'",
      ns = ns,
      withTooltip(shiny::selectInput(ns("hm_splitvar"), NULL, choices = ""),
        "Specify phenotype or gene for splitting the columns of the heatmap.",
        placement = "right", options = list(container = "body")
      ),
      withTooltip(shiny::checkboxInput(ns("hm_average_group"), "average by group", FALSE),
        "Average expression values by group."
      )
    ),
    hr(),    
    withTooltip(
      shiny::selectInput(ns("hm_samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter the relevant samples for the analysis.",
      placement = "top", options = list(container = "body")
    ),
     withTooltip(shiny::selectInput(ns("hm_features"), "Gene family:",
                                    choices = NULL, multiple = FALSE),
      "Select a gene family for filtering which genes to show in the heatmap."
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
      withTooltip(
        shiny::selectInput(ns("hm_contrast"), NULL, choices = NULL),
        "Select contrast to be used as signature.",
        placement = "right", options = list(container = "body")
      )
    )
  )

  
  bigdash::tabSettings(
    settings_items1,
    br(),
    withTooltip(shiny::actionLink( ns("hm_options"), "Advanced options",
                                  icon = icon("cog", lib = "glyphicon")),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.hm_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::radioButtons(ns("hm_clustmethod"), "Layout:",
            c("tsne", "pca", "umap"),
            inline = TRUE
          ),
          "Choose the layout method for clustering plots.",
        ),
        hr(),
        withTooltip(shiny::selectInput(ns("hm_level"), "Level:", choices = c("gene", "geneset")),
          "Specify the level analysis: gene or geneset level.",
          placement = "top", options = list(container = "body")
        ),
        hr(),        
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

  board_info = "The Clustering Board performs unsupervised clustering analysis. After having done the QC, it is probably the first way to explore your data. The main purpose is to discover patterns and subgroups in the data, show correlation with known phenotypes, detect outliers, or investigate batch effects."

  heatmap_info = HTML("The <b>Clustered Heatmap</b> is a powerful 2-way unsupervised hierarchical clustering technique that simultaneously clusters the expression matrix along rows and columns, clustering similar genes and similar samples together. The tree-like dendrogram shows the 'distance' between features and the approximate groups. The column annotations show the correlation with the phenotypes.")

  pca_info = HTML("<b>Dimensionality reduction</b> is an unsupervised clustering technique that projects the samples into a lower dimensional, here 2D, space. Samples that have similar expression profiles will cluster close together. By coloring the points by condition, we can see which phenotype best explains the clustering.")

  parallel_info = HTML("The <b>Parallel Coordinates</b> plot is great for visualizing time series or ordered experiments. By grouping samples by time points and showing them sequentially, we can see trends in the expression of groups of genes, or so-called gene modules. The figure is interactive so you can manally order the time points.")
  
  rowH  <- 350
  rowH  <- "40vh"
  fullH = "calc(100vh - 180px)"

  div(
    boardHeader(title = "Cluster Samples", info_link = ns("board_info")),
    #
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Heatmap",
        bslib::layout_column_wrap(
          width = 1,
          height = fullH,
          heights_equal = "row",          
          bs_alert(heatmap_info),
          bslib::layout_column_wrap(
            width = 1,
            height = fullH,
            style = htmltools::css(grid_template_columns = "7fr 5fr"),
            clustering_plot_splitmap_ui(
              id = ns("splitmap"),
              label = "a",
              title = "Clustered Heatmap",
              caption = "Heatmap showing gene expression sorted by 2-way hierarchical clustering.",
              info.text = "In the heatmap, red corresponds to overexpression, blue to underexpression of the gene. Gene clusters are also functionally annotated in the 'Annotate clusters' panel on the right. Hierarchical clustering can be performed on gene level or gene set level expression in which users have to specify it under the {Level} dropdown list. Under the plot settings, users can split the samples by a phenotype class (e.g., tissue, cell type, or gender) using the {split by} setting. In addition, users can specify the top N = (50, 150, 500) features to be used in the heatmap. The ordering of top features is selected under {top mode}. The criteria to select the top features are: SD - features with the highest standard deviation across all the samples,marker - features that are overexpressed in each phenotype class compared to the rest, or by PCA - by principal components. Users can also choose between 'relative' or 'absolute' expression scale. Under the {cexCol} and {cexRow} settings, it is also possible to adjust the cex for the column and row labels.",
              height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
              #
              width = c("auto","100%")
            ),
            bslib::layout_column_wrap(
              width = 1,
              #
              equal_heights = "row",
              clustering_plot_clusterannot_ui(
                id = ns("plots_clustannot"),
                title = "Functional annotation of clusters",
                info.text =  "For each cluster, functional annotation terms are ranked by correlating gene sets from more than 42 published reference databases, including well-known databases such as GO, KEGG and Gene Ontology. In the plot settings, users can specify the level and reference set to be used under the Reference level and Reference set settings, respectively.",
                caption = "Top ranked annotation features (by correlation) for each gene cluster as defined in the heatmap.",
                label = "a",
                #
                height = c("60%", TABLE_HEIGHT_MODAL),              
                width = c("100%", "100%")
              ),
              clustering_table_clustannot_ui(
                ns("tables_clustannot"),
                title = "Annotation scores",
                info.text = "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings.",
                caption = "Average correlation values of annotation terms, for each gene cluster.",
                #
                height = c("40%", TABLE_HEIGHT_MODAL),              
                width = c("auto", "100%")
              )
            )
          )
        )
      ),
      shiny::tabPanel(
        "PCA/tSNE",
        bslib::layout_column_wrap(
          width = 1,
          height = fullH,
          heights_equal = "row",
          bs_alert(HTML(pca_info)),          
          bslib::layout_column_wrap(
            width = 1,
            height = fullH,
            style = htmltools::css(grid_template_columns = "7fr 5fr"),
            clustering_plot_clustpca_ui(
              ns("PCAplot"),
              title = "PCA/tSNE plot",
              info.text = "The PCA/tSNE panel visualizes unsupervised clustering obtained by the principal components analysis (PCA), t-distributed stochastic embedding (tSNE) or the Uniform Manifold Approximation and Projection (UMAP) algorithms. This plot shows the relationship (or similarity) between the samples for visual analytics, where similarity is visualized as proximity of the points. Samples that are ‘similar’ will be placed close to each other. Users can select from three different clustering approaches (default=t-SNE).",
              caption = "Clustering plot of the dataset samples.",
              label = "",
              #
              height = c("100%", TABLE_HEIGHT_MODAL),            
              width = c("auto", "100%"),
              parent = ns
            ),
            clustering_plot_phenoplot_ui(
              id = ns("clust_phenoplot"),
              title = "Phenotype distribution",
              info.text = "This figure visualizes the distribution of the available phenotype data. The plots show the distribution of the phenotypes superposed on the t-SNE clustering. You can choose to put the group labels in the figure or as separate legend in the plot settings",
              caption = "t-SNE clustering plot of phenotype distribution for the current samples.",
              label = "",
              #
              height = c("100%", TABLE_HEIGHT_MODAL),            
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Parallel",
        bslib::layout_column_wrap(
          width = 1,
          height = fullH,
          heights_equal = "row",          
          bs_alert(parallel_info),
          bslib::layout_column_wrap(          
            width = 1,
            heights_equal = "row",
            style = htmltools::css(grid_template_columns = "8fr 4fr"),
            bslib::layout_column_wrap(
              width = 1,
              heights_equal = "row",
              clustering_plot_parcoord_ui(
                id = ns("parcoord"),
                title = "Parallel coordinates",
                info.text = "The Parallel Coordinates panel displays the expression levels of selected genes across all conditions in the analysis. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustered heatmap. The plot is interactive.",
                caption = "The interactive Parallel Coordinates plot displays the expression levels of selected genes across all conditions.",
                label = "a",
                width = c("100%", "100%"),
                #
                height = c("100%", TABLE_HEIGHT_MODAL)              
              ),
              clustering_table_parcoord_ui(
                id = ns("parcoord"),
                title = "Selected genes",
                info.text = "In this table, users can check mean expression values of features across the conditions for the selected genes.",
                caption = "Table showing the expression in each sample of the  genes displayed in the Parallel Coordinates.",
                label = "a",
                width = c("100%", "100%"),
                #
                height = c("100%", TABLE_HEIGHT_MODAL)              
              )
            ),
            clustering_plot_genemodule_ui(
              id = ns("genemodule"),
              title = "Module expression",
              info.text = "",
              caption = "",
              width = c("100%", "100%"),
              height = c("calc(100vh - 200px)", TABLE_HEIGHT_MODAL)
            )
          )
        )
      )
    )
  )
}
