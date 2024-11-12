##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ClusteringInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace

  topmodes <- c("sd", "pca", "marker")

  settings_items1 <- tagList(
    withTooltip(shiny::selectInput(ns("selected_phenotypes"), "Show phenotypes:", choices = NULL, multiple = TRUE),
      "Select phenotypes to show in heatmap and phenotype distribution plots.",
      placement = "top"
    ),
    hr(id = ns("pheno_bar")),
    withTooltip(shiny::selectInput(ns("hm_topmode"), "Top mode:", topmodes, width = "100%"),
      "Specify the criteria for selecting top features to be shown in the heatmap.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(ns("hm_ntop"), "Top N:", c(50, 150, 500), selected = 50),
      "Select the number of top features in the heatmap.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::selectInput(ns("hm_clustk"), "K modules:", 1:6, selected = 4),
      "Select the number of gene clusters.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(
        ns("hm_scale"), "Scale:",
        choices = c("relative", "absolute", "BMC"), inline = TRUE
      ),
      "Show relative (i.e. mean-centered), absolute expression values or batch-mean-centered.",
      placement = "right", options = list(container = "body")
    ),
    hr(id = ns("cluster_bar")),
    withTooltip(
      shiny::radioButtons(
        ns("hm_splitby"), "Split samples by:",
        inline = TRUE,
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
      withTooltip(
        shiny::checkboxInput(ns("hm_average_group"), "average by group", FALSE),
        "Average expression values by group."
      )
    ),
    hr(id = ns("spliby_bar")),
    withTooltip(
      shiny::selectInput(ns("hm_samplefilter"), "Filter samples:",
        choices = NULL, multiple = TRUE
      ),
      "Filter the relevant samples for the analysis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::selectInput(ns("hm_features"), tspan("Gene family:"),
        choices = NULL, multiple = FALSE
      ),
      "Select a gene family for filtering which genes to show in the heatmap."
    ),
    shiny::conditionalPanel(
      "input.hm_features == '<custom>'",
      ns = ns,
      withTooltip(
        shiny::textAreaInput(ns("hm_customfeatures"), NULL,
          value = NULL,
          height = "150px", width = "100%",
          rows = 5
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
    withTooltip(
      shiny::actionLink(ns("hm_options"), "Advanced options",
        icon = icon("cog", lib = "glyphicon")
      ),
      "Toggle advanced options.",
      placement = "top"
    ),
    shiny::br(),
    shiny::conditionalPanel(
      "input.hm_options % 2 == 1",
      ns = ns,
      shiny::tagList(
        withTooltip(
          shiny::selectInput(ns("hm_clustmethod"), "Layout:",
            choices = c("tsne", "pca", "umap", "pacmap")
          ),
          "Choose the layout method for clustering plots.",
        ),
        hr(),
        withTooltip(shiny::selectInput(ns("hm_level"), "Level:", choices = c("gene", "geneset")),
          "Specify the level analysis: gene or geneset level.",
          placement = "top", options = list(container = "body")
        ),
        hr(),
        withTooltip(shiny::checkboxInput(ns("hm_filterXY"), tspan("exclude X/Y genes"), FALSE),
          "Exclude genes on X/Y chromosomes.",
          placement = "top", options = list(container = "body")
        ),
        withTooltip(
          shiny::checkboxInput(
            ns("hm_filterMitoRibo"),
            tspan("exclude mito/ribo genes"), FALSE
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

  board_info <- "The Clustering Board performs unsupervised clustering analysis. After having done the QC, it is probably the first way to explore your data. The main purpose is to discover patterns and subgroups in the data, show correlation with known phenotypes, detect outliers, or investigate batch effects."

  heatmap_info <- HTML("The <b>Clustered Heatmap</b> is a powerful 2-way unsupervised hierarchical clustering technique that simultaneously clusters the expression matrix along rows and columns, clustering similar genes and similar samples together. The tree-like dendrogram shows the 'distance' between features and the approximate groups. The column annotations show the correlation with the phenotypes.")

  pca_info <- HTML("<b>Dimensionality reduction</b> is an unsupervised clustering technique that projects the samples into a lower dimensional, here 2D, space. Samples that have similar expression profiles will cluster close together. By coloring the points by condition, we can see which phenotype best explains the clustering.")

  parallel_info <- HTML("The <b>Parallel Coordinates</b> plot is great for visualizing time series or ordered experiments. By grouping samples by time points and showing them sequentially, we can see trends in the expression of groups of genes, or so-called gene modules. The figure is interactive so you can manually order the time points.")

  rowH <- 350
  rowH <- "40vh"
  fullH <- "calc(100vh - 180px)"

  div(
    boardHeader(title = "Cluster Samples", info_link = ns("board_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Heatmap",
        bslib::layout_columns(
          height = fullH,
          col_widths = 12,
          bs_alert(heatmap_info),
          bslib::layout_columns(
            col_widths = c(7, 5),
            height = fullH,
            clustering_plot_splitmap_ui(
              id = ns("splitmap"),
              label = "a",
              title = "Clustered Heatmap",
              caption = "Heatmap showing gene expression sorted by 2-way hierarchical clustering.",
              info.text = "Using the {cexCol} and {cexRow} options it is possible to adjust the font size for the column and row labels. Also, it is possible to select whether to display or not the legend. Gene clusters are functionally annotated in the 'Annotate clusters' panel on the right.",
              info.methods = "The heatmap is generated using the ComplexHeatmap R/Bioconductor package [1] on scaled log-expression values (z-score) using euclidean distance and Ward linkage using the fastcluster R package [2]. The available methods to select the top features are sd (standard deviation) - features with the highest standard deviation across all the samples, marker - features that are overexpressed in each phenotype class compared to the rest, or by PCA - principal component analysis (performed using the irlba R package [3]). In the heatmap, red corresponds to overexpression, blue to underexpression of the gene.",
              info.references = list(
                list(
                  "Gu Z (2016). “Complex heatmaps reveal patterns and correlations in multidimensional genomic data.” Bioinformatics.",
                  "https://doi.org/10.1093/bioinformatics/btw313"
                ),
                list(
                  "Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python.” Journal of Statistical Software, 53(9), 1–18.",
                  "https://doi.org/10.18637/jss.v053.i09"
                ),
                list(
                  "Baglama J (2022). “irlba: Fast Truncated Singular Value Decomposition and Principal Components Analysis for Large Dense and Sparse Matrices”.",
                  "https://doi.org/10.32614/CRAN.package.irlba"
                )
              ),
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
              height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            bslib::layout_columns(
              col_widths = 12,
              clustering_plot_clusterannot_ui(
                id = ns("plots_clustannot"),
                title = "Functional annotation of gene modules",
                info.text = "By setting the {Reference level} and {Reference set} it can be specified the level and reference set to be used. The gene clusters can be visualized in the Heatmap on the left.",
                info.methods = "The clusters are computed using fastcluster R package [1]. For each cluster, functional annotation terms are ranked by correlating gene sets from different public databases including: MSigDB [2], Gene Ontology [3], and Kyoto Encyclopedia of Genes and Genomes (KEGG) [4].",
                info.references = list(
                  list(
                    "Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python.” Journal of Statistical Software, 53(9), 1–18.",
                    "https://doi.org/10.18637/jss.v053.i09"
                  ),
                  list(
                    "Liberzon (2011). Molecular signatures database (MSigDB) 3.0. Bioinformatics, 27(12), 1739-1740.",
                    "https://doi.org/10.1093/bioinformatics/btr260"
                  ),
                  list(
                    "Ashburner (2000). Gene ontology: tool for the unification of biology. Nature genetics, 25(1), 25-29.",
                    "https://doi.org/10.1038/75556"
                  ),
                  list(
                    "Kanehisa (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.",
                    "https://doi.org/10.1093/nar/28.1.27"
                  )
                ),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "Top ranked annotation features (by correlation) for each gene cluster as defined in the heatmap.",
                label = "a",
                height = c("60%", TABLE_HEIGHT_MODAL),
                width = c("100%", "100%")
              ),
              clustering_table_clustannot_ui(
                ns("tables_clustannot"),
                title = "Annotation scores",
                info.text = "Mean correlation values of features in the clusters with respect to the annotation references database. The information displayed corresponds to the functional annotation of gene modules plot.",
                caption = "Average correlation values of annotation terms, for each gene cluster.",
                height = c("40%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            )
          )
        )
      ),
      shiny::tabPanel(
        "PCA/tSNE",
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          bs_alert(HTML(pca_info)),
          bslib::layout_columns(
            col_widths = c(7, 5),
            height = fullH,
            clustering_plot_clustpca_ui(
              ns("PCAplot"),
              title = "Dimensionality reduction",
              info.text = "Using the {Color/label}, {Shape} and {Label} options it is possible to control how the points are colored and shaped (acording to which available phenotypes) and it is possible to control where are the labels located respectively. There is also the option to visualize the three dimensionality reduction techniques at the same time, and the option to visualize the plot in three dimensions.",
              info.methods = "Relationship (or similarity) between the samples for visual analytics, where similarity is visualized as proximity of the points. Three clustering methods are available, t-SNE (using the Rtsne R package [1]), UMAP (using the uwot R package [2]) and PCA (using the irlba R package [3]). Samples that are ‘similar’ will be placed close to each other.",
              info.references = list(
                list(
                  "Krijthe J (2023) Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation.",
                  "https://doi.org/10.32614/CRAN.package.Rtsne"
                ),
                list(
                  "Melville J (2024) uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction.",
                  "https://doi.org/10.32614/CRAN.package.uwot"
                ),
                list(
                  "Baglama J (2022). “irlba: Fast Truncated Singular Value Decomposition and Principal Components Analysis for Large Dense and Sparse Matrices”.",
                  "https://doi.org/10.32614/CRAN.package.irlba"
                )
              ),
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
              caption = "Clustering plot of the dataset samples.",
              label = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%"),
              parent = ns
            ),
            clustering_plot_phenoplot_ui(
              id = ns("clust_phenoplot"),
              title = "Phenotype distribution",
              info.text = "Visualization of the dimensionality reduction plot coloured by the available phenotypes. The group labels can be toggled on the options.",
              info.methods = "See Dimensionality reduction",
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
              caption = "t-SNE clustering plot of phenotype distribution for the current samples.",
              label = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Parallel",
        bslib::layout_columns(
          col_widths = 12,
          height = fullH,
          bs_alert(parallel_info),
          bslib::layout_columns(
            col_widths = c(8, 4),
            bslib::layout_columns(
              col_widths = 12,
              clustering_plot_parcoord_ui(
                id = ns("parcoord"),
                title = "Parallel coordinates",
                info.text = "Control the scale of the values from the options. Also, there is the possibility of averaging by gene module. Arrange the experimental conditions by click&dragging their names on the x-axis. Highlight genes by click and dragging on any y-axis.",
                info.methods = "Expression levels of selected genes across all conditions in the analysis. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustering (performed using the fastcluster R package [1]).",
                info.references = list(
                  list(
                    "Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python.” Journal of Statistical Software, 53(9), 1–18.",
                    "https://doi.org/10.18637/jss.v053.i09"
                  )
                ),
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
                caption = "The interactive Parallel Coordinates plot displays the expression levels of selected genes across all conditions.",
                label = "a",
                width = c("100%", "100%"),
                height = c("100%", TABLE_HEIGHT_MODAL)
              ),
              clustering_table_parcoord_ui(
                id = ns("parcoord"),
                title = "Selected genes",
                info.text = "Expression levels of selected genes across all conditions in the analysis. The selection of genes is performed on the Parallel coordinates plot by click & dragging on any y-axis. If not selection is performed on the plot, all the genes are displayed.",
                caption = "Table showing the expression in each sample of the  genes displayed in the Parallel Coordinates.",
                label = "a",
                width = c("100%", "100%"),
                height = c("100%", TABLE_HEIGHT_MODAL)
              )
            ),
            clustering_plot_genemodule_ui(
              id = ns("genemodule"),
              title = "Module expression",
              info.text = "Series of histograms displaying the overall expression of each module by individual sample.",
              info.methods = "The modules are computed using the fastcluster R package [1].",
              info.references = list(
                list(
                  "Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python.” Journal of Statistical Software, 53(9), 1–18.",
                  "https://doi.org/10.18637/jss.v053.i09"
                )
              ),
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#clustering",
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
