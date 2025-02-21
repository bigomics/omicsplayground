##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WgcnaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    ## data set parameters
    shiny::selectInput(ns("selected_module"), "Select module:", choices = NULL),
    shiny::selectInput(ns("selected_trait"), "Select trait:", choices = NULL),    
    br(),
    bslib::accordion(
      id = ns("compare_accordion"),
      open = FALSE,
      bslib::accordion_panel(
        "Recompute",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::selectInput(ns("ngenes"), tspan("Number genes:"),
            choices = c(1000, 2000, 4000, 8000),
            selected = 2000
          ),
          shiny::selectInput(ns("networktype"), "Network type",
            choices = c("unsigned","signed","signed hybrid"), selected = "signed"
          ),
          shiny::sliderInput(ns("power"), "Power", 1, 20, 6),
          shiny::sliderInput(ns("cutheight"), "Merge cut height", 0.05, 0.8, 0.15, 0.05),
          shiny::selectInput(ns("minmodsize"), "Min. module size",
            choices = c(5, 10, 20, 50, 100),  selected = 20
          ),
          shiny::br(),
          shiny::actionButton(
            ns("compute"), "Recompute!",
            icon = icon("running"),
            class = "btn-outline-primary"
          )
        )
      )
    )
  )
}


WGCNA_REFS <- list(
  list("Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008).", "https://doi.org/10.1186/1471-2105-9-559")
)

WgcnaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),
      shiny::tabPanel(
        "WGCNA",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Module detection.</b> <b>(a)</b> Modules are detected using the dynamic branch cutting approach. <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold. <b>(c)</b> Topological overlap matrix visualized as heatmap. <b>(d)</b> Multi-dimensional scaling colored by WGCNA module. <b>(e)</b> Size of WGCNA modules.")),
          bslib::layout_columns(
            col_widths = c(6, 6),
            height = "35%",
            wgcna_plot_gdendogram_ui(
              ns("geneDendro"),
              label = "a",
              title = "(a) Gene dendrogram and modules",
              caption = "WGCNA gene dendrogram and modules",
              info.text = "Gene modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.",
              info.methods = "",
              info.references = WGCNA_REFS,
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/modules/mod9_CellProfiling/#wgcna",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_s_independence_ui(
              ns("topologyPlots"),
              title = "(b) Scale independence and mean connectivity",
              info.text = "Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis). WGCNA requires a user-defined soft-threshold (or 'power') based on visually observing the above graph; specifically, the scale-freeness. Combined with the negative slope indicated in the mean connectivity graph, we can assume that we have a scale-free network as required by WGCNA.",
              info.methods = "The power was automatically selected using the `pickSoftThreshold()` function with a maximum of power=20. Users can override this default by choosing the power manually on the right. The threshold suggested by the WGCNA authors is the minimum power that reaches a R² threshold of 0.80-0.95.",
              info.references = WGCNA_REFS,
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/modules/mod9_CellProfiling/#wgcna",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            height = "65%",
            wgcna_plot_TOMheatmap_ui(
              ns("TOMplot"),
              title = "(c) TOM heatmap",
              caption = "Topological Overlap Matrix (TOM) visualized as heatmap.",
              label = "c",
              info.text = "The adjacency matrix is transformed in a Topological Overlapping Matrix (TOM) considering the number of neighbors the genes share. From this matrix, a dissimilarity matrix is computed and used to determine the gene modules after dynamic branch cutting.",
              info.methods = "The adjacency matrix is computed as correlation of the standardized expression matrix. The TOM matrix is computed as described in [1]. The heatmap is plotted from a power elevated matrix TOM^7 to accentuate the clusters. ",
              info.references = WGCNA_REFS,
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/modules/mod9_CellProfiling/#wgcna",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_gclustering_ui(
              ns("umap"),
              title = "(d) MDS of features",
              caption = "Multi-dimensional scaling of features colored by WGCNA module.",
              info.text = "Multi-dimensional scaling of features colored by WGCNA module. The MDS plot shows the separation of modules and their member genes.",
              info.methods = "MDS analysis was performed by computing the first two principal components of the covariance matrix.",
              info.references = WGCNA_REFS,
              info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/modules/mod9_CellProfiling/#wgcna",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_module_barplot_ui(
              ns("moduleSize"),
              title = "(e) Module size",
              caption = "Size of modules by number of genes.",
              info.text = "The plot shows the size (number of genes) of each module. The grey module represents genes that show weak correlation and were not assigned to any primary group.",
              info.methods = "The module size is computed as the number of genes within each module.",
              info.references = WGCNA_REFS,
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      shiny::tabPanel(
        "Eigengenes",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Eigengene analysis.</b> The module eigengene of a given module is defined as the first principal component of the standardized expression profiles. <b>(a)</b> Module-trait correlation identifies modules that are significantly associated with the measured traits. <b>(b)</b> Clustering of eigengenes. <b>(c)</b> Clustering of trait vectors. <b>(d)</b> Correlation of eigengene and traits as heatmap, <b>(e)</b> as dendrogram and <b>(f)</b> as graph.")),
          bslib::layout_columns(
            col_widths = c(4,4,4,4,4,4),
            height = "calc(100vh - 181px)",
            wgcna_plot_MTrelationships_ui(
              ns("moduleTrait"),
              title = "(a) Module-Trait relationships",
              info.text = "WGCNA module and trait relationship. Network construction and module detection on the normalized, transposed, and filtered (removal of very low varying genes) gene expression matrix are performed using the WGCNA function blockwiseModules. To assess relationships between gene modules and metadata (eg. phenotypes of interest), Pearson correlation between each module eigengene (i.e., ME, the first PC component of the expression profiles of all genes mapped within the module) and the trait of interest is computed. These correlation coefficients (ranging from -1 to 1, representing fully negative and fully positive correlation, respectively), are conveniently displayed in a heatmap using the WGCNA labeledheatmap function. Positive correlation suggests that genes in the module may be upregulated as the phenotype/trait value increases. Therefore, ME-trait relationship indicates how strong a gene module is associated with a phenotype.",
              caption = "Module-trait analysis identifies modules that are significantly associated with the measured traits by quantifying the association as the correlation of the eigengenes with external traits.",
              label = "a",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_sampledendrogram_ui(
              ns("sampleDendrogram2"),
              title = "(b) Sample dendrogram + eigengenes",
              info.text = "Sample dendrogram and eigen genes",
              caption = "Sample dendrogram and eigen genes",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_sampledendrogram_ui(
              ns("sampleDendrogram"),
              title = "(c) Sample dendrogram + traits",
              info.text = "sampleClustering",
              caption = "Sample dendrogram of traits",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_eigengene_heatmap_ui(
              ns("eigenHeatmap"),
              title = "(d) Eigengene correlation heatmap",
              info.text = "Eigengene correlationheatmap",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),           
            wgcna_plot_eigengene_clustering_ui(
              ns("eigenClustering"),
              title = "(e) Eigengene dendrogram",
              info.text = "eigenClustering",
              caption = "Cluster dendrogram of Eigengenes",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_module_graph_ui(
              ns("moduleGraph"),
              title = "(f) Module graph",
              info.text = "WGCNA module graph.",
              caption = "Graph network of WGCNA modules.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      
      shiny::tabPanel(
        "Modules",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Module analysis.</b>  <b>(a)</b> Correlation of module eigengene with traits. <b>(b)</b> Partial correlation network of genes most correlated to the eigengene. </b> <b>(c)</b> Module membership (MM) as the correlation of the genes with the module eigengene. <b>(d)</b> Quantitative measures of significance such as Module Membership (MM), gene Trait Significance (GS), foldChange and network centrality. <b>(e)</b> Importance score as the geometric mean of the significance measures to identify 'driver genes' of the module.")),
          bslib::layout_columns(
            col_widths = c(3,4,5,5,7),
            wgcna_plot_module_significance_ui(
              ns("moduleSignificance"),
              title = "(a) Trait correlation",
              info.text = "For each module, we also define a trait correlation.",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),

            wgcna_plot_correlation_network_ui(
              ns("corGraph"),
              label = "c",
              title = "(b) Partial correlation network",
              info.text = "Partial correlation graph centered on module eigen-gene with top most correlated features. Green edges correspond to positive (partial) correlation, red edges to negative (partial) correlation. Width of the edges is proportional to the correlation strength of the gene pair. The regularized partial correlation matrix is computed using the 'graphical lasso' (Glasso) with BIC model selection.",
              caption = "Module enrichment plot of top most enriched genesets.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),

            wgcna_plot_module_membership_ui(
              ns("eigenCorrelation"),
              title = "(c) Module membership",
              info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
              caption = "For each module, we also define
                a quantitative measure of 'module membership' (MM) as the correlation of the module eigengene and the gene
                expression profile. This allows us to quantify the similarity of all genes to every module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),

            wgcna_plot_membership_v_trait_ui(
              ns("intraScatter"),
              title = "(d) Gene significance",
              info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
              caption = "We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            
            wgcna_table_genes_ui(
              ns("geneTable"),
              title = "(e) Significance table",
              info.text = "Genes in the selected WGCNA module.",
              caption = "Table of genes in the selected module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
            
          )
        )
      ),

      shiny::tabPanel(
        "Enrichment",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Module Enrichment.</b> <b>(a)</b> Module enrichment plot of top most enriched genesets. <b>(b)</b> Functional enrichment of the module calculated using Fisher's exact test.")),
          bslib::layout_columns(
            col_widths = c(4, 8, 12),
            height = "calc(100vh - 181px)",

            wgcna_plot_enrichment_ui(
              ns("enrichPlot"),
              title = "(a) Enrichment barplot",
              info.text = "Functional enrichment of the selected module.",
              caption = "Module enrichment plot of top most enriched genesets.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            
            wgcna_plot_geneset_heatmap_ui(
              ns("genesetHeatmap"),
              title = "(b) Geneset heatmap for module",
              info.text = "Eigengene correlationheatmap",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),           
            
            wgcna_table_enrichment_ui(
              ns("enrichTable"),
              title = "(c) Enrichment table",
              info.text = "In this table, users can check mean expression values of features across the conditions for the selected genes.",
              caption = "Functional enrichment of the module calculated using Fisher's exact test.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )

            ) ## end layout_columns (left column)
        ) ## end layout_columns (page)
      )  ## end tabPanel
    )
  )
}
