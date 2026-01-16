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
            choices = c(1000, 2000, 4000),
            selected = 2000
          ),
          shiny::selectInput(ns("power"), "Soft treshold:",
            choices = c("<auto>", 1, 3, 6, 9, 12, 20), selected = 12
          ),
          shiny::selectInput(ns("minmodsize"), "Min. module size",
            choices = c(5, 10, 20, 40, 100), selected = 20
          ),
          shiny::br(),
          shiny::actionButton(
            ns("compute"), "Recompute!",
            icon = icon("running"),
            class = "btn-outline-primary"
          )
        )
      )
    ),
    shinyjs::hidden(
      bslib::accordion(
        id = ns("report_options"),
        open = TRUE,
        bslib::accordion_panel(
          "Report options",
          icon = icon("cog", lib = "glyphicon"),
          wgcna_report_inputs(ns("wgcnaReport"))
        )
      )
    )
  )
}


WGCNA_REFS <- list(
  list("Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008).", "https://doi.org/10.1186/1471-2105-9-559"),
  list("Langfelder P, Luo R, Oldham MC, Horvath S (2011) Is My Network Module Preserved and Reproducible? PLoS Comput Biol 7(1): e1001057.", "https://doi.org/10.1371/journal.pcbi.1001057")
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
          bs_alert(HTML("<b>Module detection.</b> <b>(a)</b> Modules are detected using the dynamic branch cutting approach. <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold. <b>(c)</b> Topological overlap matrix visualized as heatmap. <b>(d)</b> Dimensionality reduction map of features colored by module. <b>(e)</b> Size of WGCNA modules.")),
          bslib::layout_columns(
            col_widths = 12,
            height = "calc(100vh - 181px)",
            row_heights = c(1, 1),
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
                title = "(d) Feature UMAP",
                caption = "UMAP dimensional reduction of features colored by WGCNA module.",
                info.text = "UMAP dimensional reduction of features colored by WGCNA module. The MDS plot shows the separation of modules and their member genes.",
                info.methods = "UMAP was performed on the non-scaled centered expression matrix. MDS analysis was performed by computing the first two principal components of the covariance matrix.",
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
        )
      ),
      shiny::tabPanel(
        "Eigengenes",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Eigengene analysis.</b> The module eigengene of a given module is defined as the first principal component of the standardized expression profiles. <b>(a)</b> Module-trait correlation identifies modules that are significantly associated with the measured traits. <b>(b)</b> Clustering of eigengenes. <b>(c)</b> Clustering of trait vectors. <b>(d)</b> Correlation of eigengene and traits as heatmap, <b>(e)</b> as dendrogram and <b>(f)</b> as graph.")),
          bslib::layout_columns(
            col_widths = c(4, 4, 4, 4, 4, 4),
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
          bs_alert(HTML("<b>Module analysis.</b>  <b>(a)</b> Summary of module. <b>(b)</b> Correlation of module eigengene with traits. <b>(c)</b> Circle network of top hub genes. </b>  <b>(d)</b> Table of importance score to identify 'driver genes' of the module. <b>(e)</b> Plot of gene significance paramters.")),
          bslib::layout_columns(
            height = "100%",
            col_widths = c(3, 9),
            wgcna_html_module_summary_ui(
              id = ns("moduleSummary"),
              title = "(a) Summary",
              info.text = "",
              caption = "Information about the Module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            bslib::layout_columns(
              col_widths = c(7, 5, 7, 5),
              height = "calc(100vh - 181px)",
              wgcna_plot_module_significance_ui(
                ns("moduleSignificance"),
                title = "(b) Trait correlation",
                info.text = "For each module, we also define a trait correlation.",
                caption = "",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              wgcna_plot_correlation_network_ui(
                ns("corGraph"),
                label = "c",
                title = "(c) Circle network of hub genes",
                caption = "Circle network of hub genes",
                info.text = "The circle plot visualizes the connection strengths between top hub genes in the module. Thickness of lines reflect the absolute correlation. Sizes of circles indicate the connectivity of the gene as quantified by its module membership. Higher connected genes are represented by larger circles.",
                info.methods = "The Pearson correlation was computed from the log expression matrix across samples. Module membership was computed using the WGCNA R package as part of the full WGCNA analysis. Plotting was performed using the igraph R package.",
                info.references = WGCNA_REFS,
                info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/modules/mod9_CellProfiling/#wgcna",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              wgcna_table_genes_ui(
                ns("geneTable"),
                title = "(d) Significance table",
                info.text = "Genes in the selected WGCNA module.",
                caption = "Table of genes in the selected module.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              wgcna_plot_membership_v_trait_ui(
                ns("memberTrait"),
                title = "(e) Gene significance",
                info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
                caption = "We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            )
          )
        )
      ),
      shiny::tabPanel(
        "Enrichment",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 181px)",
          bs_alert(HTML("<b>Module Enrichment.</b> <b>(a)</b> Enrichment heatmap of top most enriched genesets in module. <b>(b)</b> Expression heatmap of genes in selected geneset. <b>(c)</b> Functional enrichment of the module calculated using Fisher's exact test. <b>(d)</b> Top enriched genesets in module.")),
          bslib::layout_columns(
            col_widths = c(7, 5, 7, 5),
            height = "calc(100vh - 181px)",
            wgcna_plot_geneset_heatmap_ui(
              ns("genesetHeatmap"),
              title = "(a) Geneset heatmap",
              info.text = "Eigengene correlationheatmap",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_gene_heatmap_ui(
              ns("geneHeatmap"),
              title = "(b) Gene heatmap",
              info.text = "Module gene heatmap",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_table_enrichment_ui(
              ns("enrichTable"),
              title = "(c) Enrichment scores",
              info.text = "In this table, users can check mean expression values of features across the conditions for the selected genes.",
              caption = "Functional enrichment of the module calculated using Fisher's exact test.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_topgenes_ui(
              ns("topgenesPlot"),
              title = "(d) Gene frequency",
              info.text = "Gene frequency in top enriched genesets.",
              caption = "Gene frequency in top enriched genesets",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ) ## end layout_columns (left column)
        ) ## end layout_columns (page)
      ), ## end tabPanel


      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "AI Report",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("⚠️ This page contains AI-generated content. While efforts are made for accuracy, please verify information independently.")),
          bslib::layout_columns(
            col_widths = c(8,4),
            height = "calc(100vh - 180px)",            
            wgcna_html_report_ui(
              ns("wgcnaReport"),
              title = "AI Report",
              caption = "AI-generated summary report",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_report_diagram_ui(
              ns("wgcnaReport"),
              title = "AI Diagram",
              caption = "AI-generated diagram",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )
      
      
    )
  )
}
