##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

WgcnaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("selected_module"), "select module", choices = NULL),
    shiny::br(),
    shiny::br(),
    shiny::br(),    
    shiny::actionLink(ns("options"), "Recompute", icon = icon("cog", lib = "glyphicon")),
    shiny::br(), 
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        shiny::selectInput(ns("ngenes"), tspan("Number genes:"),
          choices = c(500, 1000, 2000, 4000, 8000),
          selected = 1000
        ),
        shiny::selectInput(ns("networktype"), "Network type",
          choices = c("unsigned","signed","signed hybrid"), selected = "signed"
        ),
        shiny::selectInput(ns("tomtype"), "TOM type",
          choices = c("signed","unsigned","signed Nowick","none"), selected = "signed"
        ),
        shiny::sliderInput(ns("power"), "Power", 1, 20, 6),
        shiny::sliderInput(ns("deepsplit"), "Deepsplit", 0, 4, 2),        
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
}

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
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>WGCNA module detection.</b> <b>(a)</b> Modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles. <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold. <b>(c)</b> Topological overlap matrix visualized as heatmap. <b>(d)</b> Dimensionality reduction maps colored by WGCNA module. <b>(e)</b> Graph network of WGCNA modules.")),
          bslib::layout_columns(
            col_widths = c(6, 6),
            height = "35%",
            wgcna_plot_gdendogram_ui(
              ns("geneDendro"),
              label = "a",
              caption = "WGCNA gene dendrogram and gene modules",
              info.text = "Gene modules are detected as branches of the resulting cluster tree using the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene. The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_s_independence_ui(
              ns("topologyPlots"),
              title = "Scale independence and mean connectivity",
              info.text = " Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis).",
              caption = "Scale independence and mean connectivity plots to determine the soft threshold.",
              label = "b",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            height = "65%",
            wgcna_plot_TOMheatmap_ui(
              ns("TOMplot"),
              title = "TOM heatmap",
              info.text = "WGCNA Topological Overlap Matrix (TOM) heatmap.",
              caption = "Topological overlap matrix visualized as heatmap.",
              label = "c",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_gclustering_ui(
              ns("umap"),
              title = "Gene clustering",
              info.text = "UMAP visualization of TOM correlation of genes.",
              caption = "Dimensionality reduction maps colored by WGCNA module.",
              label = "d",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_module_graph_ui(
              ns("moduleGraph"),
              title = "Module graph",
              info.text = "WGCNA module graph.",
              caption = "Graph network of WGCNA modules.",
              label = "e",
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
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>WGCNA functional analysis.</b> <b>(a)</b> Module-trait analysis identifies modules that are significantly associated with the measured traits by quantifying the association as the correlation of the eigengenes with external traits. <b>(b)</b> Partial correlation network of genes most correlated to the eigengene. <b>(c)</b> Module enrichment plot of top most enriched genesets. <b>(d)</b> Table of genes in the selected module. <b>(e)</b> Functional enrichment of the module calculated using Fisher's exact test.")),
          bslib::layout_columns(
            col_widths = c(4, 4, 4),
            height = "60%",
            wgcna_plot_MTrelationships_ui(
              ns("moduleTrait"),
              title = "Module-Trait relationships",
              info.text = "WGCNA module and trait relationship. Network construction and module detection on the normalized, transposed, and filtered (removal of very low varying genes) gene expression matrix are performed using the WGCNA function blockwiseModules. To assess relationships between gene modules and metadata (eg. phenotypes of interest), Pearson correlation between each module eigengene (i.e., ME, the first PC component of the expression profiles of all genes mapped within the module) and the trait of interest is computed. These correlation coefficients (ranging from -1 to 1, representing fully negative and fully positive correlation, respectively), are conveniently displayed in a heatmap using the WGCNA labeledheatmap function. Positive correlation suggests that genes in the module may be upregulated as the phenotype/trait value increases. Therefore, ME-trait relationship indicates how strong a gene module is associated with a phenotype.",
              caption = "Module-trait analysis identifies modules that are significantly associated with the measured traits by quantifying the association as the correlation of the eigengenes with external traits.",
              label = "a",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_correlation_network_ui(
              ns("corGraph"),
              label = "c",
              title = "Correlation network",
              info.text = "Partial correlation graph centered on module eigen-gene with top most correlated features. Green edges correspond to positive (partial) correlation, red edges to negative (partial) correlation. Width of the edges is proportional to the correlation strength of the gene pair. The regularized partial correlation matrix is computed using the 'graphical lasso' (Glasso) with BIC model selection.",
              caption = "Module enrichment plot of top most enriched genesets.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_enrichment_ui(
              ns("enrichPlot"),
              label = "c",
              title = "Module Enrichment Plot",
              info.text = "Functional enrichment of the selected module.",
              caption = "Module enrichment plot of top most enriched genesets.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_columns(
            col_widths = c(4, 8),
            height = "40%",
            wgcna_table_genes_ui(
              ns("geneTable"),
              label = "d",
              title = "Module genes",
              info.text = "Genes in the selected WGCNA module.",
              caption = "Table of genes in the selected module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_table_enrichment_ui(
              ns("enrichTable"),
              label = "e",
              title = "Module Enrichment Table",
              info.text = "In this table, users can check mean expression values of features across the conditions for the selected genes.",
              caption = "Functional enrichment of the module calculated using Fisher's exact test.",
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
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>WGCNA eigengene analysis.</b> <b>(a)</b> It is often interesting to visualizing the network of eigengenes and study the relationships among the found modules. One can use the eigengenes as represen- tative profiles and quantify module similarity by eigengene correlation. <b>(b)</b> For each module, we also define a quantitative measure of 'module membership' (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes to every module.")),
          bslib::layout_columns(
            col_widths = c(3, 9),
            height = "calc(100vh - 180px)",
            wgcna_plot_eigengene_clustering_ui(
              ns("eigenClustering"),
              title = "Eigengene clustering",
              info.text = "eigenClustering",
              label = "a",
              caption = "It is often interesting to visualizing the network of eigengenes
                and study the relationships among the found modules. One can use the eigengenes as represen- tative profiles
                and quantify module similarity by eigengene correlation.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_module_membership_ui(
              ns("eigenCorrelation"),
              label = "a",
              title = "Module membership (eigengene correlation)",
              info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
              caption = "For each module, we also define
                a quantitative measure of 'module membership' (MM) as the correlation of the module eigengene and the gene
                expression profile. This allows us to quantify the similarity of all genes to every module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      shiny::tabPanel(
        "Intramodular",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>WGCNA intramodular analysis.</b> We quantify associations of individual genes with our trait of interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module membership in interesting modules.")),
          bslib::layout_columns(
            col_widths = c(3, 9),
            height = "calc(100vh - 180px)",
            wgcna_plot_heatmap_membership_ui(
              ns("intraHeatmap"),
              label = "a",
              title = "Membership-trait heatmap",
              info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
              caption = "We quantify associations of individual genes with our trait of
                  interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between
                  the gene and the trait. For each module, we also define a quantitative measure of module membership MM
                  as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures,
                  we can identify genes that have a high significance for weight as well as high module membership in interesting modules.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            wgcna_plot_membership_v_trait_ui(
              ns("intraScatter"),
              label = "b",
              title = "Membership vs. trait correlation",
              info.text = "For each module, we also define a quantitative measure of module membership (MM) as the correlation of the module eigengene and the gene expression profile. This allows us to quantify the similarity of all genes on the array to every module.",
              caption = "We quantify associations of individual genes with our trait of
                  interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between
                  the gene and the trait. For each module, we also define a quantitative measure of module membership MM
                  as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures,
                  we can identify genes that have a high significance for weight as well as high module membership in interesting modules.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )
    )
  )
}
