##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

MultiWGCNA_Inputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::selectInput(ns("phenotype"), "Phenotype", choices = NULL),
    shiny::selectInput(ns("condition"), "Condition on phenotype:", choices = NULL),
    shiny::selectInput(ns("layers"), "Layers:", choices = NULL, multiple = TRUE),
    shiny::conditionalPanel(
      "input.layers && input.layers.indexOf('gset') > -1",
      ns = ns,
      shiny::selectInput(ns("gsfilter"), "Geneset filter:", choices = NULL)
    ),
    shiny::selectInput(ns("module"), "Module:", choices = NULL, multiple = FALSE),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("wgcna_options"),
      open = TRUE,
      bslib::accordion_panel(
        "WGCNA options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::selectInput(ns("power"), "Soft power:",
            choices = c("<auto>", 1, 3, 6, 9, 12, 20), selected = 12
          ),
          shiny::selectInput(ns("deepsplit"), "Deepsplit:", choices = c(0:4), 2),
          shiny::selectInput(ns("ngenes"), "Max. features:",
            choices = c(1000, 2000, 4000),
            2000
          ),
          shiny::selectInput(ns("minmodsize"), "Min. module size",
            choices = c(5, 10, 20, 40, 100), selected = 10
          ),
          shiny::checkboxInput(ns("consensus"), "use consensus", FALSE),
          shiny::br(),
          shiny::actionButton(ns("compute"), "Compute",
            size = "xs",
            icon = icon("refresh")
          )
        )
      )
    ),
    shiny::br(),
    shinyjs::hidden(
      bslib::accordion(
        id = ns("lasagna_options"),
        open = FALSE,
        bslib::accordion_panel(
          "Lasagna options",
          icon = icon("cog", lib = "glyphicon"),
          multiwgcna_plot_lasagna_inputs(ns("multiwgcnaLasagna"))
        )
      )
    )
  )
}


MULTIWGCNA_INFO <- "The <b>Multi-partite graph</b> shows the correlation structure between multiple sets of features. The color of the edges correspond to positive (purple) and negative (yellow) correlation. Thicker edges mean higher correlation. The sizes of the circles represent the page-rank centrality of the feature. The log2FC is indicated for the chosen comparison. The node color corresponds to up (red) and down (blue) regulation."

MultiWGCNA_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Multiomics WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Dendrograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>Multiomics WGCNA</b> is a generalization of WGCNA for integratiing multi-omics where WGCNA is performed for each layer separately. Integration is performed by computing the module correlation across layers using LASAGNA.")),
          bslib::layout_columns(
            col_widths = c(12),
            # height = "calc(100vh - 180px)",
            height = "100vh",
            multiwgcna_plot_dendrograms_ui(
              ns("multiwgcnaDendro"),
              title = "Dendrograms and module colors",
              caption = "Dendrogram of hierarchical clustering of feature co-expression patterns and modules assignment.",
              info.text = "The gene dendrogram is a highly used visualization in WGCNA. It provides a combined visual summary of the features' hierarchical clustering dendrogram and the module color assignments. It shows that co-expressed features are grouped into common biologically meaningful modules. The tree of the dendrogram reflects the co-expression similarity of the features. Each leaf (end point) of the tree corresponds to a single feature. The branching pattern shows the correlation structure: features that are highly co-expressed (strongly correlated across samples) are clustered together and connected up until lower branch points. Features that are less correlated are joined higher up in the tree. The dendrogram is typically built using hierarchical clustering on the Topological Overlap Matrix (TOM), which captures both direct and indirect co-expression relationships. The title of the dendrogram denotes the data type (eg., gx for gene expression; px for proteomics) and the power value (set in the available WGCNA options on the right). Other key WGCNA parameters can be modified and WGCNA can be recomputed by clicking at the 'Compute' button.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          ),
          bslib::layout_columns(
            col_widths = c(12),
            # height = "calc(100vh - 180px)",
            height = "100vh",
            multiwgcna_plot_power_ui(
              ns("multiwgcnaPower"),
              title = "Scale and connectivity plots",
              caption = "Scale and connectivity plots.",
              info.text = "The SFT model fit scatter plot shows the soft threshold power vs the signed scale-free topology fit index. For each datatype, the slope of the log-log connectivity plot vs. frequency is extracted and its opposite sign is computed. This value is then multiplied by the scale-free topology fit index (R^2), with higher values indicating good scale-free structure. The resulting value is plotted against the SFT power value. In WGCNA, features are connected based on their correlation across samples. Thus, the mean connectivity is a measure of overall network density, i.e., how connected the features are, on average.<br><br>The mean connectivity graph plots soft threshold power vs. mean connectivity value, for each data type. It shows the average number (or strenght) of connection per each feature, at any given value of power. A decreasing trend is an expected behavior: the mean connectivity decreases as the power increases, suggesting that the network becomes less dense (and therefore the adjacency matrix becomes sparser).<br><br>The soft threshold power for network sparsity can be changed under 'WGCNA option' and clicking at the 'Compute' button. When 'auto', our in-house IQR method is used. IQR maximizes variation of heights of the dendrogram, promoting seperability between groups. It tests a range of power values. For each power: (i) computes signed adjacency matrix; (ii) converts it to TOM similarity; (iii) clusters it to get a dendrogram; (iv) it computes the 25th, 50th, and 75th percentiles of dendrogram heights and their interquartile range (IQR). The power value with the largest IQR is picked. When clusters are well-separated, dendrogram heights are more variable and thus have higher IQR.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        ## bslib::nav_panel(
        "Module-Trait",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>Module-trait heatmaps</b> show the correlation between eigengenes and traits (i.e. phenotype conditions). Heatmaps can be created for each datatype or merged. We look for modules that are highly correlated with traits.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            multiwgcna_plot_moduletrait_ui(
              ns("multiwgcnaTrait"),
              title = "Module-trait heatmaps",
              caption = "Module-trait heatmaps",
              info.text = "The multiomics WGCNA 'module-trait' module provides information on the correlation between traits (e.g., phenotypes) and feature modules, for each data type. This is critical in multi-omics data analysis to assess the importance of distinct molecular omics in modulating the phenotype of interest. Under 'Layers', you can optionally select the data type of interest, and include gene sets ('gs'). Under plot options, you can 'Merge modules' for a combined view of all data types and modules, 'show correlation values' to add correlation coefficients, 'transpose matrix' to flip the plot. Importantly, only the top 20 most correlated WGCNA modules are displayed by default (or all if less than 20 modules are available). All other modules can be shown by clicking at the option 'show all modules'.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        ## bslib::nav_panel(
        "Module correlation",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>Module correlation heatmaps</b> show the pairwise correlation of module eigengenes across layers. Heatmaps can be shown per layer or merged for all layers.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            multiwgcna_plot_modulecorr_ui(
              ns("multiwgcnaCorr"),
              title = "Module-module correlation heatmaps",
              caption = "Module-module correlation heatmaps",
              info.text = "The multiomics WGCNA 'Module-module correlation' module provides information on the correlation between modules across layers. All possible pairwise correlations are shown. Modules are first inferred for each data type. In multi-omics data analysis, correlation between modules of distinct data types allow to assess functional convergence between data types. Under 'Layers', you can optionally select the data type of interest, and include gene sets ('gs'). Critically, the correlation heatmaps can be conditioned on a phenotype of interest. Conditioning involves scaling to better highlight the modules most heavily associated with the selected phenotypes. Under plot options, you can 'Merge modules' for a combined view of all data types and modules, 'show correlation values' to add correlation coefficients. Importantly, only the top 20 most correlated WGCNA modules are displayed by default (or all if less than 20 modules are available). All other modules can be shown by clicking at the option 'show all modules'.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        ## bslib::nav_panel(
        "WGCNA-Lasagna",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>WGCNA-LASAGNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately and then integrated using LASAGNA (Layered Approach to Simultaneous Analysis of Genomic and Network Association).")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            multiwgcna_plot_lasagna_ui(
              ns("multiwgcnaLasagna"),
              title = "Multipartite graph",
              caption = "Multipartite graph.",
              info.text = "Multipartite graph. LASAGNA (Layered Approach to Simultaneous Analysis of Genomic and Network Association) is a stacked layer model for multi-omics where each layer corresponds to a data type. Each layer (vertical bar), corresponds to a data type. Under 'Layers', you can optionally select the data type of interest, and include gene sets ('gs'). WGCNA modules, inferred for each data type, are shown as nodes within each layer. The sizes of the node represent the page-rank centrality of the module. The log2FC is indicated for the chosen comparison. The color of the edges correspond to positive (purple) and negative (yellow) correlation between WGCNA feature modules. Thicker edges mean higher correlation. Correlation coefficients can be set under 'Lasagna option' to filter out lowly correlated modules.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Feature Table",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(3, 4, 5),
            height = "100vh",
            AISummaryCardUI(
              id = ns("multiwgcnaSummary"),
              title = "Summary",
              caption = "AI-generated summary of selected WGCNA module.",
              info.text = "Summary. Description about selected WGCNA module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            bslib::layout_columns(
              col_widths = c(12),
              multiwgcna_table_modulegenes_ui(
                ns("multiwgcnaTable"),
                title = "Module members",
                caption = "Table of features mapped in the selected module.",
                info.text = "Table of features mapped in the selected module. The data type and the WGCNA module can be selected under 'Module' in the drop-down menu on the right.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              multiwgcna_table_crossgenes_ui(
                ns("multiwgcnaCrossgene"),
                title = "Highly correlated genes",
                caption = "Table of features mapped in the selected module and their correlation coefficients.",
                info.text = "Table of features mapped in the selected module and their correlation coefficient. The data type and the WGCNA module can be selected under 'Module' in the drop-down menu on the right.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            multiwgcna_table_enrichment_ui(
              ns("multiwgcnaEnrichment"),
              title = "Module enrichment",
              caption = "Table of gene sets constructed from features mapped in the selected module. Gene set enrichment score, q value, and feature overlap are reported.",
              info.text = "Table of gene sets constructed from features mapped in the selected module. The data type and the WGCNA module can be selected under 'Module' in the drop-down menu on the right. Gene set enrichment score, q value, and feature overlap are reported.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )
    ) ## end tabsetPanel
  ) ## end div
}
