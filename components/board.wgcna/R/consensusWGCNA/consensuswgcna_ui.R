##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

ConsensusWGCNA_Inputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::selectInput(ns("splitby"), "Consensus by:", choices=NULL),
    shinyjs::hidden(shiny::selectInput(ns("module"), "Module:", choices=NULL, multiple=FALSE)),
    shinyjs::hidden(shiny::selectInput(ns("trait"), "Trait:", choices=NULL)),
    shiny::br(),
    shiny::actionButton(ns("compute"), "Compute", size = "xs", icon = icon("refresh")),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("mwgcna_options"),
      # open = NULL,
      bslib::accordion_panel(
        "WGCNA options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(          
          shiny::selectInput(ns("power"),"Soft treshold:",
            choices=c("<auto>",1,3,6,9,12,20), selected=12),
          shiny::selectInput(ns("deepsplit"),"Deepsplit:", choices=0:4, 2),
          shiny::selectInput(ns("ngenes"),"Max. features:", choices=c(1000,2000,4000),
            2000),
          shiny::selectInput(ns("minmodsize"),"Min. module size:",
            choices=c(5,10,20,40,100), 10)
        )
      )
    )
  )
}


CONSENSUSWGCNA_INFO <- "The <b>Multi-partite graph</b> shows the correlation structure between multiple sets of features. The color of the edges correspond to positive (purple) and negative (yellow) correlation. Thicker edges mean higher correlation. The sizes of the circles represent the page-rank centrality of the feature. The log2FC is indicated for the chosen comparison. The node color corresponds to up (red) and down (blue) regulation."

ConsensusWGCNA_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Consensus WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Dendrograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1, 0.7),
          bs_alert(HTML("<b>Consensus WGCNA</b> is an application of WGCNA to identify modules that are conserved across two or more datasets by computing overlapping modules.")),
          bslib::layout_columns(
            col_widths = c(6, 6),
            # height = "calc(100vh - 180px)",
            height = "100vh",
            consensusWGCNA_plot_dendrograms_ui(
              ns("consensusWGCNADendro"),
              title = "Dendrograms and Module Colors",
              caption = "Dendrogram of hierarchical clustering of feature co-expression patterns and WGCNA module for consensus analysis.",
              info.text = "The gene dendrogram is a highly used visualization in WGCNA. It provides a combined visual summary of the features' hierarchical clustering dendrogram and the module color assignments. It shows that co-expressed features are grouped into common biologically meaningful modules. The tree of the dendrogram reflects the co-expression similarity of the features. Each leaf (end point) of the tree corresponds to a single feature. The branching pattern shows the correlation structure: features that are highly co-expressed (strongly correlated across samples) are clustered together and connected up until lower branch points. Features that are less correlated are joined higher up in the tree. For consensus WGCNA, it displays evidence of modules' consensus between datasets/datatypes/phenotypes. Modules that are 'in consensus' indicate shared feature correlation structure, and therefore shared regulatory networks and functional effects. Optionally, all levels for any selected phenotypes can be displayed. Annotation by each available trait can also be optionally displayed. Key consensus WGCNA algorithm parameters can be modified and WGCNA can be recomputed by clicking at the 'Compute' button.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_plot_power_ui(
              ns("consensusWGCNAPower"),
              title = "Scale and Connectivity Plots",
              caption = "Scale and connectivity plots. Scatter plots of WGCNA power values.",
              info.text = "The SFT model fit scatter plot shows the soft threshold power vs the signed scale-free topology fit index. For each dataset, datatype or phenotype, the slope of the log-log connectivity plot vs. frequency is extracted and its opposite sign is computed and This value is then multiplied by the scale-free topology fit index (R^2), with higher values indicating good scale-free structure. The resulting value is plotted against the SFT power value. In WGCNA, features are connected based on their correlation across samples. The mean connectivity is a measure of overall network density, i.e., how connected the features are, on average.<br><br>The mean connectivity graph plots soft threshold power vs. mean connectivity value, for each dataset, datatype, or phentoype shows the average number (or strenght) of connection per each feature, at any given value of power. A decreasing trend is an expected behavior: the mean connectivity decreases as the power increases, suggesting that the network becomes less dense (and therefore the adjacency matrix becomes sparser).<br><br>The soft threshold power for network sparsity can be changed under 'WGCNA option' and clicking at the 'Compute' button. Optionally, our in-house IQR method can be selected in the plot options to calculate the WGCNA power. IQR maximizes variation of heights of the dendrogram, promoting seperability between groups. It tests a range of power values. For each power: (i) computes signed adjacency matrix; (ii) converts it to TOM similarity; (iii) clusters it to get a dendrogram; (iv) it computes the 25th, 50th, and 75th percentiles of dendrogram heights and their interquartile range (IQR). The power value with the largest IQR is picked. When clusters are well-separated, dendrogram heights are more variable and thus have higher IQR.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        ## bslib::nav_panel(
        "Sample Clustering",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>Consensus WGCNA</b> is an application of WGCNA to identify modules that are conserved across two or more datasets (or datatypes, phenotypes) by clustering each dataset (or datatype, phenotype) and computing overlapping modules.<b>Sample clustering</b> shows the clustering tree (of each datasts) of their samples. The heatmap shows sample traits and module eigengenes.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            consensusWGCNA_plot_sampletree_ui(
              ns("consensusWGCNASampleTree"),
              title = "Sample Tree and Traits",
              caption = "Dendrogram and heatmap of samples and WGCNA modules for a specific dataset, datatype, or trait.",
              info.text = "The structure of the dendrogram reflects the similarity of the features. The title of each plot denotes the dataset, datatype, or phenotype for which consensus is being assessed. Each leaf (end point) of the tree corresponds to a sample. The structure of the tree, including clusters and branching pattern are indicative of the degree of dissimilarity between samples with respect to WGCNA modules. Samples that cluster together have similar WGCNA modules, i.e. similar eigengene/trait patterns. Large separations, i.e. long branches or distinct clusters, indicate samples whose expression patterns differ more strongly. On the rows of the heatmap, the WGCNA module eigengenes as inferred for a specific trait, are shown. Blue color indicates an overall lower expression of the module for a sample. Red color indicates an overall higher expression of the module for a sample. Similar correlation structures between modules and samples across phenotypes may appear for WGCNA modules exhibiting consensus between phenotypes.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      shiny::tabPanel(
        "Module-Trait",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1),
          bs_alert(HTML("<b>Consensus Module-Trait</b> analysis identifies modules that have high correlation with your phenotypes. Modules are concordant if the trait correlation have same sign in the consensus groups, i.e. always up (or down) regulated in all groups.")),
          bslib::layout_columns(
            col_widths = c(7, 5),
            height = "100vh",
            consensusWGCNA_plot_moduletrait_ui(
              ns("consensusWGCNATrait"),
              title = "Module-Trait Heatmaps",
              caption = "Heatmaps of correlation values between modules and traits.",
              info.text = "Heatmaps of correlation values between modules and traits. Heatmaps are displayed for each available groups within the selected pheonotype as well as for the consensus between groups. Correlation coefficient values can be optionally added into the heatmap. Statistical significance of the correlation coefficient values can also be optionally added into the heatmap, displayed as stars.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_plot_traitsignificance_ui(
              ns("consensusWGCNATraitSignificance"),
              title = "Trait Significance",
              caption = "Visualize the patterns of the module scores across all available groups within the selected phenotype.",
              info.text = "Scatter plots (for a continuous trait) or boxplots (for a categorical trait) of the scores for the selected module across all available groups within the selected phenotype. Modules that are concordant between phenotypes show similar patterns in score values. Modules that are discordant are expected to display opposite patterns across phenotypes. The number of plots displayed can be set using in the options.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        ),
        bslib::layout_columns(
          col_widths = c(12),
          height = "calc(100vh - 180px)",
          consensusWGCNA_plot_moduletrait_scatter_ui(
            ns("consensusWGCNATrait"),
            title = "Module-Trait Scatterplots",
            caption = "...",
            info.text = "...",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
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
          bs_alert(HTML("<b>Consensus WGCNA</b> is an application of WGCNA to identify modules that are conserved across two or more datasets (or datatypes, phenotypes) by clustering each dataset (or datatype, phenotype) and computing overlapping modules.")),
          bslib::layout_columns(
            col_widths = c(3, 4, 5),
            height = "100vh",
            wgcna_html_module_summary_ui(
              id = ns("consensusWGCNAmoduleSummary"),
              title = "Summary",
              info.text = "Summary. Description of the selected WGCNA module. The WGCNA module and phenotype/trait of interest can both be selected in the drop-down menu on the right.",
              caption = "LLM-generated information about the seleted WGCNA module and phenotype/trait of interest.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_table_modulegenes_ui(
              ns("consensusWGCNATable"),
              title = "Module Features",
              caption = "Table reporting the features mapped within the selected WGCNA module and associated with the selected phenotype/trait of interest. The table reports title, trait-specific score, score p value, and consensus status for each feature.",
              info.text = "Table reporting the features mapped within the selected WGCNA module and associated with the selected phenotype/trait of interest. For each feature, title, trait-specific score, score p value, and consensus status are displayed.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_table_enrichment_ui(
              ns("consensusWGCNAEnrichment"),
              title = "Module Enrichment",
              caption = "Table reporting the gene sets calculated from enrichment analysis of the features mapped within the selected WGCNA module. The table reports score, enrichment q value and overlap for each gene set.",
              info.text = "Table reporting the gene sets calculated from enrichment analysis of the features mapped within the selected WGCNA module. For each gene set, score, enrichment q value and overlap are displayed.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )
    ) ## end tabsetPanel
  ) ## end div
}
