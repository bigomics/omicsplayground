##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

PreservationWGCNA_Inputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::selectInput(ns("splitpheno"), "Split dataset by:", choices=NULL),
    #shiny::selectInput(ns("reference"), "Reference:", choices=NULL),    
    shiny::selectInput(ns("module"), "Module:", choices=NULL, multiple=FALSE),
    shiny::selectInput(ns("trait"), "Trait:", choices=NULL, multiple=FALSE),
    shiny::br(),
    shiny::actionButton(ns("compute"), "Compute", size = "xs", icon=icon("refresh")),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("mwgcna_options"),
      #open = NULL,
      bslib::accordion_panel(
        "WGCNA options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(          
          shiny::selectInput(ns("power"),"Soft treshold:",
            choices=c("<auto>",1,3,6,9,12,20), selected=12),
          shiny::selectInput(ns("deepsplit"),"Deepsplit:", choices=0:4, 2),
          shiny::selectInput(ns("ngenes"),"Max. features:", choices=c(1000,2000,4000),
            2000),
          shiny::selectInput(ns("minmodsize"),"Min. module size:", choices=c(5,10,20,40,100),
            10)
        )
      )
    )

  )
}


PreservationWGCNA_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Preservation WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Dendrograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1, 0.7),
          bs_alert(HTML("<b>Preservation WGCNA</b> is an application of WGCNA to test whether modules of a reference data set are preserved in others (test) datasets. A module is said to be preserved if the intramodule connectivity and density is maintained.")),
          bslib::layout_columns(
            col_widths = c(6,6),
            #height = "calc(100vh - 180px)",
            height = "100vh",
            preservationWGCNA_plot_dendrograms_ui(
              ns("preservationWGCNADendro"),
              title = "Dendrograms and Module Colors",
              caption = "Dendrogram of hierarchical clustering of feature co-expression patterns and WGCNA modules for preservation analysis.",
              info.text = "The gene dendrogram is a highly used visualization in WGCNA. It provides a combined visual summary of the features' hierarchical clustering dendrogram and the module color assignments. It shows that co-expressed features are grouped into common biologically meaningful modules. The tree of the dendrogram reflects the co-expression similarity of the features. Each leaf (end point) of the tree corresponds to a single feature. The branching pattern shows the correlation structure: features that are highly co-expressed (strongly correlated across samples) are clustered together and connected up until lower branch points. Features that are less correlated are joined higher up in the tree. For preservation WGCNA, it displays evidence of preservation of modules across datasets/datatypes/phenotypes. Modules that are 'preserved' are modules sharing feature correlation structure, possibly indicating conserved regulatory pathways and networks. Therefore, preservation WGCNA is an excellent complement to consensus WGCNA. Optionally, all levels for any selected phenotypes can be displayed. Annotation by each available trait can also be optionally displayed. Key preservation WGCNA algorithm parameters can be modified and WGCNA can be recomputed by clicking at the 'Compute' button.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            preservationWGCNA_plot_summaries_ui(
              ns("preservationWGCNASummaries"),
              title = "Preservation scores",
              caption = "Scatter plots of preservation WGCNA metrics. The following scatter plots are displayed: Zsummary vs. module size, medianRank vs. module size, MedianRank vs. Zsummary. In each scatter plot, each dot represents a module. Each dot is colored with the respective module color.",
              info.text = "Preservation WGCNA analysis relies on specific module preservation statistics, including Z-summary and medianRank. The Z-summary value combines multiple preservation metrics (density and connectivity preservation) into a single composite score that reflects how similar the module structure in the test dataset is to that in the reference dataset. A high Z-summary (e.g., >10) indicates strong preservation. A low Z-summary indicate a lowly/weakly preserved network. Notably, caution is advised when intepreting the Z-summary as it is sensitive to module size: larger modules tend to have higher Z-scores compared to small modules.<br><br>The MedianRank is a rank-based measure of module preservation, derived by ranking all modules based on the preservation measures connectivity and density. It is a composite measure that takes the average of two median rank scores: one for the module's connectivity and one for its density. Different to the Z-summary, the MedianRank is significantly more robust to variation in module size. It ranks the observed connectivity and density of a module against the distribution of all other modules in a reference network. A lower medianRank indicates a higher degree of preservation. Altogether, Z-summary and medianRank summarize how well each module’s structure and connectivity patterns are maintained between the reference and the test dataset. ",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )            
          )
        )
      ),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Module Overlap",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          #bs_alert(HTML("<b>Sample clustering</b> shows the clustering tree (of each datasts) of their samples. The heatmap shows sample traits and module eigengenes.")),
          bslib::layout_columns(
            col_widths = c(6,6),
            height = "100vh",
            preservationWGCNA_plot_overlap_ui(
              ns("preservationWGCNAOverlap"),
              title = "Module Overlap",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            preservationWGCNA_plot_eigenNetwork_ui(
              ns("preservationWGCNAEigenNetwork"),
              title = "Module Correlation and Preservation",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        ##bslib::nav_panel(      
        "Module-Trait",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          #bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(6,6),
            height = "100vh",
            bslib::layout_columns(
              col_widths = c(12),
              preservationWGCNA_plot_moduletrait_ui(
                ns("preservationWGCNAModuleTrait"),
                title = "Module Preservation and Module-Trait Consensus",
                caption = "Heatmap of WGCNA module preservation Zsummary values and WGCNA consensus Module-Trait correlation values.",
                info.text = "Heatmap of WGCNA module preservation Zsummary values. The Z-summary value combines multiple preservation metrics (density and connectivity preservation) into a single composite score that reflects how similar the module structure in the test dataset is to that in the reference dataset. A high Z-summary (e.g., >10) indicates strong preservation. A low Z-summary indicate a lowly/weakly preserved network. Caution is advised when intepreting the Z-summary as it is sensitive to module size: larger modules tend to have higher Z-scores compared to small modules. Zsummary Modules can optionally ordered by decreasing Zsummary values, hierarchical clustering order, or name.<br><br>Heatmap of WGCNA consensus Module-Trait correlation values.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              preservationWGCNA_plot_moduletrait_barplot_ui(
                ns("preservationWGCNAModuleTrait"),
                title = "Module-Trait barplots",
                caption = "Bar plots displaying correlation coefficients (y-axis) of each WGCNA module (x-axis) with the selected trait.",
                info.text = "Module-trait bar plots. Bar plots displaying correlation coefficients (y-axis) of each WGCNA module (x-axis) with the selected trait. The trait can be selected in the drop-down menu on the right. Optionally, WGCNA modules can be colored.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            preservationWGCNA_plot_traitsignificance_ui(
              ns("preservationWGCNATraitSignificance"),
              title = "Trait Significance",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),
      
      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Feature Table",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bslib::layout_columns(
            col_widths = c(7,5),
            height = "100vh",            
            bslib::layout_columns(
              col_widths = c(12),
              preservationWGCNA_table_modulegenes_ui(
                ns("preservationWGCNATable"),
                title = "Module Features",
                caption = "Table reporting the features mapped within the selected WGCNA module and associated with the selected phenotype/trait of interest. WGCNA metrics are also reported for each feature.",
                info.text = "Table reporting the features mapped within the selected WGCNA module and associated with the selected phenotype/trait of interest. For each feature, the following fields are reported: title, trait-specific score, module membership, module membership p value, log2 fold-change log2 fold-change p-value are displayed. Features are sorted by increasing score value.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              preservationWGCNA_table_enrichment_ui(
                ns("preservationWGCNAEnrichment"),
                title = "Module Enrichment",
                caption = "Table reporting the gene sets calculated from enrichment analysis of the features mapped within the selected WGCNA module. For each gene set, the following fields are reported: gene set enrichment score, enrichment meta p value (the highest p value among the single gene set enrichment methods' p values), enrichment q value, p value from each gene set enrichment method, number of featues mapped within the gene set, and features mapped within the gene set.",
                info.text = "Table reporting the gene sets calculated from enrichment analysis of the features mapped within the selected WGCNA module. For each gene set, the following fields are reported: gene set enrichment score, enrichment meta p value (the highest p value among the single gene set enrichment methods' p values), enrichment q value, p value from each gene set enrichment method, number of featues mapped within the gene set, and features mapped within the gene set.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            preservationWGCNA_plot_modulenetwork_ui(
              ns("preservationWGCNAModuleNetwork"),
              title = "Module Network",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )
      
    ) ## end tabsetPanel
  )  ## end div 
}
