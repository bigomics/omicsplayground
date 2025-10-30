##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

ConsensusWGCNA_Inputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shinyjs::hidden(shiny::selectInput(ns("splitpheno"), "Consensus by:", choices=NULL)),
    shinyjs::hidden(shiny::selectInput(ns("splitdata"), "Consensus by:", choices=NULL,
      multiple=TRUE)),              
    shinyjs::hidden(shiny::selectInput(ns("module"), "Module:", choices=NULL, multiple=FALSE)),
    shinyjs::hidden(shiny::selectInput(ns("trait"), "Trait:", choices=NULL)),
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
          shiny::selectInput(ns("minmodsize"),"Min. module size:",
            choices=c(5,10,20,40,100), 10),
          shiny::checkboxInput(ns("useLLM"),"AI summary", FALSE)
        )
      )
    )

  )
}


CONSENSUSWGCNA_INFO = "The <b>Multi-partite graph</b> shows the correlation structure between multiple sets of features. The color of the edges correspond to positive (purple) and negative (yellow) correlation. Thicker edges mean higher correlation. The sizes of the circles represent the page-rank centrality of the feature. The log2FC is indicated for the chosen comparison. The node color corresponds to up (red) and down (blue) regulation."

ConsensusWGCNA_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Consensus WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Dendrograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto", 1, 0.7),
          bs_alert(HTML("<b>Consensus WGCNA</b> is an application of WGCNA to identify modules that are conserved across two or more datasets by computing overlapping modules.")),
          bslib::layout_columns(
            col_widths = c(6,6),
            #height = "calc(100vh - 180px)",
            height = "100vh",
            consensusWGCNA_plot_dendrograms_ui(
              ns("consensusWGCNADendro"),
              title = "Dendrograms and Module Colors",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_plot_power_ui(
              ns("consensusWGCNAPower"),
              title = "Scale and Connectivity Plots",
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
        "Sample Clustering",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Sample clustering</b> shows the clustering tree (of each datasts) of their samples. The heatmap shows sample traits and module eigengenes.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            consensusWGCNA_plot_sampletree_ui(
              ns("consensusWGCNASampleTree"),
              title = "Sample Tree and Traits",
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
        "Module-Trait",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Consensus Module-Trait</b> analysis identifies modules that have high correlation with your phenotypes. Module are concordant if the trait correlation have same sign in the consensus groups, i.e. always up (or down) regulated in all groups.")),
          bslib::layout_columns(
            col_widths = c(7,5),
            height = "100vh",
            consensusWGCNA_plot_moduletrait_ui(
              ns("consensusWGCNATrait"),
              title = "Module-Trait Heatmaps",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_plot_traitsignificance_ui(
              ns("consensusWGCNATraitSignificance"),
              title = "Trait Significance",
              caption = "...",
              info.text = "...",
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
      
      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Feature Table",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(3,4,5),
            height = "100vh",            
            wgcna_html_module_summary_ui(
              id = ns("consensusWGCNAmoduleSummary"),
              title = "AI Summary",
              info.text = "",
              caption = "Information about the Module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_table_modulegenes_ui(
              ns("consensusWGCNATable"),
              title = "Module Features",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            consensusWGCNA_table_enrichment_ui(
              ns("consensusWGCNAEnrichment"),
              title = "Module Enrichment",
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
