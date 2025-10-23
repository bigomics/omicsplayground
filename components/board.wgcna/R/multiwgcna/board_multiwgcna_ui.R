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
      shiny::selectInput(ns("gsfilter"),"Geneset filter:",choices=NULL)          
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
          shiny::selectInput(ns("power"),"Soft power:",
            choices=c("<auto>",1,3,6,9,12,20), selected=12),
          shiny::selectInput(ns("deepsplit"),"Deepsplit:", choices = c(0:4), 2),
          shiny::selectInput(ns("ngenes"),"Max. features:", choices = c(1000,2000,4000),
            2000),
          shiny::selectInput(ns("minmodsize"), "Min. module size",
            choices = c(5, 10, 20, 40, 100), selected = 10
          ),
          shiny::checkboxInput(ns("consensus"),"use consensus",FALSE),
          shiny::checkboxInput(ns("useLLM"),"AI summary", FALSE),
          shiny::br(),
          shiny::actionButton(ns("compute"), "Compute", size = "xs",
            icon=icon("refresh"))
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


MULTIWGCNA_INFO = "The <b>Multi-partite graph</b> shows the correlation structure between multiple sets of features. The color of the edges correspond to positive (purple) and negative (yellow) correlation. Thicker edges mean higher correlation. The sizes of the circles represent the page-rank centrality of the feature. The log2FC is indicated for the chosen comparison. The node color corresponds to up (red) and down (blue) regulation."

MultiWGCNA_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Multiomics WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Dendrograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Multiomics WGCNA</b> is a generalization of WGCNA for integratiing multi-omics where WGCNA is performed for each layer separately. Integration is performed by computing the module correlation across layers using LASAGNA.")),
          bslib::layout_columns(
            col_widths = c(12),
            #height = "calc(100vh - 180px)",
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
            #height = "calc(100vh - 180px)",
            height = "100vh",
            multiwgcna_plot_power_ui(
              ns("multiwgcnaPower"),
              title = "Scale and connectivity plots",
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
          bs_alert(HTML("<b>Module-trait heatmaps</b> show the correlation between eigengenes and traits (i.e. phenotype conditions). Heatmaps can be created for each datatype or merged. We look for modules that are highly correlated with traits.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            multiwgcna_plot_moduletrait_ui(
              ns("multiwgcnaTrait"),
              title = "Module-trait heatmaps",
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
        "Module correlation",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Module correlation heatmaps</b> show the pairwise correlation of module eigengenes across layers. Heatmaps can be shown per layer or merged for all layers.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            multiwgcna_plot_modulecorr_ui(
              ns("multiwgcnaCorr"),
              title = "Module-module correlation heatmaps",
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
        "WGCNA-Lasagna",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(12),
            height = "100vh",
            multiwgcna_plot_lasagna_ui(
              ns("multiwgcnaLasagna"),
              title = "Module-module correlation heatmaps",
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
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(3,4,5),
            height = "100vh",            
            wgcna_html_module_summary_ui(
              id = ns("multiwgcnaSummary"),
              title = "AI Summary",
              info.text = "",
              caption = "Information about the Module.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),                      
            bslib::layout_columns(
              col_widths = c(12),
              multiwgcna_table_modulegenes_ui(
                ns("multiwgcnaTable"),
                title = "Module members",
                caption = "...",
                info.text = "...",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              multiwgcna_table_crossgenes_ui(
                ns("multiwgcnaCrossgene"),
                title = "Highly correlated genes",
                caption = "...",
                info.text = "...",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
            ),
            multiwgcna_table_enrichment_ui(
              ns("multiwgcnaEnrichment"),
              title = "Module enrichment",
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
