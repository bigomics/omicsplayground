##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##

MultiWGCNA_Inputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::selectInput(ns("phenotype"), "Phenotype", choices = NULL),
    shiny::selectInput(ns("layers"), "Layers:", choices = NULL, multiple = TRUE),
    shiny::conditionalPanel(
      "input.layers && input.layers.indexOf('gset') > -1",
      ns = ns,
      shiny::selectInput(ns("gsfilter"),"Geneset filter:",choices=NULL)          
    ),
    shiny::selectInput(ns("module"), "Module:", choices = NULL, multiple = FALSE),
    shiny::br(),
    shiny::actionButton(ns("updateplots"), "Update plots", size = "xs", icon=icon("refresh")),
    shiny::br(),
    shiny::br(),
    bslib::accordion(
      id = ns("mwgcna_options"),
      open = FALSE,
      bslib::accordion_panel(
        "WGCNA options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::selectInput(ns("clust"),"Cluster method:",c("average","complete")),
          shiny::selectInput(ns("power"),"Soft power:",
            choices=c("<auto>",1,3,6,9,12,20), "<auto>"),
          shiny::selectInput(ns("deepsplit"),"Deepsplit:", choices=0:4, 2),
          shiny::selectInput(ns("ngenes"),"Max. features:", choices=c(1000,2000,4000),
            2000)
        )
      )
    ),
    shiny::br(),
    bslib::accordion(
      id = ns("mpartite_options"),
      open = FALSE,
      bslib::accordion_panel(
        "Network options",
        icon = icon("cog", lib = "glyphicon"),
        shiny::tagList(
          shiny::checkboxInput(ns("edge_norm"),"normalize edges",TRUE),
          shiny::checkboxInput(ns("edge_pos"),"positive edges",FALSE),
          shiny::checkboxInput(ns("solveSP"),"shortest path",FALSE)          
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
    boardHeader(title = "Multi-WGCNA", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "Dendrograms",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(12),
            #height = "calc(100vh - 180px)",
            height = "100vh",
            multiwgcna_plot_dendrograms_ui(
              ns("multiwgcnaDendro"),
              title = "Dendrograms and module colors",
              caption = "...",
              info.text = "...",
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
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
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
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
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
        "Feature table",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          row_heights = c("auto",1),
          bs_alert(HTML("<b>Multi-WGCNA</b> is an application of WGCNA for multi-omics where WGCNA is performed on each layer separately.")),
          bslib::layout_columns(
            col_widths = c(6,6),
            height = "100vh",            
            multiwgcna_table_modulegenes_ui(
              ns("multiwgcnaTable"),
              title = "Module tables",
              caption = "...",
              info.text = "...",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
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
