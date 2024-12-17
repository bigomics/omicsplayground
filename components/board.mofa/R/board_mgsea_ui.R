##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MGseaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("contrast"), "Select contrast", choices = NULL),
    shiny::br(),
    shiny::br(),
    shiny::br(),    
    shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
    shiny::br(), 
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        ## shiny::selectInput(ns("ngenes"), tspan("Number genes:"),
        ##   choices = c(500, 1000, 2000, 4000, 8000),
        ##   selected = 1000
        ## ),
        ## shiny::selectInput(ns("cutheight"), "Merge cut height",
        ##   choices = c(0.05, 0.10, 0.25, 0.5, 0.9, 0.999),
        ##   selected = 0.25
        ## )
      )
    )
  )
}

MGseaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Multi-Omics GSEA",
                info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "multiGSEA",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>MultiGSEA</b> combines pathway enrichment on multiple omics layers to create a robust composite multi-omics pathway enrichment measure.")),
          bslib::layout_columns(
            col_widths = breakpoints(
              sm = c(12, 12, 12, 12),
              xl = c(7, 5, 7, 5),              
              xxxl = c(4, 3, 5, 4)
            ),
            mofa_plot_enrichment_ui(
              ns("menrichment"),
              title = "multiGSEA enrichment",
              info.text = "",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),            
            mofa_plot_mgsea_ui(
              ns("mgsea_plot"),
              title = "MultiGSEA plot",
              info.text = "MultiGSEA plot",
              caption = "The plot simultaneously visualizes the enrichment scores of two omics types in one figure. Pathway/genesets that are enriched in both modalities have a higher multi-score (see Table).",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_table_mgsea_ui(
              ns("mgsea_table"),
              title = "multiGSEA scores",
              info.text = "",
              caption = "",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            ),
            mofa_plot_pathbank_ui(
              ns("pathbank_pathway"),
              title = "Multi-omics pathway",
              caption = "Pathways that integrate proteomics and metabolomics data types in a single pathway diagram.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("100%", "100%")
            )                        
          )
        )
      )

      
    )
  )
}
