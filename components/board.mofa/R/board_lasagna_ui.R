##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LasagnaInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("contrast"), "Select comparison", choices = NULL),
    shiny::br(),
    shiny::br(),
    shiny::br(),    
    shiny::actionLink(ns("options"), "Options", icon = icon("cog", lib = "glyphicon")),
    shiny::br(), 
    shiny::conditionalPanel(
      "input.options % 2 == 1",
      ns = ns,
      shiny::tagList(
        shiny::selectInput(ns("ngenes"), tspan("Number genes:"),
          choices = c(500, 1000, 2000, 4000, 8000),
          selected = 1000
        ),
        shiny::selectInput(ns("cutheight"), "Merge cut height",
          choices = c(0.05, 0.10, 0.25, 0.5, 0.9, 0.999),
          selected = 0.25
        )
      )
    )
  )
}

my_navset_card_tab <- function(...) {
  htmltools::tagAppendAttributes(
    bslib::navset_card_tab(
      ... ,
      tags$style(HTML("@media (min-width: 1200px) {.root_navset { height: calc(100vh - 48px); }}"))
    ),
    class = "root_navset border-0"
  )
}


LasagnaUI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

##  shiny::div(
  ##boardHeader(title = "LASAGNA", info_link = ns("info")),
  ##  shiny::tabsetPanel(
  ## bslib::navset_card_tab(     
  my_navset_card_tab(
    id = ns("tabs"),
      title = boardHeader(title = "LASAGNA", info_link = ns("info")),
      
      ##----------------------------------------------------------------
      ##shiny::tabPanel(
      bslib::nav_panel(      
        "Mult-layer model",
        bslib::layout_columns(
          col_widths = 12,
#          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>LASAGNA</b> is a stacked layer model for multi-omics where each layer corresponds to a data type. The acronym stands for a <u>L</u>ayered <u>A</u>pproach to <u>S</u>imultaneous <u>A</u>nalysis of <u>G</u>enomic and <u>N</u>etwork <u>A</u>ssociations'.")),
          bslib::layout_columns(
            col_widths = c(6,6),
            height = "calc(100vh - 180px)",            
            mofa_plot_lasagna_ui(
              ns("lasagna"),
              title = "Multi-layer model",
              info.text = "LASAGNA is a acronym of 'Layered Approach to Simultaneous Analysis of Genomic and Network Associations.'",
              caption = "Scale independence and mean connectivity plots to determine the soft threshold.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_clustering_ui(
              ns("clusters"),
              title = "Feature UMAP per datatype",
              info.text = "Clustering of features",
              caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),

      
      ##----------------------------------------------------------------
      ## shiny::tabPanel(
      bslib::nav_panel(      
        "Path scoring",
        bslib::layout_columns(
          col_widths = 12,
          #height = "calc(100vh - 180px)",
          ## bs_alert(HTML("<b>LASAGNA</b> is a stacked layer model for multi-omics where each layer corresponds to a data type. The acronym stands for 'a Layered Approach to Simultaneous Analysis of Genomic and Network Associations'.")),
          bslib::layout_columns(
            col_widths = breakpoints(
              xxxl = c(6, 6),
              xl = c(12, 12),              
              sm = c(12, 12)
            ),
            bslib::layout_columns(
              col_widths = breakpoints(
                xxxl = c(6, 6, 12),
                xl = c(5, 3, 4),              
                sm = c(12, 12, 12)
              ),
              mofa_plot_lasagnaSP_ui(
                ns("lasagnaSP"),
                title = "Parallel coordinates",
                info.text = "LASAGNA is a acronym of 'Layered Approach to Simultaneous Analysis of Genomic and Network Associations.'",
                caption = "Scale independence and mean connectivity plots to determine the soft threshold.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_lasagnaSP_scores_ui(
                ns("lasagnaSP"),
                title = "Path scores",
                info.text = "Path scores",
                caption = ""
              ),
              mofa_plot_lasagnaSP_frequency_ui(
                ns("lasagnaSP"),
                title = "Term frequency",
                info.text = "Frequency of features",
                caption = "Term frequency in top ranked shortest path solutions."
              )
            ),
            mofa_table_lasagnaSP_ui(
              ns("lasagnaSP"),
              title = "Pathscore table",
              info.text = "Clustering of features",
              caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      )

    )
  ##)
}
