##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MPCSF_Inputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    shiny::hr(), shiny::br(),

    ## data set parameters
    shiny::selectInput(ns("selected_pheno"), "Select phenotype", choices = NULL),
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

MPCSF_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Multi-Omics PCSF", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "MultiPCSF",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Multi-Omics PCSF</b>")),
          bslib::layout_columns(
            col_widths = c(6,6),
            height = "calc(100vh - 180px)",            
            mofa_plot_dummy_ui(
              ns("mpcsf"),
              title = "SNF affinity matrices",
              info.text = "SNF affinity matrices",
              caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            mofa_plot_dummy_ui(
              ns("dummy"),
              title = "SNF heatmap",
              info.text = "SNF affinity matrices",
              caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            )
          )
        )
      ),


    )
  )
}
