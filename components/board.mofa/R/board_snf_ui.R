##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

SNF_Inputs <- function(id) {
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

SNF_UI <- function(id) {
  ns <- shiny::NS(id) ## namespace

  fullH <- 700 ## full height of page
  rowH1 <- 250 ## row 1 height
  rowH2 <- 440 ## row 2 height

  shiny::div(
    boardHeader(title = "Similarity Network Fusion", info_link = ns("info")),
    shiny::tabsetPanel(
      id = ns("tabs"),

      ##----------------------------------------------------------------
      shiny::tabPanel(
        "SNF Clustering",
        bslib::layout_columns(
          col_widths = 12,
          height = "calc(100vh - 180px)",
          bs_alert(HTML("<b>Similarity Network Fusion</b> (SNF) is a network-based method for multi-omics integration by taking multiple views of a network and fusing them together to construct an overall status matrix (Wang et al., 2014).")),
          bslib::layout_columns(
            ##col_widths = c(6,6),
            col_widths = breakpoints(
              xxxl = c(7, 5),
              xl = c(7, 5),              
              sm = c(12, 12)
            ),
            mofa_plot_snf_ui(
              ns("snf_affinity"),
              title = "SNF affinity matrices",
              info.text = "SNF affinity matrices",
              caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
              height = c("100%", TABLE_HEIGHT_MODAL),
              width = c("auto", "100%")
            ),
            bslib::layout_columns(
              col_widths = 12,
              mofa_plot_snf_ui(
                ns("snf_heatmap"),
                title = "SNF heatmap",
                info.text = "SNF affinity matrices",
                caption = "Each datatype affinity matrix captures the pairwise similarities between samples, highlighting high similarities among samples within the same datatype.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              ),
              mofa_plot_snfgraph_ui(
                ns("snf_cluster"),
                title = "SNF clustering of samples",
                info.text = "Partial correlation graph",
                caption = "Module enrichment plot of top most enriched genesets.",
                height = c("100%", TABLE_HEIGHT_MODAL),
                width = c("auto", "100%")
              )
              ## mofa_plot_dummy_ui(
              ##   ns("snf_dummy"),
              ##   label = "c",
              ##   title = "SNF data",
              ##   info.text = "Partial correlation graph",
              ##   caption = "Module enrichment plot of top most enriched genesets.",
              ##   height = c("100%", TABLE_HEIGHT_MODAL),
              ##   width = c("auto", "100%")
              ## )
            )
          )
        )
      )

    )
  )
}
