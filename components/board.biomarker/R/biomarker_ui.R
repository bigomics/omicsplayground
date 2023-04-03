##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

BiomarkerInputs <- function(id) {
  ns <- shiny::NS(id) ## namespace
  bigdash::tabSettings(
    class = "p-1",
    shiny::tagList(
      shiny::hr(), shiny::br(),
      withTooltip(
        shiny::selectInput(ns("pdx_predicted"), "Predicted target:",
          choices = NULL
        ),
        "Select the target variable for biomarker selection.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("pdx_filter"), "Feature filter:",
          choices = NULL
        ),
        "Select a filter for the features.",
        placement = "top"
      ),
      shiny::conditionalPanel(
        "input.pdx_filter == '<custom>'",
        ns = ns,
        withTooltip(
          shiny::div(
            class = "gene-list",
            shiny::textAreaInput(ns("pdx_select"), "Custom features:",
              value = NULL,
              height = "100px", width = "100%",
              rows = 5, placeholder = "Paste your gene list"
            )
          ),
          "Paste a custom gene list to be used as features.",
          placement = "top"
        )
      ),
      shiny::br(),
      withTooltip(
        shiny::actionButton(ns("pdx_runbutton"), 
          label = "Compute",
          class = "btn-outline-primary"
        ),
        "Click to start biomarker computation.",
        placement = "right"
      )
    )
  )
}

BiomarkerUI <- function(id) {
  ns <- shiny::NS(id)
  imgH1 <- c("30vh", "70vh") ## heights for small and fullscreen image
  imgH2 <- c("50vh", "70vh") 
  
  imgH1 <- c("280px", "70vh") ## heights for small and fullscreen image
  imgH2 <- c("400px", "70vh") 
  
  div(
    boardHeader(title = "Biomarker Selection", info_link = ns("pdx_info")),
    tagList(
      div(
        class = "row row-cols-1 row-cols-md-2 row-cols-xxxl-4",
        div(
          class = "col",
          biomarker_plot_importance_ui(ns("pdx_importance"), height = imgH1, label = "a")
        ),
        div(
          class = "col",
          biomarker_plot_boxplots_ui(ns("pdx_boxplots"), height = imgH1, label = "b")
        ),
        div(
          class = "col",
          biomarker_plot_heatmap_ui(ns("pdx_heatmap"), height = imgH2, label = "c")
        ),
        div(
          class = "col",
          biomarker_plot_decisiontree_ui(ns("pdx_decisiontree"), height = imgH2, label = "d")
        )
      )
    )
  )
}
