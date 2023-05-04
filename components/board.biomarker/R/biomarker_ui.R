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

  imgH1 <- c("calc(40vh - 120px)", "70vh") ## heights for small and fullscreen image
  imgH2 <- c("calc(60vh - 120px)", "70vh")

  div(
    boardHeader(title = "Biomarker Selection", info_link = ns("pdx_info")),
    bslib::layout_column_wrap(
      width = 1/2,
      height = "calc(100vh - 130px)",
      heights_equal = "row",
      biomarker_plot_importance_ui(
        ns("pdx_importance"),
        title = "Variable importance",
        info.text = "An importance score for each variable is calculated using multiple machine learning algorithms, including LASSO, elastic nets, random forests, and extreme gradient boosting. By combining several methods, the platform aims to select the best possible biomarkers. The top features are plotted according to cumulative ranking by the algorithms.",
        caption = "Barchart indicating the cumulative weight of a proposed biomarker based on six machine learning algorithms.",
        height = imgH1,
        width = c("auto", "100%"),
        label = "a"
      ),
      biomarker_plot_boxplots_ui(
        ns("pdx_boxplots"),
        title = "Biomarker expression",
        info.text = "These boxplots shows the expression of genes/samples of the identified features.",
        caption = "Expression boxplots of the most likely biomarkers across selected phenotypic groups.",
        height = imgH1,
        width = c("auto", "100%"),
        label = "b"
      ),
      biomarker_plot_heatmap_ui(
        ns("pdx_heatmap"),
        title = "Heatmap",
        info.text = "Expression heatmap of top gene features according to their variable importance.",
        caption = "Heatmap indicating the expression pattern across selected phenotypic groups for the most likely biomarkers.",
        height = imgH2,
        width = c("auto", "100%"),
        label = "c"
      ),
      biomarker_plot_decisiontree_ui(
        ns("pdx_decisiontree"),
        title = "Decision tree",
        info.text = "The decision tree shows a tree solution for classification based on the top most important features. The plot provides a proportion of the samples that are defined by each biomarker in the boxes.",
        caption = "Decision tree indicating expression-based biomarkers that distinguish the selected phenotypic groups.",
        height = imgH2,
        width = c("auto", "100%"),
        label = "d"
      )
    )
  )
}
