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
        shiny::selectInput(ns("pdx_samplefilter"), "Filter samples:",
          choices = NULL, multiple = TRUE
        ),
        "Filter samples for the analysis.",
        placement = "top"
      ),
      withTooltip(
        shiny::selectInput(ns("pdx_filter"), "Feature set:",
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
      uiOutput(ns("biom_button"))
    )
  )
}

BiomarkerUI <- function(id) {
  ns <- shiny::NS(id)

  imgH1 <- c("calc(40vh - 120px)", "70vh") ## heights for small and fullscreen image
  imgH2 <- c("calc(60vh - 180px)", "70vh")
  fullH <- "calc(100vh - 200px)"

  div(
    boardHeader(title = "Biomarker Selection", info_link = ns("pdx_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Feature selection",
        bslib::layout_column_wrap(
          width = 1 / 2,
          height = fullH,
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
      ), ## tabPanel 1
      shiny::tabPanel(
        "Feature-set ranking",
        bslib::layout_column_wrap(
          width = 1 / 2,
          height = fullH,
          heights_equal = "row",
          biomarker_plot_featurerank_ui(
            id = ns("featurerank"),
            title = "Feature-set ranking",
            info.text = "Ranked discriminant score for top feature sets. The plot ranks the discriminative power of the feature set (or gene family) as a cumulative discriminant score for all phenotype variables. In this way, we can find which feature set (or gene family) can explain the variance in the data the best. Correlation-based discriminative power is calculated as the average '(1-cor)' between the groups. Thus, a feature set is highly discriminative if the between-group correlation is low. P-value based scoring is computed as the average negative log p-value from the ANOVA. The 'meta' method combines the score of the former methods in a multiplicative manner.",
            caption = "Ranked discriminant score for top feature sets.",
            label = "",
            height = c("100%", TABLE_HEIGHT_MODAL),
            width = c("auto", "100%")
          )
        )
      )
    ) ## tabsetpanel
  ) ## div
}
