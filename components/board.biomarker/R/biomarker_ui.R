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
              rows = 5
            )
          ),
          "Paste a custom list to be used as features.",
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
  imgH2 <- c("calc(60vh - 180px)", "70vh")
  fullH <- "calc(100vh - 200px)"

  div(
    boardHeader(title = "Biomarker Selection", info_link = ns("pdx_info")),
    shiny::tabsetPanel(
      id = ns("tabs1"),
      shiny::tabPanel(
        "Feature selection",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = fullH,
          biomarker_plot_importance_ui(
            ns("pdx_importance"),
            title = "Variable importance",
            info.text = "Barplot displaying the most important genes for the {Predicted target}.",
            info.methods = "The importance score for each gene is calculated using multiple machine learning algorithms, including sPLS [1], elastic nets [2], random forests [3], and extreme gradient boosting [4]. By combining several methods, the best possible biomarkers are selected. The top features are plotted according to cumulative ranking by the algorithms.",
            info.references = list(
              list(
                "Le Cao (2017). “mixOmics: An R package for 'omics feature selection and multiple data integration.” PLoS computational biology, 13(11), e1005752.",
                "https://doi.org/doi:10.18129/B9.bioc.mixOmics"
              ),
              list(
                "Friedman J, Tibshirani R, Hastie T (2010). “Regularization Paths for Generalized Linear Models via Coordinate Descent.” Journal of Statistical Software, 33(1), 1–22.",
                "https://doi.org/10.18637/jss.v033.i01"
              ),
              list(
                "Liaw A, Wiener M (2002). “Classification and Regression by randomForest.” R News, 2(3), 18-22.",
                "https://doi.org/10.32614/CRAN.package.randomForest"
              ),
              list(
                "Chen T (2024) xgboost: Extreme Gradient Boosting",
                "https://doi.org/10.32614/CRAN.package.xgboost"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#biomarker-analysis",
            caption = "Barchart indicating the cumulative weight of a proposed biomarker based on six machine learning algorithms.",
            height = imgH1,
            width = c("auto", "100%"),
            label = "a"
          ),
          biomarker_plot_boxplots_ui(
            ns("pdx_boxplots"),
            title = "Biomarker expression",
            info.text = "Boxplot displaying the expression for the top genes according to their importance. The expression is grouped by the {Predicted target} selected on the computation.",
            info.methods = "See Variable importance",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#biomarker-analysis",
            caption = "Expression boxplots of the most likely biomarkers across selected phenotypic groups.",
            height = imgH1,
            width = c("auto", "100%"),
            label = "b"
          ),
          biomarker_plot_heatmap_ui(
            ns("pdx_heatmap"),
            title = "Heatmap",
            info.text = "Expression heatmap of top 40 gene features according to their importance for the {Predicted target}. By default only the samples used on the computation are displayed, using the {show all samples} plot setting all samples can be displayed.",
            info.methods = "Heatmap clustering performed with the fastcluster R package [1] using the 'ward.D2' method with euclidean distance.",
            info.references = list(
              list(
                "Müllner D (2013). “fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python.” Journal of Statistical Software, 53(9), 1–18.",
                "https://doi.org/10.18637/jss.v053.i09"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#biomarker-analysis",
            caption = "Heatmap indicating the expression pattern across selected phenotypic groups for the most likely biomarkers.",
            height = imgH2,
            width = c("auto", "100%"),
            label = "c"
          ),
          biomarker_plot_decisiontree_ui(
            ns("pdx_decisiontree"),
            title = "Decision tree",
            info.text = "Decision tree showing a solution for classification based on the top most important features for the {Predicted target}. The plot provides a proportion of the samples that are defined by each biomarker in the boxes.",
            info.methods = "See Variable importance",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#biomarker-analysis",
            caption = "Decision tree indicating expression-based biomarkers that distinguish the selected phenotypic groups.",
            height = imgH2,
            width = c("auto", "100%"),
            label = "d"
          ),
          biomarker_plot_auc_ui( ## NEW AZ
            ns("pdx_auc"),
            title = "Receiver-operating characteristic (ROC) curve and Area Under the Curve (AUC)",
            info.text = "ROC and AUC for top biomarkers selected by regression model in decision tree.",
            info.methods = "ROC and AUC are common machine learning techniques. They are generally adopted to assess the relationship between correctly classified and incorrectly classified data points. The ROC curve is generally drawn by calculating the true positive rate (TPR) and false positive rate (FPR), expressed in terms of sensitivity and specificity.",
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#biomarker-analysis",
            caption = "ROC and AUC for top biomarkers selected by regression model in decision tree.",
            height = imgH1,
            width = c("auto", "100%"),
            label = "e"
          ),
        )
      ), ## tabPanel 1
      shiny::tabPanel(
        "Feature-set ranking",
        bslib::layout_columns(
          col_widths = c(6, 6),
          height = fullH,
          biomarker_plot_featurerank_ui(
            id = ns("featurerank"),
            title = "Feature-set ranking",
            info.text = "Ranked discriminant score for top feature sets. The plot ranks the discriminative power of the feature set (or gene family) as a cumulative discriminant score for all phenotype variables. Three methods to compute the discriminant score can be selected on the plot settings under {Method}.",
            info.methods = "Correlation-based discriminative power is calculated as the average '(1-cor)' between the groups (performed using core stats R package). Thus, a feature set is highly discriminative if the between-group correlation is low. 'P-value' based scoring is computed as the average negative log p-value from the ANOVA. The 'meta' method combines the score of the former methods in a multiplicative manner. For 'meta' and 'p-value' the limma R/Bioconductor package is used [1]. In this way, is possible to find which feature set (or gene family) can explain the variance in the data the best.",
            info.references = list(
              list(
                "Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), e47.",
                "https://doi.org/10.1093/nar/gkv007"
              )
            ),
            info.extra_link = "https://omicsplayground.readthedocs.io/en/latest/methods/#biomarker-analysis",
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
