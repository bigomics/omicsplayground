##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

BiomarkerInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::tagList(
            withTooltip( shiny::actionLink(ns("pdx_info"), "Info", icon = shiny::icon("info-circle")),
                   "Show more information about this module."),
            shiny::hr(), shiny::br(),             
            withTooltip(shiny::selectInput(ns("pdx_predicted"),"Predicted target:", choices=NULL),
                   "Select the target variable for biomarker selection.",placement="top"),
            withTooltip( shiny::selectInput(ns("pdx_filter"),"Feature filter:", choices=NULL),
                   "Select a filter for the features.", placement="top"),
            shiny::conditionalPanel(
                "input.pdx_filter == '<custom>'", ns=ns,
                withTooltip(
                    shiny::div(class='gene-list',
                        shiny::textAreaInput(ns("pdx_select"), "Custom features:", value = NULL,
                                      height = "100px", width = "100%",                                       
                                      rows=5, placeholder="Paste your gene list")),
                    "Paste a custom gene list to be used as features.", placement="top")
            ),
            shiny::br(),
            withTooltip(shiny::actionButton(ns("pdx_runbutton"), label="Compute", class="run-button"),
                   "Click to start biomarker computation.", placement="right")
        )
    )
}

BiomarkerUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        flex = c(1),
        height = 760,
        shiny::tabsetPanel(
            id=ns("tabs"),
                shiny::fillCol(
                height = 800,
                flex = c(NA,0.05,1),
                tags$div(
                    HTML("<b>Biomarker selection</b>. The expression of certain genes may be used as <i>markers</i> to predict a certain
                    phenotype such as response to a therapy. Finding such <i>biomarkers</i> are of high importance in clinical applications.
                    <b>(a)</b> An importance score for each feature is calculated using multiple machine learning algorithms,
                    including LASSO, elastic nets, random forests, and extreme gradient boosting. The top features are plotted  according
                    to cumulative ranking by the algorithms. <b>(b)</b> The heatmap shows the expression distribution for the top most important features.
                    <b>(c)</b> The decision tree shows (one) tree solution for classification based on the top most important features.
                    <b>(d)</b> Boxplots show the expression of biomarker genes across the groups.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1,0.1,1),
                    height = 800,                
                    shiny::fillCol(
                        flex = c(0.5,1),
                        plotWidget(ns("pdx_importance")),
                        plotWidget(ns("pdx_heatmap"))
                    ),
                    shiny::br(), ## spacer
                    shiny::fillCol(
                        flex = c(1,0.9), 
                        plotWidget(ns("pdx_decisiontree")),
                        plotWidget(ns("pdx_boxplots"))
                    )
                )
            )
        )
    )
}