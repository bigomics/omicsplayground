##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DataViewInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>DataView.</b> Information and descriptive statistics to quickly lookup a gene, check the total counts, or view the data tables.")
        ),
        shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("data_info"), "Tutorial", icon = shiny::icon("youtube")),
                   "Show more information about this module."),
            shiny::hr(), shiny::br(),
            shinyBS::tipify( shiny::selectInput(ns("search_gene"),"Gene:", choices=NULL),
                   "Enter a gene of interest for the analysis.", placement="top"),
            shinyBS::tipify( shiny::selectInput(ns("data_samplefilter"),"Filter samples:",
                                choices=NULL, multiple=TRUE),
                   "Filter the relevant samples for the analysis.", placement="top"),
            shinyBS::tipify( shiny::selectInput(ns('data_groupby'),'Group by:', choices=NULL),
                   "Select phenotype for grouping the samples.", placement="top"),
            shiny::br(),
            shinyBS::tipify( shiny::actionLink(ns("data_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            shiny::br(),br(),
            shiny::conditionalPanel(
                "input.data_options % 2 == 1", ns=ns,
                shinyBS::tipify( shiny::radioButtons(ns('data_type'),'Data type:',
                                     choices=c("counts","logCPM"), selected="logCPM", inline=TRUE),
                       "Choose an input data type for the analysis.", placement="bottom")
            )
        )
    )
}

DataViewUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tabsetPanel(
        id = ns("tabs"),
        shiny::tabPanel("Plots",uiOutput(ns("plotsUI"))),
        shiny::tabPanel("QC",uiOutput(ns("countsUI"))),
        shiny::tabPanel("Counts",uiOutput(ns("genetableUI"))),
        shiny::tabPanel("Samples",uiOutput(ns("sampletableUI"))),
        shiny::tabPanel("Contrasts",uiOutput(ns("contrasttableUI"))),
        shiny::tabPanel("Resource info",uiOutput(ns("resourceinfoUI")))
    )
}