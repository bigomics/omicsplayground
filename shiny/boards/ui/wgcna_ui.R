##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WgcnaInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        tags$div(
         HTML("<b>WGCNA Analysis.</b> Weighted correlation network analysis (WGCNA) is a gene-level cluster analysis
        method based on pairwise correlations between genes. It allows one to define modules (clusters),
        intramodular hubs, and network nodes with regard to module membership, to study the relationships
        between co-expression modules.")
        ),
        shiny::tagList(
            shiny::actionLink(ns("info"), "Info", icon=icon("info-circle")),
            shiny::hr(), shiny::br(),             
            
            ## data set parameters
            shiny::selectInput(ns('selected_module'),'select module', choices=NULL),
            shiny::actionButton(ns("compute"),"Compute!",icon=icon("running"),
                         class="run-button"),            
            shiny::br(),
            shiny::br(),
            shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
            shiny::br(),br(),br(),
            shiny::conditionalPanel(
                "input.options % 2 == 1", ns=ns,
                shiny::tagList(
                    shiny::selectInput(ns("ngenes"),"Number genes:",
                                choices = c(500,1000,2000,4000,8000),
                                selected = 1000),
                    shiny::selectInput(ns("minmodsize"),"Min. module size",
                                choices = c(10,30,100,250),
                                selected = 30 ),
                    shiny::selectInput(ns("power"),"Power", c(2,4,6,10), selected=6),
                    shiny::selectInput(ns("deepsplit"),"deepsplit", 0:4, selected=2),
                    shiny::selectInput(ns("cutheight"),"Merge cut height",
                                choices = c(0.05, 0.10, 0.25, 0.5, 0.9, 0.999),
                                selected = 0.25)
                )
            )
        )
    )
}

WgcnaUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("WGCNA",uiOutput(ns("wgcna_UI"))),
            shiny::tabPanel("Modules",uiOutput(ns("modules_UI"))),
            shiny::tabPanel("Eigengenes",uiOutput(ns("eigen_UI"))),
            shiny::tabPanel("Intramodular",uiOutput(ns("intra_UI")))
        )
    )
    ui
}