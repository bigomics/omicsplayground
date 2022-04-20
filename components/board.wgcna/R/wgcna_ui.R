##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WgcnaInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    bigdash::tabSettings(
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
}

WgcnaUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("WGCNA",
                shiny::fillCol(
                flex = c(NA,0.02,0.75,0.08,1),
                height = 700,
                tags$div(
                HTML("<b>WGCNA module detection.</b> <b>(a)</b> Modules are detected as branches of the resulting cluster tree using
                the dynamic branch cutting approach. Genes inside a given module are summarized with the module eigengene.
                The module eigengene of a given module is defined as the first principal component of the standardized expression profiles.
                <b>(b)</b> Scale independence and mean connectivity plots to determine the soft threshold <b>(c)</b> Topological
                overlap matrix visualized as heatmap <b>(d)</b> Dimensionality reduction maps colored by WGCNA module
                <b>(e)</b> Graph network of WGCNA modules.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1,1),                
                    plotWidget(ns('geneDendro')),
                    plotWidget(ns('topologyPlots'))                
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1.3,1,0.8),
                    plotWidget(ns('TOMplot')),
                    plotWidget(ns('umap')),
                    plotWidget(ns('moduleGraph'))
                )
            )),
            shiny::tabPanel("Modules",
                shiny::fillCol(
                flex = c(NA, 0.05, 2, 0.07, 1),
                height = 700,
                tags$div(
                HTML("<b>WGCNA functional analysis.</b> <b>(a)</b> Module-trait analysis identifies modules that are significantly
                associated with the measured clinical traits by quantifying the association as the correlation of the
                eigengenes with external traits. <b>(b)</b> Partial correlation network of genes most correlated
                to the eigengene. <b>(c)</b> Module enrichment plot of top most enriched genesets. <b>(d)</b> Table of
                genes in the selected module. <b>(e)</b> Functional enrichment of the module calculated using Fisher's exact test.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1.2, 1.1, 0.9),
                    plotWidget(ns('moduleTrait')),                
                    plotWidget(ns('corGraph')),
                    plotWidget(ns('enrichPlot'))
                ),
                shiny::br(),
                shiny::fillRow(
                    flex = c(1, 0.08, 2.3),
                    tableWidget(ns("geneTable")),
                    shiny::br(),
                    tableWidget(ns("enrichTable"))
                )
            )),
            shiny::tabPanel("Eigengenes",
                shiny::fillCol(
                flex = c(NA,0.04,2),
                height = 700,
                tags$div(
                HTML("<b>WGCNA eigengene analysis.</b> <b>(a)</b> It is often interesting to visualizing the network of eigengenes
                and study the relationships among the found modules. One can use the eigengenes as represen- tative profiles
                and quantify module similarity by eigengene correlation. <b>(b)</b> For each module, we also define
                a quantitative measure of 'module membership' (MM) as the correlation of the module eigengene and the gene
                expression profile. This allows us to quantify the similarity of all genes to every module.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex=c(1,0.06,2.5),
                    plotWidget(ns('eigenClustering')),
                    shiny::br(),
                    plotWidget(ns('eigenCorrelation'))
                )
            )),
            shiny::tabPanel("Intramodular",
                shiny::fillCol(
                flex = c(NA,0.04,2,1),
                height = 700,
                tags$div(
                HTML("<b>WGCNA intramodular analysis.</b> We quantify associations of individual genes with our trait of
                interest (weight) by defining Gene Significance GS as (the absolute value of) the correlation between
                the gene and the trait. For each module, we also define a quantitative measure of module membership MM
                as the correlation of the module eigengene and the gene expression profile. Using the GS and MM measures,
                we can identify genes that have a high significance for weight as well as high module membership in interesting modules.")
                ),
                shiny::br(),
                shiny::fillRow(
                    flex=c(1,0.06,2.5),
                    plotWidget(ns('intraHeatmap')),
                    shiny::br(),
                    plotWidget(ns('intraScatter'))
                )
            ))
        )
    )
    
}
