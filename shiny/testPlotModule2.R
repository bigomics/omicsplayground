##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

library(shiny)
source("./global.R")  ## global variable
source("./modules/plotModules/dataviewTSNEPlotModule2.R", encoding = "UTF-8")

load("../data/example-data.pgx",verbose=1)
ngs = pgx.initialize(ngs)
colnames(ngs$Y)  ## possible group vectors
genes <- rownames(ngs$X)

ui = fluidPage(
    fillCol(
        flex = c(1,1),
        height = 800,
        fillRow(
            tagList(
                selectInput("gene","Gene",genes),
                selectInput("groupby","Group",colnames(ngs$Y))
            ),
            plotWidget("tSNEPlot"),
            br(),br()
        ),
        fillRow()
    )
)

server = function(input, output, session) {

    filterList <- reactive({
        list(
            search_gene = input$gene,
            data_samplefilter = NULL,
            data_type = "logCPM",
            data_groupby = input$groupby
        )
    })
    
    dataviewtSNEplotModule(
        "tSNEPlot", filterList, ngs,
        label="A", imgH=400, watermark=FALSE
    )
        
}

shiny::shinyApp(ui, server=server)



