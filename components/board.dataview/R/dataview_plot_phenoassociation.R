##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_phenoassociation_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    info_text = "<b>Phenotype clustering.</b> Clustered heatmap of sample information (i.e. phenotype data). The values corresponds to the -log10(p) value of the corresponding statistical test between two phenotype variables. A higher value corresponds to stronger 'correlated' variables. For discrete-discrete pairs the Fisher's exact test is used. For continuous-discrete pairs, the Kruskal-Wallis test is used. For continuous-continous pairs, Pearson's correlation test is used."

    opts <- shiny::tagList(
        withTooltip( shiny::checkboxInput(ns('phenoclustsamples'),'cluster samples',TRUE),
               "Cluster samples.", placement="top")
    )
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Phenotype association",
        label = label,
        info.text = info_text,
        options = opts,
        download.fmt = c("png","pdf","csv"),         
        width = c("auto","100%"),
        height = height
    )
    
}

dataview_plot_phenoassociation_server <- function(id, pgx, r.samples, watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

      plot_data  <- shiny::reactive({
            shiny::req(pgx$X,pgx$Y)
            samples <- r.samples()
            annot <- pgx$samples
            annot <- annot[samples,,drop=FALSE]
            list(annot = annot)
        })
            
        plot.RENDER <- function() {
            res <- plot_data()
            shiny::req(res)            
            pq <- pgx.testPhenoCorrelation(res$annot, plot=TRUE)
        }
        
        modal_plot.RENDER <- function() {
            plot.RENDER()
        }
        
        PlotModuleServer(
            "pltmod",
            plotlib = "base",
            plotlib2 = "base",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            ## csvFunc = plot_data,   ##  *** downloadable data as CSV
            renderFunc = shiny::renderPlot,
            renderFunc2 = shiny::renderPlot,        
            res = c(100,170)*0.85,                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


