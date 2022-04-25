##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_phenoheatmap_ui <- function(id, label='', height=c(600,800))
{

  ns <- shiny::NS(id)

  info_text = "<b>Phenotype clustering.</b> Clustered heatmap of sample information (i.e. phenotype data). Column ordering has been performed using hierarchical clustering on a one-hot encoded matrix."

  opts <- shiny::tagList(
    withTooltip( shiny::checkboxInput(ns('clustsamples'),'cluster samples',TRUE),
                "Cluster samples.", placement="top")
  )
  
  PlotModuleUI(
    ns("pltmod"),
    title = "Phenotype clustering",
    label = label,
    info.text = info_text,
    options = opts,
    download.fmt = c("png","pdf","csv"),         
    width = c("auto","100%"),
    height = height
  )
  
}

dataview_plot_phenoheatmap_server <- function(id, pgx, r.samples, watermark=FALSE)
{
  moduleServer( id, function(input, output, session) {

    ## extract data from pgx object
    plot_data  <- shiny::reactive({
      
      shiny::req(pgx$samples,r.samples())

      dbg("[dataview_plot_phenoheatmap_server:plot_data] called!")

      samples <- r.samples() 
      do.clust <- input$clustsamples
      
      annot <- pgx$samples
      annot <- annot[samples,,drop=FALSE]

      list(
        annot = annot,
        do.clust = do.clust
      )            
    })
    
    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res$annot)

      dbg("[dataview_plot_phenoheatmap_server:plot.RENDER] called!")
      dbg("[dataview_plot_phenoheatmap_server:plot.RENDER] dim(res$annot) = ",dim(res$annot))
      
      annot.ht <- ifelse(ncol(res$annot) > 10, 5, 6)
      annot.ht <- ifelse(ncol(res$annot) > 20, 4, annot.ht)
      annot.ht <- ifelse(ncol(res$annot) > 30, 3, annot.ht)

      plt <- pgx.plotPhenotypeMatrix0(
        annot = res$annot,
        annot.ht = annot.ht,
        cluster.samples = res$do.clust
      )
      ## plt <- plt %>% plotly::config(displayModeBar = FALSE)
      plt
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
      csvFunc = plot_data,   ##  *** downloadable data as CSV
      renderFunc = shiny::renderPlot,
      renderFunc2 = shiny::renderPlot,        
      res = c(60,120),                ## resolution of plots
      pdf.width = 8, pdf.height = 5,
      add.watermark = watermark
    )

  })  ## end of moduleServer
}


