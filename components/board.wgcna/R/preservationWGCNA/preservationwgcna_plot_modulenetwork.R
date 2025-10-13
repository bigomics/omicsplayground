##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_modulenetwork_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("plotgraph"), "Plot graph", TRUE),
    shiny::checkboxInput(ns("plotheatmap"), "Plot heatmap", TRUE),
    shiny::selectInput(
      inputId = ns("numgenes"),
      label = "Nr. of genes:",
      choices = c(15, 25, 50, 100, 250),
      selected = 25
    ),
    shiny::sliderInput(
      inputId = ns("rhograph"),
      label = "Graph threshold:",
      min=0.5, max=1, value=0.95, step=0.01
    ),
    shiny::sliderInput(
      inputId = ns("rhoheatmap"),
      label = "Heatmap threshold:",
      min=0.1, max=1, value=0.6, step=0.05
    ),
    shiny::sliderInput(
      inputId = ns("rgamma"),
      label = "Heatmap gamma:",
      min=1, max=12, value=6, step=1
    )
    
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_modulenetwork_server <- function(id,
                                                        rwgcna,
                                                        rmodule
                                                        ) {
  moduleServer(id, function(input, output, session) {
    
    renderPlot <- function(format=0) {
      res <- rwgcna()
      shiny::req(res)

      module <- rmodule()
      shiny::req(module)

      n <- length(res$layers) * (input$plotgraph + input$plotheatmap)
      shiny::validate(shiny::need(n>0, "Please select plot"))
      
      nr <- ceiling(sqrt(n))
      nc <- ceiling(n / nr)
      par(mfrow=c(nr,nc))
      if(format==1) par(mfrow=c(nc,nr))
      
      net <- res$layers[[1]]
      shiny::req(module %in% names(net$me.genes))
      genes <- net$me.genes[[module]]

      numgenes <- as.integer(input$numgenes)
      if(length(genes) > numgenes) {
        ##mm <- net$stats[['moduleMembership']][genes,module]
        mm <- matrixStats::colSds(net$datExpr[,genes])
        genes <- head(genes[order(-mm)], numgenes)
      }
      
      if(input$plotgraph) {
        par(mar=c(0,0,3,0))
        k=1
        for(k in 1:length(res$layers)) {
          playbase::wgcna.plotGeneNetwork(
            res$layers[[k]],
            sort(genes),
            rgamma = 1,
            min.rho = as.numeric(input$rhograph),
            edge.alpha = 0.3
          )
          title( names(res$layers)[k], line=1, cex.main=1.3)
          title( paste(module,"module"), line=0, cex.main=1, font.main=1)  
        }
      }
      
      if(input$plotheatmap) {      
        R <- cor(res$layers[[1]]$datExpr[,genes])
        ii <- hclust(as.dist(1-R), method="average")$order
        genes <- genes[ii]
        par(mar=c(3,2,3.5,7))
        k=1
        for(k in 1:length(res$layers)) {
          w <- res$layers[[k]]
          playbase::wgcna.plotModuleHeatmap(
            w,
            module = module,
            genes = genes,
            cex = 0.85,
            rgamma = as.numeric(input$rgamma),
            min.rho = as.numeric(input$rhoheatmap),
            type = "correlation",
            cluster = FALSE,
            main = ""
          )
          title( names(res$layers)[k], line=0.9, cex.main=1.3)  
        }
      }
    }

    render1 <- function() {
      renderPlot(format=0)
    }
    render2 <- function() {
      renderPlot(format=1)
    }
    
    PlotModuleServer(
      "plot",
      func = render1,
      func2 = render2,
      pdf.width = 12,
      pdf.height = 8,
      res = c(75, 100),
      add.watermark = FALSE
    )

    
  })
}



