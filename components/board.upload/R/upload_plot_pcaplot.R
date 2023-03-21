##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

upload_plot_pcaplot_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  pcaplot.opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("pcaplot.method"), "Method:",
        choices = c("pca", "tsne", "umap"),
        width = "100%"
      ), "Choose clustering method.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = "PCA/tSNE plot",
    plotlib = "plotly",
    options = pcaplot.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

upload_plot_pcaplot_server <- function(id,
                                       phenoRT,
                                       countsRT,
                                       sel.conditions,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      getCurrentWordEnrichment()
    })

    pcaplot.RENDER <- shiny::reactive({
      ## ngs <- inputData()
      ## X <- ngs$X
      pheno <- phenoRT()
      counts <- countsRT()
      if (is.null(pheno) || is.null(counts)) {
        return(NULL)
      }
      if (NCOL(pheno) == 0 || NCOL(counts) == 0) {
        return(NULL)
      }
      shiny::req(pheno)
      shiny::req(counts)

      method <- input$pcaplot.method
      X <- log2(1 + counts)
      clust <- pgx.clusterMatrix(X, dims = 2, method = method)
      names(clust)

      cond <- sel.conditions()
      if (length(cond) == 0 || is.null(cond)) {
        return(NULL)
      }
      ## par(mar=c(4,1,1,1))
      pgx.scatterPlotXY(
        clust$pos2d,
        var = cond, plotlib = "plotly",
        legend = FALSE ## , labels=TRUE
      )
    })

    PlotModuleServer(
      "plot",
      func = pcaplot.RENDER,
      plotlib = "plotly",
      pdf.width = 5, pdf.height = 5,
      res = 72,
      add.watermark = watermark
    )
  })
}
