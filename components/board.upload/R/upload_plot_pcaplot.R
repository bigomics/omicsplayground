##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


# PlotModuleUI for pcaplot
upload_plot_pcaplot_ui <- function(id,
                                   title,
                                   info.text,
                                   caption,
                                   label = "",
                                   height,
                                   width) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("pcaplot.method"), "Method:",
        choices = c("pca", "tsne", "umap"),
        width = "100%"
      ), "Choose clustering method.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    caption = caption,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    options = options,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

# PlotModuleServer for pcaplot
upload_plot_pcaplot_server <- function(id,
                                       phenoRT,
                                       countsRT,
                                       sel.conditions,
                                       watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        pheno <- phenoRT()
        counts <- countsRT()
        shiny::req(pheno)
        shiny::req(counts)
        method <- input$pcaplot.method
        X <- log2(1 + counts)
        X[is.na(X)] <- median(X, na.rm = TRUE)

        ## clust <- playbase::pgx.clusterMatrix.DEPRECATED(X, dims = 2, method = method)$pos2d
        clust <- playbase::pgx.clusterBigMatrix(X, dims = 2, method = method[1])
        clust[[1]]
      })

      plot.RENDER <- function() {
        pos2d <- plot_data()
        cond <- sel.conditions()
        shiny::req(cond)
        playbase::pgx.scatterPlotXY(
          pos2d,
          var = cond,
          plotlib = "plotly",
          legend = FALSE
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot.RENDER,
        func2 = plot.RENDER,
        csvFunc = plot_data,
        res = c(70, 140),
        pdf.width = 8, pdf.height = 8,
        add.watermark = watermark
      )
    }
  )
}
