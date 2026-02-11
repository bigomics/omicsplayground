##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_plot_power_ui <- function(
  id,
  title = "",
  info.text = "",
  caption = "",
  label = "",
  height = 400,
  width = 400
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::selectInput(
      inputId = ns("plottype"),
      label = "Plot type:",
      choices = c("sft.modelfit", "mean.k", "dendro.IQR"),
      selected = "sft.modelfit"
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    ## options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

multiwgcna_plot_power_server <- function(id,
                                         mwgcna,
                                         r_layers) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      wgcna <- mwgcna()

      layers <- r_layers()
      sel.layers <- intersect(layers, names(wgcna))
      wgcna <- wgcna[sel.layers]
      shiny::req(length(wgcna) > 0)

      exprList <- lapply(wgcna, function(w) t(w$datExpr))
      playbase::wgcna.plotPowerAnalysis_multi(
        exprList,
        maxpower = 20,
        # plots = input$plottype,
        nmax = 2000, cex = 1,
        RsquaredCut = 0.85,
        setPar = TRUE,
        cex.legend = 0.9,
        main = ""
      )
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(105, 140),
      add.watermark = FALSE
    )
  })
}
