##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_plot_power_ui <- function(
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
    shiny::checkboxInput(
      inputId = ns("showiqr"),
      label = "Show IQR",
      value = FALSE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

consensusWGCNA_plot_power_server <- function(id,
                                             mwgcna,
                                             r_layers) {
  moduleServer(id, function(input, output, session) {
    plot.RENDER <- function() {
      cons <- mwgcna()
      shiny::req(cons)

      nw <- length(cons$datExpr)
      plots <- c("sft.modelfit", "mean.k")
      if (input$showiqr) plots <- c(plots, "dendro.IQR")

      par(mfrow = c(nw, length(plots)))
      par(mar = c(5, 4, 3, 1), mgp = c(2.5, 0.8, 0))
      i <- 1
      for (i in 1:length(cons$datExpr)) {
        playbase::wgcna.plotPowerAnalysis(
          cons$datExpr[[i]],
          maxpower = 25,
          plots = plots,
          RsquaredCut = 0.85,
          setPar = FALSE,
          main = toupper(names(cons$datExpr)[i])
        )
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 12,
      res = c(90, 130),
      add.watermark = FALSE
    )
  })
}
