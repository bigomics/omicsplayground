##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_barplot_ui <- function(
  id,
  label = "",
  title = "",
  info.text = "",
  caption = "",
  height = 400,
  width = 400,
  ...
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(ns("addtraits"), "Add traits", TRUE)
  )

  PlotModuleUI(
    ns("plotmodule"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    height = height,
    caption = caption,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    ...
  )
}

wgcna_plot_module_barplot_server <- function(id,
                                             wgcna,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    csvFunc <- function() {
      res <- wgcna()
      data.frame(module = colnames(res$W), size = colSums(res$W != 0))
    }

    RENDER <- function() {
      res <- wgcna()
      W <- res$W
      msize <- colSums(W != 0)
      par(mar = c(5, 4, 1, 1))
      barplot(msize,
        ylab = "Module size (N)",
        col = res$me.colors[colnames(W)],
        width = 1, las = 3, names.arg = ""
      )
      dy <- 0.04 * max(msize)
      text(
        x = (-0.33 + 1:length(msize)) * 1.2,
        ## Move labels to just below bottom of chart.
        y = par("usr")[3] - dy,
        labels = names(msize),
        xpd = NA,
        srt = 35,
        adj = 0.965,
        cex = 0.9
      )
    }

    PlotModuleServer(
      "plotmodule",
      func = RENDER,
      csvFunc = csvFunc,
      pdf.width = 8, pdf.height = 6,
      res = c(75, 120),
      add.watermark = watermark
    )
  })
}
