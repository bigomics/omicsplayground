##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_s_independence_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = info.text,
    height,
    width,
    ... ) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    ...
  )
}

wgcna_plot_s_independence_server <- function(id,
                                             wgcna.compute,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    RENDER <- shiny::reactive({
      out <- wgcna.compute()
      shiny::req(out)

      networktype <- out$networktype
      if(is.null(networktype)) networktype <- "signed"
      
      ## Choose a set of soft-thresholding powers
      powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
      ## Call the network topology analysis function
      sft <- WGCNA::pickSoftThreshold(
        out$datExpr,
        powerVector = powers,
        networkType = networktype,
        verbose = 0
      )

      ## Plot the results:
      par(mfrow = c(1, 2), mar = c(3.3, 3, 1, 1), mgp = c(2, 0.8, 0))
      cex1 <- 0.9
      ## Scale-free topology fit index as a function of the soft-thresholding power
      base::plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        type = "n",
        xlab = "Soft threshold (power)",
        ylab = "SFT model fit (signed R^2)",
        main = NULL
      )

      text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        labels = powers, cex = cex1, col = "red"
      )
      ## this line corresponds to using an R^2 cut-off of h
      abline(h = 0.90, col = "red")

      ## Mean connectivity as a function of the soft-thresholding power
      base::plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
        type = "n",
        xlab = "Soft threshold (power)",
        ylab = "Mean connectivity",
        main = NULL
      )
      text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 120),
      add.watermark = watermark
    )
  })
}
