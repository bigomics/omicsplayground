##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_correlation_network_ui <- function(
    id,
    title,
    info.text,
    caption,
    label,
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_correlation_network_server <- function(id,
                                                  wgcna.compute,
                                                  selected_module,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    corGraph.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      k <- selected_module()
      shiny::req(k)
      genes <- out$me.genes[[k]]

      dim(out$datExpr)
      xx <- cbind(out$net$MEs[, k, drop = FALSE], out$datExpr[, genes])
      ## rho1 <- cor(xx, out$net$MEs[, k])[, 1]
      rho1 <- cor(xx, out$net$MEs[, k], use = "pairwise.complete.obs")[, 1]
      ntop <- min(nrow(xx) - 1, 20)
      topgg <- names(sort(rho1, decreasing = TRUE))

      topgg <- head(topgg, ntop)

      ## rho <- Matrix::nearPD(cor(xx[, topgg]))$mat
      rho <- Matrix::nearPD(cor(xx[, topgg], use = "pairwise.complete.obs"))$mat
      me.color <- out$me.colors[k]
      color1 <- me.color
      color1 <- c("white", me.color)[1 + 1 * (colnames(rho) == k)]
      size1 <- c(7, 10)[1 + 1 * (colnames(rho) == k)]

      qgraph::qgraph(rho,
        graph = "glasso",
        layout = "spring",
        sampleSize = nrow(xx),
        labels = rownames(rho),
        color = color1,
        tuning = 0, ## gamma for EBIClasso. 0.5=default, 0=BIC
        vsize = size1,
        cut = 0,
        maximum = .45,
        border.width = 1.5,
        threshold = TRUE
      )
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = corGraph.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 80),
      add.watermark = watermark
    )
  })
}
