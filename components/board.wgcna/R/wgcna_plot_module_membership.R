##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_plot_module_membership_ui <- function(
    id,
    label,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)

  eigenCorrelation_opts <- shiny::tagList(
    shiny::checkboxInput(ns("eigen_cov"), "covariance", FALSE)
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = eigenCorrelation_opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

wgcna_plot_module_membership_server <- function(id,
                                                wgcna.compute,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    eigenCorrelation.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      MEs <- out$net$MEs
      ## rho <- cor(MEs, out$datExpr)
      rho <- cor(MEs, out$datExpr, use = "pairwise.complete.obs")
      rho[is.na(rho) | is.infinite(rho)] <- 0

      ylab0 <- "ME correlation"
      if (input$eigen_cov) {
        sdx <- apply(out$datExpr, 2, sd, na.rm = TRUE)
        rho <- t(t(rho) * sdx**2)
        ylab0 <- "ME covariance"
      }

      n <- nrow(rho)
      nr <- ceiling(sqrt(n))
      nc <- ceiling(n / nr)

      ntop <- 15
      par(mfrow = c(nr, nc), mar = c(6, 3.1, 2.3, 1), oma = c(1, 1, 1, 1) * 0, mgp = c(2.1, 0.8, 0))
      k <- 1
      me <- names(out$me.colors) ## sorted
      for (m in me) {
        i1 <- head(order(rho[m, ]), ntop)
        i2 <- tail(order(rho[m, ]), ntop)
        barplot(sort(rho[k, c(i1, i2)]),
          ylab = ylab0, las = 3, cex.names = 0.90, main = NULL
        )
        title(m, line = 0.3)
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = eigenCorrelation.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(90, 105),
      add.watermark = watermark
    )
  })
}
