##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanomethods_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "The <strong>Volcano (methods)</strong> panel displays the volcano plots provided by different enrichment calculation methods. This provides users an quick overview of the sensitivity of the statistical methods at once. Methods showing better statistical significance will show volcano plots with 'higher' wings."

  PlotModuleUI(
    ns("plot"),
    title = "Volcano plots for all methods",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcanomethods_server <- function(id,
                                                  pgx,
                                                  gs_features,
                                                  gs_contrast,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    volcanoMethods.RENDER <- shiny::reactive({
      shiny::req(pgx, gs_features())

      cmp <- 1
      cmp <- gs_contrast()
      mx <- pgx$gset.meta$meta[[cmp]]
      fx <- unclass(mx$fc)
      qv <- unclass(mx$q)
      pv <- unclass(mx$p)

      fx[which(is.na(fx))] <- NA
      fx[which(is.infinite(fx))] <- NA
      qv[which(is.na(qv))] <- 1

      fdr <- 1
      lfc <- 0
      fdr <- 0.05
      lfc <- 1
      fdr <- as.numeric(input$gs_fdr)
      lfc <- as.numeric(input$gs_lfc)
      sel.gsets <- rownames(mx)
      sel.gsets <- COLLECTIONS[[1]]
      sel.gsets <- COLLECTIONS[[gs_features()]]

      nlq <- -log10(1e-99 + unlist(qv))
      ymax <- max(3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis

      nplots <- ncol(fx)
      nc <- max(nplots / 2, 5)
      nn <- c(2, nc)
      par(mfrow = nn, mar = c(1, 1, 1, 1) * 0.2, mgp = c(2.6, 1, 0), oma = c(1, 1, 0, 0) * 2)

      shiny::withProgress(message = "Computing volcano plots ...", value = 0, {
        for (i in 1:nplots) {
          is.sig <- (qv[, i] <= fdr & abs(fx[, i]) >= lfc)
          sig.gs <- rownames(mx)[which(is.sig)]
          sig.gs <- intersect(sel.gsets, sig.gs)

          method <- colnames(fx)[i]
          gx.volcanoPlot.XY(
            x = fx[, i], pv = qv[, i],
            use.fdr = TRUE, p.sig = fdr, lfc = lfc,
            gene = rownames(mx),
            xlab = "effect size (NES)", ylim = c(0, ymax),
            lab.cex = 0, nlab = 0, axes = FALSE,
            render = "canvas", n = 1000, highlight = sig.gs,
            cex = 1, cex.axis = 1.3, main = ""
          )

          graphics::box(lwd = 1, col = "black", lty = "solid")

          ## volcano_plot(limma, render="plotly", n=1000, cex=1, highlight=genes)
          ## draw axis if first column or last row
          is.first <- (i %% nc == 1)
          last.row <- ((i - 1) %/% nc == (nplots - 1) %/% nc)
          if (is.first) axis(2, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          if (last.row) axis(1, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          legend("top",
            legend = method, box.lty = 0,
            x.intersp = 0.3, y.intersp = 0.5,
            inset = c(0, 0.01),
            cex = 1.2, bg = "white"
          )

          shiny::incProgress(1 / nplots)
        }
      })
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = volcanoMethods.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(75, 90),
      add.watermark = watermark
    )
  })
}
