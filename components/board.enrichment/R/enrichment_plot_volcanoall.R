##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanoall_ui <- function(id, height, width) {
  ns <- shiny::NS(id)

  info_text <- "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for gene sets across all contrasts. This provides users an overview of the statistics across all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."

  PlotModuleUI(
    ns("plot"),
    title = "Volcano plots for all contrasts",
    info.text = info_text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcanoall_server <- function(id,
                                              pgx,
                                              gs_features,
                                              gs_statmethod,
                                              gs_fdr,
                                              gs_lfc,
                                              calcGsetMeta,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    volcanoAll.RENDER <- shiny::reactive({
      shiny::req(pgx)
      if (is.null(gs_features())) {
        return(NULL)
      }

      meta <- pgx$gset.meta$meta
      gsmethod <- colnames(meta[[1]]$fc)
      gsmethod <- gs_statmethod()
      if (is.null(gsmethod) || length(gsmethod) == 0) {
        return(NULL)
      }

      fdr <- 1
      lfc <- 1
      fdr <- as.numeric(gs_fdr())
      lfc <- as.numeric(gs_lfc())
      sel.gsets <- NULL
      sel.gsets <- rownames(meta[[1]])
      sel.gsets <- COLLECTIONS[[1]]
      sel.gsets <- COLLECTIONS[[gs_features()]]

      i <- 1
      mx.list <- list()
      for (i in 1:length(meta)) {
        mx.list[[i]] <- calcGsetMeta(i, gsmethod, pgx = pgx)
      }
      names(mx.list) <- names(meta)

      Q <- lapply(mx.list, function(mx) mx[, "qv"])
      names(Q) <- names(mx.list)

      ## select maximum 24 comparisons (because of space...)
      q.score <- sapply(Q, function(q) mean(tail(sort(-log10(1e-99 + q)), 100)))
      sel <- head(names(sort(q.score, decreasing = TRUE)), 20)
      Q <- Q[which(names(Q) %in% sel)]
      mx.list <- mx.list[names(Q)]
      nlq <- -log10(1e-99 + unlist(Q))
      ymax <- max(3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis

      ## ------------- layout ----------------
      nplots <- length(mx.list)
      par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0.2, mgp = c(2.6, 1, 0), oma = c(1, 1, 0, 0) * 2)
      if (nplots > 24) {
        nc <- max(ceiling(nplots / 3), 6)
        par(mfrow = c(3, nc))
      } else if (FALSE && nplots <= 3) {
        nc <- 3
        par(mfrow = c(1, nc))
      } else {
        nc <- max(ceiling(nplots / 2), 6)
        par(mfrow = c(2, nc))
      }

      shiny::withProgress(message = "Computing volcano plots ...", value = 0, {
        i <- 1
        for (i in 1:nplots) {
          mx <- mx.list[[i]]
          is.sig <- (mx[, "qv"] <= fdr & abs(mx[, "fc"]) >= lfc)
          table(is.sig)
          sig.gs <- rownames(mx)[which(is.sig)]
          if (!is.null(sel.gsets)) sig.gs <- intersect(sel.gsets, sig.gs)

          gx.volcanoPlot.XY(
            x = mx[, "fc"], pv = mx[, "qv"],
            use.fdr = TRUE, p.sig = fdr, lfc = lfc,
            gene = rownames(mx),
            xlab = "effect size (NES)", lab.cex = 0, nlab = 0,
            render = "canvas", n = 1000, highlight = sig.gs,
            cex = 1, cex.axis = 1.3, cex.main = 1.4, axes = FALSE,
            ylim = c(0, ymax), main = ""
          )

          ## draw axis if first column or last row
          graphics::box(lwd = 1, col = "black", lty = "solid")
          is.first <- (i %% nc == 1)
          last.row <- ((i - 1) %/% nc == (nplots - 1) %/% nc)
          if (is.first) axis(2, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          if (last.row) axis(1, mgp = c(2, 0.7, 0), cex.axis = 0.8)
          legend("top",
            legend = names(mx.list)[i], box.lty = 0,
            x.intersp = 0.3, y.intersp = 0.5,
            inset = c(0, 0.01), cex = 1.2, bg = "white"
          )
          shiny::incProgress(1.0 / nplots)
        }
      })
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = volcanoAll.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 85),
      add.watermark = watermark
    )
  })
}
