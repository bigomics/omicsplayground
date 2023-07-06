##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanomethods_ui <- function(
    id,
    title,
    info.text,
    caption,
    height,
    width) {
  ns <- shiny::NS(id)


  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_method"), "scale per method", FALSE),
      "Scale the volcano plots individually per method..",
      placement = "right", options = list(container = "body")
    )
  )
  
  PlotModuleUI(
    ns("plot"),
    plotlib = "grid",
    title = title,
    caption = caption,
    options = plot_options,
    info.text = info.text,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

enrichment_plot_volcanomethods_server <- function(id,
                                                  pgx,
                                                  gs_features,
                                                  gs_contrast,
                                                  gs_fdr,
                                                  gs_lfc,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx, gs_features(), gs_contrast())

      cmp <- gs_contrast()
      mx <- pgx$gset.meta$meta[[cmp]]
      F <- unclass(mx$fc)
      Q <- unclass(mx$q)
      rownames(Q) <- rownames(F) <- rownames(mx)
      F[which(is.infinite(F))] <- NA
      Q[which(is.na(Q))] <- 1

      fdr <- as.numeric(gs_fdr())
      lfc <- as.numeric(gs_lfc())
      sel.gsets <- playdata::COLLECTIONS[[gs_features()]]

      nlq <- -log10(1e-99 + unlist(Q))
      ymax <- max(3, 1.2 * quantile(nlq, probs = 0.999, na.rm = TRUE)[1]) ## y-axis
      xmax <- quantile(abs(unlist(F)), probs = 0.999, na.rm = TRUE)[1]

      pd <- list(
        F = F,
        Q = Q,
        sel.gsets = sel.gsets,
        xmax = xmax,
        ymax = ymax,
        fdr = fdr,
        lfc = lfc
      )
      pd
    })

    get_ggplots <- function(cex = 1, base_size = 12) {
      pd <- plot_data()
      shiny::req(pd)
      F <- pd$F
      Q <- pd$Q
      mx.list <- pd$mx.list
      ymax <- pd$ymax
      xmax <- pd$xmax
      nplots <- ncol(Q)
      fdr <- pd$fdr
      lfc <- pd$lfc
      sel.gsets <- pd$sel.gsets

      shiny::withProgress(message = "Computing volcano plots ...", value = 0, {
        plt <- list()
        i <- 1
        for (i in 1:nplots) {
          fc <- F[, i]
          qv <- Q[, i]
          is.sig <- (qv <= fdr & abs(fc) >= lfc)
          table(is.sig)
          sig.genes <- names(fc)[which(is.sig)]
          if (!is.null(sel.gsets)) sig.genes <- intersect(sel.gsets, sig.genes)

          xy <- cbind(x = fc, y = -log10(qv))
          tt <- colnames(Q)[i]
          #

          ymax1 <- ymax
          if(input$scale_per_method) {
            ymax1 <- 1.2*quantile(xy[,2], probs = 0.999, na.rm = TRUE)[1] ## y-axis
          }
          
          plt[[i]] <- playbase::pgx.scatterPlotXY.GGPLOT(
            xy,
            title = tt,
            cex.title = 0.75,
            var = is.sig,
            type = "factor",
            col = c("#bbbbbb", "#1e60bb"),
            legend.pos = "none", #


            hilight2 = NULL,
            xlim = xmax * c(-1, 1),
            ylim = c(0, ymax1),
            xlab = "difference  (log2FC)",
            ylab = "significance  (-log10q)",
            hilight.lwd = 0,
            hilight.col = "#1e60bb",
            hilight.cex = 1.5,
            cex = cex,
            cex.lab = 1.8 * cex,
            base_size = base_size
          ) + ggplot2::theme_bw(base_size = base_size)

          shiny::incProgress(1 / nplots)
        }
      })
      return(plt)
    }

    volcano.RENDER <- function() {
      plt <- get_ggplots(cex = 0.4, base_size = 12)
      shiny::req(plt)
      nplots <- length(plt)
      nr <- ifelse(nplots <= 6, 1, 2)
      nc <- max(ceiling(nplots / nr), 4)
      #
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }

    volcano.RENDER2 <- function() {
      plt <- get_ggplots(cex = 0.8, base_size = 18)
      shiny::req(plt)
      nplots <- length(plt)
      nr <- ifelse(nplots <= 5, 1, 2)
      nc <- max(ceiling(nplots / nr), 3)
      #
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }


    PlotModuleServer(
      "plot",
      plotlib = "grid",
      func = volcano.RENDER,
      func2 = volcano.RENDER2,
      pdf.width = 10, pdf.height = 5,
      res = c(75, 90),
      add.watermark = watermark
    )
  })
}
