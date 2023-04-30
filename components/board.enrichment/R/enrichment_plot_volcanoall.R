##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanoall_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    plotlib = "grid",
    title = title,
    caption = caption,
    info.text = info.text,
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

    plot_data <- shiny::reactive({

      shiny::req(pgx)
      if (is.null(gs_features())) {
        return(NULL)
      }

      meta     <- pgx$gset.meta$meta
      gsmethod <- colnames(meta[[1]]$fc)
      gsmethod <- gs_statmethod()
      if (is.null(gsmethod) || length(gsmethod) == 0) {
        return(NULL)
      }

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
      F <- lapply(mx.list, function(mx) mx[, "fc"])
      names(Q) <- names(mx.list)
      names(F) <- names(mx.list)

      ## select maximum 24 comparisons (because of space...)
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

    get_ggplots <- function(cex=1, base_size=12) {

      pd <- plot_data()
      shiny::req(pd)
      F <- pd$F
      Q <- pd$Q
      mx.list <- pd$mx.list
      ymax <- pd$ymax
      xmax <- pd$xmax
      nplots <- length(Q)
      fdr <- pd$fdr
      lfc <- pd$lfc
      sel.gsets <- pd$sel.gsets
      
      shiny::withProgress(message = "Computing volcano plots ...", value = 0, {
        i <- 1
        plt <- list()
        for (i in 1:nplots) {

          fc <- F[[i]]
          qv <- Q[[i]]
          is.sig1 <- (qv <= fdr & abs(fc) >= lfc)
          table(is.sig1)
          sig.genes <- names(fc)[which(is.sig1)]
          if (!is.null(sel.gsets)) sig.genes <- intersect(sel.gsets, sig.genes)
          
          xy <- cbind(x = fc, y = -log10(qv))
          tt <- names(F)[i]
          ## xmax <- max(abs(mx[,"fc"]))
          
          plt[[i]] <- playbase::pgx.scatterPlotXY.GGPLOT(
            xy,
            title = tt,
            cex.title = 0.75,
            var = is.sig1,
            type = "factor",
            col = c("#bbbbbb", "#1e60bb"),
            legend.pos = "none", ## plotlib="ggplot",
            hilight = NULL,
            # hilight2 = sig.genes,
            hilight2 = NULL,
            xlim = xmax * c(-1, 1),
            ylim = c(0, ymax),
            xlab = "difference  (log2FC)",
            ylab = "significance  (-log10q)",
            hilight.lwd = 0,
            hilight.col = "#1e60bb",
            hilight.cex = 1.5,
            cex = cex,
            cex.lab = 1.8*cex,
            base_size = base_size
          ) + ggplot2::theme_bw(base_size = base_size)
          
          shiny::incProgress(1.0 / nplots)
        }
      })
      return(plt)
    }

  
    volcano.RENDER <- function() {
      plt <- get_ggplots(cex=0.4, base_size=12)
      shiny::req(plt)    
      ## ------------- layout ----------------
      nplots <- length(plt)
      nc <- max(4,nplots)
      nr <- 1
      if (nplots > 6) {
        nc <- ceiling(nplots / 2)
        nr <- 2
      }
      if (nplots > 12 ) {
        nc <- ceiling(nplots / 3)
        nr <- 3
      }
      ## if(nr*nc > nplots) nplots <- c(nplots, rep(gridExtra::blank, nr*nc - nplots))
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }
    
    volcano.RENDER2 <- function() {
      plt <- get_ggplots(cex=0.9, base_size=16)
      shiny::req(plt)    
      ## ------------- layout ----------------
      nplots <- length(plt)
      nr = nc = 1
      nc <- 3
      nr <- 1
      if (nplots > 3) {
        nc <- ceiling(nplots / 2)
        nr <- 2
      }
      if (nplots > 8) {
        nc <- ceiling(nplots / 3)
        nr <- 3
      }    
      if (nplots > 15 ) {
        nc <- ceiling(nplots / 4)
        nr <- 4
      } 
      ##if(nr*nc > nplots) nplots <- c(nplots, rep(gridExtra::blank, nr*nc - nplots))
      gridExtra::grid.arrange(grobs = plt, nrow = nr, ncol = nc)
    }
    
    PlotModuleServer(
      "plot",
      plotlib = "grid",
      func = volcano.RENDER,
      func2 = volcano.RENDER2,
      pdf.width = 10,
      pdf.height = 5,
      res = c(72, 85),
      add.watermark = watermark
    )
    
  })  ## end module-server
} ## server
