##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_compare_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    caption = caption,
    info.text = info.text,
    cards = TRUE,
    plotlib = c("grid", "base"),
    download.fmt = c("png", "pdf"),
    card_names = c("Volcano", "Ranking")
  )
}

enrichment_plot_compare_server <- function(id,
                                           pgx,
                                           gs_contrast,
                                           gs_features,
                                           gs_statmethod,
                                           gs_fdr,
                                           gs_lfc,
                                           gset_selected,
                                           calcGsetMeta,
                                           selected_gsetmethods,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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
      sel.gsets <- playdata::COLLECTIONS[[gs_features()]]

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

    render_compare <- function() {
      shiny::req(pgx, gs_contrast())

      comp <- 1
      comp <- gs_contrast()
      if (is.null(comp)) {
        return(NULL)
      }

      gset <- rownames(pgx$gsetX)[1]
      gset <- gset_selected()
      if (is.null(gset)) {
        frame()
        text(0.5, 0.5, "Please select a geneset", col = "grey50")
        return()
      }
      gset <- gset[1]

      score <- sapply(pgx$gset.meta$meta, function(x) x[gset, "meta.fx"])

      top.up <- names(sort(score[which(score > 0)], decreasing = TRUE))
      top.dn <- names(sort(score[which(score < 0)]))
      genes <- names(which(pgx$GMT[, gset] != 0))
      genes <- toupper(sub(".*:", "", genes))
      gx.meta <- pgx$gx.meta$meta

      gsmethods <- selected_gsetmethods()

      par(mfrow = c(2, 5), mar = c(0.5, 3.2, 2.6, 0.5), mgp = c(2, 0.8, 0))
      i <- 1
      for (i in 1:5) {
        if (i > length(top.up)) {
          frame()
        } else {
          cmp <- top.up[i]
          rnk0 <- gx.meta[[cmp]]$meta.fx
          names(rnk0) <- rownames(gx.meta[[1]])
          names(rnk0) <- toupper(sub(".*:", "", names(rnk0)))

          gs.meta <- pgx$gset.meta$meta[[cmp]]
          qv0 <- max(gs.meta[gset, "q"][, gsmethods], na.rm = TRUE)

          gs1 <- playbase::breakstring(gset, 28, 50, force = FALSE)
          cmp <- paste0(gset, "\n@", cmp)
          playbase::gsea.enplot(rnk0, genes,
            names = NULL, ## main=gs,
            main = cmp, xlab = "",
            cex.main = 0.80, len.main = 72
          )
          qv1 <- formatC(qv0, format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }
      for (i in 1:5) {
        if (i > length(top.dn)) {
          frame()
        } else {
          cmp <- top.dn[i]
          rnk0 <- gx.meta[[cmp]]$meta.fx
          names(rnk0) <- rownames(gx.meta[[1]])
          names(rnk0) <- toupper(sub(".*:", "", names(rnk0)))

          gs.meta <- pgx$gset.meta$meta[[cmp]]
          qv0 <- max(gs.meta[gset, "q"][, gsmethods], na.rm = TRUE)

          gs1 <- playbase::breakstring(gset, 28, 50, force = FALSE)
          cmp <- paste0(gset, "\n@", cmp)
          playbase::gsea.enplot(rnk0, genes,
            names = NULL, ## main=gs,
            main = cmp, xlab = "",
            cex.main = 0.80, len.main = 72
          )
          qv1 <- formatC(qv0, format = "e", digits = 2)
          legend("topright", paste("q=", qv1), bty = "n", cex = 0.85)
        }
      }
    }

    compare.RENDER <- function() {
      render_compare()
    }
    
    compare.RENDER2 <- function() {
      render_compare()
    }

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
            cex.lab = 1.8*cex
          ) + ggplot2::theme_bw()
          
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
    
    plot_grid <- list(
      list(plotlib = "grid", func = volcano.RENDER, card = 1),
      list(plotlib = "base", func = compare.RENDER, card = 2)
      )
    
    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        # csvFunc = plot_data_hm1,
        res = c(80, 95), # resolution of plots
        pdf.width = 10, pdf.height = 8,
        add.watermark = watermark,
        card = x$card
      )
    })

    
    # PlotModuleServer(
    #   "plot",
    #   func = compare.RENDER,
    #   func2 = compare.RENDER,      
    #   pdf.width = 5, pdf.height = 5,
    #   res = c(95, 100),
    #   add.watermark = watermark
    # )
  })
}
