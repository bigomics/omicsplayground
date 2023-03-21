##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
connectivity_plot_cmap_FCFCplots_ui <- function(id,
                                                label = "",
                                                rowH = 660) {
  ns <- shiny::NS(id)
  info_text <- strwrap(
    "<b>FC scatter plots.</b> Scatter plots of gene expression foldchange
    values between two contrasts. Foldchanges that are similar show high
    correlation, i.e. are close to the diagonal. You can switch to enrichment
    type plots in the plot settings."
  )

  plot_opts <- shiny::tagList(
    shiny::radioButtons(
      ns("fcfc_plottype"),
      "Plot type:",
      c("scatter", "enrichment"),
      inline = TRUE
    )
  )

  PlotModuleUI(ns("plot"),
    title = "FC scatter plots",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = plot_opts,
    height = c(360, 600),
    width = c("auto", 1280)
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
connectivity_plot_cmap_FCFCplots_server <- function(id,
                                                    inputData,
                                                    cmap_contrast,
                                                    getCurrentContrast,
                                                    getTopProfiles,
                                                    getConnectivityScores,
                                                    watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      cmap_FCFCscatter <- function(fc, F, mfplots, ylab) {
        mfplots <- c(4, 5)
        ## get the foldchanges of selected comparison and neighbourhood
        F0 <- F
        F[is.na(F)] <- 0 ## really??
        names(fc) <- toupper(names(fc))
        gg <- intersect(names(fc), rownames(F)) ## uppercase for MOUSE
        fc <- fc[gg]

        nplots <- mfplots[1] * mfplots[2]
        F <- F[gg, 1:min(nplots, ncol(F)), drop = FALSE]
        F0 <- F0[gg, colnames(F), drop = FALSE]
        i <- 1
        par(
          mfrow = mfplots, mar = c(5.1, 1.6, 0.2, 0.5),
          mgp = c(2.6, 0.7, 0), oma = c(0, 3, 0, 0)
        )
        i <- 1
        for (i in 1:ncol(F)) {
          ct1 <- colnames(F)[i]
          ct1x <- sub("\\]", "]\n", ct1)
          nna <- (is.na(fc) | is.na(F0[, ct1]))
          col <- c("grey15", "grey70")[1 + nna]
          base::plot(F[, ct1], fc,
            pch = 20, cex = 0.5,
            cex.lab = 0.9, cex.axis = 0.9,
            xlab = ct1x, ylab = "", col = col
          )
          abline(v = 0, h = 0, lty = 2, lwd = 0.5)
          abline(lm(fc ~ F0[, ct1]), col = "red")
          if (i %% mfplots[2] == 1) {
            mtext(ylab, 2, line = 3, cex = 0.60)
          }
        }
      }

      cmap_FCFCenplot <- function(fc, F, mfplots, ylab, res) {
        mfplots <- c(4, 5)
        names(fc) <- toupper(names(fc))
        nplots <- mfplots[1] * mfplots[2]
        i <- 1
        par(mfrow = mfplots, mar = c(0.1, 4, 2.6, 1))
        for (i in 1:min(ncol(F), nplots)) {
          j1 <- head(order(F[, i]), 100)
          j2 <- head(order(-F[, i]), 100)
          gset.dn <- rownames(F)[j1]
          gset.up <- rownames(F)[j2]
          gset.both <- c(gset.dn, gset.up)
          rnk <- fc
          pw <- colnames(F)[i]
          gsea.enplot(abs(rnk), gset.both,
            xlab = "",
            main = pw, cex.main = 0.8, len.main = 32
          )
          R <- res[match(pw, res$pathway), , drop = FALSE]
          legend("topright",
            cex = 0.75, y.intersp = 0.85, bty = "n",
            c(
              paste("NES=", round(R$NES[1], 3)),
              paste("padj=", round(R$padj[1], 4))
            )
          )
        }
      }

      plot_data <- shiny::reactive({
        res <- list(
          pgx = inputData(),
          cmap_contrast = cmap_contrast()
        )
        return(res)
      })

      plot_RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        cmap_contrast <- res$cmap_contrast

        shiny::req(pgx, cmap_contrast)

        res1 <- getCurrentContrast()
        fc <- res1$fc
        ct <- res1$name
        F <- getTopProfiles()
        if (NCOL(F) == 0) {
          return(NULL)
        }
        F <- F[, 1:min(ncol(F), 10), drop = FALSE]

        if (input$fcfc_plottype == "scatter") {
          mfplots <- c(2, 5)
          cmap_FCFCscatter(fc, F, mfplots, ylab = ct)
        } else {
          mfplots <- c(3, 4)
          df <- getConnectivityScores()
          cmap_FCFCenplot(fc, F, mfplots, ylab, df)
        }
        p <- grDevices::recordPlot()
        p
      })

      PlotModuleServer(
        "plot",
        plotlib = "base",
        func = plot_RENDER,
        func2 = plot_RENDER,
        csvFunc = plot_data,
        res = c(90, 110),
        pdf.height = 4.5, pdf.width = 10,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
