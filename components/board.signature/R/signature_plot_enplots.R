##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics Sagl. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
signature_plot_enplots_ui <- function(id, height, width) {
  ns <- shiny::NS(id)
  info_text <- "<b>Enrichment plots.</b> Enrichment of the query signature in all constrasts. Positive enrichment means that this particular contrast shows similar expression changes as the query signature."

  PlotModuleUI(ns("plot"),
    title = "Enrichment plots",
    info.text = info_text,
    download.fmt = c("png", "pdf"),
    height = height,
    width = width
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
signature_plot_enplots_server <- function(id,
                                          inputData,
                                          sigCalculateGSEA,
                                          enrichmentContrastTable,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    enplots.RENDER <- shiny::reactive({
      ngs <- inputData()
      alertDataLoaded(session, ngs)
      if (is.null(ngs)) {
        return(NULL)
      }

      gsea <- sigCalculateGSEA()
      if (is.null(gsea)) {
        return(NULL)
      }

      ## filter with table selection/search
      ii <- enrichmentContrastTable$rows_selected()
      if (is.null(ii)) {
        ii <- enrichmentContrastTable$rows_all()
      }

      shiny::req(ii)
      ct <- rownames(gsea$output)[ii]
      F <- as.matrix(gsea$F[, ct, drop = FALSE])
      qv <- gsea$output[ct, "q"]
      gset <- gsea$gset

      if (ncol(F) == 1) {
        cex.main <- 1.1
        nc <- 1
      }
      if (ncol(F) > 1) {
        cex.main <- 1.1
        nc <- 3
        par(mfrow = c(4, 3), mar = c(0.3, 3, 3, 0.5), mgp = c(1.9, 0.7, 0), oma = c(0, 1, 0, 0))
      }
      if (ncol(F) > 12) {
        par(mfrow = c(5, 4), mar = c(0.2, 2, 3, 0.6))
        cex.main <- 0.9
        nc <- 4
      }
      for (i in 1:min(20, ncol(F))) {
        f <- colnames(F)[i]
        tt <- sub(".*\\]", "", f)
        tt <- breakstring(substring(tt, 1, 50), 28, force = TRUE)
        ylab <- ""
        if (i %% nc == 1) ylab <- "rank metric"
        gsea.enplot(F[, i], gset,
          main = tt, cex.main = cex.main,
          xlab = "", ylab = ylab
        )
        qv1 <- paste("q=", round(qv[i], digits = 3))
        legend("topright", qv1, cex = 0.9, bty = "n", adj = 0)
        if (grepl("^\\[", f)) {
          db <- sub("\\].*", "]", colnames(F)[i])
          legend("topleft", db, cex = 0.9, bty = "n", adj = 0)
        }
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = enplots.RENDER,
      res = c(90, 90), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
