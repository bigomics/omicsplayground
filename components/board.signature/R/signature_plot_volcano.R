##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
signature_plot_volcano_ui <- function(
  id,
  title,
  info.text,
  caption,
  height,
  width) {
  ns <- shiny::NS(id)
  
  PlotModuleUI(ns("plot"),
    title = title,
    info.text = info.text,
    caption = caption,
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
signature_plot_volcano_server <- function(id,
                                          pgx,
                                          sigCalculateGSEA,
                                          enrichmentContrastTable,
                                          selected_gxmethods,
                                          enrichmentGeneTable,
                                          getEnrichmentGeneTable,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    volcanoPlots.RENDER <- shiny::reactive({
      alertDataLoaded(session, pgx)
      if (is.null(pgx)) {
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

      ct <- colnames(pgx$model.parameters$contr.matrix)
      ct <- rownames(gsea$output)[ii]

      mm <- selected_gxmethods()
      meta <- playbase::pgx.getMetaMatrix(pgx, methods = mm)
      F <- meta$fc[, ct, drop = FALSE]
      qv <- meta$qv[, ct, drop = FALSE]
      score <- abs(F) * -log(qv)
      gset <- head(rownames(F), 100)
      gset <- intersect(gsea$gset, rownames(F))

      sel <- enrichmentGeneTable$rows_selected()
      sel.gene <- NULL
      if (length(sel)) {
        df <- getEnrichmentGeneTable()
        sel.gene <- df$gene[sel]
      }

      if (ncol(F) == 1) {
        cex.main <- 1.2
      } else {
        cex.main <- 1.2
        par(mfrow = c(2, 2), mar = c(2, 4, 3, 1), mgp = c(2.2, 0.8, 0))
      }
      if (ncol(F) > 4) {
        par(mfrow = c(3, 3), mar = c(1, 4, 3, 1), mgp = c(2.2, 0.8, 0))
        cex.main <- 1
      }
      if (ncol(F) > 9) {
        par(mfrow = c(4, 4), mar = c(0.2, 2, 3, 0.6))
        cex.main <- 0.9
      }
      for (i in 1:min(16, length(ct))) {
        gset2 <- head(gset[order(-score[gset, i])], 30)
        cex2 <- 0.8
        if (!is.null(sel.gene)) {
          gset2 <- sel.gene
          cex2 <- 1.3
        }
        xy <- cbind(fc = F[, i, drop = FALSE], z = -log10(qv[, i, drop = FALSE]))
        playbase::pgx.scatterPlotXY.BASE(
          xy,
          var = NULL, type = "factor", title = "",
          xlab = "differential expression (log2FC)",
          ylab = "significance (-log10q)",
          hilight = gset, hilight2 = gset2,
          cex = 0.9, cex.lab = cex2, cex.title = 1.0,
          legend = FALSE, col = c("grey80", "grey80"),
          opacity = 1
        )
        title(ct[i], cex.main = cex.main, line = 0.3)
      }
      p <- grDevices::recordPlot()
      p
    })

    PlotModuleServer(
      "plot",
      func = volcanoPlots.RENDER,
      res = c(90, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
