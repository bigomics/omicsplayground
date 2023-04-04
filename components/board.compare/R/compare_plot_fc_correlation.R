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
compare_plot_fc_correlation_ui <- function(id,
                                           height,
                                           width) {
  ns <- shiny::NS(id)
  info_text <- "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."

  PlotModuleUI(ns("plot"),
    title = "FC Correlation",
    plotlib = "base",
    label = "a",
    info.text = info_text,
    download.fmt = c("png", "pdf", "csv"),
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
compare_plot_fc_correlation_server <- function(id,
                                               pgx,
                                               dataset2,
                                               hilightgenes,
                                               input.contrast1,
                                               input.contrast2,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      pgx1 <- pgx
      pgx2 <- dataset2()
      ct1 <- head(names(pgx1$gx.meta$meta), 2)
      ct2 <- head(names(pgx2$gx.meta$meta), 2)
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()
      shiny::req(ct1)
      shiny::req(ct2)
      if (!all(ct1 %in% names(pgx1$gx.meta$meta))) {
        return(NULL)
      }
      if (!all(ct2 %in% names(pgx2$gx.meta$meta))) {
        return(NULL)
      }

      F1 <- playbase::pgx.getMetaMatrix(pgx1)$fc[, ct1, drop = FALSE]
      F2 <- playbase::pgx.getMetaMatrix(pgx2)$fc[, ct2, drop = FALSE]
      gg <- intersect(toupper(rownames(F1)), toupper(rownames(F2)))
      g1 <- rownames(F1)[match(gg, toupper(rownames(F1)))]
      g2 <- rownames(F2)[match(gg, toupper(rownames(F2)))]
      F1 <- F1[g1, , drop = FALSE]
      F2 <- F2[g2, , drop = FALSE]
      rownames(F2) <- rownames(F1) ## force same names if mouse .vs human
      colnames(F1) <- paste0("1:", colnames(F1))
      colnames(F2) <- paste0("2:", colnames(F2))
      return(cbind(F1, F2))
    })

    fcfcplot.RENDER <- function() {
      higenes <- hilightgenes()
      F <- plot_data()
      indexes <- substr(colnames(F), 1, 1)
      F1 <- F[, indexes == 1, drop = FALSE]
      F2 <- F[, indexes == 2, drop = FALSE]
      plot.SPLOM(F1, F2 = F2, cex = 0.3, cex.axis = 0.95, hilight = higenes)
      p
    }

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = fcfcplot.RENDER,
      csvFunc = plot_data,
      res = c(85, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
