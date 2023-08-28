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
    plotlib = "plotly",
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

      # Require inputs
      shiny::req(pgx$X)
      shiny::req(dataset2)
      shiny::req(input.contrast1)
      shiny::req(input.contrast2)
      pgx1 <- pgx
      pgx2 <- dataset2()
      ct1 <- input.contrast1()
      ct2 <- input.contrast2()

      # Allow only common contrats
      if (!all(ct1 %in% names(pgx1$gx.meta$meta))) {
        shiny::validate(shiny::need(all(ct1 %in% names(pgx2$gx.meta$meta)), "Warning: No common contrasts."))
        return(NULL)
      }
      if (!all(ct2 %in% names(pgx2$gx.meta$meta))) {
        shiny::validate(shiny::need(all(ct2 %in% names(pgx2$gx.meta$meta)), "Warning: No common contrasts."))
        return(NULL)
      }

      # Match matrices from both datasets
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

    plot_interactive_comp_fc <- function(plot_data, hilight = NULL, cex = 0.5, cex.axis = 1, cex.space = 0.2) {
      
      var <- plot_data()
      pos <- plot_data()

      p <- playbase::pgx.scatterPlotXY(
        pos,
        var = var,
        plotlib = "plotly", 
        cex = cex,
        hilight = hilight,
        key = colnames(var)
      )
      
      return(p)
    }

    fcfcplot.RENDER <- function() {
      higenes <- hilightgenes()
      p <- plot_interactive_comp_fc(plot_data = plot_data, cex = 0.6, cex.axis = 0.95, hilight = higenes) %>%
        plotly::layout(
          dragmode = "select",
          margin = list(l = 5, r = 5, b = 5, t = 20)
        )
      return(p)
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = fcfcplot.RENDER,
      csvFunc = plot_data,
      res = c(85, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
