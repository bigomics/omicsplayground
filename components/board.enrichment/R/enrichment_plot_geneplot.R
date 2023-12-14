##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Create UI for gene plot enrichment visualization
#'
#' @param id Id prefix for namespace
#' @param title Plot title
#' @param caption Plot caption
#' @param info.text Info text to display
#' @param height Plot height
#' @param width Plot width
#'
#' @return Shiny UI for gene plot enrichment visualization
enrichment_plot_geneplot_ui <- function(
    id,
    title,
    caption,
    info.text,
    height,
    width) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("gs_ungroup2"), "ungroup samples", FALSE),
      "Ungroup samples in the plot",
      placement = "top",
      options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "c",
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

#' Gene Expression Plot Server Function
#'
#' Server function for generating the gene expression plot module in Shiny.
#'
#' @param id Module id
#' @param pgx PGX object
#' @param gs_contrast Get the selected contrast function
#' @param gene_selected Get the selected gene function
#' @param subplot.MAR Plot margins
#' @param watermark Add watermark
#'
#' @return Shiny module server function
enrichment_plot_geneplot_server <- function(id,
                                            pgx,
                                            gs_contrast,
                                            gene_selected,
                                            subplot.MAR,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    render_subplot_geneplot <- function() {
      par(mfrow = c(1, 1), mgp = c(1.8, 0.8, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      shiny::req(pgx$X)

      comp0 <- colnames(pgx$model.parameters$contr.matrix)[1]
      comp0 <- gs_contrast()

      has.design <- !is.null(pgx$model.parameters$design)
      collapse.others <- ifelse(has.design, FALSE, TRUE)

      sel <- gene_selected()
      if (is.null(sel) || is.na(sel) || length(sel) == 0) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") %>%
          plotly::config(
            displayModeBar = FALSE
          ))
      } else {
        probe <- sel$rn
        gene <- sel$gene
        if (length(probe) > 1) {
          probe <- grep("\\[gx|\\[mrna", probe, value = TRUE)
        }
        ngrp <- length(unique(pgx$samples$group))
        grouped <- TRUE
        grouped <- !input$gs_ungroup2
        srt <- ifelse(!grouped || ngrp > 4, 30, 0)
        if (!grouped && ncol(pgx$X) > 15) srt <- 60
        playbase::pgx.plotExpression(
          pgx,
          probe,
          comp = comp0,
          logscale = TRUE,
          level = "gene",
          collapse.others = collapse.others,
          grouped = grouped,
          srt = srt,
          main = "",
          xlab = gene,
          plotlib = "plotly",
        )
      }
    }

    ## -------------------------------------------------------------------
    ## Render plots
    ## -------------------------------------------------------------------

    subplot_geneplot.RENDER1 <- function() {
      render_subplot_geneplot() %>%
        plotly_default()
    }

    subplot_geneplot.RENDER2 <- function() {
      render_subplot_geneplot() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = subplot_geneplot.RENDER1,
      func2 = subplot_geneplot.RENDER2,
      pdf.width = 5, pdf.height = 5,
      res = c(78, 100),
      add.watermark = watermark
    )
  })
}
