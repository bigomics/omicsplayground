##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' Heatmap plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
biomarker_plot_heatmap_ui <- function(id,
                                      label = "",
                                      height = c(600, 800)) {
  ns <- shiny::NS(id)
  info_text <- strwrap("<b>Biomarker heatmap.</b> Expression heatmap
                      of top gene features according to their variable
                      importance.")

  PlotModuleUI(ns("plot"),
    title = "Heatmap",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
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
biomarker_plot_heatmap_server <- function(id,
                                          calcVariableImportance,
                                          pgx,
                                          pdx_predicted,
                                          watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        shiny::req(pgx)

        res <- calcVariableImportance()

        if (is.null(res)) {
          return(NULL)
        }

        gg <- rownames(res$X)
        gg <- intersect(gg, rownames(pgx$X))
        X <- pgx$X[gg, ]

        X <- head(X[order(-apply(X, 1, sd)), ], 40) ## top50

        splitx <- NULL
        ct <- pdx_predicted()
        do.survival <- grepl("survival", ct, ignore.case = TRUE)

        splitx <- pgx$Y[colnames(X), ct]
        if (!is.categorical(splitx) || do.survival) {
          splitx <- NULL
        }

        rownames(X) <- substring(rownames(X), 1, 40)
        annot <- pgx$Y[colnames(X), ]
        sdx <- apply(X, 1, sd)

        res <- list(X = X, splitx = splitx)
      })

      plot.RENDER <- shiny::reactive({
        res <- plot_data()
        shiny::req(res)

        X <- res$X
        splitx <- res$splitx

        gx.splitmap(X,
          split = NULL, splitx = splitx, main = "  ",
          dist.method = "euclidean",
          show_colnames = FALSE, ## save space, no sample names
          show_legend = ifelse(is.null(splitx), TRUE, FALSE),
          key.offset = c(0.05, 1.03),
          show_rownames = 99,
          lab.len = 50, cexRow = 0.88, mar = c(2, 8)
        )
      })

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        csvFunc = plot_data,
        res = c(72, 435),
        pdf.width = 10, pdf.height = 10,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
