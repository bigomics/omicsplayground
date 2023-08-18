##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
biomarker_plot_heatmap_ui <- function(
    id,
    title,
    info.text,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot_options <- tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_all"), "show all samples",
        value = FALSE
      ), "Show all samples or only selected.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "base",
    caption = caption,
    info.text = info.text,
    options = plot_options,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
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
      ## return data structure for plots
      plot_data <- shiny::reactive({
        shiny::req(pgx$X)

        res <- calcVariableImportance()

        if (is.null(res)) {
          return(NULL)
        }

        gg <- rownames(res$X)
        gg <- intersect(gg, rownames(pgx$X))
        if (input$show_all) {
          kk <- colnames(pgx$X)
        } else {
          kk <- colnames(res$X)
        }
        X <- pgx$X[gg, kk]
        X <- head(X[order(-apply(X, 1, sd)), ], 40) ## top50

        splitx <- NULL
        ct <- pdx_predicted()
        do.survival <- grepl("survival", ct, ignore.case = TRUE)

        splitx <- pgx$Y[colnames(X), ct]
        if (!playbase::is.categorical(splitx) || do.survival) {
          splitx <- NULL
        }

        rownames(X) <- substring(rownames(X), 1, 40)
        annot <- pgx$Y[colnames(X), ]
        sdx <- apply(X, 1, sd)

        res <- list(X = X, splitx = splitx)
      })

      plot.RENDER <- function() {
        res <- plot_data()
        shiny::req(res)

        X <- res$X
        splitx <- res$splitx

        playbase::gx.splitmap(X,
          split = NULL, splitx = splitx, main = "  ",
          dist.method = "euclidean",
          show_colnames = FALSE, ## save space, no sample names
          show_legend = ifelse(is.null(splitx), TRUE, FALSE),
          key.offset = c(0.05, 0.98),
          show_rownames = 99,
          lab.len = 50,
          cexRow = 0.88,
          mar = c(2, 8)
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        csvFunc = plot_data,
        res = c(72, 110),
        pdf.width = 10,
        pdf.height = 10,
        add.watermark = watermark,
        filename = "biomarker_plot_heatmap"
      )
    } ## end of moduleServer
  )
}
