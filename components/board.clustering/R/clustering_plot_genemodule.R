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
#' @param width
#'
#' @export
clustering_plot_genemodule_ui <- function(
    id,
    title,
    caption,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
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
#' @param comp
#' @param pgx
#' @param res
#' @param ii
#' @param watermark
#'
#'
#'
#' @export
clustering_plot_genemodule_server <- function(id,
                                              getTopMatrix,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      res <- getTopMatrix()
      shiny::req(res)

      mat <- res$mat
      idx <- res$idx
      modx <- sapply(1:ncol(mat), function(i) tapply(mat[, i], idx, mean))
      colnames(modx) <- colnames(mat)

      return(list(
        mat = modx
      ))
    })


    render_plotly <- function(title.cex = 1, title.y = 0.95) {
      pd <- plot_data()
      shiny::req(pd)

      mat <- pd$mat
      plotly.colors <- omics_pal_d("muted_light")(8)
      plotly.colors <- rep(plotly.colors, 10)

      nplots <- nrow(mat)
      plts <- list()

      for (i in 1:nplots) {
        gx <- mat[i, ]
        names(gx) <- colnames(mat)

        anntitle <- list(
          x = 0.5, y = title.y,
          xref = "paper", yref = "paper",
          xanchor = "center", yanchor = "bottom",
          text = rownames(mat)[i],
          font = list(size = 18 * title.cex),
          align = "center", showarrow = FALSE
        )

        p <- playbase::pgx.barplot.PLOTLY(
          data = data.frame(
            gx = gx,
            xgroup = factor(names(gx), levels = names(gx))
          ),
          x = "xgroup",
          y = "gx",
          grouped = FALSE,
          fillcolor = plotly.colors[(i - 1) %% 10 + 1],
          yaxistitle = "avg expr (log2)",
          xaxistitle = "",
          annotations = anntitle
        ) %>% plotly::layout(
          plot_bgcolor = "#f2f2f2",
          margin = list(l = 20, r = 0, b = 0, t = 0),
          bargap = 0.35
        )
        plts[[i]] <- p
      }

      return(plts)
    }

    plotly.RENDER <- function() {
      ## layout in subplots
      plts <- render_plotly(title.cex = 0.92, title.y = 0.85)
      plts <- head(plts, 18)
      ncols <- ifelse(length(plts) > 4, 2, 1)
      ncols <- ifelse(length(plts) > 12, 3, ncols)
      nrows <- ceiling(length(plts) / ncols)
      plotly::subplot(
        plts,
        nrows = nrows,
        margin = c(0.02, 0.02, 0.03, 0.03), ## lrtb
        titleX = TRUE,
        titleY = TRUE,
        shareY = TRUE,
        shareX = TRUE
      ) %>%
        plotly::layout(
          font = list(size = 12),
          xaxis = list(tickfont = list(size = 14)),
          yaxis = list(tickfont = list(size = 10)),
          margin = list(l = 10, r = 10, b = 10, t = 10),
          showlegend = FALSE
        )
    }

    modal_plotly.RENDER <- function() {
      plts <- render_plotly(title.cex = 1.5, title.y = 0.88)
      plts <- head(plts, 18)
      nrows <- ifelse(length(plts) > 3, 2, 1)
      nrows <- ifelse(length(plts) > 8, 3, nrows)
      ncols <- ceiling(length(plts) / nrows)
      fig <- plotly::subplot(
        plts,
        nrows = nrows,
        margin = c(0.02, 0.02, 0.04, 0.04), ## lrtb
        titleX = TRUE,
        titleY = TRUE,
        shareY = TRUE,
        shareX = FALSE
      )

      fig <- fig %>%
        plotly::layout(
          font = list(size = 18),
          xaxis = list(tickfont = list(size = 22)),
          yaxis = list(tickfont = list(size = 14)),
          margin = list(l = 10, r = 10, b = 10, t = 10),
          showlegend = FALSE
        ) %>%
        plotly::style(
          marker.size = 20
        )
      fig
    }

    PlotModuleServer(
      "pltmod",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      plotlib = "plotly",
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 105), ## resolution of plots
      pdf.width = 14,
      pdf.height = 3.5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
