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
functional_plot_enrichmap_ui <- function(
  id,
  label = "",
  title,
  info.text,
  caption,
  info.width,
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    id = ns("plotmodule"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    info.width = info.width,
    options = NULL,
    download.fmt = c("png", "pdf", "svg"),
    height = height,
    width = width
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
functional_plot_enrichmap_server <- function(id,
                                             pgx,
                                             fa_contrast,
                                             watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      #

      ## reactive or function? that's the question...
      plot_data <- shiny::reactive({
        res <- compute_enrichmentmap(
          pgx,
          qsig = 0.99,
          ntop = 120,
          wt = 1,
          plot = FALSE
        )
        res$contrast <- fa_contrast()
        names(res)
        return(res)
      })

      get_plots <- function(cex, lwd) {
        res <- plot_data()
        ct <- head(colnames(res$F), 6)
        ct

        if (!interactive()) {
          progress <- shiny::Progress$new()
          on.exit(progress$close())
          progress$set(message = "Calculating feature-set scores", value = 0)
        }

        plist <- list()
        for (i in 1:length(ct)) {
          plist[[i]] <- plot_enrichmentmap(
            res,
            contrast = i,
            qsig = 0.05,
            cex = cex,
            lwd = lwd,
            title = ct[i],
            title.y = 0.08,
            title.x = 0.02,
            plot_bgcolor = "#cdceeb66",
            label = FALSE
          )
          if (!interactive()) shiny::incProgress(1 / length(ct))
        }
        plist
      }

      plot_RENDER <- function() {
        plist <- get_plots(cex = 0.45, lwd = 0.7)
        nr <- min(2, ceiling(length(plist) / 2))
        plotly::subplot(head(plist, 6), nrows = nr, margin = 0.025) %>%
          plotly::layout(
            #
            margin = list(l = 0, r = 0, b = 0, t = 0, pad = 10)
          )
      }

      plot_RENDER2 <- function() {
        plist <- get_plots(cex = 0.6, lwd = 0.9)
        nr <- min(2, ceiling(length(plist) / 2))
        plotly::subplot(head(plist, 6), nrows = nr, margin = 0.025) %>%
          plotly::layout(
            #
            margin = list(l = 0, r = 0, b = 0, t = 0, pad = 10)
          )
      }

      PlotModuleServer(
        "plotmodule",
        plotlib = "plotly",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        #
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
