##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
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
biomarker_plot_importance_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    caption = caption,
    label = label,
    plotlib = "base",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
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
biomarker_plot_importance_server <- function(id,
                                             calcVariableImportance,
                                             is_computed,
                                             watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        # Return NULL in case of error instead of
        # placing a shiny::req ; that way we can handle
        # the case on the plotting function and display a
        # help message.
        out <- tryCatch(
          {
            calcVariableImportance()
          },
          error = function(w) {
            NULL
          }
        )

        if (is.null(out)) {
          return(NULL)
        }
        res <- list(
          R = out$R
        )
        return(res)
      })

      plot.RENDER <- function() {
        res <- plot_data()

        ## if (is.null(res) || length(res) == 0) {
        ##   frame()
        ##   text(0.5, 0.5, "Please compute desired output on the right Settings tab", col = "grey50")
        ##   return()
        ## }
        shiny::validate(shiny::need(is_computed(), "Please select target class and run 'Compute'"))
        shiny::req(res)

        R <- res$R
        R <- R[order(-rowSums(R, na.rm = TRUE)), , drop = FALSE]
        R <- pmax(R, 0.05)

        par(mfrow = c(1, 1), oma = c(1, 1, 1, 1) * 0.2)
        par(mar = c(5, 3, 0, 1))
        R.top <- head(R, 40)
        barplot(t(R.top),
          las = 3, horiz = FALSE,
          cex.names = 0.85, ylab = "cumulative importance"
        )
        klr <- grey.colors(ncol(R))
        legend("topright",
          legend = rev(colnames(R)), fill = rev(klr),
          cex = 0.8, y.intersp = 0.8
        )
      }

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        csvFunc = plot_data,
        res = c(70, 140),
        pdf.width = 10, pdf.height = 5,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
