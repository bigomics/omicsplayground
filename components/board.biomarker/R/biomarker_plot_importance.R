##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Importance plot UI input function
#' @description A shiny Module for plotting (UI code).
#' @param id
#' @param label
#' @param height
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
  width
) {
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
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

#' Importance plot Server function
#' @description A shiny Module for plotting (server code).
#' @param id
#' @return
#' @export
biomarker_plot_importance_server <- function(id,
                                             pgx,
                                             calcVariableImportance,
                                             is_computed,
                                             watermark = FALSE) {

  moduleServer(

    id, function(input, output, session) {

      plot_data <- shiny::reactive({
        out <- tryCatch(
          {
            calcVariableImportance()
          },
          error = function(w) {
            NULL
          }
        )
        if (is.null(out)) return(NULL)
        return(list(R = out$R))
      })

      plot.RENDER <- function() {
        res <- plot_data()
        shiny::validate(shiny::need(is_computed(), "Please select target class and run 'Compute'"))
        shiny::req(res)
        R <- res$R
        R <- R[order(-rowSums(R, na.rm = TRUE)), , drop = FALSE]
        R <- pmax(R, 0.05)
        rownames(R) <- playbase::probe2symbol(rownames(R), pgx$genes, "gene_name", fill_na = TRUE)
        par(mfrow = c(1, 1), oma = c(1, 1, 1, 1) * 0.2)
        par(mar = c(8, 4, 1, 0.2), mgp = c(2.5, 0.8, 0))
        barplot(t(R), las = 3, horiz = FALSE,
          cex.names = 0.8, ylab = "cumulative rank")
        klr <- grey.colors(ncol(R))
        legend("topright",
          legend = rev(colnames(R)), fill = rev(klr),
          cex = 0.75, y.intersp = 0.75,
          inset = c(0.03, -0.03), xpd = TRUE)
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
