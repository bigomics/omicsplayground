##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Drug Connectivity plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
drugconnectivity_plot_cmap_enplot_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
  )
}

#' Drug Connectivity plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
drugconnectivity_plot_cmap_enplot_server <- function(id,
                                                     pgx,
                                                     getActiveDSEA,
                                                     cmap_table,
                                                     watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        dsea <- getActiveDSEA()

        res <- list(
          pgx = pgx,
          dsea = dsea,
          cmap_table = cmap_table
        )
        return(res)
      })

      plotly.RENDER <- shiny::reactive({
        res <- plot_data()
        pgx <- res$pgx
        dsea <- res$dsea
        cmap_table <- res$cmap_table

        ####
        dt <- dsea$table
        ii <- cmap_table$rows_selected()
        if (length(ii) == 0) {
          ii <- cmap_table$rows_all()
        }
        if (length(ii) == 0) {
          return(NULL)
        }

        ## draw enrichment plot
        d <- dt$drug[ii[1]]
        rnk <- dsea$stats
        dtype <- sub("[@_].*$", "", names(rnk))
        gmt <- names(rnk)[dtype == d]
        p1 <- playbase::gsea.enplotly(rnk, gmt, main = d)
        return(p1)
      })

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plotly.RENDER,
        func2 = plotly.RENDER,
        csvFunc = plot_data,
        res = c(80, 105),
        pdf.width = 6, pdf.height = 10,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
