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
drugconnectivity_plot_enplots_ui <- function(
  id,
  label = "",
  title,
  info.text,
  info.methods,
  info.references,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList()

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    caption = caption,
    plotlib = "base",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    options = plot_opts,
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
drugconnectivity_plot_enplots_server <- function(id,
                                                 pgx,
                                                 dsea_contrast,
                                                 dsea_method,
                                                 dsea_table,
                                                 watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        dsea_contrast <- dsea_contrast()
        dsea_method <- dsea_method()
        shiny::req(pgx$X, dsea_contrast, dsea_method)

        dt <- dsea_table$data()
        ii <- dsea_table$rows_selected()
        jj <- dsea_table$rows_all()
        shiny::req(jj) ## must have non-empty table

        shiny::validate(shiny::need(
          "drugs" %in% names(pgx),
          "missing 'drugs' in object."
        ))
        if (is.null(pgx$drugs)) {
          return(NULL)
        }

        if (is.null(dsea_contrast)) {
          return(NULL)
        }
        
        ## filter with table selection/search
        if (length(ii) > 0) {
          dt <- dt[ii, , drop = FALSE]
        }
        if (length(ii) == 0 && length(jj) > 0) {
          dt <- dt[jj, , drop = FALSE]
        }

        return(dt)
      })

      ## plot.RENDER <- shiny::reactive({
      plot.RENDER <- function() {
        dt <- plot_data()
        if (nrow(dt) == 0) {
          return(NULL)
        }
        sel.drugs <- rownames(dt)
        
        dcontrast <- dsea_contrast()
        dmethod <- dsea_method()
        nplots <- min(nrow(dt), 16)
        
        playbase::pgx.plotDrugConnectivity(
          pgx,
          contrast = dcontrast,
          db = dmethod,
          drugs = sel.drugs,
          nplots = nplots)

      }
      
      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER,
        csvFunc = plot_data,
        res = c(78, 110),
        pdf.width = 6.5, pdf.height = 12.8,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
