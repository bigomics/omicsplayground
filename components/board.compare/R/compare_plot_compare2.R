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
compare_plot_compare2_ui <- function(id,
                                     height,
                                     title,
                                     info.text,
                                     info.methods,
                                     info.extra_link,
                                     width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "base",
    label = "b",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "featuremap",
    color_selection = TRUE
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
compare_plot_compare2_server <- function(id,
                                         pgx,
                                         contrast2,
                                         ## contrast1,
                                         hilightgenes,
                                         createPlot,
                                         plottype,
                                         dataset2,
                                         getMatrices,
                                         watermark = FALSE,
                                         labeltype) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    call_createPlot <- function(get_data) {
      ct <- contrast2()
      shiny::req(ct)
      pgx1 <- pgx
      pgx2 <- dataset2()
      mat <- getMatrices()
      all.ct <- names(pgx2$gx.meta$meta)
      ## if (length(ct) == 0) ct <- all.ct[2]
      if (!any(ct %in% all.ct)) {
        return(NULL)
      }
      ct <- intersect(ct, all.ct)
      higenes <- hilightgenes()
      cex.lab <- 1.0
      ntop <- 9999
      type <- plottype()
      target_col <- mat$target_col

      ## Editor: custom gradient colors
      color_low <- if (!is.null(input$color_low)) input$color_low else "#3181de"
      color_high <- if (!is.null(input$color_high)) input$color_high else "#f23451"
      custom_col <- c(color_low, "#f8f8f8", color_high)

      ## Editor: custom labels
      if (isTRUE(input$custom_labels) && !is.null(input$label_features) && input$label_features != "") {
        custom_genes <- strsplit(input$label_features, "\\s+")[[1]]
        higenes <- custom_genes
      }

      ## Editor: color just selected
      color_sel <- isTRUE(input$color_selection)

      if (length(higenes) <= 3) cex.lab <- 1.3
      data <- createPlot(
        pgx2, pgx1, pgx2, ct, target_col, type, cex.lab,
        higenes, ntop, get_data, labeltype, col = custom_col, color_selection = color_sel
      )
      data
    }

    plot_data <- shiny::reactive({
      data <- call_createPlot(get_data = TRUE)
      return(data)
    })

    plot.RENDER <- shiny::reactive({
      call_createPlot(get_data = FALSE)
    })


    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = plot.RENDER,
      csvFunc = plot_data,
      res = c(90, 110), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
