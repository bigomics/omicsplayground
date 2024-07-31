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
compare_plot_compare1_ui <- function(id,
                                     height,
                                     title,
                                     info.text,
                                     info.methods,
                                     info.references,
                                     info.extra_link,
                                     width) {
  ns <- shiny::NS(id)
  info_text <- "<b>Expression scatter plot.</b> Scatter plots of gene expression."

  PlotModuleUI(ns("plot"),
    title =  title,
    plotlib = "base",
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = width
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
compare_plot_compare1_server <- function(id,
                                         pgx,
                                         input.contrast1,
                                         ## input.contrast2,
                                         hilightgenes,
                                         createPlot,
                                         plottype,
                                         dataset2,
                                         ## compute,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## ------- no need! trigger on input.contrast1 <- reactiveVal <- button
    ## contrast1 <- shiny::reactiveVal(FALSE)
    ## contrast2 <- shiny::reactiveVal(FALSE)
    ## shiny::observeEvent(compute(), {
    ##   contrast1(input.contrast1())
    ##   contrast2(input.contrast2())
    ## })

    call_createPlot <- function( get_data ) {

      ct1 <- input.contrast1()
      ##req(ct1)
      pgx1 <- pgx
      pgx2 <- dataset2()
      all.ct <- names(pgx1$gx.meta$meta)
      ##shiny::req(ct1)
      if(length(ct1)==0) ct1 <- all.ct[1]
      if (!any(ct1 %in% all.ct)) {
        return(NULL)
      }
      ct1 <- intersect(ct1, all.ct)
      higenes <- hilightgenes()
      cex.lab <- 1.0
      ntop <- 9999
      type <- plottype()

      if (length(higenes) <= 3) cex.lab <- 1.3
      data <- createPlot(pgx1, pgx1, pgx2, ct1, type, cex.lab, higenes, ntop, get_data )
      data
    }
    
    plot_data <- shiny::reactive({
      data <- call_createPlot( get_data = TRUE ) 
      return(data)
    })

    scatter1.RENDER <- shiny::reactive({
      dbg("[compare_plot_compare1_server:scatter1.RENDER] reacted! ")
      ## shiny::validate(shiny::need(input.contrast1(), "Please select contrasts and run 'Compute'"))
      call_createPlot( get_data = FALSE )      
    })

    PlotModuleServer(
      "plot",
      plotlib = "base",
      func = scatter1.RENDER,
      csvFunc = plot_data,
      res = c(90, 110), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
