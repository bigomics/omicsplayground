##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Decisiontree plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
biomarker_plot_decisiontree_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::selectInput(ns("plottype"), "Plot type:",
      choices = c("simple","fancy","extended"),
      selected = "fancy")
  )
  
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "base",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    options = options,
    caption = caption,
    download.fmt = c("png", "pdf"),
    width = width,
    height = height
  )
}

#' Decisiontree plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
biomarker_plot_decisiontree_server <- function(id,
                                               calcVariableImportance,
                                               is_computed,
                                               watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        res <- calcVariableImportance()
        shiny::req(res)
        shiny::req(is_computed())
        return(res)
      })

      plot.RENDER <- function() {
        res <- plot_data()
        shiny::req(res)
        
        playbase::plotDecisionTreeFromImportance(
          res, type = input$plottype )
        
        ## par(mfrow = c(1, 1), mar = c(1, 0, 2, 0))
        ## is.surv <- grepl("Surv", res$rf$call)[2]
        ## is.surv
        ## if (is.surv) {
        ##   rf <- partykit::as.party(res$rf)
        ##   partykit::plot.party(rf)
        ## } else {
        ##   ## rpart.plot::rpart.plot(res$rf)
        ##   rf <- partykit::as.party(res$rf)
        ##   is.multinomial <- length(table(res$y)) > 2
        ##   if (is.multinomial) {
        ##     ## plot(rf, type="extended")
        ##     plot(rf, type = "simple")
        ##   } else {
        ##     plot(rf, type = "simple")
        ##   }
        ## }

      }

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        res = c(60, 100),
        pdf.width = 10, pdf.height = 6,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
