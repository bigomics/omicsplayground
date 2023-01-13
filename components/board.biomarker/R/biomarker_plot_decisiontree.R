##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
biomarker_plot_decisiontree_ui <- function(id,
                                           label = "",
                                           height = c(600, 800)) {
  ns <- shiny::NS(id)
  info_text <- strwrap("The decision tree shows a tree solution
                      for classification based on the top most
                      important features.")

  PlotModuleUI(ns("plot"),
    title = "Decision tree",
    label = label,
    plotlib = "base",
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
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
                                               watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        res <- calcVariableImportance()
        shiny::req(res)

        res <- list(
          res = res
        )
        return(res)
      })

      plot.RENDER <- shiny::reactive({
        res <- plot_data()
        shiny::req(res)

        res <- res$res

        par(mfrow = c(1, 1), mar = c(1, 0, 2, 0))
        is.surv <- grepl("Surv", res$rf$call)[2]
        is.surv
        if (is.surv) {
          rf <- partykit::as.party(res$rf)
          partykit::plot.party(rf)
        } else {
          rpart.plot::rpart.plot(res$rf)
          title("Classification tree", cex = 1.2, line = 3, adj = 0.35)
        }
      })

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        csvFunc = plot_data,
        res = c(72, 315),
        pdf.width = 10, pdf.height = 6,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
