##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Decisiontree plot UI input function
#' @description A shiny Module for plotting (UI code).
#' @param id
#' @param label
#' @param height
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
  width
) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::selectInput(ns("plottype"), "Plot type:",
      choices = c("simple", "fancy", "extended"),
      selected = "fancy"
    )
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
    download.fmt = c("png", "pdf", "svg"),
    width = width,
    height = height
  )
}

#' Decisiontree plot Server function
#' @description A shiny Module for plotting (server code).
#' @param id
#' @return
#' @export
biomarker_plot_decisiontree_server <- function(id,
                                               calcVariableImportance,
                                               pgx,
                                               is_computed,
                                               watermark = FALSE) {

  moduleServer(

    id, function(input, output, session) {

      plot_data <- shiny::reactive({
        res <- calcVariableImportance()
        shiny::req(res, is_computed())
        return(res)
      })

      plot.RENDER <- function() {
        imp <- plot_data()
        shiny::req(imp)
        imp$rf$frame$var <- playbase::probe2symbol(imp$rf$frame$var, pgx$genes, "gene_name", fill_na = TRUE)
        playbase::plotDecisionTreeFromImportance(imp = NULL, rf = imp$rf, type = input$plottype)
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
