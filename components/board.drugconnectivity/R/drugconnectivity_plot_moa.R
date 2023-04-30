##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Mechanism-of-action plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
drugconnectivity_plot_moa_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height
                                         ) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("dsea_moatype"), "Plot type:", c("drug class", "target gene"), inline = TRUE),
      "Select plot type of MOA analysis: by class description or by target gene."
    )
  )
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv"),
    height = height,
    width = c("auto", "100%"),
  )
}

#' Mechanism of action plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
drugconnectivity_plot_moa_server <- function(id,
                                             getActiveDSEA,
                                             getMOA.target,
                                             getMOA.class,
                                             watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {

      getMOA <- shiny::reactive({
        moatype <- input$dsea_moatype
        if (moatype == "target gene") {
          res <- getMOA.target()
        } else if (moatype == "drug class") {
          res <- getMOA.class()
        } else {
          res <- NULL
        }
        res
      })

      plot_data <- shiny::reactive({
        res <- getMOA()
        shiny::req(res)
        return(res)
      })

      plotTopBarplot <- function(ntop) {
        res <- plot_data()
        shiny::req(res)

        jj <- unique(c(head(order(-res$NES), ntop), tail(order(-res$NES), ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        
        p <- playbase::pgx.barplot.PLOTLY(
          data = data.frame(
            x = factor(names(moa.top), levels = names(moa.top)),
            y = as.numeric(moa.top)
          ),
          x = "x",
          y = "y",
          yaxistitle = "Enrichment (NES)",
          xaxistitle = "",
          grouped = FALSE,  ## not really...
          yrange = c(-1.1, 1.1) * max(abs(as.numeric(moa.top))),
          xlen = 25
        )
        
        return(p)
      }

      plot.RENDER <- function() {
        plotTopBarplot(14)         
      }

      plot.RENDER2 <- function() {
        plotTopBarplot(24) %>%
          plotly::layout(
            font = list( size = 18 )
          )
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot.RENDER,
        func2 = plot.RENDER2,
        csvFunc = plot_data,
        res = c(70, 110),
        pdf.width = 9, pdf.height = 6,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
