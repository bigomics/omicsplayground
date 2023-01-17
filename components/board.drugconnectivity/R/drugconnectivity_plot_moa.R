##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
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
drugconnectivity_plot_moa_ui <- function(id,
                                         label = "",
                                         height = c(600, 800)) {
  ns <- shiny::NS(id)
  info_text <- strwrap("This plot visualizes the <strong>mechanism of
                       action</strong> (MOA) across the enriched drug profiles.
                       On the vertical axis, the GSEA normalized enrichment
                       score of the MOA class or gene target is plotted. You
                       can switch to visualize between MOA class or target
                       gene.")

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("dsea_moatype"), "Plot type:", c("drug class", "target gene"), inline = TRUE),
      "Select plot type of MOA analysis: by class description or by target gene."
    )
  )
  PlotModuleUI(ns("plot"),
               title = "Mechanism of action",
               label = label,
               plotlib = "base",
               info.text = info_text,
               options = plot_opts,
               download.fmt = c("png", "pdf", "csv"),
               width = c("auto", "100%"),
               height = height
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

      plot.RENDER <- shiny::reactive({
        res <- plot_data()
        shiny::req(res)

        ntop <- 16
        jj <- unique(c(head(order(-res$NES), ntop), tail(order(-res$NES), ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        par(mfrow = c(2, 1), mar = c(4, 3.5, 0.1, 0), mgp = c(1.7, 0.65, 0))
        barplot(moa.top,
                horiz = FALSE, las = 3,
                ylab = "enrichment  (NES)",
                cex.names = 0.96
        )
      })

      plot.RENDER2 <- shiny::reactive({
        res <- plot_data()
        shiny::req(res)

        ntop <- 32
        jj <- unique(c(head(order(-res$NES), ntop), tail(order(-res$NES), ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        par(mfrow = c(2, 1), mar = c(4, 3.5, 0.1, 0), mgp = c(1.7, 0.65, 0))
        barplot(moa.top,
                horiz = FALSE, las = 3,
                ylab = "enrichment  (NES)",
                cex.names = 1
        )
      })

      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER2,
        csvFunc = plot_data,
        res = c(70, 110),
        pdf.width = 10, pdf.height = 5,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
