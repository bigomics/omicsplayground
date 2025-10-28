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
    ),
    withTooltip(
      shiny::checkboxInput(ns("qweight"), "q-weighting", FALSE),
      "Apply q-value weighting for NES score values."
    )
  )
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = plot_opts,
    download.fmt = c("png", "pdf", "csv", "svg"),
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
                                             pgx,
                                             getActiveDSEA,
                                             getMOA.target,
                                             getMOA.class,
                                             watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
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

      shiny::observeEvent(pgx$X, {
        choices <- c("drug class", "target gene")
        names_choies <- c("drug class", tspan("target gene", js = FALSE))
        names(choices) <- names_choies
        shiny::updateRadioButtons(session, "dsea_moatype", choices = choices)
      })

      plotTopBarplot <- function(ntop, return_csv = FALSE) {
        res <- plot_data()
        shiny::req(res)

        res$score <- res$NES
        yaxistitle <- "score (NES)"
        if (input$qweight) {
          res$score <- res$NES * (1 - res$padj) * (1 - 1 / res$size**1)
          yaxistitle <- "score (qNES)"
        }
        jj <- unique(c(head(order(-res$score), ntop), tail(order(-res$score), ntop)))
        moa.top <- res$score[jj]
        names(moa.top) <- res$pathway[jj]

        df <- data.frame(
          x = factor(names(moa.top), levels = names(moa.top)),
          y = as.numeric(moa.top)
        )

        if (return_csv) {
          return(df)
        }

        p <- playbase::pgx.barplot.PLOTLY(
          data = df,
          x = "x",
          y = "y",
          yaxistitle = yaxistitle,
          xaxistitle = "",
          grouped = FALSE, ## not really...
          yrange = c(-1.1, 1.1) * max(abs(as.numeric(moa.top))),
          xlen = 30
        )

        return(p)
      }

      plot.RENDER <- function() {
        plotTopBarplot(14)
      }

      plot.RENDER2 <- function() {
        plotTopBarplot(24) %>%
          plotly::layout(
            font = list(size = 18)
          )
      }

      plot_data_csv <- function() {
        df <- plotTopBarplot(14, return_csv = TRUE)
        return(df)
      }

      PlotModuleServer(
        "plot",
        plotlib = "plotly",
        func = plot.RENDER,
        func2 = plot.RENDER2,
        csvFunc = plot_data_csv,
        res = c(70, 110),
        pdf.width = 9, pdf.height = 6,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
