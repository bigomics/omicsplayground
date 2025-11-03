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
#' @param width
#'
#' @export
expression_plot_fc_fc_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = "FC-FC comparison",
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
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
expression_plot_fc_fc_server <- function(id,
                                         pgx,
                                         comp, # input$gx_contrast
                                         labeltype = reactive("symbol"),
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(comp())
      comp <- comp()
      mx <- pgx$gx.meta$meta[[comp]]
      mx.features <- rownames(mx)
      mx.symbols <- pgx$genes[mx.features, "symbol"]

      if (labeltype() == "symbol") {
        names <- mx.features
        label.names <- mx.symbols
      } else {
        names <- mx.symbols
        label.names <- mx.features
      }

      pd <- list(
        mx = mx,
        comp = comp,
        names = names,
        label.names = label.names
      )

      return(pd)
    })

    plotly_plot <- function() {
      pd <- plot_data()
      shiny::req(pd)
      mx <- pd$mx
      x <- mx[, "fc", drop = FALSE]
      scatter_plot <- plotly::plot_ly(
        x = x$fc[, 1],
        y = x$fc[, 2],
        text = pd$names,
        type = "scatter",
        mode = "markers",
        marker = list(
          size = 4,
          color = omics_colors("brand_blue"),
          line = list(width = 0.5, color = "rgb(0,0,0)")
        ),
        hovertemplate = paste0(
          "<b>%{text}</b><br>",
          "baseline(log2FC): %{x:.2f}<br>",
          "custom (log2FC): %{y:.2f}<br>",
          "<extra></extra>"
        )
      ) %>%
        plotly::layout(
          xaxis = list(title = "baseline(log2FC)"),
          yaxis = list(title = "custom (log2FC)"),
          showlegend = FALSE
        ) %>%
        plotly::toWebGL()
    }

    plot_data_csv <- function() {
      pd <- plot_data()
      df <- pd$mx
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly_plot,
      csvFunc = plot_data_csv,
      res = c(80, 90), # resolution of plots
      pdf.width = 12,
      pdf.height = 5,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
