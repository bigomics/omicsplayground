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
compare_plot_fc_correlation_ui <- function(id,
                                           height,
                                           width) {
  ns <- shiny::NS(id)
  info_text <- "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."

  PlotModuleUI(ns("plot"),
    title = "FC Correlation",
    plotlib = "plotly",
    label = "a",
    info.text = info_text,
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
compare_plot_fc_correlation_server <- function(id,
                                               cum_fc,
                                               hilightgenes,
                                               input.contrast1,
                                               input.contrast2,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      # Require inputs
      shiny::req(cum_fc())
      out_data <- cum_fc()
      return(out_data)
    })

    plot_interactive_comp_fc <- function(plot_data, hilight = NULL, marker_size = 6, label_size = 6, cex.axis = 12) {
      shiny::req(plot_data())
      pos <- plot_data()

      ncol_pos <- ncol(pos)
      nrow_pos <- nrow(pos)
      sample_size <- floor(30000 / ncol_pos)
      sample_size <- ifelse(sample_size > nrow_pos, nrow_pos, sample_size)
      sample_size <- sample(rownames(pos), sample_size)
      sample_size <- c(hilight, sample_size)
      sample_size <- unique(sample_size)
      # Get data ready
      data_1 <- pos[, grep("1:", colnames(pos)), drop = FALSE]
      data_2 <- pos[, grep("2:", colnames(pos)), drop = FALSE]
      ncol_d1 <- ncol(data_1)
      ncol_d2 <- ncol(data_2)
      nplots <- ncol_d1 * ncol_d2

      # Prepare collection list
      sub_plots <- vector("list", nplots)
      counter <- 1

      # Iterate over the cols of both data sets
      for (i in seq_len(ncol_d1)) {
        for (j in seq_len(ncol_d2)) {
          # Get the data for the current plot

          pos_i <- cbind(data_1[sample_size, i, drop = FALSE], data_2[sample_size, j, drop = FALSE])
          xlab <- ifelse(i == ncol_d2, colnames(pos_i)[1], "")
          ylab <- ifelse(j == 1, colnames(pos_i)[2], "")
          text_i <- glue::glue("{rownames(pos_i)}<br> x: {pos_i[, 1]}<br> y: {pos_i[, 2]}")
          # Plot the points
          plot_i <- plotly::plot_ly(
            x = ~ pos_i[, 1],
            y = ~ pos_i[, 2],
            text = ~ rownames(pos_i),
            type = "scattergl",
            mode = "markers",
            marker = list(
              size = marker_size,
              opacity = 0.33,
              color = "#22222255",
              line = list(
                color = "#AAAAAA44",
                width = 0.2
              )
            ),
            showlegend = FALSE
          )

          # Add the text to hilighted points
          if (length(hilight) > 1) {
            hilight1 <- intersect(rownames(pos_i), hilight)
            plot_i <- plot_i %>%
              plotly::add_trace(
                x = ~ pos_i[hilight1, 1],
                y = ~ pos_i[hilight1, 2] + pos_i[hilight1, 2] * .05,
                text = ~hilight1,
                key = ~hilight1,
                type = "scattergl",
                mode = "marker+text",
                marker = list(opacity = 1, size = marker_size, color = "red"),
                textposition = "top center",
                textfont = list(color = "#464545"),
                showlegend = FALSE
              )
          }

          # Add the plot to the collection list
          suppressMessages(
          sub_plots[[counter]] <- plot_i %>%
                  plotly::layout(
                    xaxis = list(
                      title = xlab,
                      titlefont = list(size = cex.axis)
                    ),
                    yaxis = list(
                      title = ylab,
                      titlefont = list(size = cex.axis)
                    )
                  )
            )
          counter <- counter +1

        }
      }

      # Combine all plots and set plotly configs
      suppressMessages(
        all_plts <- plotly::subplot(sub_plots,
          nrows = ncol_d1,
          titleX = TRUE, titleY = TRUE,
          shareX = TRUE, shareY = TRUE
        ) %>%
          plotly::config(modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage")) %>%
          plotly::config(toImageButtonOptions = list(
            format = "svg",
            height = 800, width = 800, scale = 1.1
          )) %>%
          plotly::config(displaylogo = FALSE) %>%
          plotly::event_register("plotly_selected") %>%
          plotly::toWebGL()
      )

      return(all_plts)
    }

    fcfcplot.RENDER <- function() {
      higenes <- hilightgenes()
      p <- plot_interactive_comp_fc(plot_data = plot_data, marker_size = 6, cex.axis = 12, hilight = higenes) %>%
        plotly::layout(
          dragmode = "select",
          margin = list(l = 5, r = 5, b = 5, t = 20)
        )
      return(p)
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = fcfcplot.RENDER,
      csvFunc = plot_data,
      res = c(85, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
