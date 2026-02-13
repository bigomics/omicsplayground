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
compare_plot_fcfc_ui <- function(id,
                                 height,
                                 title,
                                 info.text,
                                 info.methods,
                                 info.extra_link,
                                 width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "scatter_highlight"
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
compare_plot_fcfc_server <- function(id,
                                     getMatrices,
                                     hilightgenes,
                                     watermark = FALSE,
                                     labeltype,
                                     pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      # Require inputs
      shiny::req(getMatrices())
      res <- getMatrices()
      F1 <- res$F1
      F2 <- res$F2
      out_data <- cbind(F1, F2)
      return(out_data)
    })

    interactive_fcfc <- function(plot_data, hilight = NULL,
                                 marker_size = 6, label_size = 6, cex.axis = 12,
                                 color_point = "#22222255", color_highlight = "red") {
      shiny::req(plot_data())
      FC <- plot_data()
      mat <- getMatrices()


      ## subsample for speed
      ncol_FC <- ncol(FC)
      nrow_FC <- nrow(FC)
      sample_size <- floor(30000 / ncol_FC)
      sample_size <- ifelse(sample_size > nrow_FC, nrow_FC, sample_size)

      genes <- sample(rownames(FC), sample_size)

      if (any(hilight %in% genes)) {
        genes <- c(intersect(hilight, genes), genes)
      }

      genes <- unique(genes)


      ## Get data ready
      data_1 <- mat$F1
      data_2 <- mat$F2
      ncol_d1 <- ncol(data_1)
      ncol_d2 <- ncol(data_2)
      nplots <- ncol_d1 * ncol_d2

      # Prepare collection list
      sub_plots <- vector("list", nplots)
      counter <- 1

      # Iterate over the cols of both data sets
      for (j in ncol_d2:1) {
        for (i in seq_len(ncol_d1)) {
          ## Get the data for the current plot

          F <- cbind(
            data_1[genes, i, drop = FALSE],
            data_2[genes, j, drop = FALSE]
          )
          xlab <- ifelse(j == 1, colnames(data_1)[i], "")
          ylab <- ifelse(i == 1, colnames(data_2)[j], "")

          ## Plot the points
          plot_i <- plotly::plot_ly(
            x = F[, 1],
            y = F[, 2],
            text = rownames(F),
            type = "scattergl",
            mode = "markers",
            marker = list(
              size = marker_size,
              opacity = 0.33,
              color = color_point,
              line = list(
                color = "#AAAAAA44",
                width = 0.2
              )
            ),
            showlegend = FALSE
          )

          # Add the text to hilighted points
          if (length(hilight) > 0) {
            hilight1 <- intersect(rownames(F), hilight)
            hilighy_label <- playbase::probe2symbol(hilight1, pgx$genes, labeltype(), fill_na = TRUE)
            plot_i <- plot_i %>%
              plotly::add_trace(
                x = F[hilight1, 1],
                y = F[hilight1, 2],
                text = hilighy_label,
                key = hilight1,
                type = "scattergl",
                mode = "marker+text",
                marker = list(opacity = 1, size = marker_size, color = color_highlight),
                textposition = "top center",
                textfont = list(color = "#464545"),
                showlegend = FALSE
              )
          }

          # Add the plot to the collection list
          ## suppressMessages(
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
          ## )
          counter <- counter + 1
        }
      }

      # Combine all plots and set plotly configs
      suppressMessages(
        all_plts <- plotly::subplot(sub_plots,
          nrows = ncol_d2,
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
      shiny::validate(shiny::need(getMatrices(), "Please select contrasts and run 'Compute'"))
      higenes <- hilightgenes()

      ## Editor: custom colors
      clr_point <- if (!is.null(input$color_point)) input$color_point else "#222222"
      clr_highlight <- if (!is.null(input$color_highlight)) input$color_highlight else "#f23451"

      ## Editor: custom labels
      if (isTRUE(input$custom_labels) && !is.null(input$label_features) && input$label_features != "") {
        custom_genes <- strsplit(input$label_features, "\\s+")[[1]]
        higenes <- custom_genes
      }

      p <- interactive_fcfc(
        plot_data = plot_data, marker_size = 6, cex.axis = 12,
        hilight = higenes,
        color_point = clr_point, color_highlight = clr_highlight
      ) %>%
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
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
