##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' @export
across_plot_boxplot_ui <- function(id,
                                   title,
                                   info.text,
                                   info.methods = NULL,
                                   info.extra_link = NULL,
                                   width,
                                   height) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        inputId = ns("show_points"),
        label = "Show points",
        value = TRUE
      ),
      "Show individual data points on boxplot.",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "b",
    options = options,
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width
  )
}

#' @export
across_plot_boxplot_server <- function(id,
                                       getPlotData,
                                       watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- reactive({
      df <- getPlotData()
      shiny::validate(shiny::need(
        !is.null(df) && nrow(df) > 0,
        "No data available. Please select genes and click 'Query'."
      ))
      return(df)
    })

    boxplot.RENDER <- function() {
      df <- plot_data()
      show_points <- isTRUE(input$show_points)

      color_by <- attr(df, "color_by")
      if (is.null(color_by)) color_by <- "dataset"
      legend_title <- tools::toTitleCase(color_by)

      genes <- unique(df$gene)
      n_genes <- length(genes)

      if (n_genes == 1) {
        p <- plotly::plot_ly(
          data = df,
          x = ~reorder(color_group, count, FUN = median),
          y = ~count,
          color = ~color_group,
          type = "box",
          colors = "Set2",
          boxpoints = if (show_points) "all" else FALSE,
          jitter = 0.3,
          pointpos = 0,
          marker = list(opacity = 0.5, size = 4),
          hovertemplate = paste(
            "<b>%{x}</b><br>",
            "Value: %{y:.2f}<br>",
            "<extra></extra>"
          )
        ) %>%
          plotly::layout(
            title = list(text = genes[1], font = list(size = 14)),
            xaxis = list(
              title = legend_title,
              tickangle = 45
            ),
            yaxis = list(title = "Expression"),
            showlegend = FALSE,
            margin = list(l = 60, r = 20, t = 40, b = 100)
          )
      } else {
        p <- plotly::plot_ly(
          data = df,
          x = ~gene,
          y = ~count,
          color = ~color_group,
          type = "box",
          colors = "Set2",
          boxpoints = if (show_points) "all" else FALSE,
          jitter = 0.3,
          pointpos = 0,
          marker = list(opacity = 0.5, size = 3),
          hovertemplate = paste(
            "<b>Gene:</b> %{x}<br>",
            "<b>Value:</b> %{y:.2f}<br>",
            "<extra></extra>"
          )
        ) %>%
          plotly::layout(
            xaxis = list(
              title = "Gene",
              tickangle = if (n_genes > 5) 45 else 0
            ),
            yaxis = list(title = "Expression"),
            boxmode = "group",
            legend = list(
              title = list(text = legend_title),
              orientation = "v",
              x = 1.02,
              y = 0.5
            ),
            margin = list(l = 60, r = 120, t = 40, b = 80)
          )
      }

      p <- p %>%
        plotly::config(
          modeBarButtonsToRemove = setdiff(all.plotly.buttons, "toImage"),
          displaylogo = FALSE
        )

      return(p)
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = boxplot.RENDER,
      csvFunc = plot_data,
      res = c(85, 100),
      pdf.width = 10, pdf.height = 8,
      add.watermark = watermark
    )
  })
}

