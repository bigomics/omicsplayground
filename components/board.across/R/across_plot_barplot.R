##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' @export
across_plot_barplot_ui <- function(id,
                                   title,
                                   info.text,
                                   info.methods = NULL,
                                   info.extra_link = NULL,
                                   width,
                                   height) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    title = title,
    plotlib = "plotly",
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width
  )
}

#' @export
across_plot_barplot_server <- function(id,
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

    barplot.RENDER <- function() {
      df <- plot_data()

      color_by <- attr(df, "color_by")
      if (is.null(color_by)) color_by <- "dataset"
      legend_title <- tools::toTitleCase(color_by)

      genes <- unique(df$gene)
      n_genes <- length(genes)

      if (n_genes == 1) {
        p <- plotly::plot_ly(
          data = df,
          x = ~reorder(sample, count),
          y = ~count,
          color = ~color_group,
          type = "bar",
          colors = "Set2",
          hovertemplate = paste(
            "<b>%{x}</b><br>",
            "Count: %{y:.2f}<br>",
            "<extra></extra>"
          )
        ) %>%
          plotly::layout(
            title = list(text = genes[1], font = list(size = 14)),
            xaxis = list(
              title = "Sample",
              showticklabels = FALSE
            ),
            yaxis = list(title = "Expression"),
            legend = list(orientation = "h", y = -0.15, title = list(text = legend_title)),
            margin = list(l = 60, r = 20, t = 40, b = 40)
          )
      } else {
        plots <- lapply(genes, function(gene) {
          gene_df <- df[df$gene == gene, ]
          plotly::plot_ly(
            data = gene_df,
            x = ~reorder(sample, count),
            y = ~count,
            color = ~color_group,
            type = "bar",
            colors = "Set2",
            showlegend = (gene == genes[1]),
            hovertemplate = paste(
              "<b>%{x}</b><br>",
              "Count: %{y:.2f}<br>",
              "<extra></extra>"
            )
          ) %>%
            plotly::layout(
              annotations = list(
                list(
                  text = gene,
                  xref = "paper", yref = "paper",
                  x = 0.5, y = 1.05,
                  showarrow = FALSE,
                  font = list(size = 12, weight = "bold")
                )
              ),
              xaxis = list(showticklabels = FALSE)
            )
        })

        p <- plotly::subplot(
          plots,
          nrows = n_genes,
          shareX = FALSE,
          shareY = FALSE,
          titleY = TRUE
        ) %>%
          plotly::layout(
            legend = list(orientation = "h", y = -0.05, title = list(text = legend_title)),
            margin = list(l = 60, r = 20, t = 20, b = 40)
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
      func = barplot.RENDER,
      csvFunc = plot_data,
      res = c(85, 100),
      pdf.width = 12, pdf.height = 6,
      add.watermark = watermark
    )
  })
}

