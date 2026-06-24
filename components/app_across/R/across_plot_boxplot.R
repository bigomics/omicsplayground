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
        label = "Show outliers",
        value = TRUE
      ),
      "Show only outlier points on the boxplot (instead of all sample points).",
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
        "No data available. Please select features and click 'Query'."
      ))
      return(df)
    })

    boxplot.RENDER <- function() {
      df <- plot_data()
      show_points <- isTRUE(input$show_points)

      has_split <- isTRUE(attr(df, "has_split"))
      color_by <- attr(df, "color_by")
      if (is.null(color_by)) color_by <- "dataset"
      legend_title <- tools::toTitleCase(color_by)

      genes <- unique(df$gene)
      n_genes <- length(genes)

      ## One box per dataset. When a split is active each dataset shows
      ## side-by-side sub-boxes (one per phenotype value), coloured by value.
      make_box <- function(gene_df, show_legend) {
        plotly::plot_ly(
          data = gene_df,
          x = ~dataset,
          y = ~count,
          color = ~color_group,
          type = "box",
          colors = "Set2",
          boxpoints = if (show_points) "outliers" else FALSE,
          jitter = 0.3,
          pointpos = 0,
          marker = list(opacity = 0.5, size = 4),
          showlegend = show_legend,
          hovertemplate = paste(
            "<b>%{x}</b><br>",
            "Value: %{y:.2f}<br>",
            "<extra></extra>"
          )
        )
      }

      if (n_genes == 1) {
        p <- make_box(df, has_split) %>%
          plotly::layout(
            title = list(text = genes[1], font = list(size = 14)),
            xaxis = list(title = "Dataset", tickangle = 45),
            yaxis = list(title = "Value"),
            boxmode = "group",
            showlegend = has_split,
            legend = list(title = list(text = legend_title), orientation = "v", x = 1.02, y = 0.5),
            margin = list(l = 60, r = if (has_split) 120 else 20, t = 40, b = 100)
          )
      } else {
        plots <- lapply(seq_along(genes), function(i) {
          gene <- genes[i]
          gene_df <- df[df$gene == gene, ]
          make_box(gene_df, has_split && i == 1) %>%
            plotly::layout(
              annotations = list(
                list(
                  text = gene,
                  xref = "paper", yref = "paper",
                  x = 0.5, y = 1.0,
                  showarrow = FALSE,
                  font = list(size = 12, weight = "bold")
                )
              ),
              xaxis = list(title = "", tickangle = 45),
              boxmode = "group"
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
            yaxis = list(title = "Value"),
            boxmode = "group",
            showlegend = has_split,
            legend = list(title = list(text = legend_title), orientation = "v", x = 1.02, y = 0.5),
            margin = list(l = 60, r = 120, t = 30, b = 80)
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

