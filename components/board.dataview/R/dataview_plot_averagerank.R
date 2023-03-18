##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_averagerank_ui <- function(id, label = "", height = c(350, 600)) {
  ns <- shiny::NS(id)
  info_text <- paste0("Ranking of the selected gene by decreasing average expression.")

  PlotModuleUI(
    ns("pltsrv"),
    title = "Average rank",
    label = label,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
  )
}

dataview_plot_averagerank_server <- function(id,
                                             pgx,
                                             r.gene = reactive(""),
                                             r.samples = reactive(""),
                                             r.data_type = reactive("counts"),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y)
      shiny::req(r.gene())

      ## dereference reactives
      gene <- r.gene()
      samples <- r.samples()
      data_type <- r.data_type()

      nsamples <- length(samples)
      if (data_type == "counts") {
        mean.fc <- sort(rowMeans(pgx$counts[, samples, drop = FALSE]), decreasing = TRUE)
        ylab <- "expression (counts)"
      }
      if (data_type == "logCPM") {
        mean.fc <- sort(rowMeans(pgx$X[, samples, drop = FALSE]), decreasing = TRUE)
        ylab <- "expression (log2CPM)"
      }

      sel <- which(sub(".*:", "", names(mean.fc)) == gene)

      pd <- list(
        df = data.frame(mean.fc = mean.fc),
        sel = sel,
        gene = gene,
        ylab = ylab
      )
      pd
    })

    plot.RENDER.save <- function() {
      pd <- plot_data()
      req(pd)

      mean.fc <- pd$df$mean.fc
      sel <- pd$sel
      gene <- pd$gene
      ylab <- pd$ylab

      par(mar = c(2.3, 3.0, 1, 0), mgp = c(2.0, 0.6, 0))
      base::plot(mean.fc,
        type = "h", lwd = 0.4,
        col = "#bbd4ee", cex.axis = 0.9,
        ylab = ylab, xlab = "ordered genes", xaxt = "n"
      )
      points(sel, mean.fc[sel], type = "h", lwd = 2, col = "black")
      text(sel, mean.fc[sel], gene, pos = 3, cex = 0.9)
    }

    plot.RENDER <- function() {
      pd <- plot_data()
      req(pd)

      mean.fc <- log2(pd$df$mean.fc)
      mean.fc <- pd$df$mean.fc
      sel <- pd$sel
      gene <- pd$gene
      ylab <- pd$ylab

      ## subsample for speed
      ii <- 1:length(mean.fc)
      if (length(ii) > 200) {
        ii <- c(1:200, seq(201, length(mean.fc), 10))
      }

      fig <-
        plotly::plot_ly(
          x = ii,
          y = mean.fc[ii],
          type = "scatter",
          mode = "lines",
          fill = "tozeroy",
          fillcolor = omics_colors("light_blue"),
          line = list(
            width = 0
          ),
          hovertemplate = ~ paste(
            ## NOTE: currently x and y cannot be shown explicitly as the format doesnt comply with the chart
            ## TODO: check if the format can be adjusted (if information is needed in tooltip)
            # "Gene rank: <b>", floor(ii), "</b><br>",
            # "Density: <b>", sprintf("%1.3f", mean.fc[ii]), "</b>",
            "<extra></extra>"
          )
        )

      fig <- fig %>%
        plotly::add_lines(
          x = sel,
          y = c(0, mean.fc[sel]),
          type = "scatter",
          mode = "lines",
          line = list(
            color = omics_colors("orange"),
            width = 5
          )
        ) %>%
        ## add a second density curve with transparent filling
        ## to add an outline overwriting annotation line
        plotly::add_lines(
          x = ii,
          y = mean.fc[ii],
          type = "scatter",
          mode = "lines",
          fill = "tozeroy",
          fillcolor = "#00000000",
          line = list(
            color = omics_colors("brand_blue"),
            width = 2.5
          )
        ) %>%
        plotly::add_annotations(
          x = sel,
          y = mean.fc[sel],
          ax = 20,
          ay = -40,
          text = gene
        )

      fig <- fig %>%
        plotly::layout(
          showlegend = FALSE,
          xaxis = list(title = "Ordered genes"),
          yaxis = list(title = stringr::str_to_sentence(ylab))
        ) %>%
        plotly_default()
      fig
    }

    modal_plot.RENDER <- function() {
      fig <- plot.RENDER() %>%
        plotly_modal_default()
      ## fig <- plotly::style(fig, marker.size = 14)
      fig
    }

    PlotModuleServer(
      "pltsrv",
      plotlib = "plotly",
      plotlib2 = "plotly",
      func = plot.RENDER,
      func2 = modal_plot.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
