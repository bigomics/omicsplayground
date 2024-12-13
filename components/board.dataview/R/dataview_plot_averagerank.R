##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_averagerank_ui <- function(
    id,
    label = "",
    height,
    width,
    title,
    info.text,
    caption) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltsrv"),
    title = title,
    label = label,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info.text,
    options = NULL,
    caption = caption,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

dataview_plot_averagerank_server <- function(id,
                                             pgx,
                                             r.gene = reactive(""),
                                             r.samples = reactive(""),
                                             r.data_type = reactive("counts"),
                                             labeltype,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$Y)
      shiny::req(r.gene())

      ## dereference reactives
      gene <- r.gene()
      samples <- r.samples()
      data_type <- r.data_type()

      if (!all(samples %in% colnames(pgx$X))) {
        return(NULL)
      }
      if (!gene %in% rownames(pgx$X)) {
        return(NULL)
      }

      nsamples <- length(samples)

      if (data_type %in% c("counts", "abundance")) {
        mean.fc <- sort(rowMeans(pgx$counts[, samples, drop = FALSE], na.rm = TRUE),
          decreasing = TRUE
        )
        ylab <- tspan("average counts", js = FALSE)
      }
      if (data_type %in% c("logCPM", "log2")) {
        mean.fc <- sort(rowMeans(pgx$X[, samples, drop = FALSE], na.rm = TRUE), decreasing = TRUE)
        ylab <- tspan("average counts (log2)", js = FALSE)
      }

      sel <- which(sub(".*:", "", names(mean.fc)) == sub(".*:", "",gene))
      #sel <- which(names(mean.fc) == gene)      

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

      mean.fc <- pd$df$mean.fc
      sel <- pd$sel
      gene <- pd$gene
      ylab <- pd$ylab
      xanchor <- "center"
      if (sel < length(mean.fc) / 5) xanchor <- "left"
      if (sel > length(mean.fc) * 4 / 5) xanchor <- "right"

      ## subsample for speed
      ii <- 1:length(mean.fc)
      if (length(ii) > 200) {
        ii <- c(1:200, seq(201, length(mean.fc), 10))
      }

      # Translate gene to labeltype
      gene <- playbase::probe2symbol(gene, pgx$genes, labeltype(), fill_na = TRUE)

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
          ax = ifelse(sel < length(mean.fc) / 2, 40, -40),
          ay = -40,
          xanchor = xanchor,
          text = gene
        )

      fig <- fig %>%
        plotly::layout(
          showlegend = FALSE,
          xaxis = list(title = tspan("Ordered genes", js = FALSE)),
          yaxis = list(title = stringr::str_to_sentence(ylab))
        ) %>%
        plotly_default()
      fig
    }

    modal_plot.RENDER <- function() {
      fig <- plot.RENDER() %>%
        plotly_modal_default()
      #
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
