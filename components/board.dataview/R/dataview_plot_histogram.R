##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_histogram_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  info_text <- paste0("Histogram of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the <code>grouped</code> setting uder the main <i>Options</i>.")

  PlotModuleUI(
    ns("pltmod"),
    title = "Counts histogram",
    label = label,
    plotlib = "plotly",
    # outputFunc = plotOutput,
    # outputFunc2 = plotOutput,
    info.text = info_text,
    options = NULL,
    download.fmt = c("png", "pdf", "csv"),
    width = c("auto", "100%"),
    height = height
  )
}

dataview_plot_histogram_server <- function(id, getCountsTable, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    gx.hist <- function(gx, n = 1000, main = "", ylim = NULL, plot = TRUE) {
      jj <- 1:nrow(gx)
      if (length(jj) > n) jj <- sample(jj, n, replace = TRUE)
      h0 <- hist(as.vector(c(gx[jj], min(gx), max(gx))),
        breaks = 120,
        plot = plot,
        main = main,
        border = FALSE,
        col = "grey",
        freq = FALSE, ## ylim=ylim,
        xlim = c(min(gx), max(gx)),
        xlab = "expression (log2)",
        cex.lab = 1
      )
      i <- 1
      H <- c()
      for (i in 1:ncol(gx)) {
        h1 <- hist(gx[jj, i], breaks = h0$breaks, plot = FALSE)
        if (plot) lines(h0$mids, h1$density, col = "black", lwd = 0.5)
        H <- cbind(H, h1$density)
      }
      colnames(H) <- colnames(gx)
      data.frame(mids = h0$mids, density = h0$density, H)
    }

    ## extract data from pgx object
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      shiny::req(res)
      hh <- gx.hist(gx = res$log2counts, n = 2000, plot = FALSE)
      pdata <- list(
        histogram = hh,
        log2counts = res$log2counts
      )
      pdata
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      par(mar = c(8, 4, 1, 2), mgp = c(2.2, 0.8, 0))
      gx.hist(gx = res$log2counts, n = 2000, plot = TRUE) # , main="histogram")
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    plotly.RENDER <- function() {
      pdata <- plot_data()
      shiny::req(pdata)

      hist <- pdata$histogram
      smoothen <- function(y) {
        loess(y ~ mid, data.frame(mid = hist$mid, y = y), span = 0.25)$fitted
      }
      y.smooth <- apply(hist[, 3:ncol(hist)], 2, smoothen)

      df <- data.frame(
        x = rep(hist$mids, ncol(hist) - 2),
        ## y = unlist(hist[,3:ncol(hist)]),
        y = as.vector(y.smooth),
        sample = as.vector(mapply(rep, colnames(hist)[-c(1, 2)], nrow(hist)))
      )

      fig <-
        plotly::plot_ly(
          data = df,
          x = ~x,
          y = ~y,
          type = "scatter",
          mode = "lines",
          split = ~sample,
          color = ~sample,
          colors = omics_pal_d(palette = "expanded")(length(unique(df$sample))) # ,
          # hovertemplate = ~paste0(
          #   "Sample: <b>", sample, "</b><br>",
          #   "Expression: <b>", x, "</b><br>",
          #   "Density: <b>", y, "</b>",
          #   "<extra></extra>"
          # )
        ) %>%
        plotly::layout(
          xaxis = list(title = "Expression"),
          yaxis = list(title = "Density"),
          ## TODO: decide if unified label or not - maybe only in zoom mode as it's that long?
          hovermode = "x unified",
          font = list(family = "Lato"),
          margin = list(l = 10, r = 10, b = 10, t = 10),
          showlegend = FALSE
        ) %>%
        plotly_default1()
      fig
    }

    modal_plotly.RENDER <- function() {
      plotly.RENDER() %>%
        plotly::layout(
          showlegend = TRUE,
          font = list(size = 18)
        )
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      # renderFunc = shiny::renderPlot,
      # renderFunc2 = shiny::renderPlot,
      res = c(90, 170) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
