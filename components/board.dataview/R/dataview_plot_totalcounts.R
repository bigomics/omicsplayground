##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_totalcounts_ui <- function(
    id,
    label = "",
    height,
    width,
    info.text,
    caption,
    title) {
  ns <- shiny::NS(id)

  options <- shiny::tagList( ## AZ
    shiny::radioButtons(
      inputId = ns("sampleqc_plottype"),
      label = "Plot type",
      choices = c(
        "Total abundance",
        "Number of detected features"
      )
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = options,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

dataview_plot_totalcounts_server <- function(id,
                                             getCountStatistics,
                                             r.data_type,
                                             r.data_groupby = reactive(""),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      data_groupby <- r.data_groupby()
      data_type <- r.data_type()

      tbl <- getCountStatistics()
      req(tbl)

      type <- tspan("counts", js = FALSE)

      logtype <- if (data_type == "log2") {
        " (log2)"
      } else if (data_type == "logCPM") {
        " (logCPM)"
      } else {
        ("")
      }

      sampleqc_plottype <- input$sampleqc_plottype

      ## ylab <- paste0("total ", type, logtype)
      ## if (data_groupby != "<ungrouped>") {
      ##   ylab <- paste0("average total ", type, logtype)
      ## }

      if (sampleqc_plottype == "Total abundance") {
        ylab <- paste0("Total ", type)
      } else if (sampleqc_plottype == "Number of detected features") {
        ylab <- "N. of detected features"
      }

      res <- list(
        df = data.frame(
          sample = names(tbl$total.counts),
          ## counts = log10(tbl$total.counts),
          counts = tbl$total.counts,
          ndetectedfeat = tbl$n.detected.features
        ),
        ylab = ylab,
        sampleqc_plottype = sampleqc_plottype
      )

      return(res)
    })


    ## plot.RENDER <- function() {
    ##  res <- plot_data()
    ##  shiny::req(res)
    ##  df <- res[[1]]
    ## ---- xlab ------ ###
    ##  names.arg <- df$sample
    ##  if (length(names.arg) > 20) { names.arg <- "" }
    ##  cex.names <- ifelse(length(names.arg) > 10, 0.8, 0.9)
    ##  par(mar = c(8, 4, 2, 0.5), mgp = c(2.2, 0.8, 0))
    ## if (res$sampleqc_plottype == "Average total abundance") { ## AZ
    ##    barplot(
    ##        df$counts / 1e6,
    ##        las = 3,
    ##        border = NA,
    ##        col = rgb(0.2, 0.5, 0.8, 0.8),
    ##        cex.names = cex.names,
    ##        cex.lab = 1,
    ##        ylab = paste(res$ylab, "(M)"),
    ##        ylim = c(0, max(df$counts) / 1e6) * 1.1,
    ##        names.arg = names.arg
    ##    )
    ## } else if (res$sampleqc_plottype == "Number of detected features") {
    ##    barplot(
    ##        df$ndetectedfeat,
    ##        las = 3,
    ##        border = NA,
    ##        col = rgb(0.2, 0.5, 0.8, 0.8),
    ##        cex.names = cex.names,
    ##        cex.lab = 1,
    ##        ylab = paste(res$ylab, "(M)"),
    ##        ylim = c(0, max(df$ndetectedfeat) / 1e6) * 1.1,
    ##        names.arg = names.arg
    ##     )
    ## }
    ## }

    ## modal_plot.RENDER <- function() { plot.RENDER() }

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      df <- res[[1]]

      if (res$sampleqc_plottype == "Total abundance") {
        fig <-
          plotly::plot_ly(
            data = df, x = ~sample, y = ~counts, type = "bar",
            marker = list(color = omics_colors("brand_blue")),
            hovertemplate = ~ paste0(
              "Sample: <b>", sample, "</b><br>",
              res$ylab, ": <b>", sprintf("%8.0f", counts), "</b>",
              "<extra></extra>"
            )
          ) %>%
          plotly_default() %>%
          plotly::layout(
            xaxis = list(title = FALSE),
            yaxis = list(title = res$ylab),
            margin = list(l = 30, r = 0, t = 0, b = 0)
          )
        fig
      } else if (res$sampleqc_plottype == "Number of detected features") {
        fig <-
          plotly::plot_ly(
            data = df, x = ~sample, y = ~ndetectedfeat, type = "bar",
            marker = list(color = omics_colors("brand_blue")),
            hovertemplate = ~ paste0(
              "Sample: <b>", sample, "</b><br>",
              res$ylab, ": <b>", sprintf("%8.0f", ndetectedfeat), "</b>",
              "<extra></extra>"
            )
          ) %>%
          plotly_default() %>%
          plotly::layout(
            xaxis = list(title = FALSE),
            yaxis = list(title = res$ylab),
            margin = list(l = 30, r = 0, t = 0, b = 0)
          )
        fig
      }
    }

    modal_plotly.RENDER <- function() {
      plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(
          margin = list(l = 35, r = 0, t = 0, b = 0)
        )
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170) * 1, ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
