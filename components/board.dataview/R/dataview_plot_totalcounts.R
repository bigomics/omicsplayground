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
  title
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::radioButtons(
      inputId = ns("sampleqc_plottype"),
      label = "Plot type",
      choices = c("Total abundance", "Number of detected features")
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
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "expression_barplot"
  )
}

dataview_plot_totalcounts_server <- function(id,
                                             getCountStatistics,
                                             r.data_type,
                                             r.samples = reactive(""),
                                             r.data_groupby = reactive(""),
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      data_groupby <- r.data_groupby()
      data_type <- r.data_type()
      samples <- r.samples()
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

      if (sampleqc_plottype == "Total abundance") {
        ylab <- paste0("Total ", type)
      } else if (sampleqc_plottype == "Number of detected features") {
        ylab <- "N. of detected features"
      }

      res <- list(
        df = data.frame(
          sample = names(tbl$total.counts),
          counts = tbl$total.counts,
          ndetectedfeat = tbl$n.detected.features
        ),
        ylab = ylab,
        sampleqc_plottype = sampleqc_plottype
      )

      return(res)
    })

    output$rank_list <- shiny::renderUI({
      res <- plot_data()
      shiny::req(res)
      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = session$ns("rank_list_basic"),
          text = NULL,
          labels = as.character(res$df$sample)
        )
      )
    })

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      df <- res[[1]]

      bar_color <- if (is.null(input$scatter_color)) get_color_theme()$secondary else input$scatter_color
      bars_order <- input$bars_order

      is_total <- res$sampleqc_plottype == "Total abundance"
      ycol <- if (is_total) "counts" else "ndetectedfeat"

      ## Apply bar ordering
      if (!is.null(bars_order)) {
        if (bars_order == "ascending") {
          df <- df[order(df[[ycol]]), ]
        } else if (bars_order == "descending") {
          df <- df[order(-df[[ycol]]), ]
        } else if (bars_order == "custom" && !is.null(input$rank_list_basic) &&
          all(input$rank_list_basic %in% as.character(df$sample))) {
          df$sample <- factor(df$sample, levels = input$rank_list_basic)
          df <- df[order(df$sample), ]
        }
      }
      df$sample <- factor(df$sample, levels = as.character(df$sample))

      if (is_total) {
        fig <-
          plotly::plot_ly(
            data = df, x = ~sample, y = ~counts, type = "bar",
            marker = list(color = bar_color),
            hovertemplate = ~ paste0(
              "Sample: <b>", sample, "</b><br>",
              res$ylab, ": <b>", sprintf("%8.0f", counts), "</b>",
              "<extra></extra>"
            )
          ) %>%
          plotly_default() %>%
          plotly::layout(
            xaxis = list(title = FALSE),
            yaxis = list(title = list(text = res$ylab, standoff = 25L)),
            margin = list(l = 30, r = 0, t = 0, b = 0)
          )
        fig
      } else {
        fig <-
          plotly::plot_ly(
            data = df, x = ~sample, y = ~ndetectedfeat, type = "bar",
            marker = list(color = bar_color),
            hovertemplate = ~ paste0(
              "Sample: <b>", sample, "</b><br>",
              res$ylab, ": <b>", sprintf("%8.0f", ndetectedfeat), "</b>",
              "<extra></extra>"
            )
          ) %>%
          plotly_default() %>%
          plotly::layout(
            xaxis = list(title = FALSE),
            yaxis = list(title = list(text = res$ylab, standoff = 25L)),
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
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
