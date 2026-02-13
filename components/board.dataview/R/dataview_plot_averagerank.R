##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_averagerank_ui <- function(id,
                                         label = "",
                                         height,
                                         width,
                                         title,
                                         info.text,
                                         caption) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_all_isoforms"),
        "Show all available feature isoforms",
        FALSE
      ),
      "Show all available isoforms (e.g., transcripts, PTMs) for each feature",
      placement = "top"
    )
  )

  PlotModuleUI(
    ns("pltsrv"),
    title = title,
    label = label,
    outputFunc = plotly::plotlyOutput,
    outputFunc2 = plotly::plotlyOutput,
    info.text = info.text,
    options = options,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "rank_plot"
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

      ann <- pgx$genes[names(mean.fc), , drop = FALSE]
      sel <- which(names(mean.fc) == gene)


      if (input$show_all_isoforms && labeltype() != "feature") {
        symbol <- playbase::probe2symbol(gene, ann, "symbol", fill_na = TRUE)
        sel2 <- which(ann$symbol == symbol) ## isoforms
        sel <- unique(c(sel2, sel))
      } else {
        sel <- sel
      }

      pd <- list(
        df = data.frame(mean.fc = mean.fc, gene = names(mean.fc)),
        sel = sel,
        gene = gene,
        ylab = ylab
      )
      pd
    })

    plot.RENDER <- function() {
      pd <- plot_data()
      req(pd)

      mean.fc <- pd$df$mean.fc
      names(mean.fc) <- pd$df$gene
      sel <- pd$sel
      gene <- pd$gene
      ylab <- pd$ylab
      xanchor <- "center"

      ## Editor: custom colors
      clr_fill <- if (!is.null(input$color_fill)) input$color_fill else "#b8d4f0"
      clr_line <- if (!is.null(input$color_line)) input$color_line else "#3181de"
      clr_highlight <- if (!is.null(input$color_highlight)) input$color_highlight else "#e3a45a"

      # subsample for speed
      ii <- 1:length(mean.fc)
      if (length(ii) > 200) {
        ii <- c(1:200, seq(201, length(mean.fc), 10))
      }

      fig <- plotly::plot_ly(
        x = ii,
        y = mean.fc[ii],
        type = "scatter",
        mode = "lines",
        fill = "tozeroy",
        fillcolor = clr_fill,
        line = list(width = 0),
        hovertemplate = ~ paste("<extra></extra>")
      )

      i <- 1
      for (i in 1:length(sel)) {
        fig <- fig %>%
          plotly::add_lines(
            x = sel[i],
            y = c(0, mean.fc[sel[i]]),
            type = "scatter",
            mode = "lines",
            line = list(color = clr_highlight, width = 5)
          )
      }

      ## add a second density curve with transparent filling
      ## to add an outline overwriting annotation line
      fig <- fig %>%
        plotly::add_lines(
          x = ii,
          y = mean.fc[ii],
          type = "scatter",
          mode = "lines",
          fill = "tozeroy",
          fillcolor = "#00000000",
          line = list(color = clr_line, width = 2.5)
        )

      i <- 1
      for (i in 1:length(sel)) {
        if (sel[i] < length(mean.fc) / 5) xanchor <- "left"
        if (sel[i] > length(mean.fc) * 4 / 5) xanchor <- "right"
        gene <- names(mean.fc)[sel[i]]
        text <- playbase::probe2symbol(gene, pgx$genes, labeltype(), fill_na = TRUE)
        fig <- fig %>%
          plotly::add_annotations(
            x = sel[i],
            y = mean.fc[sel[i]],
            ax = ifelse(sel[i] < length(mean.fc) / 2, 40, -40),
            ay = -40,
            xanchor = xanchor,
            text = text
          )
      }

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
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
