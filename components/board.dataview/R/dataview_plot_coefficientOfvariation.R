##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

dataview_plot_variationcoefficient_ui <- function(
  id,
  height,
  width,
  label = "",
  title,
  info.text,
  info.methods,
  info.extra_link,
  caption
) {
  ns <- shiny::NS(id)

  menu_grouped <- "<code>grouped</code>"

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}

dataview_plot_variationcoefficient_server <- function(id,
                                                      pgx,
                                                      r.samples,
                                                      r.groupby,
                                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X, pgx$samples)
      counts <- pgx$counts
      Y <- pgx$samples
      samples <- r.samples()
      groupby <- r.groupby()
      if (!all(samples %in% rownames(Y))) {
        return(NULL)
      }

      if (groupby == "<ungrouped>") {
        res <- as.matrix(playbase::compute_CV(counts))
        colnames(res) <- "Samples"
      } else {
        ph <- as.character(Y[samples, groupby])
        ph.groups <- split(seq_along(ph), ph)
        LL <- lapply(ph.groups, function(sel) {
          if (length(sel) <= 1) {
            return(NULL)
          }
          playbase::compute_CV(counts[, sel, drop = FALSE])
        })
        LL <- LL[!sapply(LL, is.null)]
        res <- do.call(cbind, LL)
        rm(LL)
      }

      rm(counts, samples)
      return(res)
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      boxplot(res, ylab = "CV (%)", las = 1, outcex = 0.6)
      grid()
    }

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      long.df <- reshape2::melt(res)
      colnames(long.df) <- c("gene", "sample", "value")
      fig <- playbase::pgx.boxplot.PLOTLY(
        data = long.df,
        x = "sample",
        y = "value",
        yaxistitle = "CV (%)"
      ) %>%
        plotly_default()
      fig
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    modal_plotly.RENDER <- function() {
      plotly.RENDER() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data,
      res = c(90, 170),
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark
    )
  })
}
