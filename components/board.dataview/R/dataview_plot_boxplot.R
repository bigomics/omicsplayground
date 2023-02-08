##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_plot_boxplot_ui <- function(id, label = "", height = c(600, 800)) {
  ns <- shiny::NS(id)

  menu_grouped <- "<code>grouped</code>"
  info_text <- paste0("Boxplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ", menu_grouped, " setting uder the main <i>Options</i>.")

  PlotModuleUI(
    ns("pltmod"),
    title = "Counts distribution",
    plotlib = "plotly",
    label = label,
    info.text = info_text,
    download.fmt = c("png","pdf","csv"),
    height = height
  )
}

dataview_plot_boxplot_server <- function(id, parent.input, getCountsTable, watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## extract data from pgx object
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      req(res)

      list(
        counts = res$log2counts,
        sample = colnames(res$log2counts)
      )
    })

    plot.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      par(mar = c(8, 4, 1, 2), mgp = c(2.2, 0.8, 0))
      ## ---- xlab ------ ###
      xaxt <- "l"
      names.arg <- res$sample
      if (length(names.arg) > 20) {
        names.arg <- rep("", length(names.arg))
        xaxt <- "n"
      }

      cex.names <- ifelse(length(names.arg) > 10, 0.8, 0.9)
      boxplot(
        res$counts,
        col = rgb(0.2, 0.5, 0.8, 0.4),
        names = names.arg,
        cex.axis = cex.names,
        border = rgb(0.824, 0.824, 0.824, 0.9),
        xaxt = xaxt,
        las = 3,
        cex.lab = 1,
        ylab = "counts (log2)",
        outline = FALSE,
        varwidth = FALSE
      )
    }


    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)

      df <- res$counts[, ]
      long.df <- reshape2::melt(df)
      colnames(long.df) <- c("gene", "sample", "value")

      ## boxplot
      fig <- pgx.boxplot.PLOTLY(
        data = long.df,
        x = "sample",
        y = "value",
        yaxistitle = "Counts (log2)"
        )
      fig
    }

    modal_plot.RENDER <- function() {
      plot.RENDER()
    }

    modal_plotly.RENDER <- function() {
      plotly.RENDER()
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      # renderFunc = shiny::renderPlot,
      # renderFunc2 = shiny::renderPlot,
      res = c(90, 170), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
