##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_plot_boxplot_ui <- function(
  id,
  label = "",
  height,
  title,
  caption,
  info.text
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    plotlib = "plotly",
    label = label,
    caption = caption,
    info.text = info.text,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "expression_barplot"
  )
}

dataview_plot_boxplot_server <- function(id,
                                         parent.input,
                                         getCountsTable,
                                         r.samples = reactive(""),
                                         r.data_type,
                                         watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      res <- getCountsTable()
      samples <- r.samples()
      shiny::req(res)
      list(counts = res$log2counts, sample = colnames(res$log2counts))
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

    output$rank_list <- shiny::renderUI({
      res <- plot_data()
      shiny::req(res)
      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = session$ns("rank_list_basic"),
          text = NULL,
          labels = as.character(res$sample)
        )
      )
    })

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      data_type <- r.data_type()

      if (data_type == "counts") {
        ylab <- "Counts (log2CPM)"
      } else if (data_type == "abundance") {
        ylab <- "Abundance"
      } else {
        ylab <- "Abundance (log2)"
      }

      df <- res$counts[, , drop = FALSE]
      if (nrow(df) > 1000) {
        sel <- sample(nrow(df), 1000)
        df <- df[sel, , drop = FALSE]
      }
      long.df <- reshape2::melt(df)
      colnames(long.df) <- c("gene", "sample", "value")
      long.df$sample <- as.character(long.df$sample)

      bar_color <- if (is.null(input$scatter_color)) get_color_theme()$secondary else input$scatter_color
      fill_color <- adjustcolor(bar_color, alpha.f = 0.35)
      bars_order <- input$bars_order
      samples <- res$sample

      ## Apply sample ordering
      if (!is.null(bars_order)) {
        if (bars_order == "ascending") {
          medians <- tapply(long.df$value, long.df$sample, median, na.rm = TRUE)
          samples <- names(sort(medians))
        } else if (bars_order == "descending") {
          medians <- tapply(long.df$value, long.df$sample, median, na.rm = TRUE)
          samples <- names(sort(medians, decreasing = TRUE))
        } else if (bars_order == "custom" && !is.null(input$rank_list_basic) &&
          all(input$rank_list_basic %in% res$sample)) {
          samples <- input$rank_list_basic
        }
      }
      long.df$sample <- factor(long.df$sample, levels = samples)

      ## boxplot
      fig <- playbase::pgx.boxplot.PLOTLY(
        data = long.df,
        x = "sample",
        y = "value",
        yaxistitle = ylab,
        color = bar_color,
        fillcolor = fill_color,
        linecolor = bar_color
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
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(90, 170), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
