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
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "expression_barplot"
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

    output$rank_list <- shiny::renderUI({
      res <- plot_data()
      shiny::req(res)
      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = session$ns("rank_list_basic"),
          text = NULL,
          labels = colnames(res)
        )
      )
    })

    plotly.RENDER <- function() {
      res <- plot_data()
      shiny::req(res)
      long.df <- reshape2::melt(res)
      colnames(long.df) <- c("gene", "sample", "value")
      long.df$sample <- as.character(long.df$sample)

      bar_color <- if (is.null(input$scatter_color)) get_color_theme()$secondary else input$scatter_color
      fill_color <- adjustcolor(bar_color, alpha.f = 0.35)
      bars_order <- input$bars_order
      samples <- colnames(res)

      ## Apply sample ordering
      if (!is.null(bars_order)) {
        if (bars_order == "ascending") {
          medians <- tapply(long.df$value, long.df$sample, median, na.rm = TRUE)
          samples <- names(sort(medians))
        } else if (bars_order == "descending") {
          medians <- tapply(long.df$value, long.df$sample, median, na.rm = TRUE)
          samples <- names(sort(medians, decreasing = TRUE))
        } else if (bars_order == "custom" && !is.null(input$rank_list_basic) &&
          all(input$rank_list_basic %in% colnames(res))) {
          samples <- input$rank_list_basic
        }
      }
      long.df$sample <- factor(long.df$sample, levels = samples)

      fig <- playbase::pgx.boxplot.PLOTLY(
        data = long.df,
        x = "sample",
        y = "value",
        yaxistitle = "CV (%)",
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
      csvFunc = plot_data,
      res = c(90, 170),
      pdf.width = 6,
      pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
