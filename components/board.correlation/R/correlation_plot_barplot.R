##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
correlation_plot_barplot_ui <- function(
  id,
  title,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::selectInput(ns("order_opt"), "Order by:",
        choices = c(
          "both",
          "correlation",
          "partial Correlation"
        ),
        multiple = FALSE,
        selected = "both"
      ),
      "Sort order of groups based on correlation.",
      placement = "top"
    )
  )

  PlotModuleUI(
    id = ns("plot"),
    title = title,
    label = label,
    plotlib = "plotly",
    caption = caption,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "grouped_barplot",
    palette_default = "default"
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
correlation_plot_barplot_server <- function(id,
                                            getPartialCorrelation,
                                            getGeneCorr,
                                            cor_table,
                                            watermark = FALSE,
                                            pgx,
                                            labeltype) {
  moduleServer(id, function(input, output, session) {
    # reactive function listeninng for changes in input
    plot_data <- shiny::reactive({
      df <- getPartialCorrelation()
      R <- getGeneCorr()

      sel <- cor_table$rownames_current()
      shiny::req(sel)

      ## cor_table!=R mismatch!!!
      if (length(sel) > nrow(R)) {
        return(NULL)
      }

      NTOP <- 40 ## how many genes to show in barplot
      sel <- intersect(sel, rownames(R))
      sel <- head(sel, NTOP)
      rho <- R[sel, "cor"]
      if (length(sel) == 1) names(rho) <- sel

      prho <- df$pcor
      names(prho) <- playbase::probe2symbol(rownames(df), pgx$genes, labeltype(), fill_na = TRUE)
      names(rho) <- playbase::probe2symbol(names(rho), pgx$genes, labeltype(), fill_na = TRUE)
      prho <- prho[match(names(rho), names(prho))]

      pd <- data.frame(
        "correlation" = rho,
        "partial correlation" = prho
      )

      return(pd)
    })

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      series <- c("correlation", "partial correlation")
      default_clrs <- c("#A6CEE3", "#1F78B4")
      pickers <- lapply(seq_along(series), function(i) {
        colourpicker::colourInput(
          session$ns(paste0("custom_color_", i)),
          label = series[i],
          value = default_clrs[i]
        )
      })
      shiny::tagList(pickers)
    })

    ## Editor: rank list for custom bar ordering
    output$rank_list <- shiny::renderUI({
      pd <- plot_data()
      shiny::req(pd)
      labels <- rownames(pd)
      shiny::req(length(labels) > 0)
      sortable::bucket_list(
        header = NULL,
        sortable::add_rank_list(
          text = "Drag to reorder",
          labels = labels,
          input_id = session$ns("rank_list_basic")
        )
      )
    })

    render_barplot <- function() {
      pd <- plot_data()
      shiny::req(pd)

      ## Editor: bar ordering
      bars_order <- input$bars_order
      if (!is.null(bars_order) && bars_order != "alphabetical") {
        if (bars_order == "ascending") {
          pd <- pd[order(pd[, 1]), , drop = FALSE]
        } else if (bars_order == "descending") {
          pd <- pd[order(-pd[, 1]), , drop = FALSE]
        } else if (bars_order == "custom" && !is.null(input$rank_list_basic)) {
          custom_order <- input$rank_list_basic
          valid <- custom_order[custom_order %in% rownames(pd)]
          if (length(valid) > 0) pd <- pd[valid, , drop = FALSE]
        }
      }

      fig <- playbase::pgx.stackedBarplot(
        x = pd,
        ylab = "Correlation",
        xlab = "",
        showlegend = FALSE
      )

      ## Editor: palette override
      palette <- input$palette
      if (!is.null(palette) && palette == "default") palette <- "muted_light"
      if (!is.null(palette) && palette != "original") {
        fig <- plotly::plotly_build(fig)
        n_series <- 2
        if (palette == "custom") {
          COL <- sapply(seq_len(n_series), function(j) {
            val <- input[[paste0("custom_color_", j)]]
            if (is.null(val)) c("#A6CEE3", "#1F78B4")[j] else val
          })
        } else {
          COL <- omics_pal_d(palette = palette)(8)[1:n_series]
        }
        bar_idx <- 1
        for (i in seq_along(fig$x$data)) {
          if (!is.null(fig$x$data[[i]]$type) && fig$x$data[[i]]$type == "bar") {
            fig$x$data[[i]]$marker$color <- COL[bar_idx]
            bar_idx <- min(bar_idx + 1, n_series)
          }
        }
      }

      fig
    }

    barplot.RENDER <- function() {
      render_barplot() %>%
        plotly_default()
    }

    barplot.RENDER2 <- function() {
      render_barplot() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = barplot.RENDER,
      func2 = barplot.RENDER2,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(63, 100), ## resolution of plots
      pdf.width = 6, pdf.height = 4,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
