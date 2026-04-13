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
compare_plot_cum_fc1_ui <- function(id,
                                    height,
                                    width,
                                    title,
                                    info.text,
                                    label) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    plotlib = "plotly",
    label = label,
    info.text = info.text,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    editor = TRUE,
    ns_parent = ns,
    plot_type = "barplot",
    bar_color_default = "#66C2A5"
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
compare_plot_cum_fc1_server <- function(id,
                                        pgx,
                                        labeltype,
                                        # dataset2,
                                        getMatrices,
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ## Override bar ordering default to ascending
    shiny::observeEvent(TRUE, once = TRUE, {
      shiny::updateSelectInput(
        session, "bars_order",
        selected = "ascending"
      )
    })

    plot_data <- reactive({
      res <- getMatrices()
      FC <- cbind(res$F1, res$F2)
      FC
    })

    ## Editor: rank list for custom drag-and-drop ordering
    output$rank_list <- shiny::renderUI({
      shiny::req(getMatrices())
      res <- getMatrices()
      FC <- cbind(res$F1, res$F2)
      F1 <- res$F1
      ii <- head(order(-rowMeans(FC**2, na.rm = TRUE)), 40)
      ii <- ii[order(rowMeans(FC[ii, ], na.rm = TRUE))]
      F1 <- F1[ii, , drop = FALSE]
      if (!is.null(rownames(F1))) {
        rownames(F1) <- make.names(playbase::probe2symbol(rownames(F1), pgx$genes, labeltype(), fill_na = TRUE), unique = TRUE)
      }
      labels <- rownames(F1)
      rank_list_ui(labels, ns, input_id = "rank_list_order")
    })

    plot_data_prepared <- shiny::reactive({
      shiny::req(getMatrices())

      res <- getMatrices()
      FC <- cbind(res$F1, res$F2)
      F1 <- res$F1
      ii <- head(order(-rowMeans(FC**2, na.rm = TRUE)), 40)
      ii <- ii[order(rowMeans(FC[ii, ], na.rm = TRUE))]
      F1 <- F1[ii, , drop = FALSE]

      if (!is.null(rownames(F1))) {
        rownames(F1) <- make.names(playbase::probe2symbol(rownames(F1), pgx$genes, labeltype(), fill_na = TRUE), unique = TRUE)
      }

      ## Editor: bar ordering
      bar_order <- if (!is.null(input$bars_order)) input$bars_order else "ascending"
      if (bar_order == "descending") {
        F1 <- F1[rev(seq_len(nrow(F1))), , drop = FALSE]
      } else if (bar_order == "alphabetical") {
        F1 <- F1[order(rownames(F1)), , drop = FALSE]
      } else if (bar_order == "custom" && !is.null(input$rank_list_order)) {
        custom_order <- intersect(input$rank_list_order, rownames(F1))
        if (length(custom_order) > 0) {
          F1 <- F1[custom_order, , drop = FALSE]
        }
      }

      F1
    })

    plot.RENDER <- shiny::reactive({
      F1 <- plot_data_prepared()
      shiny::req(F1)

      bar_color <- get_editor_color(input, "bar_color", "#66C2A5")
      gp <- extract_ggprism_params(input)

      if (gp$use_ggprism) {
        y_vals <- rowMeans(F1, na.rm = TRUE)
        d <- data.frame(
          x = factor(rownames(F1), levels = rownames(F1)),
          y = y_vals
        )

        p <- playbase::pgx.barplot.GGPLOT(
          data = d,
          x = "x",
          y = "y",
          fillcolor = bar_color,
          yaxistitle = "log2FC",
          xaxistitle = "",
          title = "",
          grouped = FALSE
        )

        x_map <- p$mapping$x
        if (!is.null(x_map)) {
          suppressMessages(
            p <- p +
              ggplot2::aes(fill = !!x_map) +
              ggplot2::scale_fill_manual(values = rep(bar_color, 50)) +
              ggplot2::guides(fill = "none")
          )
        }
        p <- apply_ggprism_fill(p, gp)
        p <- apply_ggprism_theme(p, gp, x_angle = 90)
        p <- apply_editor_theme(p, input)
        fig <- ggplot_as_plotly_image(p)
      } else {
        fillcolor <- rep(bar_color, ncol(F1))

        d <- data.frame(
          x = factor(rownames(F1), levels = rownames(F1)),
          y = F1
        )
        ycols <- colnames(d[, 2:ncol(d)])

        suppressWarnings(
          fig <- playbase::pgx.barplot.PLOTLY(
            data = d,
            x = "x",
            y = ycols,
            fillcolor = fillcolor,
            yaxistitle = "log2FC",
            xaxistitle = "",
            title = "",
            grouped = FALSE
          )
        )
      }
      if (!gp$use_ggprism) fig <- apply_plotly_editor_theme(fig, input)
      return(fig)
    })

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = plot.RENDER,
      csvFunc = plot_data,
      res = c(80, 98), ## resolution of plots
      pdf.width = 10, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
