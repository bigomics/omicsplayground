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
#' @param width
#'
#' @export
expression_plot_topfoldchange_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    caption = caption,
    info.text = info.text,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
    ns_parent = ns,
    editor = TRUE,
    plot_type = "barplot",
    bar_color_default = "#A6CEE3"
  )
}

#' Expression plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#' @param comp
#' @param pgx
#' @param sel
#' @param res
#' @param watermark
#'
#'
#'
#' @export
expression_plot_topfoldchange_server <- function(id,
                                                 comp,
                                                 pgx,
                                                 sel,
                                                 res,
                                                 watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## Editor: rank list for custom drag-and-drop ordering
    output$rank_list <- renderUI({
      pd <- plot_data()
      shiny::req(pd, pd[["fc.top"]])

      labels <- names(pd[["fc.top"]])
      labels <- labels[!is.na(pd[["fc.top"]])]

      sortable::bucket_list(
        header = NULL,
        class = "default-sortable custom-sortable",
        sortable::add_rank_list(
          input_id = session$ns("rank_list_basic"),
          text = NULL,
          labels = labels
        )
      )
    })

    plot_data <- shiny::reactive({
      comp <- comp()
      sel <- sel()
      res <- res()
      shiny::validate(shiny::need(!is.null(sel) && length(sel) > 0, tspan("Please select gene in the table.", js = FALSE)))

      psel <- sel
      gene <- pgx$genes[psel, "gene_name"]

      if (is.null(sel) || length(sel) == 0) {
        return(list(sel = sel))
      }
      if (is.null(comp) || length(comp) == 0) {
        return(NULL)
      }

      fc <- sapply(pgx$gx.meta$meta, function(x) x[psel, "meta.fx"])
      if (any(is.na(fc))) {
        shiny::validate(shiny::need(!is.na(fc), "Fold change for this feature is NA."))
      }

      top.up <- head(names(sort(fc[which(fc > 0)], decreasing = TRUE)), 10)
      top.dn <- head(names(sort(fc[which(fc < 0)], decreasing = FALSE)), 10)
      fc.top <- c(fc[top.up], fc[top.dn])
      fc.top <- fc.top[head(order(-abs(fc.top)), 15)]
      fc.top <- sort(fc.top)
      fc.top <- head(c(fc.top, rep(NA, 99)), 15)

      return(list(sel = sel, fc.top = fc.top, gene = gene))
    })

    plotly.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      if (is.null(pd[["sel"]]) || length(pd[["sel"]]) == 0) {
        frame()
        text(0.5, 0.5, tspan("No gene selected", js = FALSE), col = "black")
        return(NULL)
      }

      if (is.null(res) || is.null(sel)) {
        return(NULL)
      }

      fig <- playbase::pgx.barplot.PLOTLY(
        data = data.frame(x = names(pd[["fc.top"]]), y = as.numeric(pd[["fc.top"]])),
        x = "x",
        y = "y",
        title = pd[["gene"]],
        yaxistitle = "Fold change (log2)",
        xaxistitle = "",
        yrange = c(-1.1, 1.1) * max(abs(pd[["fc.top"]])),
        margin = list(l = 10, r = 10, b = 0, t = 25),
        grouped = FALSE
      )

      ## Editor: bar color (only override when user changed from default #A6CEE3)
      bar_color <- input$bar_color
      effective_color <- if (!is.null(bar_color)) bar_color else "#A6CEE3"
      if (!is.null(bar_color) && bar_color != "#A6CEE3" && !is.null(fig)) {
        fig <- plotly::plotly_build(fig)
        for (i in seq_along(fig$x$data)) {
          if (!is.null(fig$x$data[[i]]$type) && fig$x$data[[i]]$type == "bar") {
            fig$x$data[[i]]$marker$color <- bar_color
          }
        }
      }

      ## Editor: sync title color with bar color
      if (!is.null(fig)) {
        fig <- plotly::layout(fig, title = list(font = list(color = effective_color)))
      }

      ## Editor: bars order
      bars_order <- input$bars_order
      if (!is.null(bars_order) && !is.null(fig)) {
        if (bars_order == "custom" && !is.null(input$rank_list_basic)) {
          fig <- plotly::layout(fig, xaxis = list(
            categoryorder = "array",
            categoryarray = input$rank_list_basic
          ))
        } else {
          cat_order <- switch(bars_order,
            "alphabetical" = "category ascending",
            "ascending" = "total ascending",
            "descending" = "total descending",
            "trace"
          )
          fig <- plotly::layout(fig, xaxis = list(categoryorder = cat_order))
        }
      }

      fig
    }

    plot_data_csv <- function() {
      pd <- plot_data()
      df <- data.frame(name = names(pd[["fc.top"]]), fc = pd[["fc.top"]])
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      remove_margins = FALSE,
      csvFunc = plot_data_csv, ##  *** downloadable data as CSV
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
