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
expression_plot_topgenes_ui <- function(
  id,
  title,
  caption,
  info.text,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  topgenes_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("gx_logscale"), "log scale", TRUE),
      "Logarithmic scale the counts (abundance levels).",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("gx_grouped"), "group samples", TRUE),
      "Group samples by phenotype",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("gx_showothers"), "show others", FALSE),
      "Show the 'others' class (if any)",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = topgenes_opts,
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
#' @param res
#' @param ii
#' @param watermark
#'
#'
#'
#' @export
expression_plot_topgenes_server <- function(id,
                                            comp,
                                            pgx,
                                            res,
                                            rows_current,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    # #calculate required inputs for plotting ---------------------------------

    ## Editor: rank list for custom drag-and-drop ordering
    output$rank_list <- renderUI({
      pd <- plot_data()
      shiny::req(pd)

      expmat <- pd[["pgx"]]$model.parameters$exp.matrix
      ct <- expmat[, pd[["comp"]]]

      if (pd[["grouped"]]) {
        if ("contrasts" %in% names(pd[["pgx"]])) {
          contr.labels <- pd[["pgx"]]$contrasts[, pd[["comp"]]]
          labels <- unique(contr.labels[ct != 0])
        } else {
          comp1 <- sub(".*:", "", pd[["comp"]])
          labels <- rev(strsplit(comp1, split = "_vs_|_VS_")[[1]])
        }
        if (pd[["showothers"]] && any(ct == 0)) labels <- c(labels, "other")
      } else {
        labels <- rownames(expmat)
        if (!pd[["showothers"]]) labels <- rownames(expmat)[ct != 0]
      }

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
      comp <- comp() # input$gx_contrast
      shiny::req(pgx$X)

      res <- res()
      if (is.null(res) || nrow(res) == 0) {
        return(NULL)
      }

      ## filter on active rows (using search)
      res <- res[rows_current(), , drop = FALSE]
      if (nrow(res) == 0) {
        return(NULL)
      }

      grouped <- input$gx_grouped
      logscale <- input$gx_logscale
      showothers <- input$gx_showothers

      ylab <- ifelse(logscale, "Expression (log2)", "Expression")
      ny <- nrow(pgx$samples) ## ???!!
      show.names <- ifelse(!grouped & ny > 25, FALSE, TRUE)
      srt <- 35
      sumlen.grpnames <- sum(nchar(strsplit(sub(".*:", "", comp), split = "_vs_")[[1]]))
      if (show.names && sumlen.grpnames <= 20) srt <- 0

      return(list(
        res = res,
        pgx = pgx,
        comp = comp,
        grouped = grouped,
        showothers = showothers,
        ylab = ylab,
        srt = srt,
        logscale = logscale,
        show.names = show.names
      ))
    })


    render_plotly <- function(annot.y = 1.05, xaxis.fontsize = 10, title.cex = 1) {
      pd <- plot_data()
      shiny::req(pd)

      nplots <- min(8, nrow(pd[["res"]]))
      if (pd$grouped) {
        nplots <- min(18, nrow(pd[["res"]]))
      }

      ## Editor: effective bar color and title color
      bar_color <- input$bar_color
      effective_color <- if (!is.null(bar_color)) bar_color else "#A6CEE3"
      color_changed <- !is.null(bar_color) && bar_color != "#A6CEE3"

      ## Editor: bars order
      bars_order <- input$bars_order

      plts <- list()

      for (i in 1:nplots) {
        gene <- rownames(pd[["res"]])[i]

        ## manual plotly annotation for plot title (with synced color)
        annotations <- list(
          x = 0.5,
          y = annot.y,
          text = playbase::probe2symbol(gene, pgx$genes, "gene_name", fill_na = TRUE),
          font = list(size = 10 * title.cex, color = effective_color),
          xref = "paper",
          yref = "paper",
          xanchor = "bottom",
          yanchor = "top",
          showarrow = FALSE
        )

        p <- playbase::pgx.plotExpression(
          pd[["pgx"]],
          probe = gene,
          comp = pd[["comp"]],
          grouped = pd[["grouped"]],
          max.points = 200,
          logscale = pd[["logscale"]],
          collapse.others = TRUE,
          showothers = pd[["showothers"]],
          ylab = pd[["ylab"]],
          xlab = "",
          srt = pd[["srt"]],
          names = pd[["show.names"]],
          main = "",
          plotlib = "plotly",
          plotly.annotations = annotations,
          plotly.margin = list(l = 5, r = 5, b = 5, t = 20, pad = 3)
        )

        ## Editor: override bar color
        if (color_changed && !is.null(p)) {
          p <- plotly::plotly_build(p)
          for (j in seq_along(p$x$data)) {
            if (!is.null(p$x$data[[j]]$type) && p$x$data[[j]]$type == "bar") {
              p$x$data[[j]]$marker$color <- bar_color
            }
          }
        }

        ## Editor: bars order
        if (!is.null(bars_order) && !is.null(p)) {
          if (bars_order == "custom" && !is.null(input$rank_list_basic)) {
            p <- plotly::layout(p, xaxis = list(
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
            p <- plotly::layout(p, xaxis = list(categoryorder = cat_order))
          }
        }

        p <- p %>% plotly::layout(
          plot_bgcolor = "#f2f2f2",
          xaxis = list(tickfont = list(size = xaxis.fontsize))
        )

        plts[[i]] <- p
      }
      return(plts)
    }

    plotly.RENDER <- function() {
      ## layout in subplots
      plts <- render_plotly(annot.y = 1.00, xaxis.fontsize = 10, title.cex = 1)
      plts <- head(plts, 16)
      pd <- plot_data()
      ncols <- ifelse(pd[["grouped"]], 8, 4)
      nrows <- ceiling(length(plts) / ncols)
      plotly::subplot(
        plts,
        nrows = nrows,
        margin = c(0.010, 0.010, 0.04, 0.04), ## lrtb
        titleX = TRUE,
        titleY = TRUE,
        shareY = TRUE,
        shareX = TRUE
      ) %>%
        plotly::layout(
          margin = list(b = 0),
          showlegend = FALSE
        )
    }

    modal_plotly.RENDER <- function() {
      plts <- render_plotly(annot.y = 1.00, xaxis.fontsize = 14, title.cex = 1.4)
      plts <- head(plts, 18)
      pd <- plot_data()
      ncols <- ifelse(pd[["grouped"]], 6, 4)
      nrows <- ceiling(length(plts) / ncols)
      fig <- plotly::subplot(
        plts,
        nrows = nrows,
        margin = c(0.011, 0.011, 0.04, 0.03), ## lrtb
        titleX = TRUE,
        titleY = TRUE,
        shareY = TRUE,
        shareX = TRUE
      )

      fig <- fig %>%
        plotly::layout(
          font = list(size = 18),
          margin = list(b = 0),
          showlegend = FALSE
        ) %>%
        plotly::style(
          marker.size = 20
        )
      fig
    }

    plot_data_csv <- function() {
      df <- plot_data()
      df <- df$res
      df <- df[, c("AveExpr0", "AveExpr1")]
      return(df)
    }

    PlotModuleServer(
      "pltmod",
      func = plotly.RENDER,
      func2 = modal_plotly.RENDER,
      plotlib = "plotly",
      csvFunc = plot_data_csv, ##  *** downloadable data as CSV
      res = c(90, 105), ## resolution of plots
      pdf.width = 14,
      pdf.height = 3.5,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
