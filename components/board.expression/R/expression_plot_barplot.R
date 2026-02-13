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
#' @param heightt
#' @param width
#'
#' @export
expression_plot_barplot_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width
) {
  ns <- shiny::NS(id)

  plots_barplot_opts <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("barplot_grouped"), "grouped", TRUE),
      "Group expression values by conditions.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("barplot_logscale"), "log scale", TRUE),
      "Show logarithmic (log2CPM) expression values.",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("barplot_showothers"), "show others", FALSE),
      "Show the 'others' class (if any)",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(ns("pltmod"),
    title = title,
    label = label,
    info.text = info.text,
    plotlib = "plotly",
    caption = caption,
    options = plots_barplot_opts,
    download.fmt = c("png", "pdf", "svg"),
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
expression_plot_barplot_server <- function(id,
                                           comp,
                                           pgx,
                                           sel,
                                           res,
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
        ## Extract group names from contrasts or contrast name
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
      comp <- comp()
      grouped <- input$barplot_grouped
      logscale <- input$barplot_logscale
      showothers <- input$barplot_showothers
      sel <- sel()
      res <- res()

      shiny::validate(shiny::need(!is.null(sel) && length(sel) > 0, tspan("Please select gene in the table.", js = FALSE)))
      psel <- sel # rownames(res)[sel]
      gene <- psel
      srt <- ifelse(grouped, 0, 35)
      main <- pgx$genes[psel, "gene_name"]

      return(list(
        pgx = pgx,
        gene = gene,
        comp = comp,
        sel = sel,
        grouped = grouped,
        logscale = logscale,
        showothers = showothers,
        srt = srt,
        main = main
      ))
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

      fig <- playbase::pgx.plotExpression(
        pgx = pd[["pgx"]],
        probe = pd[["gene"]],
        main = pd[["main"]],
        comp = pd[["comp"]],
        grouped = pd[["grouped"]],
        showothers = pd[["showothers"]],
        max.points = 200, ## slow!!
        names = TRUE,
        logscale = pd[["logscale"]],
        srt = pd[["srt"]],
        xlab = "",
        plotlib = "plotly"
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

    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      func = plotly.RENDER,
      csvFunc = plot_data,
      res = c(80, 95), ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark,
      parent_session = session
    )
  }) ## end of moduleServer
}
