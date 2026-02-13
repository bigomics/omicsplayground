##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' Create UI for gene plot enrichment visualization
#'
#' @param id Id prefix for namespace
#' @param title Plot title
#' @param caption Plot caption
#' @param info.text Info text to display
#' @param height Plot height
#' @param width Plot width
#'
#' @return Shiny UI for gene plot enrichment visualization
enrichment_plot_geneplot_ui <- function(
  id,
  title,
  caption,
  info.text,
  height,
  width
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("ungroup"), "ungroup samples", FALSE),
      "Ungroup samples in the plot",
      placement = "top",
      options = list(container = "body")
    ),
    withTooltip(shiny::checkboxInput(ns("show_others"), "show others", FALSE),
      "Show other samples as 'others' in the plot",
      placement = "top",
      options = list(container = "body")
    ),
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "c",
    plotlib = "plotly",
    info.text = info.text,
    caption = caption,
    options = options,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "barplot",
    bar_color_default = "#A6CEE3"
  )
}

#' Gene Expression Plot Server Function
#'
#' Server function for generating the gene expression plot module in Shiny.
#'
#' @param id Module id
#' @param pgx PGX object
#' @param gs_contrast Get the selected contrast function
#' @param gene_selected Get the selected gene function
#' @param subplot.MAR Plot margins
#' @param watermark Add watermark
#'
#' @return Shiny module server function
enrichment_plot_geneplot_server <- function(id,
                                            pgx,
                                            gs_contrast,
                                            gene_selected,
                                            subplot.MAR,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## Editor: rank list for custom drag-and-drop ordering
    output$rank_list <- renderUI({
      shiny::req(pgx$X)

      comp0 <- gs_contrast()
      shiny::req(comp0)

      expmat <- pgx$model.parameters$exp.matrix
      ct <- expmat[, comp0]
      grouped <- !input$ungroup

      if (grouped) {
        if ("contrasts" %in% names(pgx)) {
          contr.labels <- pgx$contrasts[, comp0]
          labels <- unique(contr.labels[ct != 0])
        } else {
          comp1 <- sub(".*:", "", comp0)
          labels <- rev(strsplit(comp1, split = "_vs_|_VS_")[[1]])
        }
        if (input$show_others && any(ct == 0)) labels <- c(labels, "other")
      } else {
        labels <- rownames(expmat)
        if (!input$show_others) labels <- rownames(expmat)[ct != 0]
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

    render_subplot_geneplot <- function() {
      par(mfrow = c(1, 1), mgp = c(1.8, 0.8, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      shiny::req(pgx$X)

      comp0 <- gs_contrast()
      sel <- gene_selected()

      not.selected <- (is.null(sel) || length(sel) == 0)
      shiny::validate(shiny::need(
        not.selected == FALSE, tspan("Please select a gene", js = FALSE)
      ))

      if (not.selected) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") %>%
          plotly::config(
            displayModeBar = FALSE
          ))
      } else {
        probe <- sel$rn
        gene <- sel$gene
        gene <- playbase::probe2symbol(gene, pgx$genes, "gene_name", fill_na = TRUE)
        ngrp <- length(unique(pgx$samples$group))
        grouped <- !input$ungroup
        srt <- ifelse(!grouped || ngrp > 4, 30, 0)
        if (!grouped && ncol(pgx$X) > 15) srt <- 60
        has.design <- !is.null(pgx$model.parameters$design)
        collapse.others <- ifelse(has.design, FALSE, TRUE)

        fig <- playbase::pgx.plotExpression(
          pgx,
          probe,
          comp = comp0,
          logscale = TRUE,
          level = "gene",
          collapse.others = collapse.others,
          showothers = input$show_others,
          grouped = grouped,
          srt = srt,
          main = "",
          xlab = gene,
          plotlib = "plotly"
        )

        ## Editor: bar color
        bar_color <- input$bar_color
        effective_color <- if (!is.null(bar_color)) bar_color else "#A6CEE3"
        if (!is.null(bar_color) && bar_color != "#A6CEE3" && !is.null(fig)) {
          fig <- plotly::plotly_build(fig)
          for (j in seq_along(fig$x$data)) {
            if (!is.null(fig$x$data[[j]]$type) && fig$x$data[[j]]$type == "bar") {
              fig$x$data[[j]]$marker$color <- bar_color
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
    }

    ## -------------------------------------------------------------------
    ## Render plots
    ## -------------------------------------------------------------------

    subplot_geneplot.RENDER1 <- function() {
      render_subplot_geneplot() %>%
        plotly_default()
    }

    subplot_geneplot.RENDER2 <- function() {
      render_subplot_geneplot() %>%
        plotly_modal_default()
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = subplot_geneplot.RENDER1,
      func2 = subplot_geneplot.RENDER2,
      pdf.width = 5, pdf.height = 5,
      res = c(78, 100),
      add.watermark = watermark,
      parent_session = session
    )
  })
}
