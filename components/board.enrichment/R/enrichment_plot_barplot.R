##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_barplot_ui <- function(
  id,
  title,
  caption,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  height,
  width
) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("ungroup"), "ungroup samples", FALSE
      ),
      "Ungroup samples in the plot",
      placement = "top",
      options = list(container = "body")
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("show_others"), "show others", FALSE
      ),
      "Show other samples as 'others' in the plot",
      placement = "top",
      options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "b",
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
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

enrichment_plot_barplot_server <- function(id,
                                           pgx,
                                           gset_selected,
                                           gs_contrast,
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

    render_subplot_barplot <- function() {
      par(mfrow = c(1, 1), mgp = c(1.8, 0.8, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      shiny::req(pgx$X)

      gset <- rownames(pgx$gsetX)[1]
      gset <- gset_selected()
      if (is.null(gset) || length(gset) == 0) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") %>%
          plotly::config(
            displayModeBar = FALSE
          ))
      }
      gset <- gset[1]
      if (!gset %in% rownames(pgx$gsetX)) {
        return(plotly::plotly_empty(type = "scatter", mode = "markers") %>%
          plotly::config(
            displayModeBar = FALSE
          ))
      }

      comp0 <- colnames(pgx$model.parameters$contr.matrix)[1]
      comp0 <- gs_contrast()

      grouped <- TRUE
      grouped <- FALSE
      grouped <- !input$ungroup
      has.design <- !is.null(pgx$model.parameters$design)
      collapse.others <- ifelse(has.design, FALSE, TRUE)

      ngrp <- length(unique(pgx$samples$group))
      srt <- ifelse(!grouped || ngrp > 4, 30, 0)
      if (!grouped && ncol(pgx$X) > 15) srt <- 60
      fig <- playbase::pgx.plotExpression(
        pgx, gset,
        comp = comp0,
        logscale = TRUE,
        level = "geneset",
        collapse.others = collapse.others,
        showothers = input$show_others,
        grouped = grouped,
        cex = 1.1,
        srt = srt,
        main = "",
        ylab = "enrichment (avg logFC)",
        xlab = playbase::breakstring(gset, 42, 80),
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

    subplot_barplot.RENDER <- function() {
      render_subplot_barplot() %>% plotly_default()
    }

    subplot_barplot.RENDER2 <- function() {
      render_subplot_barplot() %>% plotly_modal_default()
    }

    PlotModuleServer(
      "plot",
      plotlib = "plotly",
      func = subplot_barplot.RENDER,
      func2 = subplot_barplot.RENDER2,
      pdf.width = 5, pdf.height = 5,
      res = c(72, 100),
      add.watermark = watermark,
      parent_session = session
    )
  })
}
