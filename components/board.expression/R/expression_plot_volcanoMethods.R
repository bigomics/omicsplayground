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
expression_plot_volcanoMethods_ui <- function(
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

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_plot"), "scale plots", FALSE),
      "Scale each volcano plot individually.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("pltmod"),
    title = "Volcano plots for all methods",
    label = label,
    plotlib = c("plotly", "ggplot"),
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_options,
    download.fmt = c("png", "pdf", "csv", "svg"),
    height = height,
    width = width,
    cards = TRUE,
    card_names = c("dynamic", "static"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "volcano"
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
expression_plot_volcanoMethods_server <- function(id,
                                                  pgx,
                                                  comp, # input$gx_contrast
                                                  fdr, # input$gx_fdr
                                                  lfc, # input$gx_lfc
                                                  show_pv,
                                                  genes_selected,
                                                  labeltype = reactive("symbol"),
                                                  pval_cap,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(comp())
      comp <- comp()

      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())

      mx <- pgx$gx.meta$meta[[comp]]
      mx.features <- rownames(mx)
      mx.symbols <- pgx$genes[mx.features, "symbol"]
      mx.names <- pgx$genes[mx.features, "gene_title"]
      label.names <- playbase::probe2symbol(mx.features, pgx$genes, labeltype(), fill_na = TRUE)

      pd <- list(
        mx = mx,
        pgx = pgx,
        fdr = fdr,
        lfc = lfc,
        comp = comp,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        names = mx.names,
        symbols = mx.symbols,
        features = mx.features,
        label.names = label.names
      )

      return(pd)
    })

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      sel.genes <- pd[["sel.genes"]]
      lab.genes <- pd[["lab.genes"]]
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      names <- pd[["names"]]
      label.names <- pd[["label.names"]]

      ## meta tables
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      x <- mx[, "fc", drop = FALSE]
      y <- mx[, "q", drop = FALSE]
      title_y <- "Significance (-log10q)"
      if (show_pv()) {
        y <- mx[, "p", drop = FALSE]
        title_y <- "Significance (-log10p)"
      }

      mx.features <- rownames(mx)
      mx.symbols <- pd[["symbols"]]
      mx.names <- pd[["names"]]

      pval_cap <- pval_cap()

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = x,
        Q = y,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        names = mx.features,
        label.names = label.names,
        highlight = sel.genes,
        label = lab.genes,
        title_y = title_y,
        title_x = "Effect size (log2FC)",
        share_axis = !input$scale_per_plot,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
        by_sig = FALSE,
        pval_cap = pval_cap
      )
      return(all_plts)
    }


    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(
        cex = 3, yrange = 0.05, n_rows = 2, margin_b = 50, margin_l = 70
      ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(
        yrange = 0.02, n_rows = 3, margin_b = 70, margin_l = 70
      ) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    base.plots <- function(label.cex = 4) {
      pd <- plot_data()
      shiny::req(pd)

      sel.genes <- pd[["sel.genes"]]
      lab.genes <- pd[["lab.genes"]]
      sel.genes <- playbase::probe2symbol(sel.genes, pgx$genes, labeltype(), fill_na = TRUE)
      lab.genes <- playbase::probe2symbol(lab.genes, pgx$genes, labeltype(), fill_na = TRUE)
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      names <- pd[["names"]]
      label.names <- pd[["label.names"]]

      mx <- pd[["mx"]]
      fc <- mx[, "fc", drop = FALSE][[1]] |> unclass()
      qv <- mx[, "q", drop = FALSE][[1]] |> unclass()
      rownames(fc) <- names
      rownames(qv) <- names

      gene_names <- rep(rownames(fc), each = ncol(fc))
      label.names <- rep(label.names, each = ncol(fc))
      fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(),
          names_to = "facet",
          values_to = "fc"
        )
      qv <- data.frame(qv) %>%
        tidyr::pivot_longer(
          cols = everything(),
          names_to = "facet",
          values_to = "qv"
        )
      facet <- fc$facet
      x <- fc$fc
      y <- qv$qv

      pval_cap <- pval_cap()
      y <- -log10(y + pval_cap)

      ## Editor: custom labels
      if (isTRUE(input$custom_labels)) {
        label_features <- if (is.null(input$label_features) || input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
      } else {
        label_features <- lab.genes
      }

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        sel.genes
      }
      if (!is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic") {
        highlight <- label_features
      }

      ## Editor: custom colors
      plot_colors <- c(
        up = if (!is.null(input$color_up)) input$color_up else "#f23451",
        notsig = "#707070AA",
        notsel = "#cccccc88",
        down = if (!is.null(input$color_down)) input$color_down else "#3181de"
      )

      ## Editor: label settings
      label_size <- if (!is.null(input$label_size)) input$label_size else label.cex
      marker_size <- if (!is.null(input$marker_size)) input$marker_size else 1.2
      axis_text_size <- if (!is.null(input$axis_text_size)) input$axis_text_size else 14
      box_padding <- if (is.null(input$box_padding) || is.na(input$box_padding)) 0.1 else input$box_padding
      min_segment_length <- if (is.null(input$min_segment_length) || is.na(input$min_segment_length)) 0 else input$min_segment_length
      label_box <- if (is.null(input$label_box)) TRUE else input$label_box
      segment_linetype <- if (is.null(input$segment_linetype)) 1 else as.integer(input$segment_linetype)

      ## Editor: hyperbolic cutoff settings
      use_hyperbola <- !is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic"
      hyperbola_k <- if (!is.null(input$hyperbola_k)) input$hyperbola_k else 1

      ## Editor: ggprism settings
      use_ggprism <- isTRUE(input$use_ggprism)
      ggprism_palette <- if (is.null(input$ggprism_palette)) "black_and_white" else input$ggprism_palette
      ggprism_colors <- isTRUE(input$ggprism_colors)
      ggprism_border <- isTRUE(input$ggprism_border)
      ggprism_axis_guide <- if (is.null(input$ggprism_axis_guide)) "default" else input$ggprism_axis_guide
      ggprism_show_legend <- isTRUE(input$ggprism_show_legend)
      ggprism_legend_x <- if (is.null(input$ggprism_legend_x) || is.na(input$ggprism_legend_x)) 0.95 else input$ggprism_legend_x
      ggprism_legend_y <- if (is.null(input$ggprism_legend_y) || is.na(input$ggprism_legend_y)) 0.95 else input$ggprism_legend_y
      ggprism_legend_border <- isTRUE(input$ggprism_legend_border)

      p <- playbase::ggVolcano(
        x = x,
        y = y,
        psig = fdr,
        lfc = lfc,
        facet = facet,
        names = gene_names,
        label.names = label.names,
        highlight = highlight,
        label = label_features,
        marker.size = marker_size,
        label.cex = label_size,
        ylab = "Significance (-log10q)",
        xlab = "Effect size (log2FC)",
        showlegend = FALSE,
        title = NULL,
        axis.text.size = axis_text_size,
        colors = plot_colors,
        box.padding = box_padding,
        min.segment.length = min_segment_length,
        label.box = label_box,
        segment.linetype = segment_linetype,
        use_hyperbola = use_hyperbola,
        hyperbola_k = hyperbola_k,
        use_ggprism = use_ggprism,
        ggprism_palette = ggprism_palette,
        ggprism_colors = ggprism_colors,
        ggprism_border = ggprism_border,
        ggprism_axis_guide = ggprism_axis_guide,
        ggprism_show_legend = ggprism_show_legend,
        ggprism_legend_x = ggprism_legend_x,
        ggprism_legend_y = ggprism_legend_y,
        ggprism_legend_border = ggprism_legend_border
      )

      ## Editor: margins
      if (isTRUE(input$margin_checkbox)) {
        margin_top <- ifelse(is.na(input$margin_top), 10, input$margin_top)
        margin_right <- ifelse(is.na(input$margin_right), 10, input$margin_right)
        margin_bottom <- ifelse(is.na(input$margin_bottom), 10, input$margin_bottom)
        margin_left <- ifelse(is.na(input$margin_left), 10, input$margin_left)
        p <- p + ggplot2::theme(
          plot.margin = ggplot2::margin(
            t = margin_top, r = margin_right,
            b = margin_bottom, l = margin_left, unit = "pt"
          )
        )
      }

      ## Editor: aspect ratio
      if (isTRUE(input$aspect_ratio_checkbox)) {
        ar <- if (is.na(input$aspect_ratio)) 0.5 else input$aspect_ratio
        p <- p + ggplot2::theme(aspect.ratio = ar)
      }

      p
    }

    big_base.plots <- function(label.cex = 4) {
      base.plots(label.cex = 5)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = modal_plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = big_base.plots, card = 2)
    )

    plot_data_csv <- function() {
      pd <- plot_data()
      sel.genes <- pd[["sel.genes"]]
      comp <- pd[["comp"]]
      mx <- pd[["pgx"]]$gx.meta$meta[[comp]]
      fc <- mx[which(rownames(mx) %in% sel.genes), "fc", drop = FALSE]
      qv <- mx[which(rownames(mx) %in% sel.genes), "q", drop = FALSE]
      df <- cbind(fc, qv, mx)
      return(df)
    }

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(80, 90), # resolution of plots
        pdf.width = 12,
        pdf.height = 5,
        add.watermark = watermark,
        card = x$card,
        download.contrast.name = comp,
        parent_session = session
      )
    })
  }) ## end of moduleServer
}
