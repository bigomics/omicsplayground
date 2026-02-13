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
signature_plot_volcano_ui <- function(
  id,
  title,
  info.text,
  info.methods,
  info.references,
  info.extra_link,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  plot_opts <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        inputId = ns("share_axis"),
        label = "Share X/Y axis",
        value = TRUE
      ),
      "Share X/Y axis between plots so they show the same range. Easier to compare.",
      placement = "left", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(
        inputId = ns("showlabel"),
        label = "Show labels:",
        choices = c("no", "top10", "all"),
        selected = "top10",
        inline = TRUE
      ),
      "Show labels."
    )
  )

  PlotModuleUI(ns("plot"),
    title = title,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf", "svg"),
    plotlib = c("plotly", "ggplot"),
    options = plot_opts,
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
signature_plot_volcano_server <- function(id,
                                          pgx,
                                          sigCalculateGSEA,
                                          enrichmentContrastTable,
                                          selected_gxmethods,
                                          enrichmentGeneTable,
                                          getEnrichmentGeneTable,
                                          watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ##
    plot_data <- shiny::reactive({
      shiny::req(sigCalculateGSEA())

      ## filter with table selection/search
      ii <- enrichmentContrastTable$rows_selected()
      if (is.null(ii)) {
        ii <- enrichmentContrastTable$rows_all()
      }
      shiny::req(ii)

      # Input vars
      gsea <- sigCalculateGSEA()
      ct <- rownames(gsea$output)[ii]
      mm <- selected_gxmethods()
      meta <- playbase::pgx.getMetaMatrix(pgx, methods = mm)
      fc <- meta$fc[, ct, drop = FALSE]
      qv <- meta$qv[, ct, drop = FALSE]

      features <- rownames(fc)
      symbols <- pgx$genes[rownames(fc), "symbol"]
      symbols <- playbase::probe2symbol(features, pgx$genes, "gene_name", fill_na = TRUE)

      # Get gene selected labels
      if (length(ct) == 1) {
        sel.gene <- getEnrichmentGeneTable()$symbol
        row <- enrichmentGeneTable$rows_selected()
        if (length(row)) {
          sel.gene <- sel.gene[row]
        }
      } else {
        sel.gene <- gsea$gset
      }

      pd <- list(
        fc = fc,
        qv = qv,
        features = features,
        symbols = symbols,
        gsea = gsea,
        sel.gene = sel.gene
      )

      return(pd)
    })

    plotly_plots <- function(cex = 5,
                             label.cex = 1.1,
                             yrange = 0.5,
                             n_rows = 2,
                             margin_l = 70,
                             margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)
      share_axis <- input$share_axis

      showlabel <- input$showlabel
      label <- pd[["sel.gene"]]
      if (showlabel == "no") label <- NULL
      if (showlabel == "top10") {
        if (!is.null(label)) {
          ii <- match(label, pd[["symbols"]])
          rr <- Matrix::rowMeans(pd$fc[ii, , drop = FALSE]**2) + Matrix::rowMeans(log10(pd$qv[ii, , drop = FALSE])**2)
          label <- head(label[order(-rr)], 10)
        }
      }

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = pd[["fc"]],
        Q = pd[["qv"]],
        fdr = 0.05,
        lfc = 1,
        names = pd[["features"]],
        label.names = pd[["symbols"]],
        label.cex = label.cex,
        cex = cex,
        by_sig = FALSE,
        highlight = pd[["gsea"]]$gset,
        label = label,
        share_axis = share_axis,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        interplot_margin = c(0.02, 0.02, 0.04, 0.04),
        # Remove default titles
        title_y = "",
        title_x = "",
        color_up_down = TRUE
      ) %>%
        plotly::layout(
          annotations = list(
            list(
              x = -0.0,
              xshift = -26,
              xanchor = "right",
              y = 0.5,
              text = "significance (-log10q)",
              font = list(size = 16),
              textangle = 270,
              showarrow = FALSE,
              xref = "paper",
              yref = "paper"
            ),
            list(
              x = 0.5,
              y = 0.0,
              yshift = -20,
              yanchor = "top",
              text = "effect size (log2FC)",
              font = list(size = 16),
              showarrow = FALSE,
              xref = "paper",
              yref = "paper"
            )
          )
        )

      return(all_plts)
    }

    big_plotly.RENDER <- function() {
      nr <- length(enrichmentContrastTable$rows_all())
      n_rows <- floor(sqrt(nr))
      fig <- plotly_plots(
        cex = 8,
        label.cex = 1.3,
        yrange = 0.02,
        n_rows = n_rows,
        margin_b = 45,
        margin_l = 70
      ) %>%
        plotly::style(
          marker.size = 7
        ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    base.plots <- function() {
      pd <- plot_data()
      shiny::req(pd)

      fc <- pd[["fc"]]
      qv <- pd[["qv"]]

      gene_names <- rep(rownames(fc), each = ncol(fc))
      label.names <- rep(pd[["symbols"]], each = ncol(fc))

      showlabel <- input$showlabel
      label <- pd[["sel.gene"]]
      if (showlabel == "no") label <- NULL
      if (showlabel == "top10") {
        ii <- match(label, pd[["features"]])
        rr <- Matrix::rowMeans(pd$fc[ii, , drop = FALSE]**2) + Matrix::rowMeans(log10(pd$qv[ii, , drop = FALSE])**2)
        label <- head(label[order(-rr)], 10)
      }

      ## Editor: custom labels
      if (isTRUE(input$custom_labels)) {
        label_features <- if (is.null(input$label_features) || input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
      } else {
        label_features <- label
      }

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        pd[["gsea"]]$gset
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
      label_size <- if (!is.null(input$label_size)) input$label_size else 5
      marker_size <- if (!is.null(input$marker_size)) input$marker_size else 1.2
      axis_text_size <- if (!is.null(input$axis_text_size)) input$axis_text_size else 14
      box_padding <- if (is.null(input$box_padding) || is.na(input$box_padding)) 0.1 else input$box_padding
      min_segment_length <- if (is.null(input$min_segment_length) || is.na(input$min_segment_length)) 0 else input$min_segment_length
      label_box <- if (is.null(input$label_box)) TRUE else input$label_box
      segment_linetype <- if (is.null(input$segment_linetype)) 1 else as.integer(input$segment_linetype)

      ## Editor: hyperbolic cutoff
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

      pivot.fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(),
          names_to = "facet",
          values_to = "fc"
        )
      pivot.qv <- data.frame(qv) %>%
        tidyr::pivot_longer(
          cols = everything(),
          names_to = "facet",
          values_to = "qv"
        )
      facet <- factor(pivot.fc$facet, levels = colnames(fc))
      x <- pivot.fc$fc
      y <- -log10(pivot.qv$qv + 1e-12)

      p <- playbase::ggVolcano(
        x,
        y,
        gene_names,
        lfc = 1,
        psig = 0.05,
        facet = facet,
        label = label_features,
        highlight = highlight,
        label.names = label.names,
        label.cex = label_size,
        xlab = "Effect size (log2FC)",
        ylab = "Significance (-log10p)",
        marker.size = marker_size,
        axis.text.size = axis_text_size,
        showlegend = FALSE,
        title = NULL,
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

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly_plots, func2 = big_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = base.plots, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "plot",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        # csvFunc = plot_data_csv,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 6,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  }) ## end of moduleServer
}
