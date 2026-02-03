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
expression_plot_volcano_ui <- function(id,
                                       label = "",
                                       title,
                                       info.text,
                                       info.methods,
                                       info.references,
                                       info.extra_link,
                                       caption,
                                       height,
                                       width) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("pltmod"),
    label = label,
    plotlib = c("plotly", "ggplot"),
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    title = title,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height,
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
#' @param comp1
#' @param fdr
#' @param lfc
#' @param features
#' @param res
#' @param sel1
#' @param df1
#' @param watermark
#'
#'
#'
#' @export
expression_plot_volcano_server <- function(id,
                                           comp1,
                                           fdr,
                                           lfc,
                                           show_pv,
                                           res,
                                           genes_selected,
                                           labeltype = reactive("symbol"),
                                           watermark = FALSE,
                                           pval_cap,
                                           pgx) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(res())
      shiny::req(length(comp1()) > 0)

      comp1 <- comp1()
      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      res <- res()

      symbols <- playbase::probe2symbol(
        probes = rownames(res), res, query = "symbol", fill_na = TRUE
      )

      pval_cap <- pval_cap()
      qval <- pmax(res$meta.q, pval_cap)
      pval <- pmax(res$meta.p, pval_cap)
      x <- res$logFC
      y <- -log10(qval + pval_cap)
      ylab <- "Significance (-log10q)"
      if (show_pv()) {
        y <- -log10(pval + pval_cap)
        ylab <- "Significance (-log10p)"
      }

      names <- ifelse(is.na(res$gene_title), rownames(res), res$gene_title)
      label.names <- playbase::probe2symbol(rownames(res), pgx$genes, labeltype(), fill_na = TRUE)

      shape <- rep("circle", nrow(res))
      names(shape) <- rownames(res)
      if (any(is.na(pgx$counts)) && !any(is.na(pgx$X))) {
        jj <- which(!is.na(pgx$contrasts[, comp1]))
        if (any(jj)) {
          counts <- pgx$counts[rownames(res), rownames(pgx$contrasts)[jj], drop = FALSE]
          nas <- apply(counts, 1, function(x) sum(is.na(x)))
          na.features <- names(nas)[which(nas > 0)]
          shape[match(na.features, names(shape))] <- "cross"
        }
      }

      return(list(
        x = x,
        y = y,
        ylab = ylab,
        symbols = symbols,
        features = rownames(res),
        names = names,
        shape = shape,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes,
        fdr = fdr,
        lfc = lfc,
        label.names = label.names
      ))
    })

    plotly.RENDER <- function(marker.size = 4, lab.cex = 1) {
      pd <- plot_data()
      shiny::req(pd)
      plt <- playbase::plotlyVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd$features,
        label.names = pd[["label.names"]],
        source = "plot1",
        marker.type = "scattergl",
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]],
        label.cex = lab.cex,
        shape = pd[["shape"]],
        group.names = c("group1", "group0"),
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = marker.size,
        showlegend = FALSE,
        color_up_down = TRUE
      )
      plt
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly.RENDER(marker.size = 8, lab.cex = 1.5) %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(
            font = list(size = 18)
          )
        )
      fig
    }

    base.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      names <- pd$features

      if (input$custom_labels) {
        label_features <- if (input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
        label_features <- label_features
      } else {
        label_features <- pd[["lab.genes"]]
      }

      highlight <- if (input$color_selection) {
        label_features
      } else {
        pd[["sel.genes"]]
      }
      if (input$cutoff_type == "hyperbolic") {
        highlight <- label_features
      }

      # Custom colors (used when ggprism_colors is FALSE)
      plot_colors <- c(
        up = input$color_up,
        notsig = "#707070AA",
        notsel = "#cccccc88",
        down = input$color_down
      )

      box_padding <- if (is.null(input$box_padding) || is.na(input$box_padding)) 0.1 else input$box_padding
      min_segment_length <- if (is.null(input$min_segment_length) || is.na(input$min_segment_length)) 0 else input$min_segment_length
      label_box <- if (is.null(input$label_box)) TRUE else input$label_box
      segment_linetype <- if (is.null(input$segment_linetype)) 1 else as.integer(input$segment_linetype)
      ## Hyperbolic cutoff settings
      use_hyperbola <- !is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic"
      hyperbola_k <- if (!is.null(input$hyperbola_k)) input$hyperbola_k else 1

      ## ggprism settings
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
        x = pd[["x"]],
        y = pd[["y"]],
        names = names,
        highlight = highlight,
        label = label_features,
        label.names = pd[["label.names"]],
        label.cex = input$label_size,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = input$marker_size,
        showlegend = FALSE,
        title = NULL,
        axis.text.size = input$axis_text_size,
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

      if (input$margin_checkbox) {
        margin_top <- ifelse(is.na(input$margin_top), 10, input$margin_top)
        margin_right <- ifelse(is.na(input$margin_right), 10, input$margin_right)
        margin_bottom <- ifelse(is.na(input$margin_bottom), 10, input$margin_bottom)
        margin_left <- ifelse(is.na(input$margin_left), 10, input$margin_left)
        p <- p + ggplot2::theme(
          plot.margin = ggplot2::margin(
            t = margin_top,
            r = margin_right,
            b = margin_bottom,
            l = margin_left,
            unit = "pt"
          )
        )
      }

      if (input$aspect_ratio_checkbox) {
        if (is.na(input$aspect_ratio)) {
          p <- p + ggplot2::theme(
            aspect.ratio = 0.5
          )
        } else {
          p <- p + ggplot2::theme(
            aspect.ratio = input$aspect_ratio
          )
        }
      }

      p
    }

    base.RENDER.modal <- function() {
      pd <- plot_data()
      shiny::req(pd)

      names <- pd$features
      box_padding <- if (is.null(input$box_padding) || is.na(input$box_padding)) 0.1 else input$box_padding
      min_segment_length <- if (is.null(input$min_segment_length) || is.na(input$min_segment_length)) 0 else input$min_segment_length
      label_box <- if (is.null(input$label_box)) TRUE else input$label_box
      segment_linetype <- if (is.null(input$segment_linetype)) 1 else as.integer(input$segment_linetype)

      playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = names,
        highlight = pd[["sel.genes"]],
        label = pd[["lab.genes"]],
        label.names = pd[["label.names"]],
        label.cex = 6,
        axis.text.size = 22,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = 1.8,
        showlegend = FALSE,
        title = NULL,
        box.padding = box_padding,
        min.segment.length = min_segment_length,
        label.box = label_box,
        segment.linetype = segment_linetype
      )
    }

    plot_data_csv <- function() {
      dt <- plot_data()
      df <- data.frame(dt$x, dt$y)
      colnames(df) <- c("x", "y")
      rownames(df) <- make.unique(dt$symbols)
      return(df)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.RENDER, func2 = base.RENDER.modal, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data_csv,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card,
        parent_session = session,
        download.contrast.name = comp1
      )
    })
  }) ## end of moduleServer
}
