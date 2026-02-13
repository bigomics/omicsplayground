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
expression_plot_volcanoAll_ui <- function(id,
                                          title,
                                          caption,
                                          info.text,
                                          info.methods,
                                          info.references,
                                          info.extra_link,
                                          ## labeltype,
                                          label = "",
                                          height,
                                          width) {
  ns <- shiny::NS(id)

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_plot"), "scale per plot", FALSE),
      "Scale each volcano plots individually.",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("pltmod"),
    title = title,
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
expression_plot_volcanoAll_server <- function(id,
                                              pgx,
                                              getAllContrasts,
                                              fdr,
                                              lfc,
                                              show_pv,
                                              genes_selected,
                                              labeltype = reactive("symbol"),
                                              pval_cap,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ## reactive function listening for changes in input
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)

      # Input variables
      ct <- getAllContrasts()
      F <- ct$F
      Q <- ct$Q
      P <- ct$P

      fdr <- as.numeric(fdr())
      lfc <- as.numeric(lfc())
      comp <- colnames(F)
      shiny::req(length(comp) > 0)

      ## combined matrix for output
      FQ <- data.frame(fc = F, q = Q, p = P)
      features <- rownames(FQ)
      symbols <- pgx$genes[rownames(FQ), "symbol"]
      names <- pgx$genes[rownames(FQ), "gene_title"]

      label.names <- playbase::probe2symbol(rownames(FQ), pgx$genes, labeltype(), fill_na = TRUE)

      # Input vars
      if (show_pv()) {
        ## P <- P
        title_y <- "Significance (-log10p)"
      } else {
        P <- Q
        title_y <- "Significance (-log10q)"
      }

      ## ps: FQ contains log2FC+q-value or log2FC+p-value. Depends on show_pv option.
      pd <- list(
        FQ = FQ, ## Remember: the first element is returned as downloadable CSV
        comp = comp,
        fdr = fdr,
        lfc = lfc,
        F = F,
        P = P,
        title_y = title_y,
        symbols = symbols,
        features = features,
        names = names,
        label.names = label.names,
        sel.genes = genes_selected()$sel.genes,
        lab.genes = genes_selected()$lab.genes
      )

      return(pd)
    })

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2,
                             margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      pval_cap <- pval_cap()

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = pd$F,
        Q = pd$P,
        fdr = pd$fdr,
        lfc = pd$lfc,
        cex = cex,
        names = pd$features,
        label.names = pd$label.names,
        share_axis = !input$scale_per_plot,
        title_y = pd$title_y,
        title_x = "Effect size (log2FC)",
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
        highlight = pd$sel.genes,
        label = pd$lab.genes,
        by_sig = FALSE,
        pval_cap = pval_cap
      )

      return(all_plts)
    }

    plotly.RENDER <- function() {
      fig <- plotly_plots(
        cex = 3, yrange = 0.05, n_rows = 2, margin_b = 40, margin_l = 50
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

      fc <- pd$F
      qv <- pd$P
      feature_names <- rep(rownames(fc), each = ncol(fc))
      label.names <- rep(pd$label.names, each = ncol(fc))
      pivot.fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "fc"
        )
      pivot.qv <- data.frame(qv) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "qv"
        )
      facet <- pivot.fc$facet
      x <- pivot.fc$fc

      pval_cap <- pval_cap()
      y <- -log10(pivot.qv$qv + pval_cap)

      ## Editor: custom labels
      if (isTRUE(input$custom_labels)) {
        label_features <- if (is.null(input$label_features) || input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
      } else {
        label_features <- pd[["lab.genes"]]
      }

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        pd[["sel.genes"]]
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
        x,
        y,
        names = feature_names,
        facet = facet,
        label = label_features,
        highlight = highlight,
        label.names = label.names,
        label.cex = label_size,
        xlab = "Effect size (log2FC)",
        ylab = pd$title_y,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        marker.size = marker_size,
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

    big_base.plots <- function() {
      base.plots(label.cex = 5)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly.RENDER, func2 = big_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = big_base.plots, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "pltmod",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(70, 90), # resolution of plots
        pdf.width = 12,
        pdf.height = 5,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  }) ## end of moduleServer
}
