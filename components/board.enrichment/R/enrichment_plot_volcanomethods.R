##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanomethods_ui <- function(
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


  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_method"), "Scale per method", TRUE),
      "Scale the volcano plots individually per method..",
      placement = "right", options = list(container = "body")
    )
  )

  PlotModuleUI(
    ns("plot"),
    plotlib = c("plotly", "ggplot"),
    title = title,
    caption = caption,
    options = plot_options,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    cards = TRUE,
    card_names = c("dynamic", "static"),
    download.fmt = c("png", "pdf", "svg"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "volcano"
  )
}

enrichment_plot_volcanomethods_server <- function(id,
                                                  pgx,
                                                  gs_features,
                                                  gs_contrast,
                                                  gs_fdr,
                                                  gs_lfc,
                                                  show_pv,
                                                  gset_selected,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X, gs_features(), gs_contrast())

      cmp <- gs_contrast()
      mx <- pgx$gset.meta$meta[[cmp]]
      FC <- unclass(mx$fc)
      P <- unclass(mx$p)
      Q <- unclass(mx$q)
      rownames(P) <- rownames(Q) <- rownames(FC) <- rownames(mx)
      FC[which(is.infinite(FC))] <- NA
      P[which(is.na(P))] <- 1
      Q[which(is.na(Q))] <- 1

      fdr <- as.numeric(gs_fdr())
      lfc <- as.numeric(gs_lfc())
      gset_collections <- playbase::pgx.getGeneSetCollections(
        gsets = rownames(pgx$gsetX)
      )
      sel.gsets <- gset_collections[[gs_features()]]
      nlq <- -log10(1e-99 + unlist(Q))

      S <- Q
      title_y <- "Significance (-log10q)"
      if (show_pv()) {
        S <- P
        title_y <- "Significance (-log10p)"
      }

      pd <- list(
        FC = FC,
        S = S,
        sel.gsets = sel.gsets,
        fdr = fdr,
        lfc = lfc,
        title_y = title_y,
        gset_selected = gset_selected()
      )
      pd
    })

    plotly_plots <- function(cex = 3, yrange = 0.5, n_rows = 2,
                             margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      highlight <- pd[["gset_selected"]]
      label <- pd[["gset_selected"]]
      title_y <- pd[["title_y"]]
      ## meta tables
      F <- pd$FC
      S <- pd$S
      sel.gsets <- pd$sel.gsets
      rm(pd)
      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = F,
        Q = S,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        names = rownames(F),
        title_x = "Effect size (log2FC)",
        title_y = title_y,
        share_axis = !input$scale_per_method,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
        label = label,
        highlight = highlight,
        by_sig = FALSE
      )

      return(all_plts)
    }

    # Render functions
    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 3, yrange = 0.5, n_rows = NULL, margin_b = 20, margin_l = 50) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    plotly.RENDER <- function() {
      fig <- plotly_plots(yrange = 0.02, n_rows = NULL, margin_b = 20, margin_l = 20) %>%
        plotly::style(
          marker.size = 6
        ) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    base.plots <- function() {
      pd <- plot_data()
      shiny::req(pd)

      fc <- pd[["FC"]]
      sig <- pd[["S"]] ## can be q or p value

      gene_names <- rep(rownames(fc), each = ncol(fc))
      fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "fc"
        )
      sig <- data.frame(sig) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "sig"
        )
      facet <- fc$facet
      x <- fc$fc
      y <- sig$sig
      y <- -log10(y + 1e-12)

      ## Editor: custom labels
      if (isTRUE(input$custom_labels)) {
        label_features <- if (is.null(input$label_features) || input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
      } else {
        label_features <- pd[["gset_selected"]]
      }

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        pd[["gset_selected"]]
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

      p <- playbase::ggVolcano(
        x,
        y,
        gene_names,
        facet = facet,
        label = label_features,
        highlight = highlight,
        label.cex = label_size,
        lfc = pd[["lfc"]],
        psig = pd[["fdr"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["title_y"]],
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
      list(plotlib = "plotly", func = plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = base.plots, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "plot",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(75, 90), # resolution of plots
        pdf.width = 10,
        pdf.height = 5,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  })
}
