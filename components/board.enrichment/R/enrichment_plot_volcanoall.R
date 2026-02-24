##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcanoall_ui <- function(
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
    withTooltip(
      shiny::checkboxInput(ns("scale_per_method"), "Scale per method", FALSE),
      "Scale the volcano plots individually per method..",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(ns("n_contrasts"), "Number of contrasts shown:",
        c("1", "4", "6", "all"),
        inline = TRUE, selected = "6"
      ),
      "Number of contrasts shown in Volcano plots."
    )
  )

  PlotModuleUI(
    ns("plot"),
    plotlib = c("plotly", "ggplot"),
    title = title,
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot_options,
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

enrichment_plot_volcanoall_server <- function(id,
                                              pgx,
                                              gs_contrast,
                                              gs_features,
                                              gs_statmethod,
                                              gs_fdr,
                                              gs_lfc,
                                              show_pv,
                                              calcGsetMeta,
                                              gset_selected,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(gs_features())

      meta <- pgx$gset.meta$meta
      ctx <- names(meta)[!grepl("^IA:", names(meta))]

      sel_ctx <- gs_contrast()
      if (input$n_contrasts == "all") {
        ctx1 <- ctx[1:length(ctx)]
      } else if (input$n_contrasts == "1") {
        ctx1 <- sel_ctx
      } else {
        nn <- as.numeric(input$n_contrasts) - 1
        ctx1 <- ctx[which(ctx != sel_ctx)]
        ctx1 <- ctx1[c(1:nn)]
        ctx1 <- ctx1[!is.na(ctx1)]
        ctx1 <- unique(c(sel_ctx, ctx1))
      }

      meta <- meta[ctx1]
      gsmethod <- colnames(meta[[1]]$fc)
      gsmethod <- gs_statmethod()
      if (is.null(gsmethod) || length(gsmethod) == 0) {
        return(NULL)
      }

      fdr <- as.numeric(gs_fdr())
      lfc <- as.numeric(gs_lfc())
      sel.gsets <- rownames(meta[[1]])
      gset_collections <- playbase::pgx.getGeneSetCollections(gsets = rownames(pgx$gsetX))
      sel.gsets <- gset_collections[[gs_features()]]

      # Calc. meta scores and get Q and FC
      FC <- Q <- P <- vector("list", length(meta))
      names(FC) <- names(Q) <- names(P) <- names(meta)
      i <- 1
      for (i in names(meta)) {
        mx <- calcGsetMeta(i, gsmethod, pgx = pgx)
        FC[[i]] <- mx[, "fc", drop = FALSE]
        Q[[i]] <- mx[, "qv", drop = FALSE]
        P[[i]] <- mx[, "pv", drop = FALSE]
      }

      # Prepare output matrices
      matF <- do.call(cbind, FC)
      matQ <- do.call(cbind, Q)
      matP <- do.call(cbind, P)
      colnames(matF) <- colnames(matQ) <- colnames(matP) <- names(FC)

      S <- matQ
      title_y <- "Significance (-log10q)"
      if (show_pv()) {
        S <- matP
        title_y <- "Significance (-log10p)"
      }

      pd <- list(
        FC = matF,
        S = S,
        sel.gsets = sel.gsets,
        fdr = fdr,
        lfc = lfc,
        title_y = title_y,
        gset_selected = gset_selected()
      )
      pd
    })

    plotly_plots <- function(cex = 3, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      F <- pd$FC
      S <- pd$S
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      title_y <- pd[["title_y"]]

      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = F,
        Q = S,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        names = rownames(F),
        source = "enrich_volcanoall",
        title_x = "Effect size (log2FC)",
        title_y = title_y,
        share_axis = !input$scale_per_method,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
        colors = c(
          up     = if (!is.null(input$color_up))   input$color_up   else get_color_theme()$primary,
          notsig = "#707070AA",
          down   = if (!is.null(input$color_down)) input$color_down else get_color_theme()$secondary,
          notsel = "#cccccc88"
        ),
        by_sig = FALSE,
        highlight = pd[["gset_selected"]],
        label = pd[["gset_selected"]]
      )

      return(all_plts)
    }

    modal_plotly.RENDER <- function() {
      fig <- plotly_plots(cex = 3, yrange = 0.5, n_rows = 2, margin_b = 20, margin_l = 50) %>%
        playbase::plotly_build_light(.)
      return(fig)
    }

    big_plotly.RENDER <- function() {
      fig <- plotly_plots(yrange = 0.02, n_rows = 3, margin_b = 20, margin_l = 20) %>%
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
      qv <- pd[["S"]]

      gene_names <- rep(rownames(fc), each = ncol(fc))
      fc <- data.frame(fc) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "fc"
        )
      qv <- data.frame(qv) %>%
        tidyr::pivot_longer(
          cols = everything(), # Select all columns to pivot
          names_to = "facet", # Name of the new column for timepoints
          values_to = "qv"
        )
      facet <- fc$facet
      x <- fc$fc
      y <- qv$qv
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
      list(plotlib = "plotly", func = modal_plotly.RENDER, func2 = big_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = base.plots, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "plot",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(72, 85), # resolution of plots
        pdf.width = 10,
        pdf.height = 5,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  }) ## end module-server
} ## server
