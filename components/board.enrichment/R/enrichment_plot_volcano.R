##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_volcano_ui <- function(
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

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "a",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    height = height,
    width = width,
    caption = caption,
    plotlib = c("plotly", "ggplot"),
    download.fmt = c("png", "pdf", "svg"),
    cards = TRUE,
    card_names = c("dynamic", "static"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "volcano"
  )
}

enrichment_plot_volcano_server <- function(id,
                                           pgx,
                                           gs_contrast,
                                           selected_gxmethods,
                                           gset_selected,
                                           gs_fdr,
                                           gs_lfc,
                                           show_pv,
                                           subplot.MAR,
                                           geneDetails,
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      par(mfrow = c(1, 1), mgp = c(1.2, 0.4, 0), oma = c(0, 0, 0, 0.4))
      par(mar = subplot.MAR)

      shiny::req(pgx$X)
      shiny::validate(shiny::need(!is.null(gset_selected()), tspan("Please select a geneset.", js = FALSE)))

      comp <- 1
      gs <- 1
      comp <- gs_contrast()
      shiny::req(pgx$X)

      gxmethods <- selected_gxmethods() ## from module-expression
      shiny::req(gxmethods)

      gx.meta <- pgx$gx.meta$meta[[comp]]
      meta.p <- apply(gx.meta$p[, gxmethods, drop = FALSE], 1, max) ## max p-value
      meta.q <- apply(gx.meta$q[, gxmethods, drop = FALSE], 1, max) ## max q-value
      limma1 <- data.frame(meta.fx = gx.meta$meta.fx, meta.q = meta.q, meta.p = meta.p)
      gx.annot <- pgx$genes[rownames(gx.meta), c("gene_name", "gene_title")]
      limma <- cbind(gx.annot, limma1)

      gset <- geneDetails()$feature

      gset <- playbase::probe2symbol(gset, pgx$genes, "gene_name")

      jj <- match(gset, limma$gene_name)
      sel.genes <- setdiff(limma$gene_name[jj], c(NA, "", " "))

      fdr <- 1
      fdr <- as.numeric(gs_fdr())
      fc.genes <- as.character(limma[, grep("^gene$|gene_name", colnames(limma))])
      fx <- limma[, grep("logFC|meta.fx|fc", colnames(limma))[1]]
      qval <- limma[, grep("^q|adj.P.Val|meta.q|qval|padj", colnames(limma))[1]]
      pval <- limma[, grep("^p|P.Val|Pvalue|meta.p|pval", colnames(limma))[1]]

      qval <- pmax(qval, 1e-12) ## prevent q=0
      qval[which(is.na(qval))] <- 1
      pval <- pmax(pval, 1e-12) ## prevent p=0
      pval[which(is.na(pval))] <- 1

      lfc <- 0.20
      lfc <- as.numeric(gs_lfc())

      xlim <- c(-1, 1) * max(abs(fx), na.rm = TRUE)
      ylim <- c(0, 12)

      y <- -log10(qval)
      ylab <- "Significance (-log10q)"
      ylim <- c(0, max(12, 1.1 * max(-log10(qval), na.rm = TRUE)))

      if (show_pv()) {
        y <- -log10(pval)
        ylab <- "Significance (-log10p)"
        ylim <- c(0, max(12, 1.1 * max(-log10(pval), na.rm = TRUE)))
      }

      return(list(
        x = fx,
        y = y,
        fc.genes = fc.genes,
        sel.genes = sel.genes,
        lab.cex = 1,
        fdr = fdr,
        lfc = lfc,
        ylab = ylab
      ))
    })

    plotly.RENDER <- function(marker.size = 3, lab.cex = 1) {
      pd <- plot_data()
      shiny::req(pd)
      playbase::plotlyVolcano( ## in use
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd[["fc.genes"]],
        label.names = pd[["fc.genes"]],
        source = "plot1",
        marker.type = "scattergl",
        highlight = pd[["sel.genes"]],
        label = pd[["sel.genes"]],
        label.cex = lab.cex,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = marker.size,
        displayModeBar = FALSE,
        showlegend = FALSE,
        color_up_down = TRUE
      ) %>%
        plotly::layout(margin = list(l = 0, r = 0, t = 0, b = 0))
    }

    plotly.RENDER2 <- function() {
      plotly.RENDER(marker.size = 8, lab.cex = 1.5) %>%
        plotly::layout(
          font = list(size = 18),
          legend = list(font = list(size = 18))
        )
    }

    base.RENDER <- function() {
      pd <- plot_data()
      shiny::req(pd)

      ## Editor: custom labels
      if (isTRUE(input$custom_labels)) {
        label_features <- if (is.null(input$label_features) || input$label_features == "") {
          NULL
        } else {
          strsplit(input$label_features, "\\s+")[[1]]
        }
      } else {
        label_features <- pd[["sel.genes"]]
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
      label_size <- if (!is.null(input$label_size)) input$label_size else 4
      marker_size <- if (!is.null(input$marker_size)) input$marker_size else 1
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
        x = pd[["x"]],
        y = pd[["y"]],
        title = NULL,
        names = pd[["fc.genes"]],
        label.names = pd[["fc.genes"]],
        highlight = highlight,
        label = label_features,
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = marker_size,
        label.cex = label_size,
        axis.text.size = axis_text_size,
        showlegend = FALSE,
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

    base.RENDER.modal <- function() {
      pd <- plot_data()
      shiny::req(pd)

      playbase::ggVolcano(
        x = pd[["x"]],
        y = pd[["y"]],
        names = pd[["fc.genes"]],
        label.names = pd[["fc.genes"]],
        highlight = pd[["sel.genes"]],
        label = pd[["sel.genes"]],
        psig = pd[["fdr"]],
        lfc = pd[["lfc"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["ylab"]],
        marker.size = 2,
        label.cex = 6,
        axis.text.size = 24,
        showlegend = FALSE,
        title = NULL
      )
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = plotly.RENDER, func2 = plotly.RENDER2, card = 1),
      list(plotlib = "ggplot", func = base.RENDER, func2 = base.RENDER.modal, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "plot",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card,
        parent_session = session
      )
    })
  })
}
