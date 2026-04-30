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

      ## Click-to-label data (per-facet coordinates)
      click_df <- data.frame(
        x = as.vector(matF),
        y = as.vector(-log10(S + 1e-12)),
        feature_name = rep(rownames(matF), ncol(matF))
      )

      pd <- list(
        FC = matF,
        S = S,
        sel.gsets = sel.gsets,
        fdr = fdr,
        lfc = lfc,
        title_y = title_y,
        gset_selected = gset_selected(),
        df = click_df
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
        colors = extract_volcano_colors(input),
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
      label_features <- get_custom_labels(input, rownames(pd[["FC"]]), defaults = pd[["gset_selected"]])

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        pd[["gset_selected"]]
      }
      if (!is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic") {
        highlight <- label_features
      }

      ## Editor: extract settings via helpers
      plot_colors <- extract_volcano_colors(input)
      ls <- extract_label_settings(input)
      gp <- extract_ggprism_params(input)

      ## Editor: hyperbolic cutoff
      use_hyperbola <- !is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic"

      p <- playbase::ggVolcano(
        x,
        y,
        gene_names,
        facet = facet,
        label = label_features,
        highlight = highlight,
        label.cex = ls$label_size,
        lfc = pd[["lfc"]],
        psig = pd[["fdr"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["title_y"]],
        marker.size = ls$marker_size,
        axis.text.size = ls$axis_text_size,
        showlegend = FALSE,
        title = NULL,
        colors = plot_colors,
        box.padding = ls$box_padding,
        min.segment.length = ls$min_segment_length,
        label.box = ls$label_box,
        segment.linetype = ls$segment_linetype,
        use_hyperbola = use_hyperbola,
        hyperbola_k = ls$hyperbola_k,
        use_ggprism = gp$use_ggprism,
        ggprism_border = gp$ggprism_border,
        ggprism_axis_guide = gp$ggprism_axis_guide,
        ggprism_show_legend = gp$ggprism_show_legend,
        ggprism_legend_x = gp$ggprism_legend_x,
        ggprism_legend_y = gp$ggprism_legend_y,
        ggprism_legend_border = gp$ggprism_legend_border
      )

      ## Editor: margins & aspect ratio
      p <- apply_editor_theme(p, input)

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
