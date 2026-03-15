##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#' @description A shiny Module for plotting (UI code).
#' @param id
#' @param label
#' @param height
#' @param width
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
#' @description A shiny Module for plotting (server code).
#' @param id
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

      ## Position data for click-to-label (all methods, so nearest-neighbor
      ## matches the actual plotted point regardless of which facet is clicked)
      fc_cols <- mx[, "fc", drop = FALSE]
      q_cols <- mx[, "q", drop = FALSE]
      click_df <- data.frame(
        x = as.vector(as.matrix(fc_cols)),
        y = as.vector(-log10(pmax(as.matrix(q_cols), 1e-99))),
        feature_name = rep(mx.features, ncol(fc_cols))
      )

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
        label.names = label.names,
        df = click_df
      )

      return(pd)
    })

    plotly_plots <- function(cex = 2, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      sel.genes <- pd[["sel.genes"]]
      lab.genes <- pd[["lab.genes"]]
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      names <- pd[["names"]]
      label.names <- pd[["label.names"]]

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
        colors = extract_volcano_colors(input),
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
      label_features <- get_custom_labels(input, pd[["features"]], defaults = lab.genes)

      highlight <- if (isTRUE(input$color_selection)) {
        label_features
      } else {
        sel.genes
      }
      if (!is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic") {
        highlight <- label_features
      }

      ## Editor: extract settings via helpers
      plot_colors <- extract_volcano_colors(input)
      ls <- extract_label_settings(input, defaults = list(label_size = label.cex))
      gp <- extract_ggprism_params(input)

      ## Editor: hyperbolic cutoff settings
      use_hyperbola <- !is.null(input$cutoff_type) && input$cutoff_type == "hyperbolic"

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
        marker.size = ls$marker_size,
        label.cex = ls$label_size,
        ylab = "Significance (-log10q)",
        xlab = "Effect size (log2FC)",
        showlegend = FALSE,
        title = NULL,
        axis.text.size = ls$axis_text_size,
        colors = plot_colors,
        box.padding = ls$box_padding,
        min.segment.length = ls$min_segment_length,
        label.box = ls$label_box,
        segment.linetype = ls$segment_linetype,
        use_hyperbola = use_hyperbola,
        hyperbola_k = ls$hyperbola_k,
        use_ggprism = gp$use_ggprism,
        ggprism_palette = gp$ggprism_palette,
        ggprism_colors = gp$ggprism_colors,
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

    big_base.plots <- function(label.cex = 4) {
      base.plots(label.cex = 5)
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = modal_plotly.RENDER, func2 = modal_plotly.RENDER, card = 1),
      list(plotlib = "ggplot", func = base.plots, func2 = big_base.plots, card = 2)
    )

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
  })
}
