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
    width) {
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
    download.fmt = c("png", "pdf", "svg")
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

      playbase::ggVolcano(
        x,
        y,
        gene_names,
        facet = facet,
        label = pd[["gset_selected"]],
        highlight = pd[["gset_selected"]],
        label.cex = 5,
        lfc = pd[["lfc"]],
        psig = pd[["fdr"]],
        xlab = "Effect size (log2FC)",
        ylab = pd[["title_y"]],
        marker.size = 1.2,
        showlegend = FALSE,
        title = NULL
      )
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
        card = x$card
      )
    })
  })
}
