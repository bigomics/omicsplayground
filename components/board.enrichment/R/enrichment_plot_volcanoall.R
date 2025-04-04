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
    width) {
  ns <- shiny::NS(id)

  plot_options <- shiny::tagList(
    withTooltip(shiny::checkboxInput(ns("scale_per_method"), "scale per method", FALSE),
      "Scale the volcano plots individually per method..",
      placement = "right", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(ns("enrch_volcanoall_subs"), "Subsample plot",
        c("Yes", "No"),
        inline = TRUE, selected = "No"
      ),
      "Number of top genesets to consider for counting the gene frequency."
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
    download.fmt = c("png", "pdf", "svg")
  )
}

enrichment_plot_volcanoall_server <- function(id,
                                              pgx,
                                              gs_features,
                                              gs_statmethod,
                                              gs_fdr,
                                              gs_lfc,
                                              calcGsetMeta,
                                              gset_selected,
                                              watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      shiny::req(gs_features())

      meta <- pgx$gset.meta$meta
      ctx <- names(meta)[!grepl("^IA:", names(meta))] ## rm *IA ctx.
      meta <- meta[ctx]
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
      FC <- vector("list", length(meta))
      Q <- vector("list", length(meta))
      names(FC) <- names(meta)
      names(Q) <- names(meta)
      for (i in names(meta)) {
        mx <- calcGsetMeta(i, gsmethod, pgx = pgx)
        FC[[i]] <- mx[, "fc", drop = FALSE]
        Q[[i]] <- mx[, "qv", drop = FALSE]
      }

      # Prepare output matrices
      matF <- do.call(cbind, FC)
      matQ <- do.call(cbind, Q)
      colnames(matF) <- names(FC)
      colnames(matQ) <- names(FC)

      pd <- list(
        FC = matF,
        Q = matQ,
        sel.gsets = sel.gsets,
        fdr = fdr,
        lfc = lfc,
        gset_selected = gset_selected()
      )
      pd
    })

    plotly_plots <- function(cex = 3, yrange = 0.5, n_rows = 2, margin_l = 50, margin_b = 50) {
      pd <- plot_data()
      shiny::req(pd)

      # Input vars
      F <- pd$FC
      Q <- pd$Q
      fdr <- pd[["fdr"]]
      lfc <- pd[["lfc"]]
      # Call volcano plots
      all_plts <- playbase::plotlyVolcano_multi(
        FC = F,
        Q = Q,
        fdr = fdr,
        lfc = lfc,
        cex = cex,
        names = rownames(F),
        source = "enrich_volcanoall",
        title_y = "significance (-log10q)",
        title_x = "effect size (log2FC)",
        share_axis = !input$scale_per_method,
        yrange = yrange,
        n_rows = n_rows,
        margin_l = margin_l,
        margin_b = margin_b,
        color_up_down = TRUE,
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
      qv <- pd[["Q"]]

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
        ylab = "Significance (-log10p)",
        marker.size = 1.2,
        showlegend = FALSE,
        title = NULL
      )
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
        card = x$card
      )
    })
  }) ## end module-server
} ## server
