##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

enrichment_plot_freq_top_gsets_ui <- function(
  id,
  title,
  info.text,
  info.methods,
  info.extra_link,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  topEnrichedFreq.opts <- shiny::tagList(
    withTooltip(
      shiny::radioButtons(ns("gs_enrichfreq_ntop"), "Number of top sets",
        c(5, 10, 15),
        inline = TRUE, selected = 15
      ),
      "Number of top genesets to consider for counting the gene frequency."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("gs_enrichfreq_gsetweight"),
        tspan("Weight by geneset size"), TRUE
      ),
      "Weight by (inverse) gene set size."
    ),
    withTooltip(
      shiny::checkboxInput(
        ns("gs_enrichfreq_fcweight"),
        "Weight by FC", TRUE
      ),
      "Weight by fold-change of current contrast."
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = "b",
    plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    options = topEnrichedFreq.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "csv", "svg"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "grouped_barplot",
    palette_default = "default"
  )
}

enrichment_plot_freq_top_gsets_server <- function(id,
                                                  pgx,
                                                  getFilteredGeneSetTable,
                                                  gs_contrast,
                                                  gseatable,
                                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      rpt <- getFilteredGeneSetTable()
      shiny::req(pgx$X, rpt, gs_contrast())

      comp <- gs_contrast()
      if (is.null(comp)) {
        return(NULL)
      }
      if (!(comp %in% names(pgx$gx.meta$meta))) {
        return(NULL)
      }

      ## filter on active rows (using search)
      ii <- gseatable$rows_current()
      req(ii)
      rpt <- rpt[ii, , drop = FALSE]

      ntop <- as.integer(input$gs_enrichfreq_ntop)
      gset.weight <- input$gs_enrichfreq_gsetweight
      fcweight <- input$gs_enrichfreq_fcweight
      return(
        list(
          pgx,
          rpt,
          ntop,
          gset.weight,
          fcweight
        )
      )
    })

    ## Compute bar matrix (used by render, custom_palette_ui, rank_list)
    bar_matrix <- shiny::reactive({
      dt <- plot_data()
      shiny::req(dt)
      pgx <- dt[[1]]
      rpt <- dt[[2]]
      ntop <- dt[[3]]
      gset.weight <- dt[[4]]
      fcweight <- dt[[5]]
      ngenes <- 30
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      fx <- rpt[, fx.col]
      names(fx) <- rownames(rpt)

      top <- rownames(rpt)
      top <- head(top, ntop)
      if (!all(top %in% colnames(pgx$GMT))) {
        return(NULL)
      }

      F <- 1 * (pgx$GMT[, top, drop = FALSE] > 0)
      F <- as.matrix(F)
      wt <- FALSE
      if (gset.weight) {
        F <- Matrix::t(Matrix::t(F) / Matrix::colSums(F, na.rm = TRUE))
        wt <- TRUE
      }
      F <- Matrix::t(Matrix::t(F) * sign(fx[top]))
      if (fcweight) {
        F <- Matrix::t(Matrix::t(F) * abs(fx[top]))
        wt <- TRUE
      }
      # sum duplicated rows
      F <- rowsum(F, row.names(F))

      F <- head(F[order(-Matrix::rowSums(abs(F), na.rm = TRUE)), , drop = FALSE], ngenes)
      F <- F[order(-Matrix::rowSums(F, na.rm = TRUE)), , drop = FALSE]

      sel.zero <- which(Matrix::rowSums(abs(F), na.rm = TRUE) < 1e-4)
      if (length(sel.zero)) F <- F[-sel.zero, , drop = FALSE]
      rownames(F) <- playbase::probe2symbol(rownames(F), pgx$genes, "gene_name", fill_na = TRUE)

      list(F = F, wt = wt)
    })

    ## Editor: dynamic color pickers for custom palette
    output$custom_palette_ui <- shiny::renderUI({
      shiny::req(input$palette == "custom")
      bm <- bar_matrix()
      shiny::req(bm)
      series <- colnames(bm$F)
      shiny::req(length(series) > 0)
      ## plotly default colorway as picker defaults
      plotly_colors <- c(
        "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A",
        "#19D3F3", "#FF6692", "#B6E880", "#FF97FF", "#FECB52",
        "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A"
      )
      pickers <- lapply(seq_along(series), function(i) {
        colourpicker::colourInput(
          session$ns(paste0("custom_color_", i)),
          label = series[i],
          value = plotly_colors[i]
        )
      })
      shiny::tagList(pickers)
    })

    ## Editor: rank list for custom bar ordering
    output$rank_list <- shiny::renderUI({
      bm <- bar_matrix()
      shiny::req(bm)
      labels <- rownames(bm$F)
      shiny::req(length(labels) > 0)
      sortable::bucket_list(
        header = NULL,
        sortable::add_rank_list(
          text = "Drag to reorder",
          labels = labels,
          input_id = session$ns("rank_list_basic")
        )
      )
    })

    topEnrichedFreq.RENDER <- function() {
      bm <- bar_matrix()
      shiny::req(bm)
      F <- bm$F
      wt <- bm$wt

      ## Editor: bar ordering
      bars_order <- input$bars_order
      if (!is.null(bars_order) && bars_order != "alphabetical") {
        if (bars_order == "ascending") {
          F <- F[order(rowSums(F)), , drop = FALSE]
        } else if (bars_order == "descending") {
          F <- F[order(-rowSums(F)), , drop = FALSE]
        } else if (bars_order == "custom" && !is.null(input$rank_list_basic)) {
          custom_order <- input$rank_list_basic
          valid <- custom_order[custom_order %in% rownames(F)]
          if (length(valid) > 0) F <- F[valid, , drop = FALSE]
        }
      }

      fig <- playbase::pgx.stackedBarplot(
        x = F,
        ylab = ifelse(wt, "weighted frequency", "frequency"),
        xlab = tspan("genes", js = FALSE),
        showlegend = FALSE
      )

      ## Editor: palette override
      palette <- input$palette
      if (!is.null(palette) && !palette %in% c("original", "default")) {
        fig <- plotly::plotly_build(fig)
        n_series <- ncol(F)
        if (palette == "custom") {
          COL <- sapply(seq_len(n_series), function(j) {
            val <- input[[paste0("custom_color_", j)]]
            if (is.null(val)) "#636EFA" else val
          })
        } else {
          COL <- rep(omics_pal_d(palette = palette)(8), ceiling(n_series / 8))[1:n_series]
        }
        bar_idx <- 1
        for (i in seq_along(fig$x$data)) {
          if (!is.null(fig$x$data[[i]]$type) && fig$x$data[[i]]$type == "bar") {
            fig$x$data[[i]]$marker$color <- COL[bar_idx]
            bar_idx <- min(bar_idx + 1, n_series)
          }
        }
      }

      fig
    }

    plot_data_csv <- function() {
      bm <- bar_matrix()
      shiny::req(bm)
      return(bm$F)
    }

    PlotModuleServer(
      "plot",
      func = topEnrichedFreq.RENDER,
      plotlib = "plotly",
      pdf.width = 5,
      pdf.height = 5,
      res = c(68, 100),
      csvFunc = plot_data_csv,
      add.watermark = watermark,
      parent_session = session
    )
  })
}
