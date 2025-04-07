##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_gene_map_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.references,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  plot.opts <- shiny::tagList(
    shiny::selectInput(ns("umap_nlabel"), "nr labels:",
      c(0, 10, 20, 50, 100, 1000),
      selected = 50
    ),
    shiny::sliderInput(ns("umap_gamma"), "color gamma:",
      min = 0.1, max = 1.2, value = 0.4, step = 0.1
    )
  )

  PlotModuleUI(
    ns("gene_map"),
    title = title,
    label = "a",
    plotlib = c("plotly", "ggplot"),
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    options = plot.opts,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    cards = TRUE,
    card_names = c("dynamic", "static")
  )
}

featuremap_table_gene_map_ui <- function(
    id,
    label = "",
    title,
    caption,
    info.text,
    height,
    width) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("gene_table"),
    info.text = info.text,
    height = height,
    width = width,
    caption = caption,
    title = title,
    label = "c"
  )
}

featuremap_plot_gene_map_server <- function(id,
                                            pgx,
                                            plotUMAP,
                                            sigvar,
                                            filteredProbes,
                                            watermark = FALSE,
                                            labeltype) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    getUMAP <- function() {
      pos <- pgx$cluster.genes$pos[["umap2d"]]
      colnames(pos) <- c("x", "y")
      pos
    }

    plot_data <- shiny::reactive({
      pos <- getUMAP()
      hilight <- filteredProbes()
      nlabel <- as.integer(input$umap_nlabel)

      ## select on table filter
      F <- playbase::pgx.getMetaMatrix(pgx)$fc
      F <- scale(F, center = FALSE)
      fc <- sqrt(rowMeans(F**2, na.rm = TRUE))
      names(fc) <- rownames(F)

      gg <- intersect(rownames(pos), names(fc))
      pos <- pos[gg, , drop = FALSE]
      fc <- fc[gg]
      F <- F[gg, , drop = FALSE]

      hilight.probes <- playbase::map_probes(pgx$genes, hilight)
      ##labels <- playbase::probe2symbol(rownames(pos), pgx$genes, labeltype(), fill_na = TRUE)
      labels <- playbase::probe2symbol(rownames(pos), pgx$genes, "gene_name", fill_na = TRUE)

      pd <- list(
        df = data.frame(pos, fc = fc),
        pos = pos,
        fc = fc,
        F = F,
        hilight = hilight.probes,
        labels = labels,
        nlabel = nlabel
      )
      return(pd)
    })

    render_geneUMAP <- function(cex = 1, cex.label = 1, plotlib = "plotly") {
      pd <- plot_data()
      pos <- pd$pos
      fc <- pd$fc

      hilight <- pd$hilight
      labels <- pd$labels
      nlabel <- pd$nlabel
      colgamma <- as.numeric(input$umap_gamma)

      ## dim non hilighted genes
      fc <- sign(fc) * abs(fc / max(abs(fc), na.rm = TRUE))**colgamma
      if (length(setdiff(names(fc), hilight))) {
        fc[!names(fc) %in% hilight] <- NA
      }

      if (plotlib == "plotly") {
        p <- plotUMAP(
          pos,
          fc,
          hilight,
          labels = labels,
          nlabel = nlabel,
          title = "rms(FC)",
          cex = cex,
          cex.label = cex.label,
          xlab = "UMAP-x",
          ylab = "UMAP-y",
          plotlib = "plotly",
          source = ns("gene_umap")
        ) %>%
          plotly::layout(
            dragmode = "select",
            margin = list(l = 5, r = 5, b = 5, t = 20)
          )
        p
      } else if (plotlib == "ggplot") {
        p <- plotUMAP(
          pos,
          fc,
          hilight,
          labels = labels,
          nlabel = nlabel,
          xlab = "UMAP-x",
          ylab = "UMAP-y",
          theme = ggplot2::theme_minimal(),
          cex = cex * 0.6,
          cex.label = cex.label * 0.6,
          cex.axis = cex.label * 0.6,
          bgcolor = "white",
          gridcolor = "grey90",
          plotlib = "ggplot"
        )
        p
      }
    }

    geneUMAP.RENDER <- function() {
      p <- render_geneUMAP(cex = 1, cex.label = 1) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_default()
      p
    }

    geneUMAP.RENDER_ggplot <- function() {
      p <- render_geneUMAP(cex = 1, cex.label = 1, plotlib = "ggplot")
      p
    }

    geneUMAP.RENDER2_ggplot <- function() {
      p <- render_geneUMAP(cex = 1.5, cex.label = 1.5, plotlib = "ggplot")
      p
    }

    geneUMAP.RENDER2 <- function() {
      p <- render_geneUMAP(cex = 1.2, cex.label = 1.5) %>%
        plotly::config(
          modeBarButtons = list(list("toImage", "zoom2d", "select2d", "resetScale2d"))
        ) %>%
        plotly_modal_default()
      p
    }

    plot_grid <- list(
      list(plotlib = "plotly", func = geneUMAP.RENDER, func2 = geneUMAP.RENDER2, card = 1),
      list(plotlib = "ggplot", func = geneUMAP.RENDER_ggplot, func2 = geneUMAP.RENDER2_ggplot, card = 2)
    )

    lapply(plot_grid, function(x) {
      PlotModuleServer(
        "gene_map",
        plotlib = x$plotlib,
        func = x$func,
        func2 = x$func2,
        csvFunc = plot_data,
        res = c(80, 95), # resolution of plots
        pdf.width = 10,
        pdf.height = 8,
        add.watermark = watermark,
        card = x$card
      )
    })

    # PlotModuleServer(
    #   "gene_map",
    #   plotlib = "ggplot",
    #   plotlib2 = "plotly",
    #   func = geneUMAP.RENDER_ggplot,
    #   func2 = geneUMAP.RENDER2,
    #   csvFunc = plot_data,
    #   pdf.width = 5, pdf.height = 5,
    #   add.watermark = watermark
    # )

    # ================================================================================
    # ============================= Table server module ==============================
    # ================================================================================

    geneTable.RENDER <- shiny::reactive({
      pd <- plot_data()
      shiny::req(pd, pgx$genes)

      pos <- pd$pos
      fc <- pd$fc
      F <- pd$F
      annot <- pgx$genes

      ## Retrieve gene table with rownames (symbols)
      annot_cols <- c("feature", "symbol", "human_ortholog", "gene_title")
      annot_cols <- intersect(annot_cols, colnames(annot))
      rowids <- match(rownames(F), rownames(annot))
      annot <- annot[rowids, annot_cols, drop = FALSE]
      annot <- apply(annot, MARGIN = 2, playbase::shortstring, n = 60)
      annot <- as.data.frame(annot)

      F <- cbind(rms.FC = fc, F)
      F <- round(F, digits = 3)
      df <- data.frame(annot, F, check.names = FALSE)
      df <- df[order(F[, 1], decreasing = TRUE), ]

      # Remove feature column if it is the same as symbol column
      if (mean(df$feature %in% df$symbol, na.rm = TRUE) > 0.9) {
        df$feature <- NULL
      }
      if (mean(df$symbol == df$human_ortholog, na.rm = TRUE) > 0.9 || all(is.na(df$human_ortholog))) {
        df$human_ortholog <- NULL
      }

      ## detect brush
      sel.genes <- NULL
      b <- plotly::event_data("plotly_selected", source = ns("gene_umap"))

      if (!is.null(b) & length(b)) {
        sel <- b$key
        sel.genes <- rownames(df)[rownames(df) %in% sel]
      } else {
        sel.genes <- rownames(df)
      }
      df <- df[sel.genes, , drop = FALSE]

      DT::datatable(df,
        rownames = FALSE,
        class = "compact cell-border stripe hover",
        extensions = c("Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = 240,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    geneTable.RENDER_modal <- shiny::reactive({
      dt <- geneTable.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "gene_table",
      func = geneTable.RENDER,
      func2 = geneTable.RENDER_modal,
      selector = "none"
    )
  }) ## moduleServer
}
