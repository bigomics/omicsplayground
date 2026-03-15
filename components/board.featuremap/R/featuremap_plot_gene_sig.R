##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_gene_sig_ui <- function(
  id,
  label = "",
  title,
  caption,
  info.text,
  info.methods,
  info.extra_link,
  height,
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("gene_sig"),
    title = title,
    label = "b",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg"),
    editor = TRUE,
    ns_parent = ns,
    plot_type = "featuremap",
    color_selection = TRUE,
    color_selection_default = TRUE
  )
}

featuremap_plot_gene_sig_server <- function(id,
                                            pgx,
                                            sigvar,
                                            ref_group,
                                            plotFeaturesPanel,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)
      pheno <- sigvar()
      ref <- ref_group()

      pos <- pgx$cluster.genes$pos[["umap2d"]]
      if (any(pheno %in% colnames(pgx$samples))) {
        y <- pgx$samples[, pheno]
        if (ref == "<average>") {
          refX <- rowMeans(pgx$X, na.rm = TRUE)
        } else {
          kk <- which(y == ref)
          refX <- rowMeans(pgx$X[, kk], na.rm = TRUE)
        }
        X <- pgx$X - refX
        y <- pgx$samples[, pheno]
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE], na.rm = TRUE)
        }))
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "gene")$fc
        kk <- intersect(pheno, colnames(F))
        F <- F[, kk, drop = FALSE]
      }
      if (nrow(F) == 0 || NCOL(F) == 0) {
        return(NULL)
      }
      ## Click-to-label data: all panels share the same UMAP positions
      gg <- intersect(rownames(pos), rownames(F))
      click_df <- data.frame(
        x = pos[gg, 1],
        y = pos[gg, 2],
        feature_name = gg
      )
      return(list(F, pos = pos, df = click_df))
    })

    renderPlots_ggplot <- function() {
      dt <- plot_data()
      F <- dt[[1]]
      pos <- dt$pos
      shiny::req(F, pos)

      ## Editor: custom colors
      low_color <- get_editor_color(input, "color_low", "#3181de")
      high_color <- get_editor_color(input, "color_high", "#f23451")
      custom_col <- c(low_color, "#f8f8f8", high_color)

      ## Editor: custom labels
      sel <- get_custom_labels(input, rownames(pos), defaults = NULL, pgx = pgx)

      ## Editor: color just selected
      color_sel <- is.null(input$color_selection) || isTRUE(input$color_selection)

      gg <- intersect(rownames(pos), rownames(F))
      pos1 <- pos[gg, ]
      F1 <- F[gg, , drop = FALSE]
      qq <- quantile(F1, probs = c(0.002, 0.998), na.rm = TRUE)
      zsym <- ifelse(min(F1, na.rm = TRUE) >= 0, FALSE, TRUE)

      ## Build long-format data for faceted plot
      long_pos <- pos1[rep(seq_len(nrow(pos1)), ncol(F1)), , drop = FALSE]
      long_var <- as.vector(F1)
      names(long_var) <- rownames(long_pos)
      long_facet <- rep(colnames(F1), each = nrow(F1))

      hilight_scatter <- if (color_sel) sel else NULL
      opacity <- ifelse(is.null(sel), 0.9, 0.4)
      if (!color_sel) opacity <- 0.9

      playbase::pgx.scatterPlotXY(
        long_pos,
        var = long_var,
        col = custom_col,
        zsym = zsym,
        zlim = qq,
        softmax = 1,
        cex = 0.8,
        cex.legend = 0.9,
        cex.lab = 1.2,
        hilight = hilight_scatter,
        hilight2 = sel,
        hilight.col = NULL,
        opacity = opacity,
        xlab = "UMAP-x",
        ylab = "UMAP-y",
        hilight.lwd = 0.5,
        hilight.cex = 1.3,
        set.par = FALSE,
        legend = FALSE,
        box = FALSE,
        theme = ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(size = 11, margin = ggplot2::margin(b = 4))
          ),
        facet = long_facet,
        plotlib = "ggplot"
      )
    }

    PlotModuleServer(
      "gene_sig",
      plotlib = "ggplot",
      func = renderPlots_ggplot,
      func2 = renderPlots_ggplot,
      csvFunc = plot_data,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 90),
      add.watermark = watermark,
      parent_session = session
    )
  })
}
