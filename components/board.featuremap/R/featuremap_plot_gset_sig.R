##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

featuremap_plot_gset_sig_ui <- function(
  id,
  label = "",
  title,
  info.text,
  info.methods,
  info.extra_link,
  caption,
  height,
  width
) {
  ns <- shiny::NS(id)

  info_text <- "<b>Geneset signature maps.</b> UMAP clustering of genes colored by relative log-expression of the phenotype group. The distance metric is covariance. Genes that are clustered nearby have high covariance."

  PlotModuleUI(
    ns("gset_sig"),
    title = title,
    label = "b",
    info.text = info_text,
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

featuremap_plot_gset_sig_server <- function(id,
                                            pgx,
                                            getGsetUMAP,
                                            sigvar,
                                            ref_group,
                                            plotFeaturesPanel,
                                            watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    plot_data <- shiny::reactive({
      shiny::req(pgx$X)

      pos <- pgx$cluster.gsets$pos[["umap2d"]]
      shiny::validate(shiny::need( ## Custom species has this empty
        !is.null(pgx$cluster.gsets),
        "Cluster genesets are not available."
      ))
      hilight <- NULL
      pheno <- "tissue"
      pheno <- sigvar()
      if (any(pheno %in% colnames(pgx$samples))) {
        y <- pgx$samples[, pheno]
        ref <- ref_group()
        if (ref == "<average>") {
          refX <- rowMeans(pgx$gsetX, na.rm = TRUE)
        } else {
          kk <- which(y == ref)
          refX <- rowMeans(pgx$gsetX[, kk], na.rm = TRUE)
        }
        X <- pgx$gsetX - refX
        F <- do.call(cbind, tapply(1:ncol(X), y, function(i) {
          rowMeans(X[, i, drop = FALSE], na.rm = TRUE)
        }))
      } else {
        F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
        kk <- intersect(pheno, colnames(F))
        F <- F[, kk, drop = FALSE]
      }
      if (nrow(F) == 0) {
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

      kk <- intersect(rownames(pos), rownames(F))
      F <- F[kk, , drop = FALSE]
      pos <- pos[kk, , drop = FALSE]

      ## Editor: custom colors
      low_color <- get_editor_color(input, "color_low", "#3181de")
      high_color <- get_editor_color(input, "color_high", "#f23451")
      custom_col <- c(low_color, "#f8f8f8", high_color)

      ## Editor: custom labels
      sel <- get_custom_labels(input, rownames(pos))

      ## Editor: color just selected
      color_sel <- is.null(input$color_selection) || isTRUE(input$color_selection)

      qq <- quantile(F, probs = c(0.002, 0.998), na.rm = TRUE)
      zsym <- ifelse(min(F, na.rm = TRUE) >= 0, FALSE, TRUE)

      ## Build long-format data for faceted plot
      long_pos <- pos[rep(seq_len(nrow(pos)), ncol(F)), , drop = FALSE]
      long_var <- as.vector(F)
      names(long_var) <- rownames(long_pos)
      long_facet <- rep(colnames(F), each = nrow(F))

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
      "gset_sig",
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
