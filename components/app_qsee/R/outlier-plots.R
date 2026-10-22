##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Plotly versions of the outlier QSEE plots (replaces the old base-R
## versions). Depends on shared helpers in `normalization-plots.R`.
##
## Naming convention: qsee_outlier_plot_*_plotly.
##

## ---------------------------------------------------------------------------
## Outlier board
## ---------------------------------------------------------------------------

#' Plotly version of [qsee_outlier_plot_zscores()].
qsee_outlier_plot_zscores_plotly <- function(res, ph, z.threshold = c(3, 6, 9)) {
  Z <- res$outliers$Z
  y <- as.character(factor(res$Y[, ph]))
  ymax <- max(10, max(Z, na.rm = TRUE))
  bars <- function(v, showlegend = FALSE) {
    qsee_plotly_barplot(
      x = names(v), y = as.numeric(v), fill = y, stacked = TRUE,
      ylab = "z-score", ylim = c(0, ymax), hline = z.threshold,
      showlegend = showlegend
    )
  }
  panels <- list(
    "z.outlier (mean)" = bars(res$outliers$z.outlier, showlegend = TRUE),
    "z.outlier (geom.mean)" = bars(res$outliers$z.outlier2)
  )
  for (i in seq_len(ncol(Z))) {
    panels[[colnames(Z)[i]]] <- bars(stats::setNames(Z[, i], rownames(Z)))
  }
  qsee_plotly_grid(panels,
    n_cols = 3, y_title = "z-score",
    margin = c(0.04, 0.04, 0.08, 0.12)
  )
}

#' Plotly version of [qsee_outlier_plot_heatmap()].
#'
#' @return an iheatmapr widget.
qsee_outlier_plot_heatmap_plotly <- function(res) {
  omicsplots::pgx.plot_heatmap(
    res$heatX,
    cluster_rows = TRUE, cluster_cols = TRUE,
    scale = "row.center", show_dendro = TRUE
  )
}

#' Plotly version of [qsee_outlier_plot_pca()].
qsee_outlier_plot_pca_plotly <- function(res, ph, show_labels = TRUE) {
  qsee_plotly_pca_panel(res$pca, color = factor(res$Y[, ph]), show_labels = show_labels)
}
