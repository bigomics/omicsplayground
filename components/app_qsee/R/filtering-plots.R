##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Plotly versions of the filtering QSEE plots (replaces the old base-R
## versions). Depends on shared helpers in `normalization-plots.R`.
##
## Naming convention: qsee_filtering_plot_*_plotly.
##

## ---------------------------------------------------------------------------
## Filtering board
## ---------------------------------------------------------------------------

#' Plotly version of [qsee_filtering_plot_pca_vs_topsd()].
qsee_filtering_plot_pca_vs_topsd_plotly <- function(res, colorby_var = NULL, show_labels = TRUE) {
  X <- res$X
  Y <- res$Y

  nn <- c(nrow(X), 20 * 2**seq(0, 12, 2))
  nn <- utils::head(nn[nn <= nrow(X)], 9)
  nn <- sort(nn)
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  X <- X[order(-sdx), , drop = FALSE]

  ycol <- if (!is.null(colorby_var) && colorby_var %in% colnames(Y)) {
    Y[, colorby_var]
  } else if (ncol(Y) >= 1) {
    Y[, min(2L, ncol(Y))]
  } else {
    rep(1, ncol(X))
  }
  y <- factor(ycol)

  panels <- list()
  for (n in nn) {
    X1 <- utils::head(X, n)
    pr <- stats::prcomp(t(X1))
    panels[[paste0("N= ", n)]] <- qsee_plotly_pca_panel(
      pr$x[, 1:2, drop = FALSE],
      color = y, show_labels = show_labels
    )
  }
  qsee_plotly_grid(panels, n_cols = 3, x_title = "PC1", y_title = "PC2")
}

#' Plotly version of [qsee_filtering_plot_variance_vs_topsd()].
qsee_filtering_plot_variance_vs_topsd_plotly <- function(res, topsd = NULL,
                                                         npoints = 2000) {
  X <- res$X
  nn <- c(20 * 2**seq(0, 12, 2), nrow(X))
  nn <- nn[nn <= nrow(X)]

  U <- res$u %*% diag(res$d, length(res$d))
  U <- U[order(-rowMeans(U**2)), , drop = FALSE]
  pcx <- 100 * cumsum(rowSums(U**2)) / sum(U**2)

  ## Thin the curve for the browser but always keep the marker positions.
  keep <- sort(unique(c(
    1L, length(pcx), as.integer(nn),
    unique(round(seq(1, length(pcx), length.out = min(length(pcx), npoints))))
  )))

  p <- omicsplots::pgx.plot_multiline(
    x = keep, y = pcx[keep],
    group = rep("variance explained", length(keep)),
    xlab = "Number of topSD features",
    ylab = "Variance explained (%)",
    title = "Variance vs. topSD"
  )
  p <- plotly::add_trace(
    p,
    x = nn, y = pcx[nn], type = "scatter", mode = "markers", inherit = FALSE,
    marker = list(color = qsee_plotly_pal()@chrome$foreground, size = 4),
    name = "topSD grid", showlegend = FALSE,
    hoverinfo = "text",
    text = paste0("N = ", nn, "<br>", round(pcx[nn], 1), "%")
  )
  p <- plotly::layout(p, showlegend = FALSE)

  if (!is.null(topsd) && is.finite(as.numeric(topsd))) {
    topx <- max(1L, min(as.integer(topsd), length(pcx)))
    topy <- pcx[topx]
    p <- plotly::layout(p, shapes = list(
      list(
        type = "line", x0 = topx, x1 = topx, y0 = 0, y1 = 1, yref = "paper",
        line = list(color = "#bf616a", dash = "dash", width = 1)
      ),
      list(
        type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = topy, y1 = topy,
        line = list(color = "#bf616a", dash = "dash", width = 1)
      )
    ))
  }
  p
}

#' Plotly version of [qsee_filtering_plot_sd_histogram()].
qsee_filtering_plot_sd_histogram_plotly <- function(res, topsd = NULL) {
  sdx <- matrixStats::rowSds(res$X, na.rm = TRUE)
  sdx[is.na(sdx)] <- 0

  p <- omicsplots::pgx.plot_density(
    matrix(sdx, ncol = 1, dimnames = list(NULL, "SD")),
    n_breaks = 100, xlab = "standard deviation (SD)",
    title = "Histogram of SD", showlegend = FALSE
  )
  if (!is.null(topsd) && is.finite(as.numeric(topsd))) {
    ntop <- max(1L, min(as.integer(topsd), length(sdx)))
    sdthreshold <- sort(sdx, decreasing = TRUE)[ntop]
    p <- plotly::layout(p, shapes = list(list(
      type = "line", x0 = sdthreshold, x1 = sdthreshold,
      y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#bf616a", dash = "dash", width = 1)
    )))
  }
  p
}
