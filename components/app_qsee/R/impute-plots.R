##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Plotly versions of the imputation / missing-values QSEE plots (replaces the old base-R
## versions). Depends on shared helpers in `normalization-plots.R`.
##
## Naming convention: qsee_imputation_plot_*_plotly.
##

## ---------------------------------------------------------------------------
## Missing values / imputation board
## ---------------------------------------------------------------------------

#' Plotly version of [qsee_imputation_plot_pca()].
qsee_imputation_plot_pca_plotly <- function(res, Y, ph, show_labels = TRUE) {
  pcaX <- Filter(Negate(is.null), res$pcaX)
  y <- factor(Y[, ph])
  panels <- lapply(names(pcaX), function(name) {
    qsee_plotly_pca_panel(pcaX[[name]], color = y, show_labels = show_labels)
  })
  names(panels) <- names(pcaX)
  qsee_plotly_grid(panels, n_cols = 3,
    margin = c(0.05, 0.05, 0.08, 0.10),
    x_title = "PC1", y_title = "PC2")
}

#' Plotly version of [qsee_imputation_plot_heatmaps()].
#'
#' @return named list of iheatmapr widgets, one per imputation method.
qsee_imputation_plot_heatmaps_plotly <- function(res, nmax = 200) {
  impX <- res$impX
  cor_hclust <- function(X) {
    corX <- stats::cor(X, use = "pairwise")
    corX[is.na(corX)] <- 0
    fastcluster::hclust(stats::as.dist(1 - corX), method = "ward.D2")
  }
  X <- impX[[1]]
  sel <- rowMeans(is.na(X)) < 0.50 & rowMeans(is.na(X)) > 0
  X <- X[sel, , drop = FALSE]
  X <- X[order(-matrixStats::rowSds(X, na.rm = TRUE)), , drop = FALSE]
  X <- utils::head(X, nmax)
  if (!nrow(X)) {
    return(list())
  }
  genes <- rownames(X)
  ## Row/column order is taken from the first matrix and reused for all
  ## panels, exactly like the base version.
  ii <- cor_hclust(t(X))$order
  jj <- cor_hclust(X)$order

  panels <- lapply(names(impX), function(name) {
    X1 <- impX[[name]][genes, , drop = FALSE][ii, jj, drop = FALSE]
    omicsplots::pgx.plot_heatmap(
      X1,
      cluster_rows = FALSE, cluster_cols = FALSE,
      scale = "row", show_dendro = FALSE
    )
  })
  names(panels) <- names(impX)
  panels
}

#' Plotly version of [qsee_imputation_plot_histograms()].
#'
#' Observed values vs. the values that were filled in by each method.
qsee_imputation_plot_histograms_plotly <- function(res, nmax = 1e5) {
  impX <- res$impX
  ii <- which(!is.na(impX[[1]]))
  jj <- which(is.na(impX[[1]]))
  if (length(ii) > nmax) ii <- sample(ii, nmax)
  if (length(jj) > nmax) jj <- sample(jj, nmax)
  cmap <- c(observed = "#989898", imputed = "#bf616a")

  panels <- lapply(names(impX), function(name) {
    X <- impX[[name]]
    obs <- as.numeric(X[ii])
    imp <- if (length(jj)) as.numeric(X[jj]) else numeric(0)
    imp <- imp[!is.na(imp)]
    n <- max(length(obs), length(imp))
    dat <- list(observed = c(obs, rep(NA_real_, n - length(obs))))
    if (length(imp)) dat$imputed <- c(imp, rep(NA_real_, n - length(imp)))
    p <- omicsplots::pgx.plot_density(
      dat,
      n_breaks = 80, color = cmap,
      xlab = "intensity (log2)", showlegend = TRUE
    )
    ## pgx.plot_density() builds via ggplotly(), so p$x$data is only
    ## populated once the widget is actually built -- do that before
    ## touching per-trace styling below, else the loop is a no-op.
    p <- plotly::plotly_build(p)
    ## Add filled area under each density curve
    for (i in seq_along(p$x$data)) {
      col <- p$x$data[[i]]$line$color
      rgb <- grDevices::col2rgb(col)
      p$x$data[[i]]$fill <- "tozeroy"
      p$x$data[[i]]$fillcolor <- sprintf("rgba(%d,%d,%d,0.45)", rgb[1], rgb[2], rgb[3])
    }
    p
  })
  names(panels) <- names(impX)
  qsee_plotly_grid(panels, n_cols = 3, x_title = "intensity (log2)", y_title = "Density")
}

#' Plotly version of [qsee_imputation_plot_distributions()].
qsee_imputation_plot_distributions_plotly <- function(res, Y, ph) {
  X <- res$impX[[1]]
  y <- factor(Y[, ph])

  x.avg <- rowMeans(X, na.rm = TRUE)
  x.nar <- rowMeans(is.na(X))
  x.avg2 <- cut(x.avg, breaks = 20)
  bins <- levels(x.avg2)
  ub <- suppressWarnings(as.numeric(sub(".*,", "", sub("\\]$", "", bins))))
  miss <- tapply(seq_len(nrow(X)), x.avg2, function(i) mean(is.na(X[i, , drop = FALSE])))
  miss <- as.numeric(miss[bins])
  miss[is.na(miss)] <- 0
  lab <- format(round(ub, 2), trim = TRUE)

  panels <- list()

  ## 1. stacked missing / present per average-intensity bin
  panels[["missingness vs. intensity"]] <- qsee_plotly_barplot(
    x = rep(lab, 2),
    y = c(miss, 1 - miss),
    fill = rep(c("missing", "present"), each = length(lab)),
    stacked = TRUE,
    color = c(missing = "#bf616a", present = "#d8d8d8"),
    xlab = "average intensity (log2)", ylab = "missing ratio",
    showlegend = FALSE
  )

  ## 2. missing ratio vs. average intensity, per feature
  panels[["missingness vs. intensity (features)"]] <- plotly::layout(
    omicsplots::pgx.plot_scatter(
      x = x.avg, y = x.nar,
      labels = rownames(X),
      xlab = "average intensity (log2)", ylab = "missing ratio",
      max_points = 5000
    ),
    yaxis = list(range = c(0, 1.2))
  )

  ## 3. distribution of the non-zero missing ratios
  panels[["missing ratio histogram"]] <- if (any(x.nar > 0)) {
    omicsplots::pgx.plot_density(
      matrix(x.nar[x.nar > 0], ncol = 1, dimnames = list(NULL, "missing ratio")),
      n_breaks = 40, xlab = "missing ratio", showlegend = FALSE
    )
  } else {
    qsee_plotly_empty("No missing values")
  }

  ## 4. nr. of features retained per missing-ratio threshold
  qq <- seq(0, 1, 0.05)
  qsum <- sapply(qq, function(q) sum(x.nar <= q))
  panels[["nr. features vs. threshold"]] <- omicsplots::pgx.plot_scatter(
    x = qq, y = qsum,
    labels = paste0("<= ", qq),
    xlab = "max. missing ratio threshold", ylab = "nr. features"
  )

  ## 5. missingness per sample
  panels[["missingness per sample"]] <- qsee_plotly_barplot(
    x = colnames(X), y = colMeans(is.na(X)), fill = as.character(y),
    stacked = TRUE, ylab = "missing ratio", showlegend = FALSE
  )

  qsee_plotly_grid(panels,
    n_cols = 3, axis_titles = TRUE,
    margin = c(0.05, 0.03, 0.06, 0.12)
  )
}

#' Plotly version of [qsee_imputation_plot_validation()].
qsee_imputation_plot_validation_plotly <- function(res, nmax = 1e5) {
  impX <- res$impX
  V <- res$val_set
  if (!length(V)) {
    return(qsee_plotly_empty("Please add MAR/MNAR missing values."))
  }
  jj <- V[, 1]
  actual <- V[, 2]
  if (length(jj) > nmax) {
    jj <- utils::head(jj, nmax)
    actual <- utils::head(actual, nmax)
  }
  panels <- list()
  for (i in seq_along(impX)[-1]) {
    imputed <- impX[[i]][jj]
    r <- suppressWarnings(stats::cor(imputed, actual, use = "pairwise"))
    nm <- paste0(names(impX)[i], "  (R = ", round(r, 2), ")")
    panels[[nm]] <- plotly::layout(
      omicsplots::pgx.plot_scatter(
        x = imputed, y = actual,
        xlab = "imputed value", ylab = "actual value",
        max_points = 5000
      ),
      xaxis = list(range = c(-5, 15))
    )
  }
  qsee_plotly_grid(panels,
    n_cols = 3,
    x_title = "imputed value", y_title = "actual value"
  )
}
