##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##
## Plotly versions of the BSEE (batch-effects) QSEE plots (replaces the old base-R
## versions). Depends on shared helpers in `normalization-plots.R`.
##
## Naming convention: bsee.*_plotly or bsee_plot_*_plotly.
##

## ---------------------------------------------------------------------------
## Batch-effects (bsee) board
## ---------------------------------------------------------------------------

#' Plotly version of [bsee.plot_pca_vs_methods()].
bsee.plot_pca_vs_methods_plotly <- function(res, pheno.var, show_labels = TRUE) {
  samples <- res$samples
  pheno <- if (length(pheno.var) == 1 && pheno.var %in% colnames(samples)) {
    samples[, pheno.var]
  } else {
    pheno.var
  }
  xlist <- res$methods$xlist
  pos <- res$methods$pos

  methods <- names(xlist)
  if (!is.null(pos)) methods <- intersect(methods, names(pos))
  if (!length(methods)) {
    return(qsee_plotly_empty("No clustering positions"))
  }

  y <- factor(pheno)
  panels <- lapply(methods, function(m) {
    xy <- pos[[m]]
    ax <- colnames(xy)
    if (is.null(ax)) ax <- c("dim1", "dim2")
    qsee_plotly_pca_panel(xy, color = y, xlab = ax[1], ylab = ax[2], show_labels = show_labels)
  })
  names(panels) <- methods
  qsee_plotly_grid(panels, n_cols = 3)
}

#' Plotly version of [bsee.plot_heatmap_vs_methods()].
#'
#' @return named list of iheatmapr widgets, one per correction method.
bsee.plot_heatmap_vs_methods_plotly <- function(res, nmax = 400) {
  xlist <- res$methods$xlist
  panels <- lapply(names(xlist), function(m) {
    xx <- xlist[[m]]
    xx <- utils::head(xx[order(-apply(xx, 1, stats::sd, na.rm = TRUE)), , drop = FALSE], nmax)
    xx <- xx - rowMeans(xx, na.rm = TRUE)
    xx <- abs(xx)**0.5 * sign(xx)
    qsee_heatmap_hide_legend_title(omicsplots::pgx.plot_heatmap(
      xx,
      cluster_rows = TRUE, cluster_cols = TRUE,
      scale = "none", show_dendro = FALSE
    ))
  })
  names(panels) <- names(xlist)
  panels
}

#' Plotly version of [bsee.plot_covariate_correlation_heatmap()].
bsee.plot_covariate_correlation_heatmap_plotly <- function(res) {
  B <- res$effects$covariates_plus
  if (is.null(B) || NCOL(B) < 2) {
    return(qsee_plotly_empty("No covariates"))
  }
  rho <- stats::cor(apply(B, 2, rank))
  omicsplots::pgx.plot_corrmat(rho, cluster = TRUE, zlim = c(-1, 1))
}

#' Quadrant guide reproducing the legend panel of bc.CovariateAnalysisPlot().
#' @noRd
qsee_plotly_quadrant_legend <- function() {
  labs <- data.frame(
    x = c(0.2, 0.75, 0.2, 0.75),
    y = c(0.80, 0.80, 0.20, 0.20),
    t = c(
      "strong batch-effect", "strong model-parameter",
      "nuisance parameter", "weak model-parameter"
    ),
    stringsAsFactors = FALSE
  )
  fg <- qsee_plotly_pal()@chrome$foreground
  ln <- list(color = fg, dash = "dash", width = 1)
  p <- plotly::plot_ly()
  p <- plotly::add_trace(p,
    x = c(0.5, 0.5), y = c(0, 1), type = "scatter",
    mode = "lines", line = ln, showlegend = FALSE, hoverinfo = "skip"
  )
  p <- plotly::add_trace(p,
    x = c(0, 1), y = c(0.5, 0.5), type = "scatter",
    mode = "lines", line = ln, showlegend = FALSE, hoverinfo = "skip"
  )
  p <- plotly::add_trace(
    p,
    x = labs$x, y = labs$y, type = "scatter", mode = "text",
    text = labs$t, textposition = "middle center",
    textfont = list(size = 12, color = fg),
    showlegend = FALSE, hoverinfo = "skip"
  )
  p <- plotly::layout(
    p,
    xaxis = list(
      range = c(0, 1), showticklabels = FALSE, zeroline = FALSE,
      title = list(text = "correlation with phenotype")
    ),
    yaxis = list(
      range = c(0, 1), showticklabels = FALSE, zeroline = FALSE,
      title = list(text = "correlation with PC")
    )
  )
  omicsplots::plotly_theme(p)
}

#' Plotly version of [bsee.plot_covariate_analysis()].
bsee.plot_covariate_analysis_plotly <- function(res, k = 1:3) {
  bc <- res$effects
  pp <- intersect(rownames(bc$p.pca), rownames(bc$p.values))
  pp <- grep("^pca", pp, value = TRUE, invert = TRUE)
  if (!length(pp)) {
    return(qsee_plotly_empty("No covariates"))
  }
  k <- k[which(k <= ncol(bc$p.pca))]
  pxx <- bc$p.pca[pp, k, drop = FALSE]
  py <- bc$p.values[pp, 2]
  x1 <- -log10(1e-04 + py)

  panels <- list("quadrants" = qsee_plotly_quadrant_legend())
  for (i in k) {
    y1 <- -log10(1e-04 + pxx[, i])
    pnl <- omicsplots::pgx.plot_scatter(
      x = x1, y = y1, labels = pp,
      xlab = "significance with phenotype (-log10p)",
      ylab = "significance with PC (-log10p)"
    )
    pnl <- qsee_plotly_add_labels(pnl, x1, y1, pp, nmax = 40)
    panels[[paste0("PC", i)]] <- plotly::layout(
      pnl,
      xaxis = list(range = c(-0.6, 4.6))
    )
  }
  qsee_plotly_grid(panels,
    n_cols = 2,
    x_title = "significance with phenotype (-log10p)",
    y_title = "significance with PC (-log10p)"
  )
}

#' Plotly version of [bsee.plot_scores()].
bsee.plot_scores_plotly <- function(res) {
  sel <- c("score", "genes", "gsets", "SNR", "pc1.ratio", "silhouette")
  scores <- res$methods$scores
  sel <- intersect(sel, colnames(scores))
  if (!length(sel)) {
    return(NULL)
  }
  scores <- scores[, sel, drop = FALSE]

  ## same ordering as bc.plotResults(type = "scores")
  if ("score" %in% colnames(scores)) {
    scores <- scores[order(-scores[, "score"]), , drop = FALSE]
  }

  ylabs <- c(
    "score" = "overall score",
    "genes" = "significant genes",
    "gsets" = "significant genesets",
    "avg.fc" = "average abs.logFC",
    "avg.sd" = "average SD",
    "r.genes" = "gene.coverage",
    "r.gsets" = "gset coverage",
    "SNR" = "signal-to-noise",
    "pc1.ratio" = "PC1 ratio",
    "silhouette" = "silhouette score"
  )

  panels <- lapply(colnames(scores), function(nn) {
    qsee_plotly_barplot(
      x = rownames(scores), y = as.numeric(scores[, nn]),
      ylab = if (nn %in% names(ylabs)) ylabs[[nn]] else nn
    )
  })
  names(panels) <- colnames(scores)
  ## axis_titles = TRUE keeps each panel's ylab (the human-readable name of
  ## the metric), as in the base version.
  qsee_plotly_grid(panels,
    n_cols = 3, axis_titles = TRUE,
    margin = c(0.05, 0.03, 0.06, 0.10)
  )
}

#' Shared worker for the two PVCA plots.
#'
#' Mirrors [playbase::pgx.PC_correlation_grid()] but takes the (already
#' normalized) F-statistic matrix out of `pgx.PC_correlation(plot = FALSE)`
#' and renders it with plotly instead of ggplot.
#' @noRd
bsee_plotly_pvca <- function(res, plotby = c("phenotype", "component"), nv = 3) {
  plotby <- match.arg(plotby)
  X <- res$X
  samples <- res$samples
  xlist <- res$methods$xlist

  kk <- colnames(xlist[[1]])
  X <- X[, kk, drop = FALSE]
  samples <- samples[kk, , drop = FALSE]

  B <- playbase::pgx.computeTechnicalEffects(X, nv = 1)
  bcat <- sub("[.].*", "", colnames(B))
  colnames(B) <- paste0(bcat, ":", colnames(B)) ## for collapsing
  B <- cbind(B, samples)

  panels <- lapply(names(xlist), function(m) {
    cc <- playbase::pgx.PC_correlation(
      X = xlist[[m]], Y = B, nv = nv, stat = "F",
      expand = FALSE, collapse = TRUE, plot = FALSE
    )
    R <- cc$R
    if (is.null(R) || !length(R)) {
      return(qsee_plotly_empty("No result"))
    }
    if (plotby == "phenotype") {
      ## x = covariate, stacked by principal component
      qsee_plotly_barplot(
        x = rep(rownames(R), times = ncol(R)),
        y = as.numeric(R),
        fill = rep(colnames(R), each = nrow(R)),
        stacked = TRUE, ylab = "F-statistic"
      )
    } else {
      ## x = principal component, stacked by covariate
      qsee_plotly_barplot(
        x = rep(colnames(R), each = nrow(R)),
        y = as.numeric(R),
        fill = rep(rownames(R), times = ncol(R)),
        stacked = TRUE, ylab = "F-statistic"
      )
    }
  })
  names(panels) <- names(xlist)
  qsee_plotly_grid(panels,
    n_cols = 3, y_title = "F-statistic",
    margin = c(0.04, 0.04, 0.05, 0.12)
  )
}

#' Plotly version of [bsee.plot_pvca_by_phenotype()].
bsee.plot_pvca_by_phenotype_plotly <- function(res) {
  bsee_plotly_pvca(res, plotby = "phenotype", nv = 3)
}

#' Plotly version of [bsee.plot_pvca_by_component()].
bsee.plot_pvca_by_component_plotly <- function(res) {
  bsee_plotly_pvca(res, plotby = "component", nv = 8)
}
