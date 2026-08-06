qsee_normalization_plot_boxplots <- function(res, Y, ph) {
  normX <- res$normX
  y <- factor(Y[, ph])
  cc <- 1 + as.integer(y)
  par(mfrow = c(2, 3), cex = 1, mar = c(7, 5, 3, 1), las = 2)
  for (name in names(normX)) {
    X <- normX[[name]]
    if (is.null(X) || !length(X)) {
      plot(0, type = "n", main = name); text(0, 0, "No data"); next
    }
    if (ncol(X) > 40) X <- X[, sample(ncol(X), 40)]
    boxplot(X, col = cc, main = name, ylab = "log-expression")
  }
}

qsee_normalization_plot_histograms <- function(res, Y, ph) {
  normX <- res$normX
  cc <- 1 + as.integer(factor(Y[, ph]))
  par(mfrow = c(2, 3), cex = 1, mar = c(6, 5, 3, 1), las = 2)
  for (name in names(normX)) {
    X <- normX[[name]]
    if (is.null(X) || !length(X)) {
      plot(0, type = "n", main = name); text(0, 0, "No data"); next
    }
    if (ncol(X) > 40) X <- X[, sample(ncol(X), 40)]
    playbase::gx.hist(X, main = name, breaks = 80, lines.col = cc, las = 1,
      xlab = "signal")
  }
}

qsee_normalization_plot_pca <- function(res, Y, ph) {
  pcaX <- Filter(Negate(is.null), res$pcaX)
  cc <- 1 + as.integer(factor(Y[, ph]))
  par(mfrow = c(2, 3), cex = 1, mar = c(5, 5, 3, 1))
  for (name in names(pcaX)) {
    pos <- pcaX[[name]]
    plot(pos, main = name, pch = 19, col = cc, xlab = "PC1", ylab = "PC2", cex = 2)
    text(pos, labels = rownames(pos), cex = 0.8, pos = 3)
  }
}
