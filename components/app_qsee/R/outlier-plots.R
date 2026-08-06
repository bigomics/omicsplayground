qsee_outlier_plot_zscores <- function(res, ph) {
  cc <- 1 + as.integer(factor(res$Y[, ph]))
  par(mfrow = c(2, 3), mar = c(8, 4.5, 2, 1), cex = 1, mgp = c(2.3, 0.8, 0))
  getFromNamespace("plotOutlierScores", "playbase")(res$outliers, par = FALSE, col = cc)
}

qsee_outlier_plot_heatmap <- function(res) {
  playbase::gx.heatmap(res$heatX, nmax = nrow(res$heatX), keysize = 1, mar = c(12, 8))
}

qsee_outlier_plot_pca <- function(res, ph) {
  cc <- 1 + as.integer(factor(res$Y[, ph]))
  plot(res$pca, col = cc, pch = 19, cex = 2.4, xlab = "PC1", ylab = "PC2")
}
