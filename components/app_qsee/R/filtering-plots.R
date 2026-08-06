qsee_filtering_plot_pca_vs_topsd <- function(res, colorby_var) {
  pcaexplorer.plot_pca_vs_topsd(res, colorby_var = colorby_var)
}

qsee_filtering_plot_variance_vs_topsd <- function(res, topsd) {
  pcaexplorer.plot_variance_vs_topsd(res, topsd = topsd)
}

qsee_filtering_plot_sd_histogram <- function(res, topsd) {
  pcaexplorer.plot_sd_histogram(res, topsd = topsd)
}
