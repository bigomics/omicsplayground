# Normalization comparisons for the Qsee board.
#
# This file owns simulated noise and display-ready PCA summaries.
# Normalization runs through the canonical log2-matrix family.

qsee_normalization_add_noise <- function(rawX, amount = 1) {
  if (amount <= 0) {
    return(rawX)
  }
  zeroes <- which(rawX == 0)
  n <- ncol(rawX)
  set.seed(42 + as.integer(round(amount * 10)))
  rawX <- t(t(rawX) * rnorm(n, 1, 0.2 * amount) +
    10 * amount + rnorm(n, 0, 2 * amount))
  rawX <- pmax(rawX, 0)
  rawX[zeroes] <- NA
  rawX
}

qsee_normalization_compute <- function(rawX, progress = NULL) {
  methods <- c("CPM", "CPM+quantile", "maxMedian", "maxSum", "quantile")
  normX <- list(raw = rawX)
  for (method in methods) {
    normX[[method]] <- tryCatch(
      do.call(
        playbase.preprocess::pp.normalize,
        c(
          list(X = rawX, method = method, space = "log2"),
          if (method %in% c("CPM", "CPM+quantile")) list(prior = 0) else list()
        )
      ),
      error = function(e) {
        message("pp.normalize failed for ", method, ": ", e$message)
        rawX
      }
    )
  }

  if (!is.null(progress)) progress$set(message = "Computing PCA...", value = 0.7)
  pcaX <- lapply(normX, function(X) {
    cX <- X - rowMeans(X, na.rm = TRUE)
    sel <- which(rowMeans(is.na(cX)) == 0)
    if (length(sel) < 2 || ncol(cX) < 2) {
      return(NULL)
    }
    res <- irlba::irlba(cX[sel, , drop = FALSE], nv = 2)
    pos <- res$v
    dimnames(pos) <- list(colnames(cX), c("PC1", "PC2"))
    pos
  })

  list(normX = normX, pcaX = pcaX)
}
