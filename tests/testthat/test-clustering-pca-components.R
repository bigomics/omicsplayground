## test-clustering-pca-components.R
##
## Unit tests for pgx.pcaComponents(), the on-the-fly PCA behind the X axis /
## Y axis selectors on the Clustering board. pgx$cluster$pos only stores
## PC1-PC3, so PC4/PC5 come from here.
##
## Sources the board server standalone; only function definitions live at the
## top level of that file, so no Shiny session is needed.

.board_dir <- if (dir.exists("components/board.clustering/R")) {
  "components/board.clustering/R"
} else {
  "../../components/board.clustering/R"
}

source(file.path(.board_dir, "clustering_server.R"), local = TRUE)

## features x samples, two planted groups separated along a strong axis
make_X <- function(nsamples = 12, nfeatures = 200, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(nfeatures * nsamples), nfeatures, nsamples)
  grp <- rep(c(0, 1), length.out = nsamples)
  X[1:20, ] <- X[1:20, ] + rep(grp * 8, each = 20)
  dimnames(X) <- list(paste0("gene", 1:nfeatures), paste0("s", 1:nsamples))
  list(X = X, grp = grp)
}

test_that("returns five components with sane variance percentages", {
  X <- make_X()$X
  pc <- pgx.pcaComponents(X, npc = 5)

  expect_equal(dim(pc$pos), c(ncol(X), 5L))
  expect_equal(colnames(pc$pos), paste0("PC", 1:5))
  expect_equal(rownames(pc$pos), colnames(X))

  expect_length(pc$varexp, 5)
  expect_true(all(diff(pc$varexp) <= 0)) ## non-increasing
  expect_true(all(pc$varexp > 0))
  expect_lte(sum(pc$varexp), 100)
})

test_that("PC1 separates the planted groups", {
  d <- make_X()
  pc <- pgx.pcaComponents(d$X, npc = 5)
  expect_gt(abs(cor(pc$pos[, "PC1"], d$grp)), 0.95)
})

test_that("components match prcomp up to sign", {
  X <- make_X()$X
  pc <- pgx.pcaComponents(X, npc = 5)
  ## replicate the same preprocessing: row-centered, no SD reduction at 200 rows
  ref <- prcomp(t(X - rowMeans(X)), center = FALSE)
  for (i in 1:3) {
    expect_gt(abs(cor(pc$pos[, i], ref$x[, i])), 0.999)
  }
})

test_that("sign is pinned so the plot does not mirror between runs", {
  X <- make_X()$X
  expect_equal(
    pgx.pcaComponents(X, npc = 5)$pos,
    pgx.pcaComponents(X, npc = 5)$pos,
    tolerance = 1e-4
  )
})

test_that("number of components is capped by the sample count", {
  X <- make_X(nsamples = 4)$X
  pc <- pgx.pcaComponents(X, npc = 5)
  expect_equal(ncol(pc$pos), 3L) ## ncol(X) - 1
  expect_length(pc$varexp, 3)

  ## small datasets take the exact svd path rather than irlba
  X6 <- make_X(nsamples = 6)$X
  expect_equal(ncol(pgx.pcaComponents(X6, npc = 5)$pos), 5L)
})

test_that("missing values are imputed rather than propagated", {
  X <- make_X()$X
  X[5, 3] <- NA
  X[100, 7] <- NA
  pc <- pgx.pcaComponents(X, npc = 5)
  expect_false(any(is.na(pc$pos)))
  expect_false(any(is.na(pc$varexp)))
})

test_that("reduce.sd keeps the most variable features", {
  X <- make_X(nfeatures = 300)$X
  ## the 20 planted rows carry the signal and survive the SD cut
  pc <- pgx.pcaComponents(X, npc = 5, reduce.sd = 50)
  expect_equal(dim(pc$pos), c(ncol(X), 5L))
  expect_gt(pc$varexp[1], pc$varexp[2])
})
