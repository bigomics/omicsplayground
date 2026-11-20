## test-clustering-wgcna-topmode.R
##
## Unit tests for the WGCNA top mode of the Clustering board heatmap:
## availability gating, feature selection against the filtered matrix, and
## the PCA component cap that shares the "K modules" setting.
##
## Sources the board server standalone: only function definitions live at the
## top level of that file, so no Shiny session is needed.

.board_dir <- if (dir.exists("components/board.clustering/R")) {
  "components/board.clustering/R"
} else {
  "../../components/board.clustering/R"
}

source(file.path(.board_dir, "clustering_server.R"), local = TRUE)
## clusteringTopModeChoices() lives in the board UI module
source(file.path(.board_dir, "clustering_ui.R"), local = TRUE)

## playbase, matrixStats, irlba and fastcluster are dependencies of the app,
## so they are available whenever the app itself is.

## playbase::wgcna.getTopMEfeatures() is the feature getter pgx.wgcnaTopFeatures()
## defaults to. It lives in playbase devel (the branch the app is built
## against, see docker/Dockerfile) and is not exported by the playbase
## release some development checkouts still have installed, so the feature
## selection is exercised with a stand-in: same signature and algorithm
## (top-n features of each of the k largest modules, grey excluded, ordered
## by module membership; a missing membership column ranks as NA), always
## used here so the test is deterministic.
stub_getTopMEfeatures <- function(wgcna, n = 10, k = NULL, sortby = "mm") {
  colors <- paste0("ME", wgcna$net$colors)
  names(colors) <- names(wgcna$net$colors)
  modules <- setdiff(
    names(sort(table(colors), decreasing = TRUE)),
    c("grey", "MEgrey")
  )
  if (!is.null(k)) modules <- head(modules, k)
  membership <- wgcna$stats$moduleMembership
  top <- list()
  for (module in modules) {
    genes <- names(colors)[colors == module]
    mm <- if (module %in% colnames(membership)) membership[genes, module] else NA
    top[[module]] <- head(genes[order(-mm)], n)
  }
  top
}

## pgx.wgcnaTopFeatures() with the stand-in getter; same call shape as the
## board's call, everything else stays as shipped
top_features <- function(wgcna, features, ntop, k) {
  pgx.wgcnaTopFeatures(
    wgcna, features,
    topfeatures = stub_getTopMEfeatures,
    ntop = ntop, k = k
  )
}

## ---------------------------------------------------------------------------
## fixtures
## ---------------------------------------------------------------------------

## Minimal pgx$wgcna: three colour modules plus grey, module membership
## strongly separating the modules.
make_wgcna <- function(ngenes = 30) {
  set.seed(42)
  colors <- rep(c("blue", "brown", "green", "grey"), length.out = ngenes)
  names(colors) <- paste0("g", seq_len(ngenes))
  mm <- matrix(stats::runif(ngenes, 0.4, 0.95), ngenes, 1,
    dimnames = list(names(colors), "MEblue")
  )
  mm[names(colors)[colors == "blue"], "MEblue"] <- stats::runif(
    sum(colors == "blue"), 0.7, 0.95
  )
  mm[names(colors)[colors == "brown"], "MEblue"] <- stats::runif(
    sum(colors == "brown"), 0.4, 0.6
  )
  ## the other modules get their own membership column
  mm <- cbind(mm, matrix(NA_real_, ngenes, 2,
    dimnames = list(names(colors), c("MEbrown", "MEgreen"))
  ))
  for (module in c("MEbrown", "MEgreen")) {
    in.module <- names(colors)[paste0("ME", colors) == module]
    mm[in.module, module] <- stats::runif(length(in.module), 0.7, 0.95)
    mm[setdiff(rownames(mm), in.module), module] <- stats::runif(
      length(colors) - length(in.module), 0.05, 0.4
    )
  }
  list(
    net = list(colors = colors),
    stats = list(moduleMembership = mm),
    datExpr = matrix(
      stats::rnorm(length(colors) * 5), length(colors), 5,
      dimnames = list(names(colors), paste0("s", 1:5))
    )
  )
}

make_X <- function(nfeatures = 30, nsamples = 12, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(nfeatures * nsamples), nfeatures, nsamples)
  dimnames(X) <- list(paste0("g", seq_len(nfeatures)), paste0("s", seq_len(nsamples)))
  X
}

## ---------------------------------------------------------------------------
## availability gate
## ---------------------------------------------------------------------------

test_that("WGCNA top mode is available only at gene level with modules", {
  wgcna <- make_wgcna()
  expect_true(pgx.wgcnaTopModeAvailable("gene", wgcna))
  expect_false(pgx.wgcnaTopModeAvailable("geneset", wgcna))
  expect_false(pgx.wgcnaTopModeAvailable("gene", NULL))
  expect_false(pgx.wgcnaTopModeAvailable("gene", list()))
  expect_false(pgx.wgcnaTopModeAvailable("gene", list(net = list())))
})

test_that("enhanced (hubgenes) results need gene statistics", {
  wgcna <- make_wgcna()
  no_stats <- wgcna
  no_stats$stats <- NULL
  expect_true(pgx.wgcnaTopModeAvailable("gene", wgcna, hubgenes = FALSE))
  expect_false(pgx.wgcnaTopModeAvailable("gene", no_stats, hubgenes = TRUE))
})

## ---------------------------------------------------------------------------
## feature selection
## ---------------------------------------------------------------------------

test_that("selected WGCNA features are restricted to the filtered matrix", {
  wgcna <- make_wgcna()
  all_genes <- names(wgcna$net$colors)
  retained <- setdiff(all_genes, c("g1", "g2", "g3", "g4", "g5", "g6"))

  ff <- top_features(wgcna, retained, ntop = 5, k = 3)

  expect_type(ff, "character")
  expect_true(all(ff %in% retained))
  expect_gt(length(ff), 1)
  expect_equal(anyDuplicated(ff), 0L)
  ## every surviving feature has a module colour, which the heatmap labels
  expect_false(anyNA(wgcna$net$colors[ff]))
})

test_that("no surviving features returns NULL (partial, empty, single)", {
  wgcna <- make_wgcna()
  expect_null(top_features(wgcna, character(0), ntop = 5, k = 3))

  ## only one WGCNA feature left in the filtered matrix; take the getter's
  ## own first pick so the choice does not depend on membership ties
  first_pick <- unlist(stub_getTopMEfeatures(wgcna, n = 5, k = 3),
    use.names = FALSE
  )[1]
  expect_length(
    top_features(wgcna, first_pick, ntop = 5, k = 3),
    1L
  )

  ## none of the filtered features are network genes
  expect_null(top_features(wgcna, c("other1", "other2"), ntop = 5, k = 3))
})

test_that("the heatmap labels every selected feature by its module", {
  wgcna <- make_wgcna()
  X <- make_X(nfeatures = length(wgcna$net$colors))
  ff <- top_features(wgcna, rownames(X), ntop = 5, k = 3)
  idx <- paste0("ME", wgcna$net$colors[ff])
  idx[is.na(idx)] <- "MEgrey"

  expect_length(idx, length(ff))
  expect_false(anyNA(idx))
  expect_equal(unname(idx), paste0("ME", wgcna$net$colors[ff]))
})

## ---------------------------------------------------------------------------
## PCA component cap (shares "K modules")
## ---------------------------------------------------------------------------

## request NPCA the same way getTopMatrix() does and try the irlba call
pca_request <- function(nfeatures, nsamples, k) {
  X <- make_X(nfeatures = nfeatures, nsamples = nsamples)
  clustk <- as.integer(k)
  npca <- min(max(5, clustk), ncol(X) - 1, nrow(X) - 1)
  npca <- max(1, npca)
  suppressWarnings(irlba::irlba(X - rowMeans(X), nv = npca))
  npca
}

test_that("PCA never requests as many components as the smaller dimension", {
  cases <- list(
    c(6, 13, 6), ## six filtered features, K = 6: the regression
    c(6, 13, 4),
    c(30, 8, 6), ## fewer samples than features
    c(30, 6, 4),
    c(30, 12, 1),
    c(200, 12, 6)
  )
  npca_request <- function(case) {
    npca <- min(max(5, as.integer(case[3])), case[2] - 1, case[1] - 1)
    max(1, npca)
  }
  for (case in cases) {
    ## the request itself must be accepted by irlba (regression: K = 6 with
    ## six features asked for nv = 6 and threw)
    expect_no_error(pca_request(case[1], case[2], case[3]))
    npca <- npca_request(case)
    expect_lt(npca, min(case[1], case[2]))
    expect_gte(npca, 1)
  }
  ## the historical five-component behaviour is kept when K does not need more
  expect_equal(pca_request(200, 12, 4), 5)
  expect_equal(pca_request(200, 12, 6), 6)
})

## ---------------------------------------------------------------------------
## short labels of the Top mode selector
## ---------------------------------------------------------------------------

## Regression (issue #1934): the options sit in a narrow sidebar and used to
## wrap as "sd (high variance)", "pca (top components)", ...; the criteria are
## described in the selector tooltip instead.
test_that("Top mode choices use the short names as labels", {
  for (modes in list(c("sd", "pca", "marker"), c("sd", "pca", "marker", "wgcna"))) {
    choices <- clusteringTopModeChoices(modes)
    expect_named(choices, modes)
    expect_equal(unname(choices), modes)
  }
})
