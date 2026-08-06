##
##
##

bsee.plot_pca_vs_methods <- function(res, pheno.var, cex=1) {

  X <- res$X
  samples <- res$samples
  pheno <- samples[,pheno.var]
  xlist <- res$methods$xlist

  pos <- res$methods$pos
  
  par(cex=1.2)
  playbase::bc.plotResults(
    X = X,
    xlist = xlist,
    pos = pos,
    pheno = pheno,
    samples = samples,
    type = "umap",
    scores = NULL,
    nmax = 1000,
    cex = cex
  )

}

bsee.plot_heatmap_vs_methods <- function(res) {
  X0 <- res$X
  samples <- res$samples
  xlist <- res$methods$xlist
  
  par(cex=1.2)
  playbase::bc.plotResults(
    X = X0,
    xlist = xlist,
    pos = NULL,
    pheno = NULL,
    samples = samples,
    type = "heatmap",
    scores = NULL,
    nmax = 1000,
    cex = 1
  )
}

bsee.plot_covariate_correlation_heatmap <- function(res) {
  bc <- res$effects
  par(mfrow = c(1, 1), mar = c(0, 3, 0, 3))
  playbase::bc.plotCovariateHeatmap(bc)
}

bsee.plot_covariate_analysis <- function(res) {
  bc <- res$effects
  par(mfrow = c(2, 2), mar = c(5, 5, 3.5, 2), mgp = c(2.6, 1, 0))
  playbase::bc.CovariateAnalysisPlot(bc, k = 1:3, par = FALSE, col = "#273871")
}

bsee.plot_scores <- function(res) {
  X0 <- res$X
  samples <- res$samples
  pheno <- NULL
  xlist <- res$methods$xlist
  sel <- c("score", "genes", "gsets", "SNR", "pc1.ratio", "silhouette")
  sel <- intersect(sel, colnames(res$methods$scores))
  if(length(sel)==0) return(NULL)
  scores <- res$methods$scores[, sel, drop=FALSE]
  
  playbase::bc.plotResults(
    X = X0,
    xlist = xlist,
    pos = NULL,
    pheno = pheno,
    samples = samples,
    type = "scores",
    scores = scores,
    nmax = 1000,
    ncol = 3
  )
}

bsee.plot_pvca_by_phenotype <- function(res) {
  X <- res$X
  samples <- res$samples
  xlist <- res$methods$xlist

  kk <- colnames(xlist[[1]])
  X <- X[, kk]
  samples <- samples[kk, ]

  playbase::pgx.PC_correlation_grid(
    X = X, xlist = xlist, samples = samples, 
    plotby = "phenotype",
    ncol = NULL, nv = 3, stat = "F",
    expand = FALSE, collapse = TRUE,
    plot = TRUE, horiz = FALSE,
    srt = 45, main = NULL, text.cex = 1)   
}

bsee.plot_pvca_by_component <- function(res) {
  X <- res$X
  samples <- res$samples
  xlist <- res$methods$xlist

  kk <- colnames(xlist[[1]])
  X <- X[, kk]
  samples <- samples[kk, ]
  
  playbase::pgx.PC_correlation_grid(
    X = X, xlist = xlist, samples = samples, 
    plotby = "component",
    legend.pos = c(1,1),
    ncol = NULL, nv = 8, stat = "F",
    expand = FALSE, collapse = TRUE,
    plot = TRUE, horiz = FALSE,
    srt = 0, main = NULL, text.cex = 1) 

}

if(0) {
  library(shiny)
  library(bslib)

  source("~/Playground/playbase/dev/include.R",chdir=TRUE)
  source("./qsee_visibility.R")
  source("../R/bsee-compute.R")
  source("../R/bsee-plots.R")
  pgx <- playdata::GEIGER_PGX
  pgx <- playbase::pgx.load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")
  
  X <- pgx$X
  samples <- pgx$samples
  pheno.var=1;cex=1

  Y = samples[,-1]
  nv = 3; stat = "F";
  expand = FALSE; collapse = TRUE;
  plot = TRUE; horiz = FALSE;
  main = NULL; text.cex = 1

  par(cex=1.3)
  pgx.PC_correlation(X, Y, nv=5, plot=TRUE, plotby="phenotype")
  pgx.PC_correlation(X, Y, nv=10, plot=TRUE, plotby="component")

  pgx.PC_correlation(X, Y, nv=5, plot=TRUE, plotby="phenotype", horiz=TRUE)
  pgx.PC_correlation(X, Y, nv=10, plot=TRUE, plotby="component", horiz=TRUE)
  

}
