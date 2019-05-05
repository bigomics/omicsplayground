library(tidyimpute)
library(na.tools)
setwd("~/Projects/Immunomics/tests/")

source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/pgx-functions.R")
source("../R/pgx-graph.R")
source("../R/ngs-fit.r")
source("../R/ngs-cook.r")

load("../pgx/sallusto2019b-tenx-s1000-4k.pgx",verbose=1)

head(ngs$counts)[,1:4]
mean(Matrix::colSums(ngs$counts))

x1 <- log2(0.1 + ngs$counts)
x1 <- edgeR::cpm(100*as.matrix(ngs$counts),prior.count=1,log=TRUE)

pos1 <- Rtsne( as.matrix(t(x1)), num_threads=8 )$Y
klr1 <- factor(ngs$samples$treatment)

par(mfrow=c(3,3), mar=c(4,4,2,2))
plot( pos1[,1], pos1[,2], col=klr1, pch=20, cex=0.8 )
title(main="uncorrected")

source("../R/pgx-correct.R")
require(sva)
BCMETHODS=c("MNN","limma","ComBat","BMC")
for(b in BCMETHODS) {
    batch <- ngs$samples$treatment
    cx <- pgx.removeBatchFactor( as.matrix(x1), batch=batch, method=b)
    ##cx <- scale(cx)
    pos3 <- Rtsne( as.matrix(t(cx)), num_threads=8 )$Y
    plot( pos3[,1], pos3[,2], col=klr1, pch=20, cex=0.8 )
    title(main=b)
}
