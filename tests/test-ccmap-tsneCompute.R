setwd("~/Projects/Immunomics/tests/")
source("../R/pgx-drugs.R")
source("../R/pgx-functions.R")

FILES="../files"
X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
x.drugs <- gsub("_.*$","",colnames(X))
length(table(x.drugs))
dim(X)

avgX <- readRDS(file=file.path(FILES,"l1000_es_5685drugsAVG.rds"))
## jj <- c(grep("torin",colnames(avgX)),1:100)
## avgX <- avgX[,jj]
dim(avgX)

jj <- head(order(-apply(X,2,sd)),100)
##out <- pgx.computeDrugEnrichment(obj=X[,jj], X, x.drugs, methods="cor", contrast=NULL)
out <- pgx.computeDrugEnrichment(
    obj=avgX[,], X, x.drugs, methods="cor", nprune=-1, contrast=NULL)
names(out)
out1 <- out[["cor"]]
##save(out1, file="l1000-dxrho-avgX.rda")
dim(out1$X)

##qx <- out1$X
qx <- out1$X * (1 - out1$P)
dim(qx)

library(Rtsne)
pos <- Rtsne( (qx), is_distance=FALSE, perplexity=30)$Y
dim(pos)
rownames(pos) <- rownames(qx)
write.csv(pos, file="l1000-drugs-tsne-px30T.csv")

out1$pos <- pos
save(out1, file="l1000-dxrhoT-avgX.rda")

if(0) {
    require(scatterD3)
    pos <- read.csv(file="l1000-drugs-tsne-px30T.csv", row.names=1)
    scatterD3( pos[,1], pos[,2], lab=rownames(pos))


}
