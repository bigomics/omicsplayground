
source("../R/gx-volcano.r")

##load("../pgx/GSE72056-melanoma-scRNA-vsclusters-s200-4k-LT.pgx")
load("../pgx/sallusto2019-th1star-TX-4k-LT.pgx")

fx <- ngs$gx.meta$meta[[1]]$meta.fx
qv <- ngs$gx.meta$meta[[1]]$meta.q
gene <- rownames(ngs$gx.meta$meta[[1]])
gx.volcanoPlot.XY( fx, qv, gene, render="x11")

names(ngs$gx.meta$meta[[1]]$fc)
i=1
for(i in 1:3) {
    fx <- ngs$gx.meta$meta[[1]]$fc[,i]
    qv <- ngs$gx.meta$meta[[1]]$q[,i]
    gene1 <- gene
    gene1 <- gsub(".*:","",gene)
    gx.volcanoPlot.XY( fx, qv, gene1, render="x11", nlab=40)
}

