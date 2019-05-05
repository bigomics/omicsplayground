setwd("~/Projects/Immunomics/tests/")
source("../R/pgx-drugs.R")
source("../R/pgx-functions.R")


##======================================================================
## Compute the drug enrichment of all drugs (for t-SNE)
##======================================================================

FILES="../files"
X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
x.drugs <- gsub("_.*$","",colnames(X))
length(table(x.drugs))

avgX <- readRDS(file=file.path(FILES,"l1000_es_5685drugsAVG.rds"))
## jj <- c(grep("torin",colnames(avgX)),1:100)
## avgX <- avgX[,jj]
dim(avgX)

out <- pgx.computeDrugEnrichment(
    obj=avgX[,], X, x.drugs, methods="cor", contrast=NULL)
names(out)

dxrho <- out[["cor"]]
##save(dxrho, file="l1000-dxrho-avgX.rda")

qx <- dxrho$X
wqx <- dxrho$X * (1 - dxrho$Q)
sqx <- scale(qx) / sqrt(nrow(qx)-1)

rho.sqx <- t(sqx) %*% sqx
xdist <- (1 - rho.sqx)**2
dim(xdist)

library(Rtsne)
pos <- Rtsne( t(sqx[,]), is_distance=FALSE, perplexity=30)$Y
pos <- Rtsne( xdist[,], is_distance=TRUE, perplexity=30)$Y
rownames(pos) <- colnames(qx)
dim(pos)
write.csv(pos, file="l1000-drugs-tsne-px30d2.csv")

dxrho$pos <- pos
save(dxrho, file="l1000-dxrho-avgX.rda")
load(file="l1000-dxrho-avgX.rda")


##======================================================================
## Add clustering information (mono)
##======================================================================

require(scran)
require(igraph)

descr0 <- read.csv(file.path(FILES,"L1000_repurposing_drugs.txt"),
                   sep="\t", comment.char="#")

tsne.files <- dir(".",pattern="l1000-drugs-tsne.*.csv")
f=tsne.files[1]
for(f in tsne.files) {

    cat(">>> reading",f,"...\n")
    pos <- read.csv(file=f, row.names=1)[,1:2]
    snn <- buildSNNGraph(t(pos), k=10)
    clust <- cluster_louvain(snn)
    table(clust$membership)
    colnames(pos) <- c("tsne.x","tsne.y")
    pos$membership <- clust$membership

    descr <- descr0[match(rownames(pos),descr0$pert_iname),]
    pos <- cbind(pos, descr)
    head(pos)
    
    write.csv(pos, file=f)
}


require(scatterD3)

clust <- paste0("cluster.",pos$membership)
clust.centers <- apply(pos[,1:2],2,function(x) tapply(x,clust,median))
clust.centers
clust.moa <- tapply(pos$moa,clust,function(x) paste(unique(x),collapse="|"))
clust.moa.top3 <- tapply(as.character(pos$moa),clust,function(x) {
    gg.count <- sort(table(unlist(strsplit(x,split="\\|"))),decreasing=TRUE)
    paste(names(head(gg.count,3)),collapse="|")
})
clust.target.top10 <- tapply(as.character(pos$target),clust,function(x) {
    gg.count <- sort(table(unlist(strsplit(x,split="\\|"))),decreasing=TRUE)
    paste(names(head(gg.count,10)),collapse="|")
})

posx <- rbind(pos[,1:2], clust.centers)
clustx <- c(clust, rownames(clust.centers))
lab <- rownames(posx)
jj <- match(names(clust.moa.top3),lab)
lab[jj] <- toupper(clust.moa.top3)

##lab <- gsub("\\|","\n",lab)
tooltip_text <- paste(
    paste0("<b>",c(rownames(pos), names(clust.moa)),"</b>"),
    c(as.character(pos$moa), clust.moa),
    c(as.character(pos$target), clust.target),
    sep="<br>") 

scatterD3( posx[,1], posx[,2], col_var=clustx,
          tooltip_text = tooltip_text,
          lab=lab, labels_size=10,
          point_size=36 )



if(0) {
    
    pos1 <- read.csv(file="l1000-drugs-tsne-px30.csv", row.names=1)
    pos2 <- read.csv(file="l1000-drugs-tsne-px30s.csv", row.names=1)
    pos3 <- read.csv(file="l1000-drugs-tsne-px30q.csv", row.names=1)
    pos4 <- read.csv(file="l1000-drugs-tsne-px30d.csv", row.names=1)
    pos5 <- read.csv(file="l1000-drugs-tsne-px30d2.csv", row.names=1)
    pos6 <- read.csv(file="l1000-drugs-tsne-px30d4.csv", row.names=1)
    pos7 <- read.csv(file="l1000-drugs-tsne-px30d8.csv", row.names=1)

    par(mfrow=c(3,3),mar=c(4,4,2,2))
    plot( pos1[,1], pos1[,2], pch=20, cex=0.3, main="X")
    plot( pos2[,1], pos4[,2], pch=20, cex=0.3, main="scaled X")
    plot( pos3[,1], pos2[,2], pch=20, cex=0.3, main="q weighted")
    plot( pos4[,1], pos4[,2], pch=20, cex=0.3, main="distance")
    plot( pos5[,1], pos5[,2], pch=20, cex=0.3, main="distance^2")
    plot( pos6[,1], pos6[,2], pch=20, cex=0.3, main="distance^4")
    plot( pos7[,1], pos7[,2], pch=20, cex=0.3, main="distance^8")
    
    require(scatterD3)
    scatterD3( pos1[,1], pos1[,2], lab=rownames(pos))
}
