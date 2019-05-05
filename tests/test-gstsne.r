source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")

library(Rtsne.multicore)
library(qlcMatrix)

load(file="../pgx/geiger2018-arginine-4k.pgx",verbose=1)
load(file="../pgx/rieckmann2017-immprot-4k.pgx",verbose=1)
##load(file="../pgx/GSE22886-immune-mRNA-4k.pgx",verbose=1)
load(file="gset-sparseG.rda",verbose=1)
dim(G)
names(ngs)

## -------------------- get FC matrices ----------------------------
## gene fold-changes
names(ngs$gx.meta$meta)
F <- sapply( ngs$gx.meta$meta, function(x) unclass(x$fc)[,1])
F <- F[match(rownames(G),rownames(F)),]
rownames(F) <- rownames(G)
F[is.na(F)] <- 0
F <- F / max(abs(F),na.rm=TRUE)
dim(F)

## get gene set fold-changes
names(ngs$gset.meta$meta)
mx <- ngs$gset.meta$meta
S <- sapply( mx, function(x) unclass(x$fc)[,"gsva"])
S <- S / max(abs(S),na.rm=TRUE)
dim(S)

## -------------------- create SNN graph ----------------------------
##BiocManager::install("scran", version = "3.8")
library(scran)
library(igraph)
##library(threejs)

pos <- read.csv("tsne-all-genesets.csv",row.names=1)
pos <- read.csv("tsne-all-genesets-rk4k.csv",row.names=1)
##pos <- read.csv("tsne-all-genesets-rk4k-fromX.csv",row.names=1)
##pos <- ngs$tsne2d.gset ## from gene set
##dim(pos)

pos <- as.matrix(pos)
plot(pos, cex=0.2, pch=".")

pos <- pos[which(rownames(pos) %in% colnames(G)),]
dim(pos)
pos <- pos + 1e-2*matrix(rnorm(length(pos)),nrow(pos),ncol=2)
g0 <- buildSNNGraph(t(pos), k=10)
##g <- buildSNNGraph(G1, k=20)
g0

##F1 <- F[rownames(pos),]

## --------------- create reduced graph ----------------------------

## cluster the graph
hc <- hclust_graph(g0)
dim(hc)
hc.sum <- apply(hc,2,function(x) length(table(x)))
hc.sum

## merge nodes
j <- max(which(hc.sum < 4000))
j
idx0 <- hc[,j]
pos1 <- apply(pos,2, function(x) tapply(x,idx0,mean))
dim(pos1)
ii <- match(rownames(pos),rownames(S))
S1 <- apply(S[ii,],2, function(x) tapply(x,idx0,mean,na.rm=TRUE))
S1[is.na(S1)] = 0
new.names <- tapply(rownames(pos),idx0,function(s) paste(s,collapse="\n"))
rownames(pos1) <- as.character(new.names)
##rownames(S1) <- as.character(new.names)

##g1 = g0
g1 <- buildSNNGraph(t(pos1), k=10)
g1
idx <- tapply(hc[,1], idx0, median)
summary(idx)

## reduce also the GMT matrix
tmpG <- G[,rownames(pos)]
##tmpG <- sweep( tmpG, 2, Matrix::colSums(tmpG))
tmpG <- tmpG %*% Diagonal(x =1.0/(1e-6+Matrix::colSums(tmpG))**0.5 )   ## scale columns??
G1 <- tapply(1:ncol(tmpG), idx0, function(j) Matrix::rowSums(tmpG[,j,drop=FALSE]))
G1 <- do.call(cbind, G1)
dim(G1)
colnames(G1) <- rownames(pos1)

## keep group indices
grp.idx <- tapply( rownames(pos), idx0, list)

dim(S1)
dim(G1)

if(0) {
    FILES="../files/"
    res <- createReducedGenesetGraph(ngs, max.nodes=2000)
    names(res)
    g1 <- res$graph
    idx <- V(g1)$cluster
    pos1 <- g1$layout
    S1 <- res$fc
    G1 <- res$gmt
}

## ------------- plot using base::plot ---------------------------------

klrpal = colorRampPalette(c("blue3", "grey90", "red3"))(32)
## klrpal <- jet.colors(16)
kl.col <- col2hex(rep(rainbow(32),99)[idx])

## show all contrasts
par(mfrow=c(4,4), mar=c(0.1,0.1,1,0.1))
plot(pos1[,], col=kl.col, pch=20, cex=0.6, xaxt="n", yaxt="n")
i=1
for(i in 1:ncol(S1)) {
    fc <- S1[,i]
    fc <- sign(fc) * abs(fc / max(abs(fc)))**0.33
    fc.col <- klrpal[round(16 + 15*fc)]
    plot(pos1[,], col=fc.col, ##pch=".",
         cex= (0.1+1*abs(fc))**1,
         pch=20, xaxt="n", yaxt="n",
         main=colnames(S)[i], cex.main=1.1)
}

## ------------- show single contrast ---------------------------------

## choose single FC contrast
colnames(S1)
contr=colnames(S1)[1]
fc <- S1[,contr]
head(fc)

## show marker genes
##fx <- (G %*% S)
kk <- 2:6
fx <- (G1 %*% S1)[,kk,drop=FALSE]
rho <- rowSums(fx * F[,kk,drop=FALSE])
sum(is.na(rho))

top.down <- names(head(rho[order(rho)],10))
top.up   <- names(head(rho[order(-rho)],10))
top.markers <- c(top.up, top.down)
##top.markers <- head(rownames(fx)[order(-apply(fx,1,sd))],14)
##top.markers <- head(rownames(fx)[order(-apply(F,1,sd))],14)
top.markers
##rho[top.markers]

par(mfrow=c(5,5), mar=c(0.1,0.1,1,0.1)*2)
kl.col <- col2hex(rep(rainbow(32),99)[idx])
plot(pos1[,], col=kl.col, pch=19, cex=0.5, xaxt="n", yaxt="n")

fx1 <- rowSums(S1[,kk]**2)**0.5
fx1 <- S1[,2]
fx1 <- fx1 / max(abs(fx.sd))
fx.col <- klrpal[round(16 + 15*fx1)]
plot(pos1[,], col=fx.col, pch=19, cex=0.5, xaxt="n", yaxt="n")
title(contr,cex.main=1.2)

for(m in top.markers) {
    plot(pos1[,], col=fx.col, cex=0.5+1*abs(fc), pch=20,
         xaxt="n", yaxt="n")
    wt = G1[,m]
    wt = wt / max(wt)
    jj <- which(wt>=0.1)
    wt1 = (0.1 + wt[jj])**0.5
    points( pos1[jj,], col="green2", pch="o",
           cex=0.1+0.9*wt1, font=2)
    title(m,cex.main=1.2)
}


## ------------ show top markers in gene space -----------------

par(mfrow=c(5,5), mar=c(0.1,0.1,1,0.1)*2)

tsne.genes <- ngs$tsne2d.genes
tsne.genes <- read.csv("tsne-all-genes.csv",row.names=1)
dim(tsne.genes)

if(0) {
    ## genes tSNE using GMT matrix
    G2 <- G[,order(-apply(G,2,sd))]
    G2 <- G2[,1:4000]
    G2 <- ngs$X
    tsne.genes <- Rtsne.multicore( as.matrix(G2), check_duplicates=FALSE, num_threads=24)$Y
    rownames(tsne.genes) <- rownames(G2)
}

plot( tsne.genes[,], col="grey80", pch=20, cex=0.3, xaxt="n", yaxt="n")
jj <- setdiff( match( top.markers, rownames(tsne.genes)), NA)
points( tsne.genes[jj,], col="red", pch=20, cex=0.5)


## genes tSNE using expression
tsne.gx <- Rtsne.multicore( ngs$X, check_duplicates=FALSE,
                           num_threads=24 )$Y
rownames(tsne.gx) <- rownames(ngs$X)
plot( tsne.gx[,], col="grey80", pch=20, cex=0.3, xaxt="n", yaxt="n")
jj <- setdiff( match( top.markers, rownames(tsne.gx)), NA)
points( tsne.gx[jj,], col="red", pch=20, cex=0.5)






## ------------------ plot using igraph -----------------

klr <- col2hex(c("blue3","grey80","red3")[2+sign(fc)])
klr <- substring(klr,1,7)
names(klr) <- rownames(pos1)
V(g1)$color <- klr
V(g1)$label <- NA
V(g1)$size  <- 2.4*(0.1+abs(fc))**0.66
ee <- get.edges(g1, E(g1))
E(g1)$color <- paste0(klr[ee[,1]],"22")
##E(g)$color <- NA

## plot
par(mfrow=c(1,1), mar=c(0.1,0.1,1,0.1))
plot.igraph(g1, layout=pos1[,], vertex.frame.color=NA,
            ##xlim=range(pos1[,1]), ylim=range(pos1[,2]),
            edge.width=0.2, rescale=TRUE)

## ------------------ plot using visNetwork (zoomable) -----------------
library(visNetwork)

vx <- rownames(pos1)
V(g1)$name <- vx
V(g1)$label <- vx
V(g1)$title <- gsub("\n","<br>",vx)
V(g1)$label.size <- 0.01*V(g1)$size**2
V(g1)$size <- 25*(abs(fc)**0.5)
V(g1)$color <- paste0(col2hex(c("blue3","grey80","red3")[2+sign(fc)]),"44")
visdata <- toVisNetworkData(g1, idToLabel=FALSE)
pos2 = cbind(pos1[,1], -pos1[,2])

graph <- visNetwork(nodes = visdata$nodes, edges = visdata$edges,
                    height="1200px", width="1600px") %>%
    visNodes(font=list(size=4))  %>%
    visEdges(hidden=FALSE, width=0.02, color=list(opacity=0.1))  %>%
    visInteraction(hideEdgesOnDrag = TRUE) %>%
    visIgraphLayout(layout="layout.norm", layoutMatrix=pos2)
graph

