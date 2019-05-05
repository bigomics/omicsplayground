source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")
library(Rtsne.multicore)
library(Rtsne)
library(qlcMatrix)
library(corpora)
##devtools::install_github("bwlewis/rthreejs")

load(file="../files/allMA.rda", verbose=1)
dim(allM)
load("../pgx/geiger2018-arginine-8k-LT.pgx", verbose=1)


GMT <- readRDS("../files/gset-sparseG-XL.rds")
GMT <- t(GMT)
dim(GMT)
head(rownames(GMT))


names(ngs$gx.meta$meta)
query <- ngs$gx.meta$meta[[3]]$meta.fx
names(query) <- rownames(ngs$gx.meta$meta[[3]])

gg <- intersect(names(query), rownames(allM))
length(gg)
gg <- intersect(gg, rownames(GMT))
length(gg)

require(qlcMatrix)
qr <- cbind(rank(query[gg]))
rho1 <- corSparse(GMT[gg,], qr)[,1]
rho2 <- cor(allM[gg,], qr, use="pairwise")[,1]
names(rho1) <- colnames(GMT)
names(rho2) <- colnames(allM)

top1 <- head(names(sort(-abs(rho1))),1000)
top2 <- head(names(sort(-abs(rho2))),1000)

X <- cbind( "QUERY:signature"=query[gg], GMT[gg,top1])
X <- cbind( X, allM[gg,top2] )

X1 <- head(X[order(-apply(X,1,sd)),],1000)
X1 <- X1[,which(apply(X1,2,sd)>0)]

require(Rtsne.multicore)
R <- cor(as.matrix(X1), use="pairwise")
pos <- Rtsne.multicore( (1-abs(R)), is_distance=TRUE)$Y
rownames(pos) <- rownames(R)
colnames(pos) <- c("x","y")

head(sort(R[1,]))
tail(sort(R[1,]))

lab <- rep("",nrow(pos))
rdist <- colSums((t(pos) - pos[1,])**2)
jj <- head(order(rdist),100)
lab[jj] <- rownames(pos)[jj]

grp <- sub(":.*","",rownames(pos))
grp <- sub("].*","]",grp)
table(grp)
klr <- c(2, rep(1,nrow(pos)-1))
qq <- c("QUERY", rep("signature",nrow(pos)-1))
rr <- abs(R[rownames(pos),1])**2

require(scatterD3)
scatterD3( pos[,1], pos[,2], tooltip_text=rownames(pos),
          lab=lab, labels_size=15,
          size_var=rr, size_range = c(20,300),
          symbol_var=qq, col_var=grp)

require(plotly)
plot_ly( data=data.frame(pos), x=~x, y=~y,
       name=rownames(pos),
       lab=lab, labels_size=15,
       color=grp)

