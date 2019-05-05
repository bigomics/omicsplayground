##BiocManager::install("crossmeta", version = "3.8")
#BiocManager::install("ccmap", version = "3.8")

setwd("~/Projects/Immunomics/tests/")
##library(crossmeta)
library(ccmap)
library(ccdata)
source("../R/pgx-drugs.R")

data(cmap_es)
data(l1000_es)
dim(cmap_es)
dim(l1000_es)

rankScore <- function(query_sig, ref_set, n=200) {
    gg <- intersect(rownames(ref_set), names(query_sig))
    gg <- setdiff(gg, c("",NA,"NA"))
    if(n>0) {
        q <- query_sig[gg]
        gg <- gg[unique(c(head(order(q),n),head(order(-q),n)))]
    }
    suppressWarnings( rho <- cor( as.matrix(ref_set[gg,]), query_sig[gg],
                                 use="pairwise"))
    ##rho <- rho[order(-rowMeans(rho**2,na.rm=TRUE)),,drop=FALSE]
    rho <- rho[order(-rowMeans(rho,na.rm=TRUE)),,drop=FALSE]
    if(NCOL(rho)==1) rho <- rho[,1]
    return(rho)
}

load(file="../files/allMA.rda", verbose=1)
head(colnames(allM),20)
query_sig <- allM[,"[guarda2018]IL15_WT_Torin_vs_IL15_WT_nd"]
str(query_sig)

##------------------------------------------------------------
## Query best scoring drugs
##------------------------------------------------------------

## just uses cosine distance
system.time( top_cmap  <- query_drugs(query_sig, cmap_es, ngenes=200) )
system.time( top_l1000 <- query_drugs(query_sig, l1000_es, ngenes=200) )

head(top_cmap,10)
head(top_l1000,20)

## with my 'method'...
system.time( rnk_cmap  <- rankScore(query_sig, cmap_es, n=200))
system.time( rnk_l1000 <- rankScore(query_sig, l1000_es, n=200))
head( rnk_cmap,10)
head( rnk_l1000,10)

X <- readRDS("../files/l1000_es_5685drugs.rds")
x.drugs <- sub("_.*$","",colnames(X))
tail(sort(table(x.drugs)))
F <- cbind(query_sig,query_sig)
F <- cbind(query_sig)

obj=F;methods="GSEA";contrast=NULL;nprune=200
out <- pgx.computeDrugEnrichment(F, X, x.drugs, methods="cor")
str(out[[1]])
head(out[[1]]$X)


gg <- intersect(names(query_sig),rownames(X))
x <- query_sig[gg]
sel <- which(x.drugs=="torin-2")
y <- X[gg,sel]
y.avg <- rowMeans(y)
y <- cbind(y.avg, y)
rho <- cor(y,x,use="pairwise")[,1]
dim(y)

par(mfrow=c(10,7), mar=c(1,1,1,1)*0)
i=1
ii <- order(-rho)
for(i in ii) {
    plot(x, y[,i], pch=20, cex=0.5)
    abline(lm(y[,i] ~x), col="red")
    legend("topleft",colnames(y)[i],cex=0.5)
}




##------------------------------------------------------------
## Query combinations
##------------------------------------------------------------

# query all 856086 combinations (takes ~2 minutes on Intel Core i7-6700)
top_combos <- query_combos(query_sig, cmap_es)
head(top_combos)

# query only combinations with LY-294002
top_combos <- query_combos(query_sig, cmap_es, include='LY-294002', ncores=1)
head(top_combos)

##------------------------------------------------------------
## Using ML
##------------------------------------------------------------

# Times on Intel Core i7-6700 with MRO+MKL
# requires ~8-10GB of RAM
method  <- 'ml'
include <- names(head(top_cmap))
include

# query all 856086 combinations (~2 hours)
# top_combos <- query_combos(query_sig, 'cmap', method)

# query combinations with top single drugs (~1 minute)
top_combos <- query_combos(query_sig, 'cmap', method, include)



##================================================================================
##================================================================================
##================================================================================

source("../R/gset-gsea.r")

data(cmap_es)
data(l1000_es)
dim(cmap_es)
dim(l1000_es)

colnames(l1000_es) <- paste0("[L1000] ",colnames(l1000_es))
colnames(cmap_es) <- paste0("[CMAP] ",colnames(cmap_es))
gg <- intersect(rownames(l1000_es),rownames(cmap_es))
X <- cbind( l1000_es[gg,], cmap_es[gg,])

x.drugs <- gsub(".*\\][ ]|_.*$","",colnames(X))
head(sort(table(x.drugs),decreasing=TRUE))

d="vorinostat"
for(d in dd) {
    jj <- which(x.drugs==d)
    pairs( as.matrix(X[,head(jj,12)]), pch=".", cex=3)
}

dd <- names(which(table(x.drugs)>120))
dd
pdf("test-ccmap-heatmap.pdf",w=16,h=9)
for(d in dd) {
    jj <- head(which(x.drugs==d),120)
    dx <- tanh(1* (X[,jj]) )
    ##dx <- tanh(1* scale(X[,jj]) )
    par(mfrow=c(1,1))
    gx.heatmap(dx, mar=c(15,8), scale="none",
               key=FALSE, keysize=0.5)
    title(d)
}
dev.off()

d="geldanamycin"
d="tozasertib"
jj <- head(which(x.drugs==d),120)
dx <- X[,jj]


##================================================================================
##================================================================================
##================================================================================

require(fgsea)
X <- l1000_es
x.drugs <- gsub("_.*$","",colnames(X))
length(table(x.drugs))

gmt <- tapply(colnames(X), x.drugs, list)
gmt <- gmt[which(sapply(gmt,length)>=15)]
length(gmt)
gg <- intersect(rownames(l1000_es),rownames(cmap_es))

sel <- which(x.drugs %in% names(gmt))
length(sel)
saveRDS( X[,sel], file="../files/l1000_es_5685drugs.rds")

common.drugs <- intersect(names(gmt), colnames(cmap_es))
length(common.drugs)
gmt <- gmt[common.drugs]

DX <- gmt2mat(gmt)
dim(DX)

best.matches <- list()
i=1
for(i in 1:length(gmt)) {
    dx <- names(gmt)[i]
    qx <- cmap_es[gg,dx]
    rho <- cor( l1000_es[gg,rownames(DX)], qx)[,1]
    if(1) {
        res <- fgsea(gmt, stats=rho, nperm=1000)
        ##res <- res[order(res$padj,-(res$NES)),]
        res <- res[order(-res$NES),]
        head(res)[,1:7]
        best <- res$pathway
    } else {
        score <- t(DX) %*% rho
        best <- names(sort(score[,1],decreasing=TRUE))
    }
    best.matches[[i]] <- head(best,100)
}
names(best.matches) <- names(gmt)

qc <- sapply(names(best.matches), function(b) (b %in% head(best.matches[[b]],50)))
table(qc)
names(which(qc))
top.gmt <- names(which(qc))
sort(top.gmt)
jj <- which(x.drugs %in% top.gmt)
length(jj)
length(top.gmt)
dim(X[,jj])
saveRDS( X[,jj], file="../files/l1000_es_qc50_166drugs.rds")


pdf("test-ccmap-DSEA-gsrnkcor-TOP.pdf",w=16,h=9)
i=1
for(i in 1:length(top.gmt)) {
    dx <- top.gmt[i]
    qx <- cmap_es[gg,dx]
    rho <- cor( l1000_es[gg,rownames(DX)], qx)[,1]
    if(0) {
        res <- fgsea(gmt, stats=rho, nperm=1000)
        ##res <- res[order(res$padj,-(res$NES)),]
        res <- res[order(-res$NES),]
        head(res)[,1:7]
        bestmatch <- res$pathway
    } else {
        rnk <- scale(rank(rho))[,1]
        score <- t(DX) %*% rnk
        bestmatch <- names(sort(score[,1],decreasing=TRUE))
    }

    par(mfrow=c(5,6), mar=c(2,4,2,1))
    gsea.enplot( rho, gmt[[dx]], cex.main=1.3,
                main=paste0(dx,"   @",dx) )
    for(i in 1:29) {
        d1 <- bestmatch[i]
        cx <- ifelse(d1==dx,1.3,1)
        gsea.enplot( rho, gmt[[d1]], cex.main=cx,
                    main=paste0(d1,"   @",dx) )
    }

}
dev.off()


##================================================================================
##================================================================================
##================================================================================
source("../R/pgx-drugs.R")
source("../R/pgx-functions.R")

load("../pgx/guarda2019-myc-12k-LT.pgx")

FILES="../files"
X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
x.drugs <- gsub("_.*$","",colnames(X))
length(table(x.drugs))
avgX <- readRDS(file=file.path(FILES,"l1000_es_5685drugsAVG.rds"))

contrast <- c("IL15_Torin_vs_IL15_nd_1","IL15_Torin_vs_NT_nd_1")

out <- pgx.computeDrugEnrichment(
    ngs, X, x.drugs, methods=c("cor","GSEA"),
    avgX=avgX, contrast=contrast)
names(out)

F <- pgx.getMetaFoldChangeMatrix(ngs)$fc
dim(F)
rownames(F) <- toupper(rownames(F))
out <- pgx.computeDrugEnrichment(
    F, X, x.drugs, methods=c("cor","GSEA"),
    avgX=avgX, contrast=contrast)

names(out)
head(out[[1]]$X)
head(out[[2]]$X)

dim(avgX)
dx.rho <- pgx.computeDrugEnrichment(
    obj=avgX, X, x.drugs, methods="cor", contrast=NULL)

dim(dx.rho[[1]]$X)
xdist <- 1- dx.rho[[1]]$X * (1 - dx.rho[[1]]$Q)
pos <- Rtsne( xdist, is_distance=TRUE)$Y
rownames(pos) <- rownames(xdist)
write.csv(pos, file="l1000-drugs-tsne.csv")

require(scatterD3)
scatterD3( pos[,1], pos[,2], lab=rownames(xdist))
