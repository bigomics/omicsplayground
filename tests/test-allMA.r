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

calcAllPairsMA <- function(X, y) {

    pairs <- t(combn(unique(y),2))
    M <- c()
    A <- c()
    Q <- c()
    i=1
    for(i in 1:nrow(pairs)) {
        j0 <- which( y == pairs[i,1])
        j1 <- which( y == pairs[i,2])

        ## compute mean, diff, q-value
        f1 <- Matrix::rowMeans(X[,j1,drop=FALSE]) - Matrix::rowMeans(X[,j0,drop=FALSE])
        m1 <- Matrix::rowMeans(X[,c(j0,j1),drop=FALSE])

        ## add reverse
        f2 <- cbind(f1, -f1)
        m2 <- cbind(m1, m1)
        colnames(f2) <- c( paste(pairs[i,2:1],collapse="_vs_"),
                          paste(pairs[i,1:2],collapse="_vs_") )
        colnames(m2) <- colnames(f2)

        M <- cbind(M, f2)
        A <- cbind(A, m2)
    }
    return(list(M=M, A=A, Q=Q))
}

PGX.FILES <- c("geiger2018-arginine-8k-LT.pgx",
               "geiger2018b-liver-fltSC-8k-LT.pgx",
               "gonzalez2018-nkvaccination.pgx",
               "GSE22886-immune-mRNA-8k-LT.pgx",
               "GSE10846-dlbcl-mRNA-8k.pgx",
               "GSE72056-melanoma-scRNA-8k-LT.pgx",
               "GSE98638-liver-scRNA-8k-LT.pgx",
               "guarda2018-IL15-8k-LT.pgx",
               "rieckmann2017-immprot-8k-LT.pgx",
               "sallusto2018-rorc-SC-8k-LT.pgx",
               "sallusto2019-th1star-SC-12x-LT.pgx",
               "schmiedel2018-DICE-mRNA-8k-LT.pgx")

## public
PGX.FILES <- c("geiger2018-arginine-8k-LT.pgx",
               "GSE22886-immune-mRNA-8k-LT.pgx",
               "GSE10846-dlbcl-mRNA-8k.pgx",
               "GSE98638-liver-scRNA-8k-LT.pgx",
               "GSE72056-melanoma-scRNA-vsphenotype-TC8k-LT.pgx",
               "rieckmann2017-immprot-8k-LT.pgx",
               "schmiedel2018-DICE-mRNA-8k-LT.pgx")

allX <- list()
allY <- list()
pgx=PGX.FILES[2]
for(pgx in PGX.FILES) {
    load(file.path("../pgx",pgx),verbose=1)
    rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
    allX[[pgx]] <- ngs$X
    ##allY[[pgx]] <- gsub("[-_.]","",ngs$samples$group)
    allY[[pgx]] <- ngs$samples$group
}
names(allX)

## find common genes
##gg <- Reduce(intersect, lapply(allX, rownames))
gg.tbl <- table(unlist(sapply(allX, rownames)))
##gg <- names(which( gg.tbl >= length(allX)))
table(gg.tbl)
gg <- head(names(sort(-gg.tbl)),2000)
length(gg)
allX <- lapply(allX, function(x) {
    x <- x[match(gg,rownames(x)),]
    rownames(x) <- gg
    return(x)
})

## compute all pairs MA
allM <- c()
allA <- c()
i=1
for(i in 1:length(allX)) {
    X=allX[[i]]
    y=allY[[i]]
    res <- calcAllPairsMA(X, as.character(y))
    bx <- strsplit(PGX.FILES[i], split="[-_]")[[1]][1]
    colnames(res$A) <- paste0("[",bx,"]",colnames(res$A))
    colnames(res$M) <- paste0("[",bx,"]",colnames(res$M))
    allM <- cbind(allM, res$M)
    allA <- cbind(allA, res$A)
}
dim(allM)
dim(allA)
head(colnames(allA))
head(colnames(allM))

save(allM, allA, file="../files/allMA-pub.rda")

load(file="../files/allMA.rda", verbose=1)


batch <- factor(sub("\\].*","]",colnames(allM)))
table(batch)
nbatch <- length(unique(batch))
M0 <- scale(apply(allM,2,rank,na.last="keep"))
M0 <- M0[order(-apply(M0,1,sd)),]
A0 <- scale(apply(allA,2,rank,na.last="keep"))

M1 <- head(M0,1000)
R <- cor(M1, use="pairwise")

grep("star", rownames(R),value=TRUE)
query="[geiger2018b]CD8_HCC_vs_CD8_Blood"
query="[sallusto2019]Th1star_ut_vs_Th1_ut"
##query="[GSE10846]ABC_vs_GCB"
head(sort(R[query,],decreasing=TRUE),20)

nb <- names(head(sort(R[query,],decreasing=TRUE),40))
gx.heatmap( R[nb,nb], scale="none", mar=c(1,1)*20 )

nb.sorted <- names(sort(R[query,],decreasing=TRUE))
NTOP=10
nb1 <- head(grep("GSE98",nb.sorted,value=TRUE),NTOP)
nb1 <- head(nb.sorted,NTOP)
nb <- c(query, nb1)
nb
length(nb)
m1 <- rowMeans(allM[,nb],na.rm=TRUE)
f1 <- rowSums(abs(allM[,nb]) > 1, na.rm=TRUE)
a1 <- rowMeans(allA[,nb], na.rm=TRUE)
f1x <- f1 + rnorm(length(f1))*diff(range(f1,na.rm=TRUE))*0.02

## Volcano-like plot
par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(m1, f1x, pch=20, col="grey70",
     xlim=c(-1,1)*max(abs(m1),na.rm=TRUE),
     xlab="average FC", ylab="frequency (#sets)")
abline(v=0)
gg <- head(order(-abs(m1)),20)
text(m1[gg], f1x[gg], names(m1)[gg], col="blue", cex=1.1)
names(m1)[gg]

## MA plot
cex1 <- 2*((1+f1)/max(f1,na.rm=TRUE))
plot(a1, m1, pch=20, cex=cex1, col="grey70",
     ylab="fold-change (log2R)",
     xlab="average signal (log2)")
abline(h=0)
text(a1[gg], m1[gg], names(m1)[gg], col="blue", cex=1.1)
names(m1)[gg]


require(corrplot)
mx <- allM[gg,nb]
ii <- hclust(dist(mx))$order
jj <- hclust(dist(t(mx)))$order
mx <- mx[ii,jj]
##gx.heatmap( mx, mar=c(15,10), cexRow=1, cexCol=1, scale="none")
par(mfrow=c(1,1), mar=c(4,4,2,2))
mx1 <- mx * (abs(mx)>1)
corrplot( mx1, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
         tl.cex=0.7, tl.col="grey33", tl.srt=30,
         mar=c(1,0,1,0)*0.5 )



pdf("test-fcsne-vsN-rnk.pdf")
n=1000
for(n in c(10,50,200,500,1000)) {

    M1 <- head(M0,n)
    R <- cor(M1, use="pairwise")
    if(0) {
        A1 <- A0[rownames(M1),]
        R1 <- cor(A1, use="pairwise")
        R <- (R * R1)
    }
    pos <- Rtsne( (1-R), is_distance=TRUE, check_duplicates=FALSE)$Y
    rownames(pos) <- colnames(M1)
    pal1 <- rainbow(nbatch)
    klr1 <- pal1[as.integer(batch)]
    plot( pos, pch=20, cex=0.2, col=klr1, main=paste("N=",n) )
    legend("topright", legend=levels(batch), fill=pal1,
           cex=0.6, y.intersp=0.8)

}
dev.off()


