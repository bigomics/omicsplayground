
if(0) {
    load(file="../pgx/rieckmann2017-immprot-4k.pgx",verbose=1)
    imm.sig <- t(apply(ngs$count, 1, function(x) tapply(x, ngs$samples$group, mean)))
    imm.sig <- head(imm.sig[order(-apply(log(100+imm.sig),1,sd)),],500)
    ##imm.sig <- head(imm.sig[order(-apply(imm.sig,1,sd)),],500)
    colnames(imm.sig) <- sub("_S$",".resting",colnames(imm.sig))
    colnames(imm.sig) <- sub("_A$",".activated",colnames(imm.sig))
    write.csv(imm.sig, file="../files/immprot-signature1000.csv")
    remove(ngs)
}

require(limma)
require(corrplot)
require(FARDEEP)
source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")

FILES="../files"
load("../pgx/geiger2018b-liver-fltSC-8k-LT.pgx")
load("../pgx/sallusto2018-rorc-SC-8k-LT.pgx")
load("../pgx/sallusto2019-th1star-SC-12x-LT.pgx")
y <- ngs$samples$group
STAR <- t(apply(2**ngs$X,1,function(x) tapply(x,y,mean)))
head(STAR)                                             
load("../pgx/schmiedel2018-DICE-mRNA-8k-LT.pgx")
load("../pgx/tenx-pbmc1k-8k-LT.pgx")

source("../R/pgx-deconv.R")
##load(file=rda.file,verbose=1)

IMMPROT <- read.csv(file.path(FILES,"immprot-signature1000.csv"),row.names=1)
DICE    <- read.csv(file.path(FILES,"DICE-signature1000.csv"),row.names=1)
dim(IMMPROT)

## create large compilation of reference profiles
gg <- sort(unique(c(rownames(DICE),rownames(LM22),rownames(IMMPROT))))
META <- cbind( LM22[match(gg,rownames(LM22)),],
              IMMPROT[match(gg,rownames(IMMPROT)),],
              DICE[match(gg,rownames(DICE)),],
              STAR[match(gg,rownames(STAR)),])
rownames(META) <- gg
META <- t(t(META)/colSums(META,na.rm=TRUE)) * 1e6
sdx <- apply(log(1+META),1,sd,na.rm=TRUE)*rowMeans(!is.na(META))
META <- head(META[order(-sdx),],1000)
dim(META)
sum(is.na(META))
META <- t(exp(impute(t(log(1+META)),what="mean"))-1)

ref <- META
ref <- DICE
ref <- STAR
ref <- IMMPROT
ref <- TISSUE
ref <- as.matrix(LM22)

X <- 2**ngs$X
dim(X)
if(ncol(X)<60) X <- X[,sample(ncol(X),60,replace=TRUE)]
if(ncol(X)>60) X <- X[,sample(ncol(X),60,replace=FALSE)]
dim(X)

jj <- intersect(rownames(X),rownames(ref))
X <- X[jj,]
ref <- ref[jj,]
dim(ref)

methods = setdiff(DECONV.METHODS,c("CIBERSORT","EPIC","Abbas","FARDEEP"))
methods
dres <- pgx.deconvolution(X, ref=ref, methods=methods)
scores <- dres$results

mm <- names(scores)
mm
grp <- ngs$samples$cluster
for(i in 1:length(mm)) {
    scores[[i]] <- apply(scores[[i]],2,function(x)
        tapply(x,grp,mean))
}


Y <- t(scores[["meta"]])
Y[which(is.na(Y))] <- 0
Y <- head(Y[order(-rowMeans(Y**2,na.rm=TRUE)),],24)
ii <- rownames(Y)[hclust(dist(Y))$order]
jj <- colnames(Y)[hclust(dist(t(Y)))$order]

cex=1
par(mfrow=c(3,3), mar=c(11,1,2,12))
for(m in mm) {
    Y <- t(scores[[m]])
    Y1 <- pmax(Y[ii,jj],0)
    if(m=="cor") Y1 <- sign(Y1) * Y1**4
    ##corrplot(Y1)
    Y1 <- t(t(Y1)/ (1e-8 + colSums(Y1,na.rm=TRUE)))
    ##gx.imagemap(Y1, clust=FALSE)
    image( 1:ncol(Y1), 1:nrow(Y1), t(Y1),
          xaxt="n", yaxt="n")
    mtext(colnames(Y1), side=1, at=1:ncol(Y1), las=3,
          cex=0.75*cex, line=0.5 )
    mtext(rownames(Y1), side=4, at=1:nrow(Y1), las=1,
          cex=0.85*cex, line=0.5 )
    title(main=m,cex=cex,line=0.5)
}
