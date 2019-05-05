source("../R/gx-heatmap.r")
source("../R/gx-combat.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")

load("../pgx/geiger2018-arginine-4k.pgx",verbose=TRUE)
##load("../pgx/GSE72056-melanoma-scRNA-vsCLUST-TCELL-8k.pgx",verbose=TRUE)


genes <- head( rownames(ngs$X))

X <- ngs$X[1:5,]

makeGeneContrasts <- function(X) {
    calls = apply(X, 1, function(x) c(-1,1)[1 + 1*(x > mean(x,na.rm=TRUE))])
    dim(calls)
    rownames(calls) <- colnames(X)
    i=1
    for(i in 1:ncol(calls)) {
        jj <- which(calls[,i]!=0)
        design1 <- cbind(1, 1*(calls[jj,i]==1) )
        vfit <- lmFit(X[,jj], design=design1)
        efit <- eBayes(vfit, trend=TRUE)
        top3 = topTable(efit, coef=2, sort.by="none",number=Inf, adjust.method="BH")
        head(top3)
    }

    res <- gx.limma(X, calls[,1])


    
}




genes <- head(rownames(ngs$counts))
genes
call.matrix = sapply(genes, function(g) c("low","high")[1 + (X[g,] > mean(X[g,],na.rm=TRUE))])
rownames(call.matrix) <- colnames(ngs$counts)
call.matrix


## ----------------- chck my old function
y <- exp.matrix[,1]
y[which(y==0)] <- NA
top2 <- gx.limma(X, y, fdr=1, lfc=0, trend=TRUE)[rownames(X),]
head(top2)
