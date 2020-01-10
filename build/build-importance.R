##install.packages("xlsx")
##install.packages("readxl")
##install.packages("gdata")
library(xlsx)
library(readxl)
library(gdata)
library(tidyverse)

if(0) {
    b1 <- read.xlsx("OncoLnc_TableS1.xlsx",1)
    b2 <- read.xls("OncoLnc_TableS1.xls",1)
    b3 <- readxl::readxl("OncoLnc_TableS1.xls",1)
    b4 <- tidyverse::readxl("OncoLnc_TableS1.xls",1)
    b5 <- read_excel("OncoLnc_TableS1.xls",1)

    b2 <- read.xls("../pub/oncolnc/OncoRank_TableS1.xls",1)
    ##b3 <- readxl::readxl("../pub/oncolnc/OncoRank_TableS1.xls",1)
    b4 <- tidyverse::readxl("../pub/oncolnc/OncoRank_TableS1.xls",1)
    b5 <- read_excel("../pub/oncolnc/OncoRank_TableS1.xls",1)
}

dir("../pub/oncolnc")
path = "../pub/oncolnc/OncoRank_TableS1.xls"
oncorank <- read.xls(path,1)

rnk <- oncorank$REC.Score
names(rnk) <- oncorank$TCGA.Name
par(mfrow=c(2,1))
barplot( sort(rnk[head(order(-abs(rnk)),60)]),
        las=3, cex.names=0.5)

path = "../pub/oncolnc/OncoLnc_TableS1.xls"
xdata <- lapply(1:21, function(i) read_excel(path,i))
names(xdata) <- excel_sheets(path)

X <- lapply(xdata, function(x) {
    f <- unclass(x[,"Cox coefficient"][[1]])
    names(f) <- unclass(x["TCGA Name"][[1]])
    return(f)
})

genes <- Reduce(intersect, sapply(X,names))
X <- sapply(X, function(x) x[genes])
colnames(X) <- names(xdata)

X <- scale(X, center=FALSE)
pos.X <- pmax(X,0)
neg.X <- pmin(X,0)

jj <- head(order(-rowSums(abs(X))),60)
jj <- jj[order(-rowSums(X[jj,]))]

heatmap(X[jj,k], mar=c(4,18), scale="none",
        cexCol=0.8, cexRow=0.7)

par(mfrow=c(2,1), mar=c(4,4,2,2))
N=100
N=50
barplot( t(pos.X[jj,]), beside=FALSE, las=3,
        ylim=c(-1,1)*30,
        main = "cumulative risk",
        ylab="cumulative hazard rate", cex.names=0.55)
barplot( t(neg.X[jj,]), beside=FALSE, add=TRUE,
        cex.names=0.001)


library(GSVA)
load("../playground/lib/gmt-all.rda",verbose=1)
gmt <- gmt.all
gmt <- gmt.all[grep("HALLMARK",names(gmt.all))]

sX <- X - rowMeans(X)
Z <- gsva(sX, gmt)
##Z <- gsva(sX, gmt, method="ssgsea")
 dim(Z)

##rownames(Z) <- sub(".*HALLMARK_","",rownames(Z))
##Z <- scale(Z, center=FALSE)
pos.Z <- pmax(Z,0)
neg.Z <- pmin(Z,0)

k <- 1:ncol(Z)
k <- 1:2
k <- 2
jj <- 1:nrow(Z)
jj <- grep("HALLMARK",rownames(Z))
jj <- grep("^DRUG",rownames(Z))
jj <- head(jj[order(-rowSums(abs(Z[jj,k,drop=FALSE])))],60)
jj <- jj[order(-rowSums(Z[jj,k,drop=FALSE]))]
maxz <- max(rowSums(pos.Z[jj,k]),rowSums(neg.Z[jj,k]))

heatmap(Z[jj,k], mar=c(4,18), scale="none",
        cexCol=0.8, cexRow=0.7)

par(mfrow=c(2,1), mar=c(4,4,2,2))
par(mfrow=c(1,2), mar=c(4,4,2,2));frame()
N=100
N=60
barplot( t(pos.Z[jj,k]), beside=FALSE, horiz=TRUE, las=1,
        xlim=c(-1.3,1.3)*maxz,
        main = "increased risk (shorter survival)",
        xlab="cumulative hazard rate", cex.names=0.55)
barplot( t(neg.Z[jj,k]), beside=FALSE, add=TRUE,
        horiz=TRUE, las=1, cex.names=0.001)
cc1 <- grey.colors(length(k))
legend("topright", colnames(Z)[k], fill=cc1,
       cex=0.6, y.intersp=0.8)



library(data.table)
A <- fread("~/Downloads/Achilles_gene_effect.csv")
dim(A)





