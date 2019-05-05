source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/gset-gsea.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")

library(Rtsne.multicore)
library(qlcMatrix)

##load(file="../pgx/geiger2018-arginine-4k.pgx",verbose=1)
load(file="../pgx/GSE10846-dlbcl-mRNA-8k.pgx",verbose=1)
load(file="../files/gset-sparseG-XL.rda",verbose=1)
dim(G)
names(ngs)


## -------------------- get FC matrices ----------------------------
## gene fold-changes
names(ngs$gx.meta$meta)
F <- sapply( ngs$gx.meta$meta, function(x) unclass(x$fc)[,1])
F <- F[match(colnames(G),rownames(F)),,drop=FALSE]
rownames(F) <- colnames(G)
F[is.na(F)] <- 0
F <- F / max(abs(F),na.rm=TRUE)
dim(F)
head(F)

## get gene set fold-changes
names(ngs$gset.meta$meta)
mx <- ngs$gset.meta$meta
S <- sapply( mx, function(x) unclass(x$fc)[,"gsva"])
S <- S / max(abs(S),na.rm=TRUE)
dim(S)
S <- S[match(rownames(G),rownames(S)),,drop=FALSE]
rownames(S) <- rownames(G)
S[is.na(S)] <- 0
S <- S / max(abs(S),na.rm=TRUE)
dim(S)


## ------- create rank vector
colnames(F)
k=1
fc <- F[,k]
rnk <- S[,k]
rnk <- rnk[which(rnk!=0)]

top.genes <- head(names(fc)[order(-abs(fc))],4000)
head(top.genes)
member.sets <- apply(G[,top.genes],2,function(x) rownames(G)[which(x!=0)])


## ------- run GSEA
library(fgsea)
res <- fgsea(member.sets, rnk, nperm=10000, minSize = 10, maxSize = Inf)
res <- as.data.frame(res)
rownames(res) <- res$pathway

res <- res[order(-abs(res$NES)),]
head(res[,1:5],20)

res <- res[order(res$padj),]
head(res[,1:5],20)

top.genes2 <- head(res$pathway,50)
head(top.genes2)

combined.score <- fc[rownames(res)] * abs(res$NES)
##combined.score <- fc[rownames(res)] * res$NES
top.genes2 <- names(combined.score)[order(-abs(combined.score))]
head(top.genes2)

par(mfrow=c(5,5), mar=c(2,2,3,1))
for(gene in head(top.genes2,25)) {
    gset = rownames(G)[which(G[,gene]!=0)]
    gsea.enplot(rnk, gset, main=gene, cex.main=1.2)
    old.rank <- rank(-abs(fc))[gene]
    new.rank <- rank(-abs(combined.score))[gene]
    qv <- res[gene,"padj"]
    tt <- c(paste("new.rank=",new.rank),
            paste("old.rank=",old.rank),
            paste("q=",qv) )
    legend("topright",tt)
}


## ------- run GSEA

xnes <- res$NES
names(xnes) <- rownames(res)
xfc <- F[names(xnes),1]

par(mfrow=c(1,1), mar=c(4,4,2,2))
plot( xfc, xnes, pch=".", cex=2)
abline(h=0, v=0, lty=1, lwd=0.4)

top.fc <- head(names(fc[order(-abs(fc))]),20)
points( xfc[top.fc], xnes[top.fc], pch=20, cex=0.3, col="red")
text( xfc[top.fc], xnes[top.fc], labels=top.fc, cex=0.8, col="red")

##======================================================================
##======================================================================
##======================================================================

dim(S)
term="immune"
term="growth"
term="metabol"
term="virus"
my.gset <- rownames(S)[grep(term,rownames(S),ignore.case=TRUE)]
table(my.gset %in% names(rnk))
gsea.enplot(rnk, my.gset, main=term, cex.main=1.2)

head(rownames(S))

title.words <- sapply(rownames(S), strsplit, split="[_: -]")
title.words <- sapply(title.words, function(s) tolower(gsub("\\(|\\)","",s)))
word.freq <- table(unlist(title.words))
tail(sort(word.freq),50)
length(word.freq)

getWordPairs <- function(w) apply(rbind( w, c(w[-1],"")),2,paste,collapse=" ")
title.pairs <- lapply(title.words, getWordPairs)
pair.freq <- table(unlist(title.pairs))
tail(sort(pair.freq),50)
length(pair.freq)

##ww <- names(word.freq[word.freq>=15])
ww <- names(pair.freq[pair.freq>=15])
length(ww)
##term.sets <- mclapply( ww[], function(w) names(title.words)[which(sapply(title.words,function(a)(w %in% a)))])
term.sets <- mclapply( ww[], function(w) names(title.pairs)[which(sapply(title.pairs,function(a)(w %in% a)))])
names(term.sets) <- ww[1:length(term.sets)]
names(term.sets)

head(term.sets)
head(rnk)

res <- fgsea(term.sets, rnk, nperm=10000, minSize = 10, maxSize = Inf)
res <- as.data.frame(res)
rownames(res) <- res$pathway

res <- res[order(-abs(res$NES)),]
head(res[,1:5],40)

term="apoptotic signaling"
term="rna polymerase"
term="cell cycle"
my.gset <- term.sets[[term]]
gsea.enplot(rnk, my.gset, main=term, cex.main=1.2)


res <- res[order(res$pval),]
head(res[,1:5],20)



top.genes2 <- head(res$pathway,50)
head(top.genes2)
