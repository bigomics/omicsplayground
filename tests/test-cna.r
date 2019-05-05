source("../R/gx-util.r")
source("../R/gset-gsea.r")
source("../R/gset-fisher.r")

load(file="../pgx/GSE10846-dlbcl-mRNA-8k.pgx",verbose=1)
load("../files/gmt-all.rda",verbose=1)

gset <- gmt.all[["H:HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]
gset
jj <- grep("HALLMARK",names(gmt.all))
jj <- sample(length(gmt.all),600)
gmt <- gmt.all[jj]
G <- sapply(gmt.all[jj], function(s) 1*(rownames(ngs$X) %in% s))
dim(G)

X <- ngs$X[,1:500]
cX <- X - rowMeans(X)
system.time( z1 <- GSVA::gsva(X, gmt, parallel.sz=20 ) )
system.time( z1c <- GSVA::gsva(cX , gmt, parallel.sz=20 ) )
system.time( z2 <- GSVA::gsva(X, gmt, parallel.sz=20, method="ssgsea" ))
system.time( z2c <- GSVA::gsva( cX, gmt, parallel.sz=20, method="ssgsea" ))
dim(z1)

ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)[,1]
system.time( z3 <- corSparse( G, X))
system.time( z4 <- cosSparse( G, X))
system.time( z5 <- corSparse( G, apply(X,2,rank) ) )
system.time( z5c <- corSparse( G, apply(cX,2,rank) ) )

z5[is.na(z5)] <- 0
z5c[is.na(z5c)] <- 0

dimnames(z3) <- list(colnames(G),colnames(X))
dimnames(z4) <- dimnames(z5) <- dimnames(z5c) <- dimnames(z3)

require(affy)
z6  <- normalizeQuantiles(z5)
z7  <- t(normalizeQuantiles(t(z5)))
z7d  <- t(scale(t(z5)))
dimnames(z6) <- dimnames(z7) <- dimnames(z7d) <- dimnames(z3)

system.time({
    cells_rankings <- AUCell_buildRankings(X, nCores=24, plotStats=TRUE)
    z8 <- assay(AUCell_calcAUC(gmt, cells_rankings))
})
cells_rankings <- AUCell_buildRankings(cX, nCores=24, plotStats=TRUE)
z8c <- assay(AUCell_calcAUC(gmt, cells_rankings))


pdf("test-SSEmethods-xSAMPLES.pdf",w=10,h=10)
par(mfrow=c(1,1))
i=1
kk = intersect(intersect(rownames(z1),rownames(z3)),rownames(z8))
for(i in 1:10) {
    k <- kk[i]
    Z <- cbind( GSVA=z1[k,], cGSVA=z1c[k,],
               SSGSEA=z2[k,], cSSGSEA=z2c[k,],
               COR=z3[k,], COS=z4[k,],
               RNKCOR=z5[k,], cRNKCOR=z5c[k,],
               qn.RNKCOR=z6[k,], tqn.RNKCOR=z7[k,], ts.RNKCOR=z7d[k,],
               AUCell=(z8[k,]), logAUCell=log(z8[k,]) )
    pairs(Z, main=k, pch=19, cex=0.1)
}
dev.off()


pdf("test-SSEmethods-xGSETS.pdf",w=10,h=10)
s <- intersect(intersect(rownames(z1),rownames(z3)),rownames(z8))
kk = intersect(intersect(colnames(z1),colnames(z3)),colnames(z8))
par(mfrow=c(1,1))
i=10
for(i in 1:10) {
    k=kk[i]
    Z <- cbind( GSVA=z1[s,k], cGSVA=z1c[s,k],
               SSGSEA=z2[s,k], cSSGSEA=z2c[s,k],
               COR=z3[s,k], COS=z4[s,k],
               RNKCOR=z5[s,k], cRNKCOR=z5c[s,k],
               qn.RNKCOR=z6[s,k], tqn.cRNKCOR=z7[s,k], ts.cRNKCOR=z7d[s,k],
               AUCell=z8[s,k], logAUCell=log(z8[s,k]) )
    pairs(Z, main=k, pch=19, cex=0.1)
}
dev.off()


##--------- AUCell
X <- ngs$X[,]
dim(X)

system.time( z5 <- corSparse( G, apply(X-rowMeans(X),2,rank) ) )

system.time({
    cells_rankings <- AUCell_buildRankings(X, nCores=24, plotStats=TRUE)
    cells_AUC <- AUCell_calcAUC(gmt, cells_rankings)
})

z8 <- assay(cells_AUC)[1,]
Z <- cbind( RCOR=z5[1,], AUCell=z8, log.AUCell=log(z8) )
pairs(Z)


Z <- cbind( RCOR=z5[,1], AUCell=assay(cells_AUC)[,1])
pairs(Z)


