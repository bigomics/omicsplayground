source("../R/gx-heatmap.r")
source("../R/gx-combat.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/gx-b3plot.r")


load("../pgx/GSE72056-melanoma-scRNA-vsclusters-s200-4k-LT.pgx")

X <- head(ngs$X,200)
ct <- ngs$Y$cell.type
gx.splitmap(X, splitx=ct, show_colnames=FALSE,
            annot.ht=2.5, column_title_rot=90,
            col.annot=ngs$Y)


install.packages("ContourFunctions", version = "3.8")
##BiocManager::install("openCyto", version = "3.8")
library(openCyto)
library(ContourFunctions)
library(MASS)
library(RColorBrewer)


## some pretty colors
k <- 11
my.cols <- rev(brewer.pal(k, "RdYlBu"))
gene1 = "CD8A"
gene2 = "CD4"
x1 <- ngs$X[gene1,]
x2 <- ngs$X[gene2,]
m1 <- mean(x1)
m2 <- mean(x2)

z0 <- kde2d( x1[], x2[], n=50)
j1 <- which(x1 < m1 & x2 > m2)
z1 <- kde2d( x1[j1], x2[j1], n=50)
j2 <- which(x1 > m1 & x2 < m2)
z2 <- kde2d( x1[j2], x2[j2], n=50)
j3 <- which(x1 > m1 & x2 > m2)
z3 <- kde2d( x1[j3], x2[j3], n=50)

par(mfrow=c(1,1))
plot(x1, x2, xlab=gene1, ylab=gene2, col="grey70", pch=19, cex=1)
abline(h=mean(x1), v=mean(x2), lwd=1)
##filled.contour(z0, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
contour(z1, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
contour(z2, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)
contour(z3, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=2)

N = length(x1)
d1 = 0.02*max(x1)
d2 = 0.04*max(x2)
text( max(x1), max(x2), paste(round(100*length(j3)/N,2),"%"), pos=2, cex=2, col="gray60")
text( m1-d1, max(x2), paste(round(100*length(j1)/N,2),"%"), pos=2, cex=2, col="gray60")
text( max(x1), m2 - d2, paste(round(100*length(j2)/N,2),"%"), pos=2, cex=2, col="gray60")
text( m1-d1, m2 - d2, paste(round(100*length(j2)/N,2),"%"), pos=2, cex=2, col="gray60")

pgx.cytoPlot(ngs, gene1="CD4", gene2="CD8A")
pgx.cytoPlot(ngs, gene1="CD8A", gene2="CD4")

require(lattice)
contourplot(z3)


legend("topright", paste("R=", round(cor(x1,x2),2)), bty="n")








x <- ngs$X[1,]
y <- ngs$samples$group
par(mar=c(8,4,4,2))
gx.b3plot( x, y, srt=20, xlab="hello" )
mtext("xlab", side=1, line=5)

##install.packages("sinaplot")
library(sinaplot)
par(mfrow=c(3,1))
boxplot( x ~ y )
jj <- which( y %in% names(which(table(y)>3)))
sinaplot( x[jj] ~ y[jj] )
beeswarm( x ~ y, cex=0.3, pch="o" )
vioplot( x[jj], y[jj] )

