source("../R/gx-heatmap.r")
source("../R/gx-combat.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")

load("../pgx/GSE72056-melanoma-scRNA-vsphenotype-TC8k.pgx")

kk <- which(!(ngs$samples$patient %in% c("CY59","CY78")))
Y <- ngs$samples[kk,]
p2rx7.ratio <- tapply(Y$P2RX7, as.character(Y$patient), function(s) sum(s=="pos") / length(s))
p2rx7.avg   <- tapply(ngs$X["P2RX7",kk], as.character(Y$patient), mean)

sort(p2rx7.ratio)
length(p2rx7.ratio)

dim(ngs$gsetX)
zX <- t(apply( ngs$gsetX[,kk], 1, function(x) tapply(x, as.character(Y$patient), mean, na.rm=TRUE)))
dim(zX)
sum(is.na(zX))

rho <- cor(t(zX), p2rx7.ratio, method="spearman")[,1]
head(sort(rho,decreasing=TRUE))
head(sort(rho))

par(mfrow=c(5,5), mar=c(4,4,2,2)*0.5)
jj <- head(order(-abs(rho)),25)
##jj <- head(order(-rho),12)
##jj <- head(order(+rho),25)
for(i in jj) {
    plot( p2rx7.ratio, zX[i,], pch=19, cex=0.5)
    tt <- substring(colnames(zX)[i],1,40)
    title(tt,cex.main=0.8)

    plot( p2rx7.avg, zX[i,], pch=19, cex=0.5)
    tt <- substring(colnames(zX)[i],1,40)
    title(tt,cex.main=0.8)
}

par(mfrow=c(3,3))
i=jj[1]
plot( p2rx7.ratio, zX[i,], pch=19, cex=0.9)

gx.p2rx7 <- ngs$X["P2RX7",]
klr <- rep(rainbow(16),99)[ngs$samples$patient]
plot( gx.p2rx7, ngs$gsetX[i,], pch=19, cex=0.9, col=klr)




gmt <- readRDS("../files/gmt-all.rds")
grep("0007049", names(gmt), value=TRUE)
grep("cell cycle", names(gmt), value=TRUE,ignore.case=TRUE)
gmt.cc <- gmt[grep("cell cycle", names(gmt),ignore.case=TRUE)]
genes.cc <- head(names(sort(table(unlist(gmt.cc)),decreasing=TRUE)),200)

X <- ngs$X
dim(X)
y <- sampleTable[colnames(X),]$P2RX7
genes.cc1 <- intersect(genes.cc, rownames(X))
ccX <- X[genes.cc1,]
annot <- ngs$samples[colnames(X),]

par(mfrow=c(1,1))
gx.heatmap(ccX, col.annot=annot)

jj <- intersect(names(gmt.cc),rownames(ngs$gsetX))
ccZ <- ngs$gsetX[jj,]
dim(ccZ)

rho0 <- cor( t(ccZ), X["P2RX7",] )[,1]
head(sort(rho0,decreasing=TRUE),10)
head(sort(rho0,decreasing=FALSE),10)

rho0 <- cor( t(ngs$gsetX), X["P2RX7",] )[,1]
head(sort(rho0,decreasing=TRUE),20)
head(sort(rho0,decreasing=FALSE),20)


##-------------------------------------------------------------------
## NNM magic...
##-------------------------------------------------------------------

X <- log2(10 + ngs$counts)
dim(X)
y <- sampleTable[colnames(X),]$P2RX7
cX <- gx.nnmcorrect(X[,], y[], k=5)

annot <- sampleTable[colnames(cX),]
gx.heatmap(X, col.annot=annot)
gx.heatmap(cX, col.annot=annot)

rho1 <- cor( t(X), X["P2RX7",] )[,1]
head(sort(rho1,decreasing=TRUE),40)

rho2 <- cor( t(cX), cX["P2RX7",] )[,1]
head(sort(rho2,decreasing=TRUE),100)


##================================================================================
## ccRemover
##
## See vignette: https://cran.r-project.org/web/packages/ccRemover/vignettes/ccRemover_tutorial.html
##
##================================================================================

##install.packages("ccRemover")
library(ccRemover)
X <- log2(10 + ngs$counts)
meanx <- rowMeans(X)
X <- X - meanx

gene_names <- rownames(X)
cellcycle_idx <- gene_indexer(
    gene_names, species = "human", name_type = "symbol" )
gene_names[cellcycle_idx]

if_cc <- (1:length(gene_names) %in% cellcycle_idx)
dat <- list(x=X, if_cc=if_cc)
xhat <- ccRemover(dat, bar=FALSE)
xhat <- xhat + meanx


##================================================================================
## scLVM
##
## See vignette: scLVM/R/tutorials/scLVM_vignette.html
##
##================================================================================

library(genefilter)
library(statmod)
require(ggplot2)
library(gplots)
require(DESeq2)
library(scLVM)

##limix_path = '/Users/flo/software/limix-master/build/release.darwin/interfaces/python'
##configLimix(limix_path)

data(data_Tcells)
help(data_Tcells)

dataMouse[ 1:5, 1:4 ]
geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(dataMouse), 1, 4 ) ] )

#2. calculate normalisation for counts
countsMmus <- dataMouse[ which( geneTypes=="ENSM" ), ]
countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), ]
lengthsMmus <- dataMouse[ which( geneTypes=="ENSM" ), 1 ]
lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]

sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMmus <- sfERCC #also use ERCC size factor for endogenous genes

#normalise read counts
eps = 1e-3
eps = 0
nCountsERCC <- t( t(eps + countsERCC) / sfERCC )
nCountsMmus <- t( t(eps + countsMmus) / sfERCC )
techNoise = fitTechnicalNoise(nCountsMmus, nCountsERCC=nCountsERCC, fit_type = 'counts')
techNoise

##techNoiseLogFit = fitTechnicalNoise(nCountsMmus, fit_type = 'log', use_ERCC = FALSE, plot=FALSE)
##techNoiseLogVarFit = fitTechnicalNoise(nCountsMmus, fit_type = 'logvar', use_ERCC = FALSE, plot=FALSE)

is_het = getVariableGenes(nCountsMmus, techNoise$fit, method = "fdr",
                          threshold = 0.1, fit_type="counts",
                          sfEndo=sfMmus, sfERCC=sfERCC)

##is_hetLog = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, plot=TRUE)
##is_hetLogVar = getVariableGenes(nCountsMmus, techNoiseLogVarFit$fit, plot=TRUE)


##In order to fit the latent cell cycle factor use genes annotated in GO (term GO:0007049).
ens_ids_cc <- getEnsembl('GO:0007049')

Y = t(log10(nCountsMmus+1)) #normalised trandformed read counts
genes_het_bool = as.vector(is_het) #variable genes
geneID = rownames(nCountsMmus) #gene IDs
tech_noise = as.vector(techNoise$techNoiseLog) #technical noise


sclvm = new("scLVM")
sclvm = init(sclvm, Y=Y, tech_noise = tech_noise)
ens_ids_cc <- getEnsembl('GO:0007049')
CellCycleARD = fitFactor(sclvm, geneSet=ens_ids_cc, k=20,use_ard = TRUE)

plot(seq(1, length(CellCycleARD$X_ard)), CellCycleARD$X_ard, xlab = '# Factor', ylab = 'Variance explained')
title('Variance explained by latent factors')


CellCycle = fitFactor(sclvm,geneSet = ens_ids_cc,k=1)
#Get cell-cycle factor
Kcc = CellCycle$K
Xcc = CellCycle$X

#Plot inferred similarity matrix
image(Kcc,xaxt = "n", yaxt = "n", xlab = 'cells', ylab = 'cells')
title('Similarity matrix based on cell cycle')




