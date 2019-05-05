library(BiocManager)
##install("mixOmics")
require(mixOmics)

{
    load("../files/gmt-all.rda",verbose=1)
    load(file="../pgx/GSE10846-dlbcl-mRNA-8k-LT.pgx",verbose=1)
    Y = ngs$samples$dlbcl.type
    cs <- gmt.all[["COMPARTMENTS:Cell_surface"]]
    pk <- gmt.all[["MSIGDB:PROTEIN_KINASES"]]
    tf <- gmt.all[["CUSTOM:Transcription factors"]]
    gg <- rownames(ngs$X)
    cs.X <- (ngs$X[intersect(cs, gg),])
    pk.X <- (ngs$X[intersect(pk, gg),])
    tf.X <- (ngs$X[intersect(tf, gg),])
    data = list( CS=cs.X, PK=pk.X, TF=tf.X)
    length(data)
}

{
    require(mixOmics)
    data('breast.TCGA')
    Y = breast.TCGA$data.train$subtype    
    table(Y)
    data = list(mrna =  t(breast.TCGA$data.train$mrna),
                mirna =  t(breast.TCGA$data.train$mirna),
                prot =  t(breast.TCGA$data.train$protein))
    length(data)
}

{    
    library(Biobase)
    ##BiocManager::install("PharmacoGx")
    library(PharmacoGx)
    data("CCLEsmall")
    
    availablePSets()
    ##CCLE <- downloadPSet("CCLE") ## 245.5 MB!!    
    pset = CCLEsmall
    pset = CCLE
    
    genes <- fNames(pset, "rna")
    slotNames(pset)
    names(pset@molecularProfiles)
    ic50 <- summarizeSensitivityProfiles(
        pSet=pset, sensitivity.measure='ic50_published',
        summary.stat="median", verbose=FALSE)
    auc <- summarizeSensitivityProfiles(
        pSet=pset, sensitivity.measure='auc_published',
        summary.stat="median", verbose=FALSE)

    cnv <- summarizeMolecularProfiles(
        pset, cellNames(pset),
        mDataType="cnv", verbose=FALSE)
    mut <- summarizeMolecularProfiles(
        pset, cellNames(pset),
        mDataType="mutation", summary.stat="and",
        verbose=FALSE)
    rna <- summarizeMolecularProfiles(
        pset, cellNames(pset),
        mDataType="rnaseq", verbose=FALSE)
    table( colnames(ic50)==colnames(cnv))
    cnv <- as.matrix(cnv)
    mut <- 1*(as.matrix(mut)=="1")
    rna <- as.matrix(rna)
    mut[is.na(mut)] <- 0
    mut <- mut[rowSums(mut)>=3,]
    ##mut <- mut + 1e-3*matrix(rnorm(length(mut)),nrow(mut),ncol(mut))
    dim(mut)
    
    head(rownames(rna))
    if( grepl("ENSG",rownames(rna)[1])) {
        library(org.Hs.eg.db)
        symbols <- mapIds(org.Hs.eg.db, keys=rownames(rna),
                          keytype="ENSEMBL", column="SYMBOL")
        jj <- which(!is.na(symbols))
        rna <- rna[jj,]
        rownames(rna) <- symbols[jj]
        rna <- rna[order(-apply(rna,1,sd)),]
        rna <- rna[!duplicated(rownames(rna)),]
    }
    
    rownames(ic50)
    ##drug = "Crizotinib"
    drug = "Erlotinib"
    drug = "lapatinib"
    drug = "PLX4720"
    Y <- ic50[drug,]
    Y <- auc[drug,]   
    names(Y) <- colnames(ic50)
    
    has.all <- ( colMeans(is.na(cnv)) < 0.05 &
                 colMeans(is.na(mut)) < 0.05 &
                 colMeans(is.na(rna)) < 0.05 )
    is.tissue <- grepl("skin", CCLE@cell$Site.Primary)
    is.tissue <- grepl("breast", CCLE@cell$Site.Primary)
    is.tissue <- grepl("lung", CCLE@cell$Site.Primary)    
    kk <- which(has.all & !is.na(Y) & is.tissue)
    kk <- which(has.all & !is.na(Y) )
    length(kk)
    data <- list( MUT=mut[,kk], CNV=cnv[,kk], RNA=rna[,kk])
    Y <- Y[colnames(data[[1]])]
    Y <- factor(c("NOT","RES")[1 + 1*(Y<median(Y))])
    lapply(data, dim)

    library(ggplot2)
    par(mfrow=c(2,2))
    boxplot( data$CNV["ERBB2",] ~ Y)
    points( as.integer(Y), data$CNV["ERBB2",], pch=20, cex=4)
    boxplot( data$RNA["ERBB2",] ~ Y)
    points( as.integer(Y), data$RNA["ERBB2",], pch=20, cex=4)
    
    medianImpute <- function(X) {
        mx <- apply(X,1,median,na.rm=TRUE)
        X0=X
        X0[is.na(X0)]=0
        X0 + is.na(X)*mx
    }
    data1 <- lapply(data, medianImpute)
    lapply(data, function(x) sum(is.na(x)))
    lapply(data1, function(x) sum(is.na(x)))
    
    boxplot( data$CNV["ERBB2",] ~ Y)
    
}

source("../R/gx-limma.r")
source("../R/pgx-proteomics.R")
source("../R/pgx-functions.R")
source("../R/pgx-predict.R")

R <- list()
methods=c("glmnet","randomforest","xgboost","pls")
i=1
par(mfrow=c(3,1), mar=c(12,4,2,2))
for(i in 1:3) {

    ##X = data2[[i]]
    X = data1[[i]]
    dim(X)
    
    ## prioritize with LIMMA
    res <- gx.limma(X, Y, fdr=1, lfc=0)
    pp <- head(rownames(res)[order(-abs(res[,"logFC"]))],100)
    y1 <- 1*(Y=="RES")
    length(pp)
    
    ## prioritize with Feature Importance
    xx <- pgx.variableImportance(X[pp,], y1, methods=methods)    
    nx <- t(t(xx) / (1e-8+apply(xx,2,max)))
    nx <- apply(xx,2,rank) / nrow(xx)
    ##jj <- order(-exp(rowMeans(log(1e-8+nx))))
    jj <- order(-rowMeans(nx))
    nx <- nx[jj,]
    if(1) {
        barplot( t(head(nx,60)), beside=FALSE, las=3,
                cex.names=2)
        title(names(data1)[i])
    }
    R[[i]] <- head(nx,100)
}
names(R) <- names(data1)
lapply(R,dim)

sel <- lapply(R, function(x) head(rownames(x),100))
data2 <- lapply(1:3, function(i) data1[[i]][sel[[i]],])
names(data2) <- names(data1)
lapply(data2, function(x) head(rownames(x)))

res <- pgx.makeTriSystemGraph(data2, Y, nfeat=50, numedge=100)    
dim(res$X[[1]])
dim(res$W[[1]])

plotDiablo(res, ncomp = 1)
cimDiablo(res, margins=c(10,20), transpose=TRUE)

par(mfrow=c(1,1))
##plotLoadings(res, comp = 1, contrib = 'max', method = 'median',
##             layout=NULL, size.legend=1.5, ndisplay=20,
##             size.name = 1.8, size.title = 1.2)	
mixPlotLoadings(res, cex=2, showloops=FALSE)
mixPlotLoadings(res, cex=2, showloops=TRUE)

mixHivePlot(res, ngs=NULL, ct=NULL, showloops=TRUE,
            cex=1.5, numlab=6)     
legend("topleft",cex=4, 
       colnames(res$W[[1]]), fill=brewer.pal(8,"Set2")[1:3],
       inset=c(-1.1,0),xpd=NA)

## plot samples
klrpal <- c("blue2","orange2")
klrpal <- brewer.pal(n=3, name = "RdBu")[c(3,1)]
par(mfrow=c(1,3),mar=c(4,4,4,2)*2)
klr1 <- klrpal[factor(Y)]
klr1[is.na(klr1)] <- "grey85"
pairs(sapply(res$variates,function(x) x[,1]),
      cex=4, pch=19,col=klr1)

source("../R/gx-heatmap.r")
##source("../R/gx-limma.r")
NFEAT=25
top.features <- lapply( res$loadings[1:3], function(x)
    head(names(sort(-abs(x[,1]))),NFEAT))
X <- do.call(cbind, res$X)[,unlist(top.features)]
dim(X)
X <- X[,head(order(-apply(X,2,sd)),60)]
ptype <- names(res$X)[as.integer(sub(":.*","",colnames(X)))]
##ptype1 <- c(cs="CELL.SURFACE", pk="PROTEIN.KINASE", tf="TRANSCRIPTION.FACTOR")[ptype]
gx.splitmap( t(X), split=ptype, splitx=NULL,
            scale="row", dist.method="pearson",
            cexRow=1.5, cexCol=1.5,
            col.dist.method="pearson")





##======================================================================
##======================================================================
##======================================================================

require(mixOmics)
data('breast.TCGA')
Y = breast.TCGA$data.train$subtype    
data = list(mrna =  breast.TCGA$data.train$mrna,
            mirna =  breast.TCGA$data.train$mirna,
            prot =  breast.TCGA$data.train$protein)
length(data)

## set number of component per data set
NCOMP = 3
list.keepX = list(mrna = rep(20, 3), mirna = rep(10,3), prot = rep(10,3))

## set up a full design where every block is connected 
design = matrix(1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) =  0
design      
res <- block.splsda(
    X = data, Y = Y, ncomp = NCOMP, keepX = list.keepX, design = design)


plotDiablo(res, ncomp = 1)
plotLoadings(res, comp = 2, contrib = 'max', method = 'median',
             size.name = 1, size.title = rel(0.7))	
cimDiablo(res)



plotIndiv(res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
plotVar(res, var.names = FALSE, ##
        ##style = '3d', pch = rep("sphere",3),
        style = 'graphics', pch = c(16, 17, 15),
        cex = c(1,1,1)*1.2, legend = TRUE, 
        col = c('darkorchid', 'brown1', 'lightgreen'))    
plotArrow(res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
circosPlot(res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

##  network(res, blocks = c(1,2,3), color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)




##======================================================================
##======================================================================
##======================================================================

library(mixOmics)
data(breast.TCGA)
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
summary(Y)


list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)
plotVar(MyResult.diablo)

MyResult.diablo2 <- block.plsda(X, Y)
plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2,3),
          title = 'BRCA with DIABLO')

plotVar(MyResult.diablo, var.names = c(FALSE, FALSE, TRUE),
        legend=TRUE, pch=c(16,16,1))

plotDiablo(MyResult.diablo, ncomp = 1)
circosPlot(MyResult.diablo, cutoff=0.7)

cimDiablo(MyResult.diablo,
          color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8,20), legend.position = "right")

plotLoadings(MyResult.diablo, comp = 2, contrib = "max")
network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')




set.seed(123) # for reproducibility in this vignette
MyPerf.diablo <- perf(MyResult.diablo, validation = 'Mfold', folds = 5, 
                   nrepeat = 10, 
                   dist = 'centroids.dist')







