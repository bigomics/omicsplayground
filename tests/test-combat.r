library(tidyimpute)
library(na.tools)
setwd("~/Projects/Immunomics/tests/")

source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/pgx-functions.R")
source("../R/pgx-graph.R")
source("../R/ngs-fit.r")
source("../R/ngs-cook.r")

load("../pgx/GSE10846-dlbcl-mRNA-8k.pgx",verbose=1)
ngs$samples <- tidy.dataframe(ngs$samples)
ngs$samples <- impute(ngs$samples, .na=na.mean)
ngs$samples$gender[is.na(ngs$samples$gender)] <- "male"
    
##load("../pgx/hpbergen2019-fabry-ASIS-4k.pgx")

X=ngs$X
contr="dlbcl.type:ABC_vs_GCB"
computeNumSig <- function(ngs, X, contr=NULL) {
    ##group <- ngs$samples$group
    ##design <- model.matrix(~ 0+group)
    design <- ngs$model.parameters$design
    contr.matrix <- ngs$model.parameters$contr.matrix
    if(is.null(contr)) contr <- colnames(contr.matrix)
    contr.matrix <- contr.matrix[,contr,drop=FALSE]
    res <- ngs.fitContrastsWithLIMMA(
        X, is.logcpm=TRUE,
        contr.matrix, design,
        quantile.normalize=FALSE,
        method="limma", trend=TRUE,
        conform.output=FALSE, plot=FALSE)
    names(res)
    names(res$tables)
    fc <- sapply(res$tables,function(x) x$logFC)
    qv <- sapply(res$tables,function(x) x$adj.P.Val)
    ##numsig <- sum( abs(fc) > 1 & qv < 0.05)
    numsig <- mean(colSums(qv < 0.05))
    numsig
    return(numsig)
}


library(limma)
library(preprocessCore)
##qX <- normalize.quantiles(ngs$X)
qX <- normalize.quantiles(ngs$X)
dimnames(qX) <- dimnames(ngs$X)
bX <- removeBatchEffect(ngs$X, batch=ngs$samples$Chemotherapy)

computeNumSig(ngs, X=ngs$X, contr="dlbcl.type:ABC_vs_GCB")
computeNumSig(ngs, X=qX, contr="dlbcl.type:ABC_vs_GCB")
computeNumSig(ngs, X=bX, contr="dlbcl.type:ABC_vs_GCB")

computeNumSig(ngs, X=ngs$X, contr=NULL)
computeNumSig(ngs, X=qX, contr=NULL)

gx.heatmap(ngs$X, nmax=500, softmax=TRUE, col.annot=ngs$Y)
gx.heatmap(qX, nmax=500, softmax=TRUE, col.annot=ngs$Y)


require(sva)
nuispar1 <- c("patientID.211218","Volum.ul","Storage",
             "Batch.number","Sequence.facility",
             "numtreshhold.210119.cat","percent.mapped.210119",
             ".libsize",".cell_cycle")
nuispar1 <- c("gender","age","LDH.ratio","Chemotherapy","OS.status","OS.years",
              ".cell_cycle")
nuispar2 <- unlist(apply(combn(nuispar1,2),2,function(x) list(x)),recursive=FALSE)
nuispar3 <- unlist(apply(combn(nuispar1,3),2,function(x) list(x)),recursive=FALSE)
nuispar4 <- unlist(apply(combn(nuispar1,4),2,function(x) list(x)),recursive=FALSE)
nuispar5 <- unlist(apply(combn(nuispar1,5),2,function(x) list(x)),recursive=FALSE)
nuispar6 <- unlist(apply(combn(nuispar1,6),2,function(x) list(x)),recursive=FALSE)
nuispar <- c(sapply(nuispar1, list), nuispar2, nuispar3, nuispar4)
##nuispar <- c( "", sapply(nuispar1, list), nuispar2, nuispar3, nuispar4, nuispar5, nuispar6)
length(nuispar)

nuispar <- c( "no_correction", nuispar)
##nuispar <- nuispar[1:4]

aX <- ngs$X
##aX <- normalizeQuantiles(ngs$X)

normz <- c("nonorm","cpm","TMM","RLE","quantile","SVA")
contr=NULL
runComputeNumSig <- function(contr) {
    k="cpm"
    numsig <- c()
    for(k in normz) {
        aX <- NULL
        if(k=="nonorm") aX <- log(1+ngs$counts)
        if(k=="cpm") aX <- edgeR::cpm(ngs$counts, log=TRUE)
        if(k=="TMM") aX <- log2(1+normalizeTMM(ngs$counts))
        if(k=="RLE") aX <- log2(1+normalizeRLE(ngs$counts))
        if(k=="quantile") aX <- normalizeQuantiles(log2(1+ngs$counts))
        if(k=="SVA") {
            require(sva)
            mod1 = model.matrix( ~ group, data=ngs$samples)
            mod0 = cbind(mod1[,1])
            logcpm <- edgeR::cpm(ngs$counts, log=TRUE)
            sv = sva( logcpm, mod1, mod0, n.sv=NULL)$sv
            aX <- removeBatchEffect( logcpm, covariates=sv, design=mod1)
            dim(aX)
        }

        ## solve for all combinations
        pp <- nuispar[[10]]
        for(pp in nuispar) {
            bX = aX
            Y <- ngs$samples[colnames(bX),]
            if(pp[1]!="no_correction") {
                model.formula <- formula(paste("~ 0 + ",paste(pp,collapse=" + ")))
                bvar <- model.matrix(model.formula, data=Y)
                design <- model.matrix( ~ ngs$samples$group)                
                bX <- removeBatchEffect(bX, covariates=bvar, design=design)
            }
            pp1 <- paste0(k,":",paste(pp,collapse="+"))
            pp1
            numsig[pp1] <- computeNumSig(ngs, X=bX, contr=contr)
        }
    }
    return(numsig)
}

require(parallel)
contr=NULL
contr="dlbcl.type:ABC_vs_GCB"
##numsig <- mclapply(1:3, function(i) runComputeNumSig(contr))
numsig <- runComputeNumSig(contr)
tail(sort(unlist(numsig)),200)

save(numsig, file="numsig-results-dlbcl-it10.rda")

load("numsig-results-dlbcl-it10.rda")
tail(sort(unlist(numsig)),20)


source("../R/pgx-correct.R")
batch <- c("gender","age","LDH.ratio","Chemotherapy","OS.status")
contrast="dlbcl.type:ABC_vs_GCB"
res <- pgx.optimizeBatchCorrection(ngs, contrast, batch=batch,
                                   nparam=NULL, niter=1)
tail(sort(res),200)

##================================================================================
##============================ VISUALIZE =========================================
##================================================================================

sv <- function(k) {
    res <- svd(X)
    if(length(k)>1) {
        x1 <- res$u[,k,drop=FALSE] %*% diag(res$d[k]) %*% t(res$v[,k,drop=FALSE])
    } else {
        x1 <- res$u[,k,drop=FALSE] %*% t(res$v[,k,drop=FALSE])
    }
    dim(x1)
    rownames(x1) <- rownames(X)
    colnames(x1) <- colnames(X)
    x1 <- x1[order(-apply(x1,1,sd)),]
    x1
}

xheatmap <- function(x, main="") {
    gx.heatmap(x, col.annot=annot, annot.ht=0.8, softmax=1,
               nmax=50, keysize=1.0, key=FALSE, main=main)
}

res <- svd(X)
dim(res$u)
idx <- max.col(abs(res$u[,1:10]))
table(idx)

a1 <- annot

## ----------- show PC components with annotation
sv.top <- lapply(1:5,function(i) head(sv(i),15))
for(i in 1:length(sv.top)) {
    rownames(sv.top[[i]]) <- paste0("PC",i,":",rownames(sv.top[[i]]))
}

sv.top1 <- do.call(rbind, sv.top)
svx <- sub(":.*","",rownames(sv.top1))
table(svx)
gx.splitmap(sv.top1, split=svx, col.annot=a1)

sv.gene <- rownames(sv.top1)
sv.gene <- sub(".*:","",sv.gene)
sv.top2 <- X[sv.gene,]
gx.splitmap( sv.top2, split=svx, col.annot=a1, softmax=1)

## -------- show correlation between PC and annotation
colnames(res$v) <- paste0("PC",1:ncol(res$v))
m1 <- apply(a1, 2, function(x) model.matrix(~0+x))
m1 <- lapply(m1, function(m) {colnames(m)=sub("^x","",colnames(m));m})

rho <- sapply(m1, function(m) rowMeans(abs(cor(res$v[,1:5], m))))
heatmap( t(rho), scale="none", mar=c(10,10))

rho1 <- sapply(m1, function(m) cor(res$v[,1:5], m))
rho1 <- do.call(cbind, rho1)
##heatmap( t(rho1), scale="none", mar=c(10,10))
heatmap( t(abs(rho1)), scale="none", mar=c(10,10))


par(mfrow=c(1,2),mar=c(8,4,4,2))
klr1=grey.colors(5)
barplot( t(rho), scale="none", mar=c(10,10), col=klr1, las=3)
legend("topright", legend=colnames(rho), cex=0.8, y.intersp=0.8,fill=klr1)

barplot( rho, scale="none", mar=c(10,10), las=3, col=klr1)
legend("topright", legend=rownames(rho), cex=0.8, y.intersp=0.8,fill=klr1)

## calculate batch effect p-values
pc.rho <- function(m) apply(res$v[,1:5], 2, function(v) cor.test(m,v)$p.value)
pv.rho <- sapply(m1, function(x) apply(x,2,pc.rho))
pv.rho

pv.rho1 <- do.call(cbind, pv.rho)
heatmap( t(pv.rho1), scale="none", mar=c(10,10))
heatmap( t(-log10(pv.rho1))**0.2, scale="none", mar=c(10,10))


##================================================================================
##=========================== CORRECTION =========================================
##================================================================================


require(sva)
group <- as.character(ngs$samples$dlbcl.type)
batch <- as.character(ngs$samples$Chemotherapy)

bX <- ComBat(ngs$X, batch=batch, mod=NULL)
xheatmap(bX)

cX <- ComBat(ngs$X, batch=batch, mod=design)
xheatmap(cX)

require(limma)
bX <- limma::removeBatchEffect(ngs$X, batch=batch)
xheatmap(bX)

cX <- limma::removeBatchEffect(ngs$X, batch=batch, design=design)
xheatmap(cX)


require(sva)
mod1 = model.matrix( ~ group, data=ngs$samples)
mod0 = cbind(mod1[,1])
##mod0 = model.matrix( ~ 1, data=ngs$samples)

sv = sva(ngs$X, mod1, mod0, n.sv=NULL)$sv
dX <- removeBatchEffect( ngs$X, covariates=svseq, design=mod1)
xheatmap(dX)

sv = sva(ngs$X, mod1, mod0=NULL, n.sv=NULL)$sv
dX <- removeBatchEffect( ngs$X, covariates=svseq, design=mod1)
xheatmap(dX)


par(mfrow=c(2,2))
for(i in 1:4) plot(svseq[,i],pch=19,col="blue")



library(BatchQC)
library(bladderbatch)
data(bladderdata)
pheno <- pData(bladderEset)
edata <- exprs(bladderEset)
batch <- pheno$batch
condition <- pheno$cancer
batchQC(edata, batch=batch, condition=condition,
        report_file="batchqc_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)


library(BatchQC)
data(example_batchqc_data)
batch <- batch_indicator$V1
condition <- batch_indicator$V2
batchQC(signature_data, batch=batch, condition=condition,
        report_file="batchqc_signature_data_report.html", report_dir=".",
        report_option_binary="111111111",
        view_report=FALSE, interactive=TRUE)

