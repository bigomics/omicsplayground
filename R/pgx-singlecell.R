##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


pgx.poolCells <- function(counts, pos, npools, sample.id, cluster) {
    
    sample.id <- as.character(sample.id)
    cluster <- as.character(cluster)
    
    ## create subclustering for pooling
    pool.id <- rep(NA,ncol(counts))
    names(pool.id) <- colnames(counts)
    i=1
    pool.groups <- paste0("s",sample.id,"-c",cluster)
    ##pool.groups <- as.character(clust$Cluster)
    unique(pool.groups)
    table(pool.groups)
    
    nsamples <- length(unique(sample.id))
    ##npools   <- 2000
    clustsizes <- table(cluster)
    rel.clustsizes <- clustsizes / sum(clustsizes)
    totalcounts <- tapply(Matrix::colSums(counts),sample.id,sum)
    rel.totalcounts <- totalcounts / sum(totalcounts)
    
    message("[pgx.poolCells] detecting pools...")
    u <- pool.groups[1]
    for(u in unique(pool.groups)) {
        jj <- which(pool.groups==u)
        i <- as.integer(gsub("s|-c.*","",u))
        j <- as.character(gsub(".*-c","",u))        
        nc <- npools * rel.totalcounts[i] * rel.clustsizes[j] 
        nc <- as.integer(ceiling(nc))
        nc
        if(length(jj) > nc && nc>1) {
            pos.jj <- pos[jj,] + 1e-6*matrix(rnorm(length(pos[jj,])),nrow=length(jj))
            K <- kmeans( pos.jj, centers=nc, iter.max=1000 )
            table(K$cluster)
            pool.id[jj] <- paste0(u,"-p",K$cluster)
        } else {
            pool.id[jj] <- paste0(u,"-p",1:length(jj))
        }
    }
    
    length(pool.id)
    length(table(pool.id))
    
    message("[pgx.poolCells] pooling cells...")
    pool.counts <- tapply(1:ncol(counts), pool.id, function(k) Matrix::rowSums(counts[,k,drop=FALSE]))
    pool.counts <- do.call(cbind, pool.counts)
    dim(pool.counts)
    
    res <- list()
    res$pool.id <- pool.id
    res$pool.counts <- pool.counts
    
    return(res)
}


## qc.filter=FALSE;filter.a=2.5;sct=FALSE
## qc.filter=FALSE;filter.a=2.5;sct=TRUE
pgx.SeuratBatchCorrect <- function(counts, batch, qc.filter=FALSE, sct=FALSE) {
    
    ##
    ## From Seurat vignette: Integration/batch correction. Note there
    ## is no QC filtering for samples on ribo/mito content. You need
    ## to do that before.
    ##
    library(Seurat)
    nbatch <- length(unique(batch))
    message("[pgx.SeuratIntegration] Processing ",nbatch," batches...")        
    obj.list <- list()
    i=1
    b=batch[1]
    batches <- unique(batch)
    batches
    for(i in 1:length(batches)) {
        sel <- which(batch==batches[i])
        ##sel <- head(sel,200)
        counts1 <- counts[,sel]
        obj <- CreateSeuratObject(counts1)
        if(sct) {
            obj <- SCTransform(obj, vars.to.regress = NULL, verbose = FALSE)
            ##obj <- SCTransform(obj)
        } else {
            obj <- NormalizeData(obj, normalization.method="LogNormalize",
                                 scale.factor=10000, verbose=FALSE)
            obj <- FindVariableFeatures(obj, selection.method="vst", verbose=FALSE)
        }
        obj$batch <- b
        obj.list[[i]] <- obj
    }
    
    anchor.features=NULL
    sct
    if(sct) {
        ## See: https://satijalab.org/seurat/v3.0/integration.html
        anchor.features <- SelectIntegrationFeatures(
            object.list = obj.list, nfeatures = 3000)
        obj.list <- PrepSCTIntegration(
            object.list = obj.list,
            anchor.features = anchor.features, 
            verbose = FALSE)        
    } else {
        anchor.features=2000
    }
    
    message("[pgx.SeuratIntegration] Finding anchors...")
    options(future.globals.maxSize = 8*1024^3) ## set to 8GB
    ##NUM.CC=10;k.filter=10
    NUM.CC = max(min(20, min(table(batch))-1),1)
    NUM.CC
    ndims <- sapply(obj.list,ncol)
    mindim <- max(min(ndims)-1,1)
    mindim
    normalization.method = ifelse(sct,"SCT","LogNormalize")
    message("[pgx.SeuratIntegration] NUM.CC = ",NUM.CC)
    message("[pgx.SeuratIntegration] normalization.method = ",normalization.method)
    anchors <- FindIntegrationAnchors(
        obj.list, dims = 1:NUM.CC,
        k.filter = min(200,mindim),
        k.anchor = min(5,mindim),
        k.score = min(30,mindim),
        anchor.features = anchor.features,
        normalization.method = normalization.method,
        verbose = FALSE)
    
    message("[pgx.SeuratIntegration] Integrating data...")    
    integrated <- IntegrateData(
        anchorset = anchors,
        ##k.weight = min(100,mindim),
        dims = 1:NUM.CC,
        normalization.method = normalization.method,        
        verbose = FALSE)
    dim(integrated)
    
    key <- ifelse(sct, "SCT", "integrated")
    integrated.matrix <- as.matrix(integrated[[key]]@data)
    integrated.matrix <- exp(integrated.matrix)-1 ## natural log
    dim(integrated.matrix)
    
    return(integrated.matrix)    
}


pgx.SeuratFiltering <- function(counts, a=2.5)
{
    ## QC filter of (single) cells like Seurat
    ## See https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

    ##--------------------------------------------------------------
    ## calculate percentages
    ##--------------------------------------------------------------        
    mt.genes <- grep("^MT-",rownames(counts),ignore.case=TRUE,value=TRUE)
    rb.genes <- grep("^RP[SL]",rownames(counts),ignore.case=TRUE,value=TRUE)
    percent.mito <- colSums(counts[mt.genes,])/Matrix::colSums(counts)*100
    percent.ribo <- colSums(counts[rb.genes,])/Matrix::colSums(counts)*100
    nfeature <- colSums(counts>0)
    ncounts  <- colSums(counts)

    if(0) {
        a=2.5
        nfeature.th <- mean(nfeature) + a * sd(nfeature)
        ncounts.th <- mean(ncounts) + a * sd(ncounts)    
        mito.th <- mean(percent.mito) + a * sd(percent.mito)    
        ribo.th <- mean(percent.ribo) + a * sd(percent.ribo)
        mito.th
        ribo.th
        
        par(mfrow=c(2,2))
        hist(nfeature, breaks=100)
        abline(v=nfeature.th, col="red")
        hist(ncounts, breaks=100)
        abline(v=ncounts.th, col="red")
        
        hist(percent.mito, breaks=100)
        abline(v=mito.th, col="red")    
        hist(percent.ribo, breaks=100)
        abline(v=ribo.th, col="red")
    }

    selectInlier <- function(x, a=2.5) {
        xmin <- mean(x) - a * sd(x)
        xmin <- max(xmin, 0.01*mean(x))
        xmax <- mean(x) + a * sd(x)
        x > xmin & x < xmax
    }

    ## sel <- nfeature > 200 & nfeature < 7500 & percent.mito < 10 & percent.ribo < 50
    ## sel <- nfeature < nfeature.th & ncounts < ncounts.th &
    ## percent.mito < mito.th & percent.ribo < ribo.th
    sel <- selectInlier(nfeature,a) & selectInlier(ncounts,a) &
        selectInlier(percent.mito,a) & selectInlier(percent.ribo,a)
    table(sel)

    counts <- counts[,sel]
    counts
}


if(0) {

    ##
    ## From Seurat vignette: Batch correction
    ##
    library(Seurat)

    batch = c("BioReplicate1", "BioReplicate2", "BioReplicate3","BioReplicate4",
              "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6",
              "Tattoo1","Tattoo2","Tattoo3","Tattoo4")
    treatment = c("control","control","control","control",
                  "treated","treated","treated","treated","control","control",
                  "neg_control","neg_control","neg_control","neg_control")
    
    obj.list <- list()
    i=1
    for(i in 1:length(outputs)) {
        data.10x <- Read10X(data.dir=file.path(outputs[i],"/outs/filtered_feature_bc_matrix"))
        dim(data.10x)
        ## data.10x <- data.10x[,1:1000]  ## just subsample FTM...
        celseq <- CreateSeuratObject(data.10x, min.cells=5)
        ##celseq <- FilterCells(celseq, subset.names = "nGene", low.thresholds = 800)
        celseq <- subset(celseq, subset = nFeature_RNA > 500)
        ##-----
        celseq[["percent.mt"]] <- PercentageFeatureSet(celseq, pattern = "^mt-")
        celseq <- subset(celseq, subset = percent.mt < 5)
        ##head(celseq@meta.data, 5)
        ##VlnPlot(celseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.3)
        ##-----
        ##celseq <- NormalizeData(celseq)
        celseq <- NormalizeData(object = celseq, normalization.method = "LogNormalize", scale.factor = 10000)
        celseq <- FindVariableFeatures(celseq, selection.method="vst",
                                       do.plot = F, display.progress = F)
        ##celseq$batch <- gsub(".*/|_scRNA","",outputs[i])
        celseq$batch <- batch[i]
        celseq$treatment <- treatment[i]
        obj.list[[i]] <- celseq
    }
    
    NUM.CC = 20
    anchors    <- FindIntegrationAnchors(obj.list, dims = 1:NUM.CC,
                                         ##anchor.features=9999,
                                         verbose=FALSE)
    integrated <- IntegrateData(anchorset = anchors,
                                dims = 1:NUM.CC,
                                verbose=FALSE)
    dim(integrated)
    DefaultAssay(integrated) <- "integrated"
    dim(integrated)
    
    ## obj <- CreateSeuratObject(counts)    
    ## obj <- AddMetaData(obj, sample.id, col.name = "Sample.id")    
    ## obj
    ## slotNames(obj[["RNA"]])
    ## table(Idents(obj))    

}


if(0) {

    ##
    ## From Seurat vignette: standard QC filtering
    ##
    ##
    library(Seurat)
    obj <- CreateSeuratObject(counts)    
    obj <- AddMetaData(obj, sample.id, col.name = "Sample.id")    
    obj
    slotNames(obj[["RNA"]])
    table(Idents(obj))    
    
    head(obj@meta.data)
    mt.genes <- rownames(obj)[grep("^MT-",rownames(obj),ignore.case=TRUE)]
    C <- GetAssayData(object = obj, slot = "counts")    
    percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
    hist(percent.mito, breaks=100)
    obj <- AddMetaData(obj, percent.mito, col.name = "percent.mito")
    
    rb.genes <- rownames(obj)[grep("^RP[SL]",rownames(obj),ignore.case=TRUE)]
    percent.ribo <- colSums(C[rb.genes,])/Matrix::colSums(C)*100
    obj <- AddMetaData(obj, percent.ribo, col.name = "percent.ribo")
    
    VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo"),
            group.by = "Sample.id",pt.size = 0.1, ncol=4) + NoLegend()
    FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mito")
    FeatureScatter(obj, feature1 = "percent.ribo", feature2 = "nFeature_RNA")
    
    ## QC select cells
    obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 &
                           percent.mito < 10 & percent.ribo < 50)
    dim(obj)
    
    ## total count normalization
    hist(colSums(exp(obj[["RNA"]]@data[,1:1000])),breaks=100)
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)    
    hist(colSums(exp(obj[["RNA"]]@data[,1:1000])),breaks=100)
    
    ## Find highly variable top 2000 genes
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
    VariableFeaturePlot(obj)
    
    ## scale (=standardize????) really??? you will loose fold-change!
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features=all.genes)
    
    obj <- RunPCA(obj, features=VariableFeatures(object=obj))
    ElbowPlot(obj)    
}

