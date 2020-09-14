##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

pgx.seuratSingleCellQCFilter <- function(counts)
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

    nfeature.th <- mean(nfeature) + 3 * sd(nfeature)
    ncounts.th <- mean(ncounts) + 3 * sd(ncounts)    
    mito.th <- mean(percent.mito) + 3 * sd(percent.mito)    
    ribo.th <- mean(percent.ribo) + 3 * sd(percent.ribo)
    mito.th
    ribo.th

    if(0) {
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

    selectInlier <- function(x, a=3) {
        xmin <- mean(x) - a * sd(x)
        xmin <- max(xmin, 0.01*mean(x))
        xmax <- mean(x) + a * sd(x)
        x > xmin & x < xmax
    }

    ## sel <- nfeature > 200 & nfeature < 7500 & percent.mito < 10 & percent.ribo < 50
    ## sel <- nfeature < nfeature.th & ncounts < ncounts.th &
    ## percent.mito < mito.th & percent.ribo < ribo.th
    sel <- selectInlier(nfeature,3) & selectInlier(ncounts,3) &
        selectInlier(percent.mito,3) & selectInlier(percent.ribo,3)
    table(sel)

    counts <- counts[,sel]
    counts
}
