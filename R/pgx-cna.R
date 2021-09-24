##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

if(0) {
    
    devtools::install_github("broadinstitute/inferCNV")
    

    CreateInfercnvObject(raw_counts_matrix, gene_order_file, annotations_file,
                         ref_group_names, delim = "\t")


}

##nsmooth=80;downsample=10
##load("../pgx/tcga-prad-gx.pgx")
refgroup=NULL
pgx.inferCNV <- function(ngs, refgroup=NULL, progress=NULL ) {

    ## InferCNV: Inferring copy number alterations from tumor single
    ## cell RNA-Seq data
    ##
    ## https://github.com/broadinstitute/inferCNV/wiki
    ##
    ## BiocManager::install("infercnv")
    ##devtools::install_github("broadinstitute/infercnv", ref="RELEASE_3_9")
    require(org.Hs.eg.db)
    symbol <- as.vector(as.list(org.Hs.egSYMBOL))
    chrloc <- as.list(org.Hs.egCHRLOC)
    chr <- as.vector(sapply(chrloc,function(x) names(x)[1]))
    pos <- abs(as.integer(as.vector(sapply(chrloc,function(x) x[1]))))
    chr[sapply(chr,is.null)] <- NA
    chr <- as.character(unlist(chr))
    chr <- chr[match(ngs$genes$gene_name,symbol)]
    pos <- pos[match(ngs$genes$gene_name,symbol)]
    genes <- data.frame(chr = paste0("chr",chr),
                        start=pos-1000,
                        stop=pos-1000) ## fake start/stop
    rownames(genes) <- ngs$genes$gene_name
    Matrix::head(genes)

    ## filter known genes
    jj <- which( genes$chr %in% paste0("chr",c(1:22,"X","Y")) &
                 !is.na(genes$start) & !is.na(genes$stop) )
    length(jj)
    genes <- genes[jj,]
    
    ## prepare data objects
    gg <- intersect(rownames(genes),rownames(ngs$counts))
    length(gg)
    data = ngs$counts[gg,]
    genes = genes[gg,]
    annots = ngs$samples[,"group",drop=FALSE]

    if(FALSE && is.null(refgroup)) {
        ## if no reference group is given, we create a reference by
        ## random sampling of genes.
        ##
        ##        
        ref = t(apply(data, 1, function(x) sample(x,50,replace=TRUE)))
        dim(ref)
        colnames(ref) <- paste0("random.",1:ncol(ref))
        data = cbind(data, ref)
        annots = matrix(c(annots[,1], rep("random",ncol(ref))),ncol=1)
        rownames(annots) = colnames(data)
        colnames(annots) = "group"
        refgroup = c("random")
    }

    ## take out tiny groups
    selgrp <- names(which(table(annots[,1]) >= 2))
    kk <- which(annots[,1] %in% selgrp)
    data <- data[,kk]
    annots <- annots[colnames(data),,drop=FALSE]
    
    ## From inferCNV vignette
    infercnv_obj <- infercnv::CreateInfercnvObject(
                                  raw_counts_matrix = data, 
                                  gene_order_file = genes,
                                  annotations_file = annots,
                                  ref_group_names = refgroup)
    
    out_dir = "/tmp/Rtmpn8rPtL/file19b68b27f09/"
    out_dir = tempfile()
    ##system(paste("mkdir -p",outdir))
    ##unlink(out_dir,recursive=TRUE)   
    cat("DBG pgx.inferCNV:: setting out_dir=",out_dir,"\n")
    
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff = 1, ## 
                                  out_dir = out_dir, 
                                  cluster_by_groups = TRUE, 
                                  ##denoise = TRUE,
                                  ##HMM = TRUE,  ## slow...
                                  num_threads = 4,
                                  no_plot = FALSE)

    
    ##dir(out_dir)
    img.file <- paste0(out_dir,"/infercnv.png")
    ##cnv <- read.table(file.path(out_dir,"expr.infercnv.dat"),check.names=FALSE)
    suppressWarnings(cnv <- data.table::fread(file.path(out_dir,"expr.infercnv.dat"),check.names=FALSE))
    symbol = cnv[[1]]
    cnv <- as.data.frame(cnv, check.names=FALSE)[2:ncol(cnv)]
    cnv <- as.matrix(cnv)
    rownames(cnv) <- symbol

    dim(cnv)
    dim(genes)
    
    genes = genes[rownames(cnv),]
    pos = (genes$start + genes$stop)/2
    ichr <- as.integer(sub("X",23,sub("Y",24,sub("chr","",genes$chr))))
    jj = order(ichr, pos)
    pos = pos[jj]
    chr = as.character(genes$chr)[jj]
    logcnv = log2(cnv[jj,]/mean(cnv,na.rm=TRUE))  ## logarithmic


    
    img <- png::readPNG(img.file)

    res <- list(cna=logcnv, chr=chr, pos=pos, png=img)

    if(0) {
        
        par(mfrow=c(1,1))
        grid::grid.raster(img)

        x11()
        pgx.plotCNAHeatmap(ngs, res, pca.filter=-1, clip=0, annot="group",
                           lwd=2, lab.cex=1)
        ##pgx.plotCNAHeatmap(ngs, res, pca.filter=10, clip=0)
        pgx.plotCNAHeatmap(ngs, res, pca.filter=100, clip=0.15)
        pgx.plotCNAHeatmap(ngs, res, pca.filter=40, clip=0.15)

        grp = (ngs$samples$group)
        annot = ngs$samples[,1:2]
        jj = seq(1,nrow(res$cna),50)
        
        gx.splitmap(
            res$cna[jj,],
            split=res$chr[jj], splitx=grp,
            scale="row.center", col.annot=annot,
            cluster.rows=FALSE)
        
        X = res$cna[jj,]
        annot = annot
        idx = as.character(res$chr)[jj]
        splitx=grp
        xtips=ytips=NULL
        scale="row.center"

        source("../R/pgx-plotting.R")
        pgx.splitHeatmapX(
            X = res$cna[jj,], lmar=200, 
            annot = annot, row_clust=FALSE,
            idx=res$chr[jj], splitx=grp,
            xtips=NULL, ytips=NULL,
            row_annot_width=0.03, scale="row.center",
            colors=NULL, label_size=11 )
    }

    ## clean up folder??
    unlink(out_dir,recursive=TRUE)

    return(res)    
}

##nsmooth=80;downsample=10
pgx.CNAfromExpression <- function(ngs, nsmooth=40)
{
    ## This estimates CNV by local smoothing of relative expression
    ## values.
    ##
    ##
    require(org.Hs.eg.db)    
    symbol <- as.vector(as.list(org.Hs.egSYMBOL))
    chrloc <- as.list(org.Hs.egCHRLOC)
    chr <- as.vector(sapply(chrloc,function(x) names(x)[1]))
    pos <- abs(as.integer(as.vector(sapply(chrloc,function(x) x[1]))))
    chr[sapply(chr,is.null)] <- NA
    chr <- as.character(unlist(chr))
    chr <- chr[match(ngs$genes$gene_name,symbol)]
    pos <- pos[match(ngs$genes$gene_name,symbol)]
    genes <- data.frame(chr=chr, pos=pos)
    rownames(genes) <- ngs$genes$gene_name

    sel <- which(!is.na(genes$chr) & !is.na(genes$pos))
    genes <- genes[sel,]
    Matrix::head(genes)

    if(!is.null(ngs$counts)) {
        cna <- log2(100 + edgeR::cpm(ngs$counts))  ## moderated log2
    } else {
        cna <- ngs$X
    }
    gg <- intersect(rownames(genes),rownames(cna))
    cna <- cna[gg,]
    genes <- genes[gg,]
    dim(cna)

    ##---------------------------------------------------------------------
    ## order genes and matrix according genomic position
    ##---------------------------------------------------------------------
    jj <- which(genes$chr %in% c(1:22,"X","Y"))
    genes <- genes[jj,]
    genes$chr <- factor(genes$chr,levels=c(1:22,"X","Y"))
    jj <- order(genes$chr, genes$pos)
    genes <- genes[jj,]
    cna <- cna[rownames(genes),]
    cna0 <- cna
    dim(cna0)
    
    ##---------------------------------------------------------------------
    ## apply 'crude' moving average filter (THIS SHOULD BE IMPROVED!)
    ##---------------------------------------------------------------------
    ##cna <- cna - rowMeans(cna,na.rm=TRUE)
    mavg <- function(x,n=nsmooth){ stats::filter(x,rep(1/n,n), sides=2, circular=TRUE)}
    cna <- t(scale(t(cna),center=FALSE))  ## z-score
    cna <- apply(cna,2,mavg)
    cna <- cna - apply(cna,1,median,na.rm=TRUE)
    rownames(cna) <- rownames(cna0)
    dim(cna)
        
    res <- list(cna=cna, chr=genes$chr, pos=genes$pos)
    return(res)
}


pgx.plotCNAHeatmap <- function(ngs, res, annot=NA, pca.filter=-1, lwd=1,
                               downsample=10,
                               order.by="clust", clip=0, lab.cex=0.6 )
{
    
    cna <- res$cna
    chr <- res$chr
    chr <- as.character(chr)
    pos <- res$pos
    table(chr)

    ##---------------------------------------------------------------------
    ## Downsample if needed
    ##---------------------------------------------------------------------
    if(downsample > 1) {
        ## Downsample
        cat("downsampling CNA matrix...\n")
        n <- downsample
        jj <- as.vector(sapply(1:((nrow(cna)+n)/n),rep,n))[1:nrow(cna)]
        cna <- apply(cna, 2, function(x) tapply(x,jj,mean))
        ##g1 <- tapply(rownames(cna0),jj,function(x) x[1])
        gg <- tapply(rownames(res$cna),jj,paste,collapse=",")
        rownames(cna) <- gg
        j1 <- which(!duplicated(jj))
        chr <- chr[j1]
        pos <- tapply(pos,jj,mean)
    }

    ##---------------------------------------------------------------------
    ## take out small groups/chromsomes
    ##---------------------------------------------------------------------
    ii <- which(chr %in% names(which(table(chr)>3)) )
    cna <- cna[ii,]
    pos <- pos[ii]
    chr <- chr[ii]
    table(chr)
    
    ## ensure order on chrpos
    ichr <- as.integer(sub("X",23,sub("Y",24,sub("chr","",chr))))
    jj <- order(ichr, pos)
    cna <- cna[jj,]
    chr <- chr[jj]
    pos <- pos[jj]
    
    ## center/scale
    cna <- cna - rowMeans(cna,na.rm=TRUE)
    cna <- cna / max(abs(cna),na.rm=TRUE)
    cna <- tanh(1.3*cna)
    cna <- t(t(cna) - apply(cna,2,median))
    dim(cna)

    if(pca.filter>0) {
        k=20
        k=pca.filter
        ##plot(sv$d)
        k <- ceiling(min(0.33*ncol(cna),k))
        sv <- irlba::irlba(cna, nv=k)
        cna2 <- sv$u[,1:k] %*% diag(sv$d[1:k]) %*% t(sv$v[,1:k])
        colnames(cna2) <- colnames(cna)
        rownames(cna2) <- rownames(cna)
        cna <- cna2
    }
    
    ## sort/order
    hc <- NULL
    sv1 <- NULL
    if(order.by=="pc1") {
        ## by default order on SV1
        sv1 <- irlba::irlba(cna,nv=1)$v[,1]
        jj <- order(sv1)
        cna <- cna[,jj]
        sv1 <- sv1[jj]
    } else {
        ## order by hierarchical clustering
        jj <- Matrix::head(order(-apply(cna,1,sd)),1000)
        hc <- fastcluster::hclust(dist(t(cna[jj,])), method="ward.D2")
        cna <- cna[,hc$order]
    }

    ## create annation matrix
    ann.mat <- NULL
    if(!is.null(annot)) {
        if(is.na(annot)) {
            k <- c(grep("cell.type|tissue|cluster|group",
                        colnames(ngs$samples),ignore.case=TRUE),1)[1]
        } else {
            k <- match(annot, colnames(ngs$samples))
        }
        k
        y <- as.character(ngs$samples[colnames(cna),k])
        table(y)
        ny <- length(setdiff(unique(y),NA))
        if(ny>=2) {
            y[is.na(y)] <- "_"
            ann.mat <- model.matrix( ~ 0 + y)
            colnames(ann.mat) <- sub("^y","",colnames(ann.mat))
            rownames(ann.mat) <- colnames(cna)
        }
    }
    dim(ann.mat)

    BLUERED2 <- colorRampPalette(c("blue3","white","red3"))

    ## ---------- do plotting ------------

    par(mgp=c(0.8,0.4,0))
    wa <- 0.1
    if(!is.null(ann.mat)) wa <- 0.05 + 0.016*ncol(ann.mat)
    plotly::layout( matrix(1:3,1,3), widths=c(0.2,0.7,wa))

    if(!is.null(hc)) {
        par(mar=c(8,2,12,0))
        plot(as.dendrogram(hc),horiz=TRUE,leaflab="none",
             yaxs="i", xaxt="n", yaxt="n" )
    } else if(!is.null(sv1)) {
        par(mar=c(8,3,12,0.3))
        barplot(sv1, horiz=TRUE, border=NA, col="grey50", width=0.1,
                space=0, yaxs="i", xaxt="n")
        mtext("PC1",side=2, cex=0.8)
    } else {
        frame()
    }

    ## main heatmap
    par(mar=c(8,0.2,12,0))
    cna0 <- cna
    cna0 <- tanh( 3*cna0 )
    cna0[which(abs(cna0) < clip)] <- NA
    Matrix::image( 1:nrow(cna), 1:ncol(cna), cna0[,], col=BLUERED2(16),
          ylab="samples", xlab="DNA copy number  (log2R)",
          yaxt="n", yaxs="i", xaxt="n", xaxs="i",
          zlim=c(-1,1)*1.0 )

    ichr <- as.integer(sub("X",23,sub("Y",24,sub("chr","",chr))))
    chrbrk <- which(diff(ichr)!=0)
    chrmid <- c(0,chrbrk) + diff(c(0,chrbrk,nrow(cna)))/2
    abline(v=chrbrk, col="grey50", lty=1, lwd=lwd)
    chrlen <- length(unique(chr))
    j0 <- seq(1,chrlen,2)
    j1 <- seq(2,chrlen,2)
    mtext(unique(chr)[j0], side=3, at=chrmid[j0], cex=lab.cex, line=0.25 )
    mtext(unique(chr)[j1], side=3, at=chrmid[j1], cex=lab.cex, line=0.9 )

    if(!is.null(ann.mat)) {
        dim(ann.mat)
        par(mar=c(8,0.5,12,2))
        Matrix::image( 1:ncol(ann.mat), 1:nrow(ann.mat), t(ann.mat),
              col = rev(grey.colors(2)), xlab="", ylab="",
              yaxt="n", yaxs="i", xaxt="n", xaxs="i")
        mtext( colnames(ann.mat), side=3, at=1:ncol(ann.mat),
              las=3, cex=lab.cex, line=0.25 )
    } else {
        frame()
    }

    ## done plotting

}
