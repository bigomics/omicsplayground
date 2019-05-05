if(0) {
    library("devtools")
    install_github("broadinstitute/inferCNV")
    library(infercnv)

    CreateInfercnvObject(raw_counts_matrix, gene_order_file, annotations_file,
                         ref_group_names, delim = "\t")


}


##nsmooth=80;downsample=10
pgx.CNAfromExpression <- function(ngs, nsmooth=40, downsample=10) {

    require(org.Hs.eg.db)
    head(ngs$genes)

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
    head(genes)

    if(!is.null(ngs$counts)) {
        cna <- log2(100 + ngs$counts)
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
    cna <- t(scale(t(cna),center=FALSE))
    cna <- apply(cna,2,mavg)
    cna <- cna - apply(cna,1,median,na.rm=TRUE)
    rownames(cna) <- rownames(cna0)
    dim(cna)

    ##---------------------------------------------------------------------
    ## Statistics
    ##---------------------------------------------------------------------
    cna.var <- apply(cna,1,var)
    names(cna.var) <- rownames(cna)
    head(sort(cna.var,decreasing=TRUE),40)
    
    ##---------------------------------------------------------------------
    ## Downsample if needed
    ##---------------------------------------------------------------------
    if(downsample > 1) {
        ## Downsample
        cat("downsampling CNA matrix...\n")
        n <- downsample
        jj <- as.vector(sapply(1:((nrow(cna)+n)/n),rep,n))[1:nrow(cna)]
        cna <- apply(cna, 2, function(x) tapply(x,jj,mean))
        g1 <- tapply(rownames(cna0),jj,function(x) x[1])
        gg <- tapply(rownames(cna0),jj,paste,collapse=",")
        rownames(cna) <- gg
        genes <- genes[g1,]
        rownames(genes) <- gg
        dim(cna)
    }

    ##---------------------------------------------------------------------
    ## take out small groups/chromsomes
    ##---------------------------------------------------------------------
    ii <- which(genes$chr %in% names(which(table(genes$chr)>3)) )
    cna <- cna[ii,]
    genes <- genes[ii,]
    genes$chr <- as.character(genes$chr)
    table(genes$chr)
    
    res <- list(cna=cna, chr=genes$chr, pos=genes$pos, cna.var=cna.var)
    return(res)
}


pgx.plotCNAHeatmap <- function(ngs, res, annot=NA, pca.filter=20,
                               order.by="clust" )
{
    require(irlba)
    ##source("../R/gx-heatmap.r")
    cna <- res$cna
    chr <- res$chr
    table(chr)

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
        sv <- irlba(cna, nv=k)
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
        sv1 <- irlba(cna,nv=1)$v[,1]
        jj <- order(sv1)
        cna <- cna[,jj]
        sv1 <- sv1[jj]
    } else {
        ## order by hierarchical clustering
        jj <- head(order(-apply(cna,1,sd)),1000)
        hc <- hclust(dist(t(cna[jj,])), method="ward.D2")
        cna <- cna[,hc$order]
    }

    ## create annation matrix
    ann.mat <- NULL
    if(!is.null(annot)) {
        if(is.na(annot)) {
            k <- c(grep("cell.type|tissue|cluster",colnames(ngs$Y)),1)[1]
        } else {
            k <- match(annot, colnames(ngs$Y))
        }
        k
        y <- as.character(ngs$Y[colnames(cna),k])
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

    ##source("../R/gx-heatmap.r")
    par(mgp=c(0.8,0.4,0))
    wa <- 0.1
    if(!is.null(ann.mat)) wa <- 0.05 + 0.016*ncol(ann.mat)
    layout( matrix(1:3,1,3), widths=c(0.2,0.7,wa))

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
    cna0[which(abs(cna0)<0.15)] <- NA
    cna0 <- tanh( 3*cna0 )
    image( 1:nrow(cna), 1:ncol(cna), cna0[,], col=BLUERED2(16),
          ylab="samples", xlab="DNA copy number  (log2R)",
          yaxt="n", yaxs="i", xaxt="n", xaxs="i",
          zlim=c(-1,1)*1.0 )
    nchr <- as.integer(sub("X",23,sub("Y",24,chr)))
    chrbrk <- which(diff(nchr)!=0)
    chrmid <- c(0,chrbrk) + diff(c(0,chrbrk,nrow(cna)))/2
    abline(v=chrbrk, col="grey50", lty=1, lwd=1)
    chrlen <- length(unique(chr))
    j0 <- seq(1,chrlen,2)
    j1 <- seq(2,chrlen,2)
    mtext(unique(chr)[j0], side=3, at=chrmid[j0], cex=0.55, line=0.25 )
    mtext(unique(chr)[j1], side=3, at=chrmid[j1], cex=0.55, line=0.9 )

    if(!is.null(ann.mat)) {
        dim(ann.mat)
        par(mar=c(8,0.5,12,2))
        image( 1:ncol(ann.mat), 1:nrow(ann.mat), t(ann.mat),
              col = rev(grey.colors(2)), xlab="", ylab="",
              yaxt="n", yaxs="i", xaxt="n", xaxs="i")
        mtext( colnames(ann.mat), side=3, at=1:ncol(ann.mat),
              las=3, cex=0.5, line=0.25 )
    } else {
        frame()
    }

    ## done plotting

}
