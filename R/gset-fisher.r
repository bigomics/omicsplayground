##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##========================================================================
##======================== Fisher test based =============================
##========================================================================

##fdr=1.0;mc=TRUE;sort.by="zratio";nmin=3;min.genes=10;max.genes=500;method="fisher";common.genes=TRUE;check.background=TRUE;verbose=TRUE

gset.fisher2 <- function(genes.up, genes.dn, genesets, background=NULL,
                         fdr=0.05, mc=TRUE, sort.by="zratio", nmin=3, verbose=1,
                         min.genes=15, max.genes=500, method="fast.fisher",
                         check.background=TRUE, common.genes=TRUE )
{    
    ft.up = gset.fisher(genes=genes.up, genesets=genesets, background=background,
                        fdr=1, mc=mc, sort.by=sort.by, nmin=nmin, verbose=verbose,
                        min.genes=min.genes, max.genes=max.genes, method=method,
                        check.background=check.background, common.genes=common.genes)
    ft.dn = gset.fisher(genes=genes.dn, genesets=genesets, background=background,
                        fdr=1, mc=mc, sort.by=sort.by, nmin=nmin, verbose=verbose,
                        min.genes=min.genes, max.genes=max.genes, method=method,
                        check.background=check.background, common.genes=common.genes)
    ft.up = ft.up[rownames(ft.dn),]
    ft.sign = c(-1,1)[1+1*(ft.up$p.value < ft.dn$p.value)]
    ft1 = cbind(sign=ft.sign, ft.up)[which(ft.sign>0),,drop=FALSE]
    ft2 = cbind(sign=ft.sign, ft.dn)[which(ft.sign<0),,drop=FALSE]
    ft.res = rbind(ft1, ft2)
    ft.res$sign = ft.res$sign * ft.res$odd.ratio
    ft.res = ft.res[ which(ft.res$q.value <= fdr),,drop=FALSE]
    return(ft.res)
}

gset.fisher <- function(genes, genesets, background=NULL, 
                        fdr=0.05, mc=TRUE, sort.by="zratio", nmin=3,
                        min.genes=15, max.genes=500, method="fast.fisher",
                        check.background=TRUE,  common.genes=TRUE, verbose=1)
{
    ## NEED RETHINK
    
    if(0) {
        min.genes=-10;max.genes=500;background=NULL;fdr=1;mc=TRUE;sort.by="zratio";nmin=3;verbose=1
        genesets=readRDS("../files/gmt-all.rds")
        genes=head(rownames(ngs$X),300);genesets=ngs$gmt.all;background=rownames(ngs$X)
    }
    if( is.null(background) ) {
        ##background <- unique(c(unlist(genesets),genes))
        background <- unique(unlist(genesets))
        if(verbose>0)
            cat("setting background to ",length(background),"genes covered\n")
    }

    if(check.background) {
        ## restrict on background
        genes <- intersect(genes, background)
        genesets <- parallel::mclapply( genesets, function(s) intersect(s, background))
    }
    
    ## select
    if(!is.null(min.genes) && min.genes>0) {
        genesets.len <- sapply(genesets,length)
        genesets <- genesets[order(-genesets.len)]
        if(sum(duplicated(names(genesets)))>0) {
            if(verbose>0) cat("warning: duplicated gene set names. taking largest.\n")
            genesets <- genesets[which(!duplicated(names(genesets)))]
        }
        genesets.len <- sapply(genesets,length)
        genesets <- genesets[ genesets.len >= min.genes & genesets.len <= max.genes]
        if(verbose>0) cat("testing",length(genesets),"genesets with",length(genes),
                          "genes (background",length(background),"genes)\n")
        length(genesets)
        if(length(genesets)==0) {
            cat("warning: no gene sets passed size filter\n")
            rr <- data.frame(p.value=NA, q.value=NA, ratio0=NA, ratio1=NA,
                             zratio=NA, n.size=NA, n.overlap=NA, genes=NA)
            rownames(rr) <- NULL
            return( rr[0,] )
        }
    }

    ## odd ratio
    ## see http://jura.wi.mit.edu/bio/education/hot_topics/enrichment/Gene_list_enrichment_Mar10.pdf
    n.size = sapply(genesets,length)
    bg0 <- setdiff(background,genes)
    nbackground0 = length(background)
    nbackground1 = length(bg0)

    ## this can become slow... PLEASE OPTIMIZE!!!!
    a <- unlist(parallel::mclapply(genesets, function(x) sum(x %in% genes)))
    ##b <- unlist(parallel::mclapply(genesets, function(x) sum(!(x %in% genes))))
    b <- (n.size - a)
    c <- unlist(parallel::mclapply(genesets, function(x) sum(!(genes %in% x))))
    ##d <- unlist(parallel::mclapply(genesets, function(x) sum(!(bg0 %in% x))))
    d <- (nbackground1 - b)
    odd.ratio <- (a/c) / (b/d)    ## note: not exactly same as from fishertest

    ## intersection genes (slow..)
    commongenes <- NULL
    if(common.genes) {
        commongenes <- unlist(lapply(genesets, function(x) paste(sort(intersect(genes,x)),collapse="|")))
    }

    ## compute fisher-test (should be one-sided?)
    ##gg <- unique(c(unlist(genesets),genes))
    test.fisher <- function(gs) {
        a0 <- table(background %in% gs, background %in% genes)
        if(NCOL(a0)==1 || colSums(a0)[2]==0) return(NA)
        fisher.test(a0, alternative="greater")$p.value
    }
    test.chisq <- function(gs) {
        a0 <- table(background %in% gs, background %in% genes)
        if(NCOL(a0)==1 || colSums(a0)[2]==0) return(NA)        
        chisq.test(a0)$p.value
    }
    pv <- rep(NA, length(genesets))
    names(pv) <- names(genesets)
    if(method=="fast.fisher") {
        ## this is really fast...
        pv <- rep(NA,length(a))
        ii <- 1:length(a)
        ii <- which((a+c)>0)
        pv1 <- corpora::fisher.pval( a[ii], (a+b)[ii], c[ii], (c+d)[ii], alternative="greater")
        pv[ii] <- pv1
    } else if(method=="fisher") {
        if(mc) {
            ##cat("multicore testing\n")
            pv <- unlist(parallel::mclapply( genesets, test.fisher ))
        } else {
            i=1
            for(i in 1:length(genesets)) {
                pv[i] <- test.fisher( genesets[[i]] )
            }
        }
    } else if(method=="chisq") {
        if(mc) {
            ##cat("multicore testing\n")
            pv <- unlist(parallel::mclapply( genesets, test.chisq ))
        } else {
            for(i in 1:length(genesets)) {
                pv[i] <- test.chisq( genesets[[i]] )
            }
        }
    } else {
        stop("unknown method")
    }

    ## compute q-value
    qv <- rep(NA,length(pv))
    qv <- p.adjust(pv, method="fdr")

    ## results
    v1 = as.character(paste0(a,"/",n.size))
    rr <- data.frame(p.value=pv, q.value=qv, odd.ratio=odd.ratio, overlap=v1)
    if(!is.null(commongenes)) {
        rr <- cbind(rr, genes=commongenes)
    }
    rownames(rr) <- names(genesets)

    ## sort
    if(nrow(rr)>0) {
        ## filter
        jj <- which(rr$q.value <= fdr & n.size >= nmin )
        rr <- rr[jj,]
        ## sort
        if(sort.by %in% c("pvalue","p.value","p")) {
            rr <- rr[ order(rr$p.value),]
        } else {
            rr <- rr[ order(rr$odd.ratio,decreasing=TRUE),]
        }
    }
    dim(rr)
    rr
}

##========================================================================
##======================= end of file ====================================
##========================================================================

