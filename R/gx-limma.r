##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


########################################################################
## Compute LIMMA from matrices or from GCT/CLS file
##
########################################################################

require(limma)
REF.CLASS = c("ctrl","ctr","control","dmso","nt","0","non","no","ref",
              "wt","wildtype","untreated","normal","false","healthy")

gx.limma <- function(gep, pheno, fdr=0.05, compute.means=TRUE, lfc=0.20,
                     max.na=0.20, ref=REF.CLASS, trend=FALSE, verbose=1 )
{
    require(limma)
    if(0) {
        fdr=0.05;compute.means=TRUE;lfc=0.20;ref=REF.CLASS
        max.na=0.2;ref=REF.CLASS;trend=FALSE;verbose=1
    }
    if(sum(duplicated(rownames(gep)))>0) {
        cat("WARNING:: matrix has duplicated rownames\n")
    }
    ## detect single sample case
    is.single = (max(table(pheno))==1)
    if(is.single) {
        cat("WARNING:: no replicates, no stats...\n")
        gep <- cbind(gep,gep)
        pheno <- c(pheno,pheno)
    }

    ## filter probes and samples
    ii <- which( rowMeans(is.na(gep)) <= max.na )
    jj <- which(!is.na(pheno) )
    if(verbose>0) cat(sum(is.na(pheno)>0),"with missing phenotype\n")
    gep0 <- gep[ii,jj]
    pheno0 <- as.character(pheno[jj])
    gep0 <- gep0[!(rownames(gep0) %in% c(NA,"","NA")),]
    if(verbose>0) {
        cat("analyzing",ncol(gep0),"samples\n")
        cat("testing",nrow(gep0),"features\n")
    }

    ## auto-detect reference
    pheno.ref <- c()
    ref.detected <- FALSE
    ref <- toupper(ref)
    ## is.ref <- grepl(paste(ref,collapse="|"),pheno0)
    is.ref <- (toupper(pheno0) %in% toupper(ref))
    ref.detected <- (sum(is.ref)>0 && sum(!is.ref)>0)    

    ##if(!is.null(ref) && sum( toupper(pheno0) %in% ref)>0 ) {
    if(ref.detected) {
        pheno.ref <- unique(pheno0[which(toupper(pheno0) %in% toupper(ref))])
        if(verbose>0) cat("setting reference to y=",pheno.ref,"\n")
        bb <- c(pheno.ref, sort(setdiff(unique(pheno0),pheno.ref)) )
    } else {
        if(verbose>0) cat("WARNING: could not auto-detect reference\n")
        bb <- as.character(sort(unique(pheno0)))
        if(verbose>0) cat("setting reference to first class",bb[1],"\n")
    }
    if(length(bb)!=2) {
        stop("gx.limma::fatal error:only two class comparisons")
        return
    }
    design <- cbind(1, pheno0==bb[2])
    colnames(design) <- c( "WT", "2vs1" )
    d1 <- colnames(design)[1]
    d2 <- colnames(design)[2]
    fit <- lmFit( gep0, design)
    fit <- eBayes(fit, trend=trend)
    top <- topTable(fit, coef=d2, number=nrow(gep0))
    if("ID" %in% colnames(top)) {
        rownames(top) <- top$ID
        top$ID <- NULL
    }
    top <- top[rownames(gep0),]
    head(top)
    
    ## only significant
    top <- top[ which(top$adj.P.Val <= fdr & abs(top$logFC)>=lfc ), ]
    if(verbose>0) cat("found",nrow(top),"significant at fdr=",fdr,"and minimal FC=",lfc,"\n")

    if(compute.means && nrow(top)>0 ) {
        avg <- t(apply(gep0[rownames(top),], 1,
                         function(x) tapply(x, pheno0, mean, na.rm=TRUE)))
        avg <- avg[,as.character(bb),drop=FALSE]
        colnames(avg) <- paste0("AveExpr.",colnames(avg))
        top <- cbind(top, avg)
    }
    top$B <- NULL

    if(is.single) {
        top$P.Value <- NA
        top$adj.P.Val <- NA
        top$t <- NA
    }

    ## reorder on fold change
    top <- top[ order(abs(top$logFC),decreasing=TRUE),]
    ##colnames(top) <-   sub("logFC","logR",colnames(top))

    ## unlist???
    ##top = do.call(cbind, top)
    return(top)
}

gx.meanFstats <- function(gep, pheno) {
    getF <- function(x,y) {
        ii <- which(!is.na(y))
        y1 <- y[ii]
        if(class(y1)=="factor") y1 <- factor(as.character(y1))
        design <- model.matrix(~y1)
        fit <- lmFit( x[,ii], design)
        fit <- eBayes(fit, trend=TRUE)
        top <- topTableF(fit, number=nrow(x))
        mean(top$F)
    }
    fstat <- c()
    px <- tidy.dataframe(pheno)  ## get variable types correct
    for(p in c("random",colnames(px))) {
        if(p=="random") {
            y <- sample(c("a","b"), ncol(gep), replace=TRUE)
        } else {
            y <- px[,p]
            p
        }
        fstat[p] <- getF(gep,y)
    }
    fstat
}

gx.limmaF <- function(gep, pheno, fdr=0.05, compute.means=TRUE, lfc=0.20,
                      max.na=0.20, ref=REF.CLASS, trend=FALSE, verbose=1 )
{
    require(limma)
    if(0) {
        fdr=0.05;compute.means=TRUE;lfc=0.20;ref=REF.CLASS;max.na=0.20;
        trend=TRUE;verbose=1
    }
    if(sum(duplicated(rownames(gep)))>0) {
        cat("matrix has duplicated rownames. please remove.\n")
    }
    ## detect single sample case
    is.single = (max(table(pheno))==1)
    if(is.single) {
        cat("warning: no replicates, no stats...\n")
        gep <- cbind(gep,gep)
        pheno <- c(pheno,pheno)
    }

    ## filter probes and samples
    ii <- which( rowMeans(is.na(gep)) <= max.na )
    jj <- which(!is.na(pheno) )
    if(verbose>0) cat(sum(is.na(pheno)>0),"with missing phenotype\n")
    gep0 <- gep[ii,jj]
    pheno0 <- as.character(pheno[jj])
    gep0 <- gep0[!(rownames(gep0) %in% c(NA,"","NA")),]
    if(verbose>0) {
        cat("analyzing",ncol(gep0),"samples\n")
        cat("testing",nrow(gep0),"features\n")
        cat("in",length(unique(pheno0)),"groups\n")
    }

    ## auto-detect reference
    pheno.ref <- c()
    ref.detected <- FALSE
    ref <- toupper(ref)
    ## is.ref <- grepl(paste(ref,collapse="|"),pheno0)
    is.ref <- (toupper(pheno0) %in% toupper(ref))
    ref.detected <- (sum(is.ref)>0 && sum(!is.ref)>0)
    ref.detected

    ##if(!is.null(ref) && sum( toupper(pheno0) %in% ref)>0 ) {
    if(ref.detected) {
        pheno.ref <- unique(pheno0[which(toupper(pheno0) %in% toupper(ref))])
        if(verbose>0) cat("setting reference to y=",pheno.ref,"\n")
        bb <- c(pheno.ref, sort(setdiff(unique(pheno0),pheno.ref)) )
        pheno1 <- relevel(factor(pheno0), ref=bb[1])
    } else {
        if(verbose>0) cat("WARNING: could not auto-detect reference\n")
        bb <- as.character(sort(unique(pheno0)))
        if(verbose>0) cat("setting reference to first class",bb[1],"\n")
        pheno1 <- relevel(factor(pheno0), ref=bb[1])        
    }
    if(0 && length(bb)!=2) {
        stop("gx.limma::fatal error:only two class comparisons")
        return
    }

    ##design <- cbind(1, pheno0==bb[2])
    design <- model.matrix(~ pheno1)
    colnames(design)
    colnames(design)[2:ncol(design)] <- paste0(levels(pheno1)[-1],"_vs_",levels(pheno1)[1])
    colnames(design) <- gsub("\\(|\\)","",colnames(design))
    fit <- lmFit( gep0, design)
    fit <- eBayes(fit, trend=trend)
    ##top <- topTable(fit, coef=NULL, number=nrow(gep0))
    top <- topTableF(fit, number=nrow(gep0))
    head(top)
    top$B <- NULL
    if("ID" %in% colnames(top)) {
        rownames(top) <- top$ID
        top$ID <- NULL
    }
    top <- top[,setdiff(colnames(top),colnames(design)),drop=FALSE]
    top <- top[rownames(gep0),]
    
    ## compute average
    avg <- do.call(cbind,tapply(1:ncol(gep0), pheno1, function(i)
        rowMeans(gep0[,i,drop=FALSE])))
    ##top <- cbind( top, avg[rownames(top),])
    if(!"logFC" %in% colnames(top)) {
        maxFC <- apply(avg,1,max,na.rm=TRUE) - apply(avg,1,min,na.rm=TRUE)
        top$logFC <- NULL
        top <- cbind(logFC=maxFC, top)
        rownames(top) <- rownames(gep0)
    }
    
    ## only significant
    top <- top[ which(top$adj.P.Val <= fdr & abs(top$logFC)>=lfc ), ]
    if(verbose>0) cat("found",nrow(top),"significant at fdr=",fdr,"and minimal FC=",lfc,"\n")

    if(compute.means && nrow(top)>0 ) {
        avg1 <- avg[rownames(top),]
        colnames(avg1) <- paste0("AveExpr.",colnames(avg1))
        top <- cbind(top, avg1)
    }
    top$B <- NULL

    if(is.single) {
        top$P.Value <- NA
        top$adj.P.Val <- NA
        top$t <- NA
        top$F <- NA
    }

    ## reorder on fold change
    top <- top[order(abs(top$logFC),decreasing=TRUE),]
    ##colnames(top) <-  sub("logFC","logR",colnames(top))

    ## unlist
    ##top = do.call(cbind, top)
    return(top)
}




## two-factorial design, no interaction
gx.limma.paired <- function(gep, pheno, pair, fdr=0.05, lfc=0.20, ref=REF.CLASS,
                            compute.means=FALSE, trend=FALSE )
{
    ##fdr=0.05;lfc=0.20;ref=REF.CLASS;compute.means=TRUE
    ## LIMMA
    require(limma)
    cat("Paired LIMMA\n")
    cat("analyzing",ncol(gep),"samples\n")
    gep <- gep[!(rownames(gep) %in% c(NA,"","NA")),]
    pheno <- as.character(pheno)
    pair  <- as.character(pair)

    ## check
    a <- as.character(pair)
    b <- as.character(pheno)
    v1 <- sort(unique(a))
    v2 <- sort(unique(b))
    if( !is.null(ref) && sum(v1 %in% ref)>0) {
        r1 <- sort(intersect(v1, ref))
        cat("setting reference to",r1,"\n")
        v1 <- c( r1, setdiff(v1,r1))
    }
    if( !is.null(ref) && sum(v2 %in% ref)>0) {
        r2 <- sort(intersect(v2, ref))
        cat("setting reference to",r2,"\n")
        v2 <- c( r2, setdiff(v2,r2))
    }
    a <- factor(a, levels=v1)
    b <- factor(b, levels=v2)
    v1 <- levels(a)
    v2 <- levels(b)
    cat("pairs:",paste(v1),"\n")
    cat("phenotypes:",paste(v2),"\n")
    if( length(v2)>2  ) {
        cat("gx.limma.paired:: fatal error: only two-groups implemented\n")
        return(NULL)
    }
    if( length(v1)<2 ) {
        cat("gx.limma.paired:: fatal error: no pairs\n")
        return(NULL)
    }

    ## setup LIMMA (paired t-test)
##  design <- model.matrix( ~ a * b)
    design <- model.matrix( ~ a + b)

    ## perform fitting
    fit0 <- lmFit( gep, design)
    fit2 <- eBayes(fit0, trend=trend)

    ## extract toptable
    bcoef <- grep("^b",colnames(fit2$coefficients))
    top <- topTable(fit2, coef=bcoef, number=nrow(gep))
    if(colnames(top)[1]=="ID") {  ## old style
        rownames(top) <- top[,"ID"]
        top$ID <- NULL
    }
    top <- top[rownames(gep),]

    ## only significant
    if(fdr < 1) {
        kk <- which( top$adj.P.Val <= fdr & abs(top$logFC) >= lfc)
        top <- top[kk,]
    }
    cat("found",nrow(top),"significant at fdr=",fdr,"and minimum logFC=",lfc,"\n")

    ## compute means if requested
    gep.m <- NULL
    if(compute.means && nrow(top)>0 ) {
        ff <- paste( pair, pheno, sep="." )
        gep.m <- t(apply( gep[rownames(top),],1,
                         function(x) tapply(x, ff, mean, na.rm=TRUE)))
        top$AveExpr <- NULL
        top$AveExpr <- gep.m
    } else {
        ff <- pheno
        gep.m <- t(apply( gep[rownames(top),],1,
                         function(x) tapply(x, ff, mean, na.rm=TRUE)))
        top$AveExpr <- NULL
        top$AveExpr <- gep.m
    }

    ## fold-change
    jj <- order( -abs(top$logFC) )
    top <- top[jj,]
    if(!is.null(gep.m)) gep.m <- gep.m[rownames(top),]

    ## results
    top$B <- NULL
    return(top)
}

## two-factorial design
##ref=c("CTRL","DMSO","NT")
gx.limma.two.factorial <- function(gep, factors, fdr=0.05, lfc=0.20, trend=FALSE,
                                   ref=REF.CLASS, compute.means=TRUE)
{

    ## LIMMA
    require(limma)
    cat("Two-factorial LIMMA\n")
    cat("analyzing",ncol(gep),"samples\n")
    gep <- gep[!(rownames(gep) %in% c(NA,"","NA")),]

    ## create factors
    fct <- vector("list",ncol(factors))
    names(fct) <- colnames(factors)
    j=1
    for(j in 1:ncol(factors)) {
        vv <- sort(unique(as.character(factors[,j])))
        vv
        if(sum(vv %in%ref)>0) {
            vv0 <- vv[which(vv %in% ref)]
            vv <- c(sort(vv0),setdiff(vv,vv0))
        }
        fct[[j]] <- factor( factors[,j], levels=vv )
    }

    ## check
    a <- fct[[1]]
    b <- fct[[2]]
    v1 <- levels(a)
    v2 <- levels(b)
    cat("factors 1:",paste(v1),"\n")
    cat("factors 2:",paste(v2),"\n")
    if(! (length(fct)==2 && length(v1)==2 && length(v2)==2) ) {
        cat("gx.limma2:: fatal error: only 2-factorial with 2 levels implemented\n")
        return(NULL)
    }

    ## setup LIMMA
    design <- model.matrix( ~ a * b)
##    design <- model.matrix( ~ a + b)
    colnames(design)[2] <- paste(rev(v1),collapse="vs")
    colnames(design)[3] <- paste(rev(v2),collapse="vs")
    colnames(design)[4] <- paste(colnames(factors)[1],colnames(factors)[2],sep="*")

    ## perform fitting
    fit0 <- lmFit( gep, design)
    cc0 <- paste(v1[1],".",paste(rev(levels(fct[[2]])),collapse="vs"),sep="")
    cc1 <- paste(v1[2],".",paste(rev(levels(fct[[2]])),collapse="vs"),sep="")
    cont.matrix <- cbind( c(0,0,1,0), c(0,0,1,1), diff=c(0,0,0,1) )
    colnames(cont.matrix) <- c(cc0, cc1, "Diff")
    rownames(cont.matrix) <- colnames(fit0$coefficients)
    fit1 <- contrasts.fit( fit0, cont.matrix )
    fit2 <- eBayes(fit1, trend=trend)

    ## extract toptable
    top1 <- topTable(fit2, coef=colnames(cont.matrix)[1], number=nrow(gep))
    top2 <- topTable(fit2, coef=colnames(cont.matrix)[2], number=nrow(gep))
    top3 <- topTable(fit2, coef=colnames(cont.matrix)[3], number=nrow(gep))
    top1 <- top1[rownames(gep),]
    top2 <- top2[rownames(gep),]
    top3 <- top3[rownames(gep),]

    ## only significant
    kk <- rownames(gep)
    sig <- decideTests(fit2, p.value=fdr, lfc=lfc )
    vennDiagram(sig, cex=0.8)
    title(sub=paste("fdr=",fdr,sep=""))
    if(fdr < 1) {
        sig.up   <- which( sig[,1] > 0 & sig[,2] > 0 & sig[,3]==0)
        sig.down <- which( sig[,1] < 0 & sig[,2] < 0 & sig[,3]==0)
        length(sig.up)
        length(sig.down)
        kk <- rownames(sig)[c(sig.up, sig.down)]
        top1 <- top1[kk,]
        top2 <- top2[kk,]
        top3 <- top3[kk,]
        sig <- sig[kk,]
    }
    cat("found",nrow(sig),"significant at fdr=",fdr,"and minimal FC=",lfc,"\n")

    ## compute means if requested
    gep.m <- NULL
    if(compute.means && nrow(sig)>0 ) {
        ff <- paste( factors[,1], factors[,2], sep="." )
        gep.m <- t(apply( gep[rownames(sig),],1,function(x) tapply(x, ff, mean)))
    }

    ## fold-change
    logfc <- cbind( top1$logFC, top2$logFC, top3$logFC )
    colnames(logfc) <- colnames(cont.matrix)
    rownames(logfc) <- rownames(top1)
    jj <- order( -abs(rowMeans(logfc[,1:2])) )
    top1 <- top1[jj,]
    top2 <- top2[jj,]
    top3 <- top3[jj,]
    sig <- sig[jj,]
    logfc <- logfc[jj,]
    if(!is.null(gep.m)) gep.m <- gep.m[jj,]

    ## pq-values
    pv <- cbind( top1$P.Value, top2$P.Value, top3$P.Value)
    qv <- cbind( top1$adj.P.Val, top2$adj.P.Val, top3$adj.P.Val)
    tt <- cbind( top1$t, top2$t, top3$t)
    rownames(pv) <- rownames(qv) <- rownames(top1)
    colnames(pv) <- colnames(qv) <- colnames(sig)

    ## simple summary
    ss0 <- data.frame( logFC = rowMeans(logfc[,1:2]),
                      t = rowMeans(tt[,1:2]),
                      P.Value = rowMeans(pv[,1:2]),
                      adj.P.Val = rowMeans(qv[,1:2]),
                      AveExpr = gep.m
                      )

    ## results
    res <- c()
    res$fdr <- fdr
    res$means <- gep.m
    res$limma <- list( top1, top2, top3)
    names(res$limma) <- colnames(cont.matrix)
    res$sig <- sig
    res$logFC <- logfc
    res$p.value <- pv
    res$q.value <- qv
    res$summary <- ss0

    return(res)
}

gx.test.groups <- function(sig, class.label, fdr=0.20,
                           test.method=c("wilcox","limma","ttest","fisher"),
                           running.name=NULL,
                           output.to.file=TRUE)
{
    bb <- sort(unique(as.character(class.label)))
    if(length(bb)!=2) {
        stop("meti.limma:: currently only 2 classes")
    }
    if(length(test.method)>1) test.method <- test.method[1]

    ## cleanup NA
    kk <- which(!is.na(class.label))
    sig <- sig[,kk]
    class.label <- class.label[kk]

    ## do LIMMA
    pv <- rep(NA, nrow(sig))
    names(pv) <- rownames(sig)
    if(test.method=="limma") {
        cat("performing Limma test\n")
        require(limma)
        design <- cbind( 1, class.label==bb[2])
        if(is.null(running.name)) running.name <-  paste(bb[2],bb[1],sep="vs")
        colnames(design) <- c( "WT", running.name )
        d1 <- colnames(design)[1]
        d2 <- colnames(design)[2]
        fit <- lmFit( sig, design)
        fit <- eBayes(fit)
        tt <- topTable(fit, coef=d2, number=nrow(sig))
        pv <- tt$P.Value[match(names(pv), tt$ID)]
    } else if(test.method=="wilcox") {
        cat("performing Wilcox rank test\n")
        for(i in 1:nrow(sig)) {
            pv[i] <- wilcox.test( sig[i,] ~ class.label )$p.value
        }
    } else if(test.method=="ttest") {
        cat("performing T-test\n")
        for(i in 1:nrow(sig)) {
            pv[i] <- t.test( sig[i,] ~ class.label )$p.value
        }
    } else if(test.method=="fisher") {
        cat("performing Fisher test\n")
        cbeta <- matrix(cut(sig, breaks=c(-1,median(sig),99), label=c(0,1)),
                        nrow=nrow(sig))
        ii <- which( apply(cbeta,1,function(x) length(setdiff(unique(x),c("NA",NA)))>1 ))
        for(i in ii) {
            pv[i] <- fisher.test( cbeta[i,], class.label )$p.value
        }
    } else {
        stop("unknown test method")
    }

    ## qvalue
    library(qvalue)
    qv <- rep(NA, nrow(sig))
    kk <- which(!is.na(pv))
    qv[kk] <- qvalue(pv[kk])$qvalue

    ## return object
    rr <- data.frame( ID=rownames(sig) )
    rr$Sig0 <- rowMeans( sig[,which(class.label==bb[1])],na.rm=TRUE )
    rr$Sig1 <- rowMeans( sig[,which(class.label==bb[2])],na.rm=TRUE )
    rr$DiffSig <- (rr$Sig1 - rr$Sig0)
    rr$P.Value <- pv
    rr$Q.Value <- qv
    colnames(rr) <- sub("Sig0",paste("AveSig.",bb[1],sep=""),colnames(rr))
    colnames(rr) <- sub("Sig1",paste("AveSig.",bb[2],sep=""),colnames(rr))

    ## order on absolute difference
    rr <- rr[ which(rr$Q.Value < fdr), ]
    rr <- rr[ order(-abs(rr$DiffSig)), ]
    ##  rr <- rr[order(rr$pv),]

    return(rr)
}


##ref.class="CTRL"
gx.snrtest <- function(X,y,ref.class,nperm=200) {
    ## http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
    this.X <- X
    this.y <- y
    calc.snr <- function(this.X,this.y) {
        j1 <- which(this.y!=ref.class)
        j0 <- which(this.y==ref.class)
        ma <- rowMeans(this.X[,j1],na.rm=TRUE)
        mb <- rowMeans(this.X[,j0],na.rm=TRUE)
        sa <- rowSums((this.X[,j1] - ma)**2/(length(j1)-1))**0.5
        sb <- rowSums((this.X[,j0] - mb)**2/(length(j0)-1))**0.5
        sa <- max(0.2*abs(ma),sa,0.2)
        sb <- max(0.2*abs(mb),sb,0.2)
        (ma - mb) / (sa + sb)
    }
    x0 <- calc.snr(X,y)
    pv <- x0*0
    i=1
    for(i in 1:nperm) {
        x1 <- calc.snr(X, sample(y))
        pv <- pv + 1/nperm * (x1 > x0*sign(x0))
    }
    pv <- pmax(pv,1.0/nperm)
    cat("\n")
    head(sort(pv))
    pos.class <- setdiff(y,ref.class)[1]
    logFC <- ma - mb
    qv <- p.adjust(pv, method="fdr")
    res <- data.frame(logFC, snr, pv, qv, ma, mb)
    colnames(res) <- c("logFC","SNR","P.Value","adj.P.Val",
                       paste0("AveExpr.",pos.class),
                       paste0("AveExpr.",ref.class))
    res
}

seq.limma <- function(countdata, y, method="edgeR") {
    ## https://bioinformatics-core-shared-training.github.io/RNAseq-R/rna-seq-de.nb.html
    library(edgeR)
    library(limma)
    library(Glimma)
    library(gplots)

    if( min(countdata)<0 || !all(countdata%%1==0) ) {
        cat("WARNING:: input X should be integer counts! Proceed on own risk\n")
    }

    ## Identify genes with at least 0.5 CPM in at least 2 samples
    myCPM <- cpm(countdata)
    thresh <- myCPM > 0.5
    keep <- rowSums(thresh) >= 2
    ## Subset the rows of countdata to keep the more highly expressed genes
    counts.keep <- countdata[keep,]

    ## Convert to an edgeR object
    group = factor(y)
    dgeObj <- DGEList(counts.keep, group=group)

    ## Perform TMM normalisation
    dgeObj <- calcNormFactors(dgeObj, method="TMM")

    ## Define design matrix
    design <- model.matrix( ~ group)

    ## Estimating the dispersion
    ##dgeObj <- estimateCommonDisp(dgeObj)
    dgeObj <- estimateGLMCommonDisp(dgeObj, design)
    dgeObj <- estimateGLMTrendedDisp(dgeObj)
    dgeObj <- estimateTagwiseDisp(dgeObj)

    ## Fit the linear model
    fit <- glmFit(dgeObj, design)
    names(fit)

    ##Conduct likelihood ratio tests
    res.test <- glmLRT(fit, coef=2)
    toptable = topTags(res.test, n=nrow(countdata))@.Data[[1]]
    head(toptable)



    ## calculate means
    logcpm <- cpm(dgeObj,log=TRUE)
    xmean = c()
    for(y0 in unique(y)) {
        m1 = rowMeans( logcpm[,which(y==y0)], na.rm=TRUE)
        xmean = cbind(xmean, m1)
    }
    colnames(xmean) = paste0("mean.",unique(y))
    head(xmean)
    xmean = cbind( mean=rowMeans(xmean), xmean )

    if(0) {
        PvsV <- makeContrasts(statuspregnant - statusvirgin, levels=design)
        lrt.pVsV <- glmLRT(fit, contrast=PvsV)
        topTags(lrt.pVsV)
    }



}
