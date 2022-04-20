##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


######################################################################## 
##
## PROTEOMICS AND SILAC FUNCTIONS
##
######################################################################## 

## sep="\t";collapse.gene=TRUE;use.LFQ=FALSE;filter.contaminants=TRUE
prot.readProteinGroups <- function(file, meta=NULL, sep="\t", collapse.gene=TRUE,
                                   is.log=FALSE, use.LFQ=FALSE,
                                   filter.contaminants=TRUE)
{
    ## Split data file
    message("reading proteinGroups file ",file)
    ## D = read.csv(file, sep=sep, check.names=FALSE)
    ## D = read.csv(file, sep="\t", check.names=FALSE)
    D = data.table::fread(file, check.names=FALSE)
    D = data.frame(D, check.names=FALSE)
    ## D = data.frame(D, check.names=TRUE) ## need dots????
    dim(D)
    colnames(D)
    Matrix::head(D)[,1:10]

    ##col.required <- c("Gene.names","Protein.names")
    ## Filter contaminants
    contaminant.cols <- c("Reverse","Only identified by site","Potential contaminant")
    contaminant.cols <- intersect(contaminant.cols,colnames(D))
    contaminant.cols
    D$is.contaminant <- (rowSums(D[,contaminant.cols,drop=FALSE]=='+',na.rm=TRUE)>=1)
    table(D$is.contaminant)
    if(filter.contaminants) {
        ##D$Gene.names[which(D$is.contaminant)]
        D <- D[which(!D$is.contaminant),]
    }
    dim(D)
    
    ## parse gene annotation
    genes = D[,c("Majority protein IDs","Gene names","Protein names")]
    Matrix::head(genes)
    colnames(genes) = c("protein_id","gene_name","gene_title")
    gg = as.character(genes$gene_name)
    gg = sapply(gg, function(x) strsplit(x,split=";")[[1]][1]) ## take just FIRST gene
    genes$gene_alias <- genes$gene_name 
    genes$gene_name <- gg
    genes$is_contaminant <- D$is.contaminant
    
    ## give unique rownames
    rownames(D) = paste0("tag",1:nrow(D),":",genes$gene_name)  ## add one gene to tags
    rownames(genes) = rownames(D)
    dim(D)
    dim(genes)
    
    ## extract data blocks (use LFQ as intensity)
    counts = D[,grep("^Intensity ",colnames(D),value=TRUE)]
    dim(counts)
    if(use.LFQ) {
        sel <- grep("^LFQ Intensity",colnames(D))
        if(length(sel)==0) {
            stop("FATAL:: could not find LFQ Intensity data")
        }
        counts = D[,sel,drop=FALSE]
    }
    if(is.log) {
        counts <- 2**as.matrix(counts)
    } else {
        counts <- as.matrix(sapply(counts,as.numeric))  ## from integer64
    }
    colnames(counts) <- sub("Intensity ","",colnames(counts))
    colnames(counts) <- sub("LFQ intensity ","",colnames(counts))
    sum(is.na(counts))
    summary(Matrix::colSums(counts,na.rm=TRUE))

    ## collapse by gene
    if(collapse.gene) {
        cat("collapsing to FIRST gene\n")
        gene <- genes$gene_name
        counts <- apply(counts, 2, function(x) tapply(x, gene, sum, na.rm=TRUE))
        prots <- tapply(as.character(genes$protein_id), gene, function(s)
            paste(unique(unlist(strsplit(s,split=";"))),collapse=";"))
        genes <- genes[match(rownames(counts),genes$gene_name),]
        genes$protein_id <- NULL
        genes$protein_ids <- prots[match(rownames(counts),names(prots))]
        rownames(genes) <- genes$gene_name
    }
    
    if(!is.null(meta) && is.character(meta)) {
        if(grepl("csv$",meta)) {
            meta <- read.csv(meta, header=TRUE)
        } else {
            meta <- read.delim(meta, header=TRUE, sep="\t")
        }
    }
    if(is.null(meta)) {
        meta <- data.frame(sample.name=colnames(counts))
    }

    res <- c()
    res$samples <- meta
    res$genes <- genes
    res$counts <- as.matrix(counts)

    message("sample info template created but please complete ngs$samples")
    
    return(res)
}


##file="proteinGroups.txt";meta="meta.txt";use.LFQ=FALSE;unit="intensity";is.log2=FALSE
##meta="samples.csv";
proteus.readProteinGroups <- function(file="proteinGroups.txt", meta="meta.txt",
                                      unit="intensity", use.LFQ=FALSE, is.log2=FALSE,
                                      collapse.gene=FALSE, na.zero=FALSE) 
{
    
    

    ##------------------------------------------------------------
    ## Read protein data
    ##------------------------------------------------------------

    if(is.character(meta)) {
        if(grepl("csv$",meta)) {
            meta <- read.csv(meta, header=TRUE)
        } else {
            meta <- read.delim(meta, header=TRUE, sep="\t")
        }
    } 

    if(!"sample" %in% colnames(meta)) {
        stop("metadata file must have 'sample' column")
    }
    if(!"condition" %in% colnames(meta)) {
        stop("metadata file must have 'condition' column")
    }

    samples_with_data <- sub("Intensity ","",grep("^Intensity ",colnames(data.table::fread(file,nrow=5)),value=TRUE))
    samples_in_meta <- meta$sample

    if(any(!samples_in_meta %in% samples_with_data)) {
        sel_nodata <- setdiff(samples_in_meta, samples_with_data)
        message("WARNING: no data for samples: ",sel_nodata," (discarding)")
        meta <- meta[which(samples_in_meta %in% samples_with_data),]
    }
    if(any(!samples_with_data %in% samples_in_meta)) {
        sel_nometa <- setdiff(samples_with_data, samples_in_meta)
        message("WARNING: no meta-information for samples: ",sel_nometa," (discarding)")
    }
    
    measure.cols <- NULL
    measure.cols <- paste("Intensity",meta$sample)
    if(use.LFQ) {
        measure.cols <- paste("LFQ intensity",meta$sample)
    }
    names(measure.cols) <- meta$sample
    pdat <- proteus::readProteinGroups(file, meta, measure.cols=measure.cols)

    summary(pdat)
    ##summary(Matrix::colSums(pdat$tab, na.rm=TRUE))

    if(is.log2) {
        pdat$tab <- 2**pdat$tab
    }
    
    ## set NA to zero if needed
    if(na.zero) {
        pdat$tab[is.na(pdat$tab)] <- 0
    }
    
    ##------------------------------------------------------------
    ## create protein2gene translation vector
    ##------------------------------------------------------------
    pfile <- as.data.frame(data.table::fread(file))
    rownames(pfile) <- pfile[,"Majority protein IDs"]
    table(rownames(pdat$tab) %in% rownames(pfile))
    pfile <- pfile[rownames(pdat$tab),]
    dim(pfile)
    
    ##------------------------------------------------------------
    ## convert to copynumber 
    ##------------------------------------------------------------    
    if(unit=="intensity") {
        ## nothing
    } else if(unit=="copynumber") {
        cat("converting to copy number\n")
        ##
        ## actually all the constant mulitiplications do not matter
        ## because the scaling/normalization will cancel it out....
        ##
        ##source("R_FUNCTIONS_ROGER.R")
        if(!("Mol. weight [kDa]" %in% colnames(pfile))) {
            stop("missing column 'Mol. weight [kDa]' in data file")
        }
        mol.kda <- pfile[,"Mol. weight [kDa]"]
        cell.mass <- 1000  ## picogram guestimate
        pdat$tab <- prot.calcCopyNumber(pdat$tab, mol=mol.kda, y=cell.mass)       
    } else {
        stop("invalid unit")
    }
    
    ## get genes
    px <- sapply(rownames(pdat$tab),strsplit,split=";")
    if(!("Gene names") %in% colnames(pfile)) stop("missing 'Gene names' in file")
    pdat$gene.names <- as.character(pfile[,"Gene names"])
    ## take just FIRST gene
    pdat$gene <- sapply(pdat$gene.names, function(x) strsplit(x,split=";")[[1]][1]) 
    
    ## collapse to FIRST gene
    if(collapse.gene) {

        cat("collapsing to FIRST gene\n")
        gene <- as.character(pdat$gene)
        pdat$tab <- apply(pdat$tab, 2, function(x) tapply(x, gene, sum, na.rm=TRUE))
        pdat$detect <- apply(pdat$detect, 2, function(x) tapply(x, gene, sum, na.rm=TRUE))            

        ## collapse stats
        S <- tapply(1:nrow(pdat$stats), pdat$stats$condition, function(i) pdat$stats[i,,drop=FALSE])
        S <- lapply(1:length(S), function(k) {
            s1 <- apply(S[[k]][,3:5], 2, function(x) tapply(x, gene, mean, na.rm=TRUE))
            data.frame(id=rownames(s1), condition=names(S)[k], s1)
        })            
        S <- do.call(rbind, S)
        rownames(S) <- NULL        

        ## merge multiple proteins/genes names
        prots2 <- tapply(as.character(pdat$proteins), gene, function(s)
            paste(unique(unlist(strsplit(s,split=";"))),collapse=";"))
        genes2 <- tapply(as.character(pdat$gene.names), gene, function(s)
            paste(unique(unlist(strsplit(s,split=";"))),collapse=";"))

        jj <- which(!(rownames(pdat$tab) %in% c("",NA," ")))
        pdat$tab <- pdat$tab[jj, ]
        pdat$proteins <- prots2[match(rownames(pdat$tab),names(prots2))]
        pdat$gene.names <- genes2[match(rownames(pdat$tab),names(genes2))]
        pdat$gene <- rownames(pdat$tab)
        pdat$detect <- pdat$detect[rownames(pdat$tab),]>0        
        ##pdat$stats <- S[S$id %in% pdat$proteins,]
        pdat$stats <- S[S$id %in% rownames(pdat$tab),]        
    }
    
    ## return extended proteus object
    remove(pfile)
    return(pdat)
}

##counts=prot$tab;plot=TRUE;qnormalize=TRUE;prior.count=1;zero.thr=0.25;scaling=1e6
prot.normalizeCounts <- function(counts, scale=1e6, scaling="by.column",
                                 qnormalize=TRUE, zero.offset=0.01, 
                                 zero.thr=0.25 )

{    

    ##------------------------------------------------------------
    ## start with original counts from MaxQuant
    ##------------------------------------------------------------
    X <- counts
    sum(is.na(X))
    sum(X==0)
    sum(X==0,na.rm=TRUE)
    min(X,na.rm=TRUE)
    X[ X <= zero.thr ] <- NA  ## treat zero as NA for the moment
    which.zeros <- Matrix::which(X==0 | is.na(X), arr.ind=TRUE)
    length(which.zeros)
    
    ##------------------------------------------------------------
    ## normalize to CPM (or otherwise)
    ##------------------------------------------------------------
    ##scale=1e6
    if(is.null(scale)) {
        scale <- mean(Matrix::colSums(X,na.rm=TRUE))
        message("set scale parameter to = ",scale)
    }

    if(scaling=="by.column") {
        message("scaling by columns to ",scale," total counts")
        X <- t(t(X) / Matrix::colSums(X,na.rm=TRUE)) * scale  ## counts-per-million
    } else if(scaling=="global") {
        message("global average scaling to ",scale," counts")
        X <- X / mean(Matrix::colSums(X,na.rm=TRUE)) * scale  ## counts-per-million
    } else {
        stop("FATAL:: unknown scaling method. scaling = ",scaling)
    }
    Matrix::colSums(X, na.rm=TRUE)
    
    ##------------------------------------------------------------
    ## set zero offset
    ##------------------------------------------------------------
    ##zero.offset=0.01;zero.thr=0.25
    if(zero.offset > 0 && zero.offset < 1) {
        q0 <- quantile(X,probs=zero.offset,na.rm=TRUE)
        log2(q0)
        message("set log2(x)=0 to ",log2(q0))
        X <- X / q0 ## scale to unit
    }
       
    ##------------------------------------------------------------
    ## choose quantile or median normalization (on not-missing values)
    ##------------------------------------------------------------
    if(qnormalize) {
        ##gx.hist(log2(X))
        jj <- sample(ncol(X),10)
        X <- limma::normalizeQuantiles(X)
        ##X <- proteus::normalizeMedian(X)
        X <- pmax(X,0)
    }

    ##------------------------------------------------------------
    ## put back 'missing' values to zero
    ##------------------------------------------------------------
    X[which.zeros] <- 0
    
    return(X)
}

prot.testTwoGroups <- function(X, group1, group2, method="limma",
                               labels=NULL, gene=NULL) {

    ##X <- silac$LFQ.ratio
    ##method = "geiger"
    out1 <- NULL ## important to avoid old results
    p1=p2=p3=NULL  ## just for testing

    if(!is.null(labels) && all(c(group1,group2) %in% labels)) {
        group1 <- colnames(X)[which(labels %in% group1)]
        group2 <- colnames(X)[which(labels %in% group2)]
    }
    
    if(method == "geiger") {
        ## Your code
        out1 <- volcano(X[,group1], X[,group2], rownames(X))
        out1$Gene <- NULL
        colnames(out1) <- c("logFC","P.Value")
        dim(out1)
        Matrix::head(out1)

    } else if(method %in% c("t.welch","t.equalvar")) {
        ## faster t-test
        
        j1 <- which( colnames(X) %in% group1)
        j2 <- which( colnames(X) %in% group2)
        if(method=="t.welch")
            out0 <- matrixTests::row_t_welch( X[,j1], X[,j2] )
        if(method=="t.equalvar")
            out0 <- matrixTests::row_t_equalvar( X[,j1], X[,j2] )
        out0$pvalue[is.na(out0$pvalue)] <- 1
        out0$qvalue <- p.adjust(out0$pvalue, method="fdr")
        out1 <- out0[,c("mean.diff","pvalue","qvalue")]
        colnames(out1) <- c("logFC","P.Value","adj.P.Val")
        Matrix::head(out1)

    } else if(method == "limma" && is.null(labels) ) {
        ## See e.g. https://bioconductor.org/help/course-materials/2010/BioC2010/limma2.pdf
        
        jj <- which( colnames(X) %in% c(group1,group2))
        y <- 1*(colnames(X)[jj] %in% group1)
        design <- model.matrix( ~ y )
        colnames(design) <- c("Intercept","group1_vs_group2")
        fit <- limma::eBayes(limma::lmFit(X[,jj], design), trend=TRUE, robust=TRUE)
        out1 <- limma::topTable(fit, coef=2, sort.by="none", number=Inf)
        out1 <- out1[,c("logFC","P.Value","adj.P.Val")]

    } else if(method == "limma" && !is.null(labels) ) {
        ## Same as above but we retain all samples in the model
        
        design <- model.matrix( ~ 0 + labels)
        colnames(design) <- sub("labels","", colnames(design))
        rownames(design) <- colnames(X)        
        
        level1 <- unique(labels[match(group1,colnames(X))])
        level2 <- unique(labels[match(group2,colnames(X))])
        ##aa0 <- paste0(level1," - ",level2)
        ##contr <- limma::makeContrasts( aa0, levels=design)
        levels <- colnames(design)
        contr <- matrix( 1*(levels %in% level1) - 1*(levels %in% level2), ncol=1)
        rownames(contr) <- levels
        contr        
        contr <- contr[match(colnames(design),rownames(contr)),,drop=FALSE]
        rownames(contr) <- colnames(design)
        colnames(contr)[1] <- "contrast1"
        contr[is.na(contr)] <- 0
        fit1 <- limma::lmFit(X[,], design)
        fit2 <- limma::contrasts.fit(fit1, contr)
        fit2 <- limma::eBayes(fit2, trend=TRUE, robust=TRUE)        
        out1 <- limma::topTable(fit2, coef=1, sort.by="none", number=Inf)
        out1 <- out1[,c("logFC","P.Value","adj.P.Val")]
        
    } else if(method == "msms.edgeR") {
        
        s <- colnames(X)
        grp <- c(NA,"group1","group2")[1 + 1*(s %in% group1) + 2*(s %in% group2) ]
        pd <- data.frame( sample=colnames(X), group=grp)
        rownames(pd) <- colnames(X)
        fd <- data.frame(proteins=rownames(X), gene=rownames(X))
        rownames(fd) <- rownames(X)
        jj <- which(!is.na(grp))
        e <- MSnSet( 2**X[,jj], fd, pd[jj,])
        null.f <- "y~1"
        alt.f <- "y~group"
        div <- apply(exprs(e),2,sum,na.rm=TRUE)
        out1 <- msms.edgeR(e, alt.f, null.f, div=div, fnm="group")
        out1$adj.P.Val <- p.adjust( out1$p.value, method="fdr")
        out1$LogFC <- -out1$LogFC  ## seems switched
        out1 <- out1[,c("LogFC","p.value","adj.P.Val")]
        colnames(out1) <- c("logFC","P.Value","adj.P.Val")        

    } else {
        stop("ERROR:: unknown method")
    }

    ## set NaN p.values to 1??
    out1$P.Value[which(is.nan(out1$P.Value))] <- 1

    ## add FDR q.value
    if(!("adj.P.Val" %in% names(out1))) {
        out1$adj.P.Val <- p.adjust( out1$P.Value, method="fdr")
    }

    ## good practice to add group means
    avg <- cbind( rowMeans(X[,group1],na.rm=TRUE), rowMeans(X[,group2],na.rm=TRUE))
    ##colnames(avg) <- paste0("avg.",ct)
    colnames(avg) <- paste0("avg.",c("group1","group2"))
    out1 <- cbind(out1, avg)
    if(!is.null(gene)) {
        out1 <- cbind(gene, out1)
    }
    
    return(out1)
}

##pcex=0.7;cex=0.7;lfc=1;psig=0.05
prot.plotVolcano <- function(res, use.q=TRUE, psig=0.05, lfc=1, pmin=1e-12,
                             pcex=0.7, cex=0.7, plot.type="default", ...)
{
    
    pdata <- cbind( res[,c("logFC","P.Value")], res$gene)
    colnames(pdata) <- c("EFFECTSIZE", "P", "Gene")
    if(use.q) pdata$P <- res$adj.P.Val
    sig <- (pdata$P < psig & abs(pdata$EFFECTSIZE) > lfc)
    highlight <- pdata$Gene[which(sig)]
    pdata$P <- pmin + pdata$P    

    if(plot.type=="volcanoly") {
        
        volcanoly( pdata[,], snp="Gene", highlight=highlight,
                  effect_size_line = c(-1,1)*lfc,
                  genomewideline = -log10(pig) )
    } else {

        klr <- c("black","red")[1 + 1*sig]
        ylab="significance (-log10p)"
        if(use.q) ylab="significance (-log10q)"
        plot( pdata$EFFECTSIZE, -log10(pdata$P),
             xlab="effect size (log2FC)", ylab=ylab, ...,
             pch=20, cex=pcex, col=klr)
        abline(v=c(-1,1)*lfc, lty=2, col="grey70")
        abline(h=-log10(psig), lty=2, col="grey70")
        jj <- which(sig)
        if(length(jj)) {
            text( pdata$EFFECTSIZE[jj], -log10(pdata$P)[jj],
                 pdata$Gene[jj], cex=cex, adj=0.5,
                 offset=0.4, pos=3)
        }
    }        
}

##======================================================================
##=================== ROGER FUNCTIONS ==================================
##======================================================================
##
##

prot.filterReverseContaminantOnly <- function(x, rev=T, con=T, only=T){
  if(rev==T){x <- x[x$Reverse!="+",]}
  if(con==T & sum(colnames(x)=="Potential.contaminant")==1){x <- x[x$Potential.contaminant!="+",]}
  if(con==T & sum(colnames(x)=="Contaminant")==1){x <- x[x$Contaminant!="+",]}
  if(only==T){x <- x[x$Only.identified.by.site!="+",]}
  return(x)
}

prot.imputeMissing <- function(X, method, downshift=1.8, width=0.3,
                               k=10, q=0.01, groups=NULL, zero.na=TRUE)
{
    ## in count space!!!

    ## treat zeros as missing values
    if(zero.na) X[X==0] <- NA

    if(method=="zero") {
        X[is.na(X)] <- 0
    }
    if(method=="min") {
        X[is.na(X)] <- min(X[!is.na(X) & X>0])
    }
    if(method=="quantile") {
        X[is.na(X)] <- quantile(X[!is.na(X) & X>0], probs=q)[[1]]
    }
    if(method=="gaussian") {
        logX <- log(1+X)
        impX <- .RGimputation(logX, width=width, downshift=downshift, bycolumn=TRUE)
        X <- exp(pmax(impX,0))-1
    }
    if(method=="group.median") {
        if(is.null(groups)) stop("'group.median' method needs groups")
        X <- prot.medianImpute(X, groups) 
    }
    if(method=="nmf") {
        if(is.null(groups)) stop("'nmf' method needs groups")
        X <- prot.nmfImpute(X, groups, k=k) 
    }
    if(method=="none") {
        ## nothing
    }
    return(X)
}

k=10
prot.nmfImpute <- function(X, groups, k=10, r=0.5) {
    setZERO <- function(x,y) {
        b=tapply(x, y, function(a) {
            if(mean(is.na(a))>=r) a[is.na(a)]=0;a}) 
        names(b) <- NULL
        unlist(b)[names(x)]
    }    
    X[X==0] <- NA
    sum(is.na(X))
    impX <- t(apply(X,1,setZERO,y=groups))
    sum(is.na(impX))
    out <- NNLM::nnmf(impX[,],k=k,check.k=0)
    hatX <- (out$W %*% out$H)
    impX[is.na(impX)] <- hatX[is.na(impX)]
    return(impX)    
}

prot.medianImpute <- function(X, groups) {
    logmedian <- function(x) {exp(median(log(x)))}
    medianImputeZERO <- function(x,y) {
        ##b=tapply(x, y, function(a) {a[a==0]=median(a);a})  ## no NA
        b=tapply(x, y, function(a) {a[a==0]=logmedian(a);a})  ## no NA
        names(b) <- NULL
        unlist(b)[names(x)]
    }    
    X[which(is.na(X))] <- 0
    impX <- t(apply(X,1,medianImputeZERO,y=groups))
    return(impX)
}

###IMPUTATION by downshifted Gaussian
.RGimputation <- function(x, width=0.3, downshift=1.8, bycolumn=T){
    if(bycolumn==T){
        for(i in 1:ncol(x)){
            x[,i][is.na(x[,i])] <- rnorm(sum(is.na(x[,i])), mean(as.numeric(x[,i]), na.rm=T) - downshift*sd(as.numeric(x[,i]), na.rm=T), width*sd(as.numeric(x[,i]), na.rm=T))
        }
    } else {
        x[is.na(x)] <- rnorm(sum(is.na(x)), mean(as.numeric(x), na.rm=T) - downshift*sd(as.numeric(x), na.rm=T), width*sd(as.numeric(x), na.rm=T))
    }
    return(x)
}

prot.calcCopyNumber <- function(data, mol, y){
    ## data should be a numeric matrix
    ## mol is the molecular weight
    ## y is the mass of the cell in PICOGRAM (pg)
    TotalIntenisty <- apply(data, 2, sum)
    ## mass in gram
    zz <- vector(length=0)
    Mass <- vector(length=0)
    MassProtein <-	for(i in 1:length(TotalIntenisty)){
                            zz <- (data[,i] * y) / TotalIntenisty[i] * 10^-12
                            Mass <- cbind(Mass, zz)
                        }
    colnames(Mass) = colnames(data)
                                        # calculate moles
    Mol <- Mass/(mol*1000)
    Copy <- Mol * 6.022140857*10^23
    return(Copy)
}


##======================================================================
##=================== SPECIAL SILAC FUNCTIONS ==========================
##======================================================================
##
##

##datafile="proteinGroups.txt"
silac.readDataFile <- function(datafile, remove.outliers=TRUE) {

    cat(">>> reading datafile from",datafile,"\n")
    df <- read.delim(datafile)
    rownames(df) <- paste0("tag",1:nrow(df))  ## give proper rownames

    ## ------------ Filter contaminants (important)
    keep <- (df$Reverse !="+" & df$Only.identified.by.site !="+" & df$Potential.contaminant !="+")
    table(keep)
    df1 <- df[keep,]

    ## use LFQ values that are normalized across the samples
    ##dfLFQ <- df1[,grep("LFQ", colnames(df1))]
    LFQ.L <- as.matrix( df1[,grep("LFQ.intensity.L.", colnames(df1))] )
    LFQ.H <- as.matrix( df1[,grep("LFQ.intensity.H.", colnames(df1))] )
    table( colnames(LFQ.L) ==  sub(".H.",".L.",colnames(LFQ.H) ))  ## check!

    ## --------- rename samples
    ## names3 <- colnames(LFQ.L)
    ## names3.1 <- sub("LFQ.intensity.L.170905_MB_RG_", "", names3)
    ## names3.2 <- sub("LFQ.intensity.L.180222_MB017_", "", names3.1)
    ## names3.3 <- sub("LFQ.intensity.L.[[:digit:]]+_MB009_", "", names3.2)
    ## LFQ.names <- names3.3
    LFQ.names <- gsub(".*MB[0-9]*_|.*MB_RG_|\"","",colnames(LFQ.L))

    ## Use new sample nomenclature
    names.new.df <- read.delim("samples_renamed.txt", stringsAsFactors=FALSE)
    ##names.new <- names.new.df$Renamed
    sum(duplicated(names.new.df$Renamed))

    ## match names
    org.names <- gsub(".*MB[0-9]*_|.*MB_RG_|\"","",names.new.df$Original)
    org.names <- sub("Rest","Resting",org.names)
    names.new.df2 <- names.new.df[match(LFQ.names, org.names),]

    ## clean names.... (OMG)
    clean.names <-  names.new.df2$Renamed
    clean.names <- gsub("_SILAC","-SILAC",clean.names)
    clean.names <- gsub("_Naive","-Naive",clean.names)
    clean.names <- gsub("_Rest","-Rest",clean.names)
    clean.names <- gsub("_Tcell","-Tcell",clean.names)
    clean.names <- gsub("3-Methyladenosine","3Methyladenosine",clean.names)
    ##clean.names <- gsub("^DoXYY","DoXY",clean.names)
    clean.names <- gsub("^DorE","DoE",clean.names)
    clean.names <- gsub("_","",clean.names)

    ## write.csv( cbind(LFQ.names=LFQ.names,
    ##                  new.names=names.new.df$Renamed,
    ##                  matched.names=names.new.df2$Renamed,
    ##                  clean.names=clean.names ), file="checknames.csv")

    ## check duplicated
    ndup <- sum(duplicated(clean.names))
    ndup
    if(ndup >0) {
        stop("FATAL. Duplicated clean names.")
    }

    ## update names
    colnames(LFQ.H) <- clean.names
    colnames(LFQ.L) <- clean.names

    ## ------------ annotation table (keep separate) ---------
    ## Include proteinID, gene names and MW in dataframe A
    proteins <- df1[, c("Protein.IDs", "Gene.names", "Mol..weight..kDa.")]
    dim(proteins)

    ## ------------ filter samples  ---------------------------------------
    ## Remove wrong sample and bad donors
    if(remove.outliers) {
        bad.samples <- grep("Do109|Do4221-Tcell-Mem",colnames(LFQ.H),value=TRUE)
        keep.samples <- setdiff( colnames(LFQ.H), bad.samples )
        LFQ.L <- LFQ.L[,keep.samples]
        LFQ.H <- LFQ.H[,keep.samples]
    }

    ## ------------ filter probes ---------------------------------------
    ## Delete Factors from the list with severe miss identifications in the control group
    LFQ.ratio <- LFQ.H / (LFQ.H + LFQ.L + 1e-8)  ## only temporary, later we calculate again
    Ctrls <- LFQ.ratio[,grep("Treat=NO-SILAC=0h", colnames(LFQ.ratio))]
    dim(Ctrls)
    ##Ctrls$tot <- apply(Ctrls, 1, sum)
    ##Ctrls$Gene.names <- P$Gene.names
    ##Ctrls.order = Ctrls[order(Ctrls$tot, decreasing=T),]
    ##P.cut <- P[which(Ctrls$tot<1.9),]
    keep <- (rowSums(Ctrls) < 1.9)
    ##keep <- (rowMeans(Ctrls < 0.2) > 0.8)
    table(keep)
    LFQ.L <- LFQ.L[keep,]
    LFQ.H <- LFQ.H[keep,]
    proteins <- proteins[keep,]

    ## ------------ collapse by gene??
    gene1 <- as.character(proteins$Gene.names)
    ##gene1 <- sapply(strsplit(gene1,split=";"),"[",1)
    sum(duplicated(gene1))
    LFQ.L <- apply( LFQ.L, 2, function(x) tapply(x, gene1, sum, na.rm=TRUE))
    LFQ.H <- apply( LFQ.H, 2, function(x) tapply(x, gene1, sum, na.rm=TRUE))
    jj <- match(rownames(LFQ.L),gene1)
    proteins <- proteins[jj,]
    rownames(proteins) <- rownames(LFQ.L)
    dim(LFQ.L)

    ## ------------ create sample annotation
    ##head( colnames(LFQ.H) )
    samples <- sapply( colnames(LFQ.H), strsplit, split="-")
    samples <- data.frame( do.call(rbind, samples) )
    colnames(samples) <- c("donor","cell.type","subtype","state","treatment","SILAC")
    rownames(samples) <- colnames(LFQ.H)  ## should be already

    samples$treatment <- gsub("Treat=","",samples$treatment)
    samples$SILAC <- gsub("SILAC=","",samples$SILAC)
    rownames(samples)

    ## add mass in picograms (PLEASE CHECK!!!)
    samples$mass.pg <- 25
    samples$mass.pg[which(samples$state %in% c("Act23h","Act48h") )] <- 75

    dim(samples)
    Matrix::head(samples)
    apply(samples,2,table)

    ## ------------ define groups
    groups <- list()
    nav <- ( samples$subtype == "Naive" & samples$state == "Rest" & samples$treatment == "NO" )
    groups[["NaiveRest_0h"]]  <- which( nav & samples$SILAC == "0h")
    groups[["NaiveRest_6h"]]  <- which( nav & samples$SILAC == "6h")
    groups[["NaiveRest_12h"]] <- which( nav & samples$SILAC == "12h")
    groups[["NaiveRest_24h"]] <- which( nav & samples$SILAC == "24h")
    groups[["NaiveRest_48h"]] <- which( nav & samples$SILAC == "48h")

    mem <- ( samples$subtype == "Mem" & samples$state == "Rest" & samples$treatment == "NO" )
    groups[["MemRest_0h"]]  <- which( mem & samples$SILAC == "0h")
    groups[["MemRest_6h"]]  <- which( mem & samples$SILAC == "6h")
    groups[["MemRest_12h"]] <- which( mem & samples$SILAC == "12h")
    groups[["MemRest_24h"]] <- which( mem & samples$SILAC == "24h")
    groups[["MemRest_48h"]] <- which( mem & samples$SILAC == "48h")

    groups[["NaiveRest"]] <- which(nav)
    groups[["MemRest"]]   <- which(mem)

    ## special groups
    groups[["NaiveRest_24h_Baf"]] <- which( samples$subtype == "Naive" & samples$state == "Rest" &
                                            samples$treatment == "Bafilomycin24h" )
    groups[["NaiveRestWashout"]] <- which( samples$subtype == "Naive" & samples$state == "Rest" &
                                           samples$treatment == "CHX24hwashout24h" )
    groups[["NaiveRestCHX"]] <- which( samples$subtype == "Naive" & samples$state == "Rest" &
                                       samples$treatment == "CHX24h" )
    groups[["MemoryRest_24h_IL2"]] <- which( samples$subtype == "Mem" & samples$state == "Rest" &
                                             samples$treatment == "IL2" )
    groups[["MemRestWashout"]] <- which( samples$subtype == "Mem" & samples$state == "Rest" &
                                         samples$treatment == "CHX24hwashout24h" )
    groups[["MemRestCHX"]] <- which( samples$subtype == "Mem" & samples$state == "Rest" &
                                     samples$treatment == "CHX24h" )

    ## IMPORTANT:: immediately convert indices to names. much safer!!!
    groups <- lapply(groups, function(ii) rownames(samples)[ii])
    sapply(groups, length)  ## ALWAYS CHECK!!!


    ## ------------ add derived quantities: ratio
    LFQ.total <- (LFQ.L + LFQ.H)   ## total expression
    LFQ.ratio <- LFQ.H / (LFQ.L + LFQ.H + 1e-8)  ## notice adding small  number to avoid NaN

    molwt <- proteins$Mol..weight..kDa.
    copy.number <- silac.calcCopyNumber( data=LFQ.total, mol.weight=molwt, y=samples$mass.pg )

    dim(copy.number)

    sum(is.nan(LFQ.ratio))
    dim(LFQ.ratio)

    ## prepare output object
    output <- list( samples=samples, groups=groups, proteins=proteins,
                   LFQ.ratio=LFQ.ratio, LFQ.H=LFQ.H, LFQ.L=LFQ.L,
                   copy.number=copy.number)

    return(output)
}

silac.ttest <- function(X, group1, group2, method="limma") {

    ##X <- silac$LFQ.ratio
    ##method = "geiger"
    out1 <- NULL ## important to avoid old results
    p1=p2=p3=NULL  ## just for testing
    if(method == "geiger") {
        ## Your code
        out1 <- volcano(X[,group1], X[,group2], rownames(X))
        out1$Gene <- NULL
        colnames(out1) <- c("logFC","P.Value")
        dim(out1)
        Matrix::head(out1)
        p1 <- out1
    } else if(method == "genefilter") {
        ## faster t-test
        
        jj <- which( colnames(X) %in% c(group1,group2))
        y <- 1*(colnames(X)[jj] %in% group1) + 2*(colnames(X)[jj] %in% group2)
        out1 <- genefilter::rowttests( X[,jj], factor(y) )
        out1$statistic <- NULL
        colnames(out1) <- c("logFC","P.Value")
        Matrix::head(out1)
        p2 <- out1
    } else if(method == "limma") {
        ## See e.g. https://bioconductor.org/help/course-materials/2010/BioC2010/limma2.pdf
        
        jj <- which( colnames(X) %in% c(group1,group2))
        y <- 1*(colnames(X)[jj] %in% group1)
        design <- model.matrix( ~ y )
        colnames(design) <- c("Intercept","group1_vs_group2")
        fit <- limma::eBayes(limma::lmFit(X[,jj], design))
        out1 <- limma::topTable(fit, coef=2, sort.by="none", number=Inf)
        out1 <- out1[,c("logFC","P.Value","adj.P.Val")]
        Matrix::head(out1)
        p3 <- out1
    } else {
        stop("ERROR:: unknown method")
    }

    if(0) {
        ## testing area, Show the differences in p-value of methods.
        pv <- cbind( t.test=p1$P.Value, t.genefilter=p2$P.Value, limma=p3$P.Value )
        jj <- which( !is.nan(p1$P.Value))
        pairs(pv[jj,])
    }

    ## set NaN p.values to 1??
    out1$P.Value[which(is.nan(out1$P.Value))] <- 1

    ## add FDR q.value
    if(!("adj.P.Val" %in% names(out1))) {
        out1$adj.P.Val <- p.adjust( out1$P.Value, method="fdr")
    }

    ## good practice to add group means
    avg <- cbind( rowMeans(X[,group1],na.rm=TRUE), rowMeans(X[,group2],na.rm=TRUE))
    ##colnames(avg) <- paste0("avg.",ct)
    colnames(avg) <- paste0("avg.",c("group1","group2"))
    out1 <- cbind(out1, avg)

    return(out1)
}

##obj=silac
silac.volcanoPlot <- function(obj, group1, group2, psig=0.05, lfc=0.2) {
    if(class(group1[1])=="character" && group1[1] %in% names(obj$groups)) {
        group1 <- obj$groups[[group1]]
    }
    if(class(group2[1])=="character" && group2[1] %in% names(obj$groups)) {
        group2 <- obj$groups[[group2]]
    }

    ##Test differences naive memory
    P <- as.matrix(obj$LFQ.ratio)
    group1 <- intersect(group1, colnames(P))
    group2 <- intersect(group2, colnames(P))
    genes <- as.character(obj$proteins$Gene.names)
    pdata <- volcano( P[, group1], P[, group2], genes)
    jj <- which( !is.na(pdata$P) & !is.nan(pdata$P))
    sig <- which(pdata$P < psig & abs(pdata$EFFECTSIZE) > lfc)
    highlight <- pdata$Gene[sig]
    ##pdata$GENE <- rownames(pdata)
    volcanoly( pdata[jj,], snp="Gene", highlight=highlight,
              effect_size_line = c(-1,1)*lfc,
              genomewideline = -log10(0.05) )
    invisible(pdata)
}


silac.calcCopyNumber <- function(data, mol.weight, y){
    ## data should be a numeric matrix
    ## mol is the molecular weight
    ## y is the mass of the cell in PICOGRAM (pg)
    total.intensity <- Matrix::colSums(data,na.rm=TRUE)
    mass <- t( t(data) * y / total.intensity) * 1e-12
    ##colnames(mass) = colnames(data)
    ## calculate moles
    cn <- mass / (mol.weight *1000) * 6.022140857e23
    return(cn)
}


###
# FIT WIBULL
###
##x=timeN;y=P[j,naive]
fit.weibull2 <- function(x,y) {
    y = as.numeric(y)
    x = as.numeric(x)
    if(all(y==0)) {
        xfit = seq(0,48,1)
        return(list(x=xfit, y=xfit*0, t50=Inf))
    }
    jj <- which( x==0 | y>0)  ## skip zero measurements (outlier?) if not x==0
    x = x[jj]
    y = y[jj]
    cdf.weibull <- function(x, lambda, k) (1 - exp(-(x/lambda)^k))
    fit = nls( y ~ cdf.weibull(x,lambda,k),
              start=list(lambda=50,k=1),
              lower=list(lambda=0.01,k=0.001),
              algorithm="port")
    summary(fit)
    coef(fit)
    xfit = seq(0,48,1)
    lambda0 = coef(fit)["lambda"]
    k0 = coef(fit)["k"]
    yfit = cdf.weibull(xfit, lambda0, k0)
    t50 = lambda0*log(2)**(1/k0)  ## time at 50%
    t50
    list(x=xfit, y=yfit, t50=t50)
}

##obj=silac;samples=NULL;protein=p
silac.fitWeibull <- function(obj, protein, samples=NULL, samples2=NULL)
{
    ##samples <- obj$groups[["NaiveRest"]]
    ##samples2 <- obj$groups[["MemRest"]]
    if(!is.null(samples) && class(samples[1])=="character" &&
       samples[1] %in% names(obj$groups)) {
        samples <- obj$groups[[samples]]
    }
    if(!is.null(samples2) && class(samples2[1])=="character" &&
       samples2[1] %in% names(obj$groups)) {
        samples2 <- obj$groups[[samples2]]
    }
    if(is.null(samples)) samples <- rownames(obj$samples)
    timeN <- obj$samples[samples,"SILAC"]
    timeN <- as.integer(sub("h$","",timeN))
    timeN

    ##protein="LGALS3"
    ## Nice examples are TCF7, ENO1, GAPDH, FOXO1, B2M, CD5, SQSTM1, CD3E, CXCR4, FOXP1, RPS9 as examples
    ##p=grep("^LGALS3$", P$Gene.names)
    tag <- rownames(obj$proteins)[which(obj$proteins$Gene.names == protein)]
    if(length(tag)>1) cat("warning:: protein has multiple tags!\n")
    tag <- tag[1]  ## really???
    q1 <- obj$LFQ.ratio[tag,samples]
    ##q1 <- colMeans(obj$LFQ.ratio[tag,samples,drop=FALSE])
    plot(timeN, q1, pch=17, col="blue",
         ylim=c(0,1), ylab="", xlab="Time [h]", main=protein,
         cex=1.4, xaxt="n", yaxt="n", cex.lab=1.5, cex.main=1.4)
    axis(1,at=c(0,6,12,24,48),  las=0, cex.axis=1.5)
    axis(2,at=c(0,0.25,0.5,0.75,1),  las=2, cex.axis=1.5)
    ##par(new=T)
    abline(h=0.5, lty=3)
    abline(h=0.25, lty=3)
    abline(h=0.75, lty=3)
    fit1 = fit.weibull2(timeN, q1)
    lines(fit1$x, fit1$y, lty=1, col="blue")
    tt = paste0("t50= ",round(fit1$t50,digits=2),"h")

    fit2=NULL
    if(!is.null(samples2)) {
        timeM <- obj$samples[samples2,"SILAC"]
        timeM <- as.integer(sub("h$","",timeM))
        timeM
        ##q2 <- colMeans(obj$LFQ.ratio[tag,samples2,drop=FALSE])
        q2 <- obj$LFQ.ratio[tag,samples2]
        points(timeM, q2, pch=19, col="#EF4136")
        fit2 = fit.weibull2(timeM, obj$LFQ.ratio[tag,samples2])
        lines(fit2$x, fit2$y, lty=1, col="#EF4136")
        tt = c(tt, paste0("t50= ",round(fit2$t50,digits=2),"h"))
    }

    legend("bottomright",legend=tt, pch=c(17,19), col=c("blue","#EF4136"), cex=1, bty="n")

    res <- list( fit=fit1, fit2=fit2)
    return(res)
}

##obj=data;samples=NULL
silac.plotProteins <- function(obj, proteins, samples=NULL) {
    if(!is.null(samples) && class(samples[1])=="character" &&
       samples[1] %in% names(obj$groups)) {
        samples <- obj$groups[[samples]]
    }
    ##proteins <- c("CD3E", "CD3G", "CD3D", "CD247","TRBV28;TRBC2", "TRBC1")
    ss <- match(proteins, obj$proteins$Gene.names)
    ss <- setdiff(ss, NA)
    pp <- rownames(obj$proteins)[ss]
    Q <- obj$LFQ.ratio
    if(is.null(samples)) samples <- colnames(Q)
    plot.df <- data.frame(Q[pp, samples ])
    row.names(plot.df) <- proteins
    PL <- t(plot.df)
    boxplot(t(plot.df), ylim=c(0,1), col="darkgrey", las=3)
    stripchart(list(PL[,1],PL[,2],PL[,3],PL[,4]), vertical=T, add=T, pch=1, method="jitter", cex=1.5)
}

##obj=silac;samples="NaiveRest_6h";main=NULL;minq=3
silac.plotDistribution <- function(obj, samples, minq=3, main=NULL)
{
    if(!is.null(samples) && class(samples[1])=="character" &&
       samples[1] %in% names(obj$groups)) {
        if(is.null(main)) main <- samples
        samples <- obj$groups[[samples]]
    }
    ##Ranking according to the average at 24h
    Q <- as.matrix(obj$LFQ.ratio[,samples])
    valid <- (rowSums(Q>0) >= minq)
    Q <- Q[valid,]
    dim(Q)
    meanQ <- rowMeans(Q, na.rm=TRUE)
    plot(sort(meanQ,decreasing=TRUE), main=main)
    hist( meanQ, breaks=20, col="darkgrey", main=main)
}
