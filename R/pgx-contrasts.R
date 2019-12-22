
##-----------------------------------------------------------------------------
## LIMMA  helper functions
##-----------------------------------------------------------------------------

normalizeTMM <- function(counts, log=FALSE, method="TMM") {
    require(edgeR)
    dge <- DGEList(as.matrix(counts), group=NULL)
    dge <- calcNormFactors(dge, method=method)
    edgeR::cpm(dge, log=log)
}
normalizeRLE <- function(counts, log=FALSE) {
    require(edgeR)
    dge <- DGEList(as.matrix(counts), group=NULL)
    dge <- calcNormFactors(dge, method="RLE")
    edgeR::cpm(dge, log=log)
}

normalizeContrasts <- function(contr.matrix) {
    ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
    for(i in 1:ncol(contr.matrix)) {
        m <- contr.matrix[,i]
        m[is.na(m)] <- 0
        contr.matrix[,i] <- 1*(m>0)/sum(m>0) - 1*(m<0)/sum(m<0)
    }
    contr.matrix
}

makeContrastsFromPairs <- function(main.group, ref.group, groups=NULL, comparisons=NULL) {

    if(is.null(comparisons)) comparisons <- paste0(main.group,"_vs_",ref.group)
    ##main.group <- sapply(strsplit(comparisons, split=split),"[",1)
    ##ref.group  <- sapply(strsplit(comparisons, split=split),"[",2)
    main.group <- as.character(main.group)
    ref.group <- as.character(ref.group)
    groups1 <- unlist(strsplit(main.group, split="[+]"))
    groups2 <- unlist(strsplit(ref.group, split="[+]"))
    if(is.null(groups)) groups <- sort(unique(c(groups1, groups2)))
    contr.matrix <- matrix(0, nrow=length(groups), ncol=length(comparisons))
    i=1
    for(i in 1:length(comparisons)) {
        main <- strsplit(main.group[i], split="[+]")[[1]]
        ref  <- strsplit(ref.group[i], split="[+]")[[1]]
        ct <- 1*(groups %in% main) - 1*(groups  %in% ref)
        ct
        contr.matrix[,i] <- ct
    }
    colnames(contr.matrix) <- comparisons
    rownames(contr.matrix) <- groups

    ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
    for(i in 1:ncol(contr.matrix)) {
        m <- contr.matrix[,i]
        m[is.na(m)] <- 0
        contr.matrix[,i] <- 1*(m>0)/sum(m>0) - 1*(m<0)/sum(m<0)
    }
    contr.matrix
}

makeDirectContrasts <- function(Y, ref, na.rm=TRUE) {
    contr.matrix <- c()
    i=2
    for(i in 1:ncol(Y)) {
        m1 <- NULL
        if(!is.na(ref[i]) && ref[i]!="all" ) {
            x <- as.character(Y[,i])
            x[is.na(x)|x=="NA"] <- "_"
            m1 <- model.matrix( ~ 0 + x)
            colnames(m1) <- sub("^x","",colnames(m1))
            m1 <- m1 - m1[,ref[i]]  ## +1/-1 encoding
            m1 <- m1[,which(colnames(m1)!=ref[i]),drop=FALSE]  ## remove refvsref...
            m1 <- m1[,!colnames(m1) %in% c("NA","_"),drop=FALSE]
            colnames(m1) <- paste0(colnames(m1),"_vs_",ref[i])
        } else if(!is.na(ref[i]) && ref[i]=="all" ) {
            x <- as.character(Y[,i])
            x[is.na(x)|x=="NA"] <- "_"            
            m1 <- model.matrix( ~ 0 + x)
            colnames(m1) <- sub("^x","",colnames(m1))
            ##m1 <- m1 - m1[,ref[i]]  ## +1/-1 encoding
            m1 <- t(t(m1==1) / colSums(m1==1) - t(m1==0) / colSums(m1==0))
            ##m1 <- m1[,which(colnames(m1)!=ref[i]),drop=FALSE]
            m1 <- m1[,!colnames(m1) %in% c("NA","_"),drop=FALSE]            
            colnames(m1) <- paste0(colnames(m1),"_vs_all")
        } else {
            levels <- names(table(Y[,i]))
            levels <- setdiff(levels, c(NA,"NA"))
            levels
            if(length(levels)>1) {
                cc <- makeFullContrasts(levels)
                x <- as.character(Y[,i])
                x[is.na(x)] <- "_"
                mm <- model.matrix( ~ 0 + x)
                colnames(mm) <- sub("^x","",colnames(mm))
                m1 <- mm[,rownames(cc)] %*% cc
            }
        }
        if(!is.null(m1)) {
            colnames(m1) <- paste0(colnames(Y)[i],":",colnames(m1))
            contr.matrix <- cbind(contr.matrix, m1)
        }
    }
    ##colnames(contr.matrix) <- colnames(Y)
    ##colnames(contr.matrix) <- paste0(colnames(Y),":",colnames(contr.matrix))

    ## take out any empty comparisons
    contr.matrix <- contr.matrix[,which(colSums(contr.matrix!=0)>0),drop=FALSE]

    ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
    for(i in 1:ncol(contr.matrix)) {
        m <- contr.matrix[,i]
        m[is.na(m)] <- 0
        contr.matrix[,i] <- 1*(m>0)/sum(m>0) - 1*(m<0)/sum(m<0)
    }
    rownames(contr.matrix) <- rownames(Y)
    contr.matrix
}

makeFullContrasts <- function(levels, by.sample=FALSE) {
    levels <- sort(unique(as.character(levels)))
    cc <- t(combn(levels,2))
    contr.matrix <- c()
    for(i in nrow(cc):1) {
        ctr <- 1*(levels==cc[i,1]) - 1*(levels==cc[i,2])
        contr.matrix <- cbind(ctr,contr.matrix)
        colnames(contr.matrix)[1] <- paste(cc[i,],collapse="_vs_")
    }
    rownames(contr.matrix) <- levels
    if(by.sample) {
        design <- model.matrix( ~ 0 + levels )
        colnames(design) <- sub("^levels","",colnames(design))
        rownames(design) <- names(clusters)
        design <- design[,rownames(contr.matrix)]
        contr.matrix <- design %*% contr.matrix
    }
    return(contr.matrix)
}

makeClusterContrasts <- function(clusters, min.freq=0.01, full=FALSE,
                                 by.sample=FALSE ) {
    ## make model matrix for cluster_i vs. rest
    if(0) {
        min.size <- pmax(3, 0.01*length(clusters))
        small.clusters <- names(which(table(clusters) < min.size))
        clusters[ which(clusters %in% small.clusters)] <- "cluster0"
        sort(table(clusters))
    }
    idx <- sort(unique(as.character(clusters)))
    m1 <- model.matrix( ~ 0 + idx)
    colnames(m1) <- sub("^idx","",colnames(m1))
    rownames(m1) <- colnames(m1)
    colnames(m1) <- paste0(colnames(m1),"_vs_all")
    ##m1 <- m1 - 1/(nrow(m1)-1)*(m1==0)
    m1 <- t(t(m1==1) / colSums(m1==1) - t(m1==0) / colSums(m1==0))

    diag(m1) <- 1
    if(full==TRUE) {
        m2 <- makeFullContrasts(unique(clusters))
        m2 <- m2[rownames(m1),]
        m1 <- cbind(m1, m2)
    }
    if(by.sample) {
        design <- model.matrix( ~ 0 + clusters )
        colnames(design) <- sub("^clusters","",colnames(design))
        rownames(design) <- names(clusters)
        design <- design[,rownames(m1)]
        m1 <- design %*% m1
    }
    return(m1)
}


##mingrp=3;slen=8;ref=NULL
pgx.makeAutoContrast <- function(df, mingrp=3, slen=8, ref=NULL) {

    shortestunique <- function(xx,slen=3) {
        k <- min(which(!sapply(1:max(nchar(xx)),function(i) any(duplicated(substring(xx,1,i))))))
        substring(xx,1,max(k,slen))
    }
    
    autoContrast1 <- function(x, ref1, slen, mingrp) {
        if(is.null(ref1)) ref1 <- NA
        x <- as.character(x)
        nx <- table(x)
        too.small <- names(which(nx < mingrp))
        if(length(too.small)) x[which(x %in% too.small)] <- NA
        nx <- table(x)
        if(length(nx)<2) return(NULL)
        x <- factor(x)
        if(!is.na(ref1)) x <- relevel(x, ref=ref1)
        xlevels <- gsub("[^[:alnum:]]","",levels(x))
        levels(x) <- shortestunique(xlevels,slen=slen)
        xref <- gsub("[^[:alnum:]]","",levels(x)[1])
        nn <- length(nx)
        nn
        if(nn<2) {
            return(NULL)
        } else if(nn == 2) {
            ct <- model.matrix(~x)[,2,drop=FALSE]
            colnames(ct) <- paste0(levels(x)[2],"_vs_",levels(x)[1])
        } else if(nn >= 3) {
            if(is.na(ref1)) {
                ct <- model.matrix(~0 + x)
                colnames(ct) <- paste0(levels(x),"_vs_rest")
            } else {
                ct <- model.matrix(~0 + x)
                colnames(ct) <- paste0(levels(x),"_vs_",xref)
                i=1
                for(i in 1:ncol(ct)) {
                    j <- which(!(x %in% levels(x)[c(1,i)]))
                    j <- intersect(as.character(j),rownames(ct))
                    ct[j,i] <- NA
                }
                ct <- ct[,2:ncol(ct),drop=FALSE] ## remove REFvsREF                
            }
        }
        ## NA can make ct smaller than full
        ct <- ct[match(1:length(x),rownames(ct)),,drop=FALSE]
        rownames(ct) <- 1:length(x)
        ct[is.na(ct)] <- 0
        ct
    }

    if(!is.null(ref) && length(ref)!=ncol(df)) ref <- head(rep(ref,99),ncol(df))
    
    K <- c()
    for(i in 1:ncol(df)) {
        ref1 <- NA
        if(!is.null(ref)) ref1 <- ref[i]
        x <- df[,i]
        if(!(ref1 %in% x)) ref1 <- NA
        ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$"
        detect.ref <- any(grepl(ref.pattern,x,ignore.case=TRUE))
        if(is.na(ref1) & detect.ref) {
            ref1 <- grep(ref.pattern,x,ignore.case=TRUE,value=TRUE)
            ref1 <- sort(ref1)[1]
            cat("reference auto-detected:",ref1,"\n")
        }
        ct <- autoContrast1(x, ref=ref1, slen=slen, mingrp=mingrp)
        dim(ct)
        if(!is.null(ct)) {
            colnames(ct) <- paste0(colnames(df)[i],":",colnames(ct))
            K <- cbind(K,ct)
        }
    }
    dim(K)
    
    rownames(K) <- rownames(df)
    head(K)
    
    kcode <- apply(K,1,paste,collapse="-")
    xc <- factor(kcode, levels=unique(kcode))  ## experimental condition
    levels(xc) <- paste0("group",1:length(levels(xc)))

    jj <- which(!duplicated(kcode))
    K2 <- K[jj,,drop=FALSE]
    rownames(K2) <- xc[jj]
    head(K2)

    ## Translate coding 0/NA/1 to -1/0/+1 coding of contrast
    K[K==0] <- -1
    K[is.na(K)] <- 0
    K2[K2==0] <- -1
    K2[is.na(K2)] <- 0
        
    list(group = xc, contr.matrix = K2, exp.matrix=K)
}

##=====================================================================================
##=========================== END OF FILE =============================================
##=====================================================================================
