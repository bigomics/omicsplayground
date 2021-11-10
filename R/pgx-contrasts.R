##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## source(file.path(RDIR,"pgx-functions.R"))

##----------------------------------------------------------------------
## Auto-detection helper functions
##----------------------------------------------------------------------

pgx.getContrastGroups <- function(pgx, contrast, as.factor=TRUE) {
    exp.matrix <- pgx$model.parameters$exp.matrix
    grp <- contrastAsLabels(exp.matrix[,contrast,drop=FALSE],as.factor=as.factor)     
    if(NCOL(grp)==1) {
        grp <- grp[,1]
        names(grp) <- rownames(exp.matrix)
    }
    grp    
}

pgx.detect_batch_params <- function(Y) {
    grep("batch|cell.line|patient|mouse|sample|repl|strain", colnames(Y), value=TRUE)
}

pgx.detect_timevar <- function(Y) {
    has.time <- any(grepl("time|hour|hr|day", colnames(Y)))
    has.time
    ttvar <- grep("time|hour|hr|day", colnames(Y),value=TRUE)
    ttvar
    ##ttvar <- ttvar[sapply(ttvar, function(v) all(table(Y[,v])>1))]
    ttvar
}

pgx.getConditions <- function(exp.matrix, nmax=3) {
    ##
    ##
    ##
    group <- apply(exp.matrix,1,paste,collapse="_")
    table(group)
    group <- factor(group)
    if(ncol(exp.matrix) > nmax) {
        ngroup <- length(unique(group))
        levels(group) <- paste0("group",1:ngroup)
    }
    as.character(group)
}

pgx.expMatrix <- function(pheno, contr.matrix) {

    dbg("[pgx.expMatrix] called")
    
    ctx <- rownames(contr.matrix)
    ## already an experiment contrast
    if(length(ctx)==nrow(pheno) && all(ctx %in% rownames(pheno))) {
        contr.matrix <- contr.matrix[match(rownames(pheno),ctx),,drop=FALSE]
        return(contr.matrix)
    }    
    group.col <- names(which(apply(pheno, 2, function(x) all(ctx %in% x))))[1]
    if(is.null(group.col) || length(group.col)==0 || is.na(group.col) ) {
        message("[pgx.expMatrix] WARNING: could not resolve group column. No exp.matrix\n")
        return(NULL)
    }

    dbg("[pgx.expMatrix] group.col =",group.col)
    grp <- pheno[,group.col]

    dbg("[pgx.expMatrix] dim(contr.matrix) =",paste(dim(contr.matrix)))
    dbg("[pgx.expMatrix] grp =",paste(grp))
    dbg("[pgx.expMatrix] ctx =",paste(ctx))
    
    exp.matrix <- contr.matrix[match(grp,ctx),,drop=FALSE]
    dbg("[pgx.expMatrix] dim(exp.matrix) =",paste(dim(exp.matrix)))

    rownames(exp.matrix) <- rownames(pheno)
    exp.matrix
}

##-----------------------------------------------------------------------------
## Contrast creation functions
##-----------------------------------------------------------------------------

pgx.makeStratifiedContrastsDF <- function(data, vars, strata, ref) {

    dstrata <- data[,strata]
    S <- model.matrix( ~ 0 + dstrata)
    colnames(S) <- sub("^dstrata","",colnames(S))

    df <- cbind(data[,vars,drop=FALSE], strata=dstrata)
    md <- makeDirectContrasts(df,c(ref,"others"))

    M <- md$contr.matrix
    C0 <- M[,grep("strata",colnames(M),invert=TRUE),drop=FALSE]
    C1 <- M[,grep("strata",colnames(M)),drop=FALSE]
    colnames(C1)
    colnames(C1) <- paste0(gsub("strata:|_vs_others","",colnames(C1)))
    contr.matrix <- c()
    i=1
    for(i in 1:ncol(C0)) {
        m1 <- (C1==1) * C0[,i]
        ##colnames(m1) <- paste0(colnames(C1),"::",colnames(C0)[i])
        colnames(m1) <- paste0(colnames(C0)[i],"@",colnames(C1))
        contr.matrix <- cbind(contr.matrix, m1)
    }

    ## construct exp matrix
    grp <- md$group
    G <- model.matrix( ~ 0 + grp)
    colnames(G) <- sub("^grp","",colnames(G))
    G <- G[,rownames(contr.matrix)]
    exp.matrix <- G %*% contr.matrix
    rownames(exp.matrix) <- rownames(data)
    
    ## check levels
    sel <- ( Matrix::colSums(contr.matrix==-1)>0 &
             Matrix::colSums(contr.matrix==+1)>0 )
    table(sel)
    contr.matrix <- contr.matrix[,sel]
    exp.matrix <- exp.matrix[,sel]
    
    res <- list(contr.matrix = contr.matrix,
                group = md$group,
                exp.matrix = exp.matrix)
    res
}

pgx.makeStratifiedContrasts <- function(Y, strata, ref) {
    
    S <- model.matrix( ~ 0 + strata)
    colnames(S) <- sub("^strata","",colnames(S))

    df <- data.frame(Y,strata)
    md <- makeDirectContrasts(df,c(ref,"others"))

    M <- md$contr.matrix
    C0 <- M[,grep("strata",colnames(M),invert=TRUE),drop=FALSE]
    C1 <- M[,grep("strata",colnames(M)),drop=FALSE]
    colnames(C1)
    colnames(C1) <- paste0(gsub("strata:|_vs_others","",colnames(C1)))
    contr.matrix <- c()
    i=1
    for(i in 1:ncol(C0)) {
        m1 <- (C1==1) * C0[,i]
        ##colnames(m1) <- paste0(colnames(C1),"::",colnames(C0)[i])
        colnames(m1) <- paste0(colnames(C0)[i],"@",colnames(C1))
        contr.matrix <- cbind(contr.matrix, m1)
    }

    ## construct exp matrix
    grp <- md$group
    G <- model.matrix( ~ 0 + grp)
    colnames(G) <- sub("^grp","",colnames(G))
    G <- G[,rownames(contr.matrix)]
    exp.matrix <- G %*% contr.matrix
    rownames(exp.matrix) <- rownames(Y)
    
    ## check levels
    sel <- ( Matrix::colSums(contr.matrix==-1)>0 &
             Matrix::colSums(contr.matrix==+1)>0 )
    table(sel)
    contr.matrix <- contr.matrix[,sel]
    exp.matrix <- exp.matrix[,sel]
    
    
    res <- list(contr.matrix = contr.matrix,
                group = md$group,
                exp.matrix = exp.matrix)
    res
}

## for compatibility...
makeDirectContrasts2 <- function(Y, ref, na.rm=TRUE) {
    makeDirectContrasts(Y=Y, ref=ref, na.rm=na.rm)
}

expmat2contrast <- function(exp.matrix) {
    group <- apply(exp.matrix,1,paste,collapse="_")
    table(group)
    n.group <- length(unique(group))
    group <- factor(group)
    if(ncol(exp.matrix)>3) {
        levels(group) <- paste0("group",1:n.group)
    }
    sel <- !duplicated(group)
    contr.matrix <- exp.matrix[sel,]
    rownames(contr.matrix) <- group[sel]
    group
    list(contr.matrix=contr.matrix, group=group, exp.matrix=exp.matrix)
}

##Y=ngs$samples;na.rm=TRUE
makeDirectContrasts <- function(Y, ref, na.rm=TRUE)
{

    ## check enough levels
    nlevel <- apply(Y,2,function(y) length(unique(y)))
    if(any(nlevel<2)) {
        notlevely <- colnames(Y)[which(nlevel<2)]
        message("warning:: not enough levels: ",notlevely)
    }
    ii <- which(nlevel>1)
    ii
    if(length(ii)==0) {
        message("warning:: no valid phenotypes")
        return(NULL)
    }

    
    Y <- Y[,ii,drop=FALSE]
    ref <- ref[ii]

    ## make contrast
    exp.matrix <- makeDirectContrasts000(Y=Y, ref=ref, na.rm=na.rm, warn=FALSE) 
    exp.matrix <- sign(exp.matrix)    
    no.vs <- grep("_vs_|_VS_",colnames(exp.matrix),invert=TRUE)
    no.vs
    if(length(no.vs)>0) {
        colnames(exp.matrix)[no.vs] <- paste0(colnames(exp.matrix)[no.vs],":Y_vs_N")
    }    
    exp.matrix0 <- exp.matrix
    if(all(grepl("_vs_|_VS_",colnames(exp.matrix0)))) {
        exp.matrix0 <- contrastAsLabels(exp.matrix0)
    }
    group <- pgx.getConditions(exp.matrix0)
    table(group)
    if(length(levels(group)) > 0.5*nrow(exp.matrix)) {
        cat("WARNING:: contrast matrix looks degenerate. consider removing a contrast.\n")
    }

    contr.matrix <- exp.matrix[which(!duplicated(group)),,drop=FALSE]
    rownames(contr.matrix) <- group[which(!duplicated(group))]
    
    list(contr.matrix=contr.matrix, group=group, exp.matrix=exp.matrix)
}

##Y=pgx$samples$celltype;ref='*'
makeDirectContrasts000 <- function(Y, ref, na.rm=TRUE, warn=FALSE) {
    ## if(warn) warning("makeDirectContrasts is deprectated. please use makeDirectContrasts2()")
    if(NCOL(Y)==1) Y <- data.frame(Y=Y)

    ## check
    all <- c("all","other","others","rest")
    full <- c('*','full')
    has.ref <- rep(NA,ncol(Y))
    for(i in 1:ncol(Y)) has.ref[i] <- (ref[i] %in% Y[,i] || ref[i] %in% c(all,full))
    has.ref
    if(!all(has.ref)) {
        stop("ERROR:: reference ", which(!has.ref), " not in phenotype matrix\n")
        return(NULL)
    }

    contr.matrix <- c()
    if(length(ref)<ncol(Y)) ref <- Matrix::head(rep(ref,99),ncol(Y))
    ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$|scrambl|none|dmso|vehicle"
    i=1    
    for(i in 1:ncol(Y)) {
        m1 <- NULL
        ref1 <- ref[i]
        ## if(ref1 %in% c("*","full")) ref1 <- NA  ##???
        x <- as.character(Y[,i])
        x[is.na(x)|x=="NA"] <- "_"
        detect.ref <- any(grepl(ref.pattern,x,ignore.case=TRUE))
        detect.ref
        if(is.na(ref1) && detect.ref) {
            ref1 <- grep(ref.pattern,x,ignore.case=TRUE,value=TRUE)
            ref1 <- sort(ref1)[1]
            cat("reference auto-detected:",ref1,"\n")
        }
        cref <- as.character(ref1)
        m1 <- model.matrix( ~ 0 + x)
        colnames(m1) <- sub("^x","",colnames(m1))
        if(ref1 %in% full) {
            levels <- names(table(x))
            levels <- setdiff(levels, c(NA,"NA"))
            levels
            if(length(levels)>1) {
                cc <- makeFullContrasts(levels)
                m1 <- m1[,rownames(cc)] %*% cc
            }
        } else if(!is.na(ref1) && !(ref1 %in% all) ) {
            m1 <- m1 - m1[,cref]  ## +1/-1 encoding
            m1 <- m1[,which(colnames(m1)!=cref),drop=FALSE]  ## remove refvsref...
            m1 <- m1[,!colnames(m1) %in% c("NA","_"),drop=FALSE]
            colnames(m1) <- paste0(colnames(m1),"_vs_",ref1)
        } else if(!is.na(ref1) && (ref1 %in% all) ) {
            ##m1 <- m1 - m1[,cref]  ## +1/-1 encoding
            m1 <- t(t(m1==1) / Matrix::colSums(m1==1) - t(m1==0) / Matrix::colSums(m1==0))
            ##m1 <- m1[,which(colnames(m1)!=cref),drop=FALSE]
            m1 <- m1[,!colnames(m1) %in% c("NA","_"),drop=FALSE]            
            colnames(m1) <- paste0(colnames(m1),"_vs_others")
        } else {
            stop("[makeDirectContrasts000] FATAL")
        }
        if(!is.null(m1)) {
            mm <- gsub("[: ]","_",colnames(Y)[i])
            colnames(m1) <- paste0(mm,":",colnames(m1))
            contr.matrix <- cbind(contr.matrix, m1)
        }
    }

    ## take out any empty comparisons
    contr.matrix <- contr.matrix[,which(Matrix::colSums(contr.matrix!=0)>0),drop=FALSE]

    ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
    for(i in 1:ncol(contr.matrix)) {
        m <- contr.matrix[,i]
        m[is.na(m)] <- 0
        contr.matrix[,i] <- 1*(m>0)/sum(m>0) - 1*(m<0)/sum(m<0)
    }
    rownames(contr.matrix) <- rownames(Y)
    sign(contr.matrix)
}

makeFullContrasts <- function(labels, by.sample=FALSE) {
    levels <- sort(unique(as.character(labels)))
    cc <- t(combn(levels,2))
    contr.matrix <- c()
    for(i in nrow(cc):1) {
        ctr <- 1*(levels==cc[i,1]) - 1*(levels==cc[i,2])
        contr.matrix <- cbind(ctr,contr.matrix)
        colnames(contr.matrix)[1] <- paste(cc[i,],collapse="_vs_")
    }
    rownames(contr.matrix) <- levels
    if(by.sample) {
        design <- model.matrix( ~ 0 + labels )
        colnames(design) <- sub("^labels","",colnames(design))
        rownames(design) <- names(labels)
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
    colnames(m1) <- paste0(colnames(m1),"_vs_others")
    ##m1 <- m1 - 1/(nrow(m1)-1)*(m1==0)
    m1 <- t(t(m1==1) / Matrix::colSums(m1==1) - t(m1==0) / Matrix::colSums(m1==0))

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

##mingrp=3;contrasts=c("genotype:mut_vs_WT")
pgx.makeSpecificContrasts <- function(df, contrasts, mingrp=3)
{

    K <- c()
    i=1
    for(i in 1:length(contrasts)) {
        ct <- contrasts[i]
        if(!(grepl("[:]",ct) && grepl("_vs_",ct))) next()
        ph <- sub("[:].*","",ct)
        if(!ph %in% colnames(df)) next()
        
        groups <- strsplit(sub(".*[:]","",ct),split="_vs_")[[1]]
        y <- -1*(df[,ph] == groups[2])  + 1*(df[,ph] == groups[1])
        K <- cbind(K, y)
        colnames(K)[ncol(K)] <- ct
    }
    rownames(K) <- rownames(df)

    K0 <- contrastAsLabels(K)
    group <- pgx.getConditions(K0)
    table(group)
    if(length(levels(group)) > 0.5*nrow(K)) {
        cat("WARNING:: contrast matrix looks degenerate. consider removing a contrast.\n")
    }

    contr.matrix <- K[which(!duplicated(group)),,drop=FALSE]
    rownames(contr.matrix) <- group[which(!duplicated(group))]    
    res <- list(contr.matrix=contr.matrix, group=group)
    
    return(res)
}

## mingrp=3;slen=20;ref=NULL;fix.degenerate=FALSE;skip.hidden=TRUE
pgx.makeAutoContrastsStratified <- function(df, strata.var, mingrp=3, slen=20, ref=NULL, 
                                           fix.degenerate=FALSE, skip.hidden=TRUE )
{

    df1 <- df[,-match(strata.var,colnames(df)),drop=FALSE]
    strata <- df[,strata.var]
    strata.levels <- unique(strata)
    s=strata.levels[1]
    ct.all <- NULL
    for(s in strata.levels) {
        sel <- which(strata == s)
        if(length(sel) < mingrp) next
        suppressWarnings(
            ct1 <- pgx.makeAutoContrasts( df1[sel,,drop=FALSE], mingrp=mingrp, slen=slen,
                                         ref=ref, fix.degenerate=fix.degenerate,
                                         skip.hidden=skip.hidden)
        )
        if(is.null(ct1)) next
        ct1x <- ct1$exp.matrix
        colnames(ct1x) <- paste0(colnames(ct1x),"@",s)
        ss <- rownames(df1)[sel]
        if(is.null(ct.all)) {
            ct.all <- data.frame(sample=ss, ct1x, check.names=FALSE)
        } else {
            df2 <- data.frame(sample=ss, ct1x, check.names=FALSE)
            ct.all <- dplyr::full_join(ct.all, df2, by="sample")
        }
    }
    
    dim(ct.all)
    if(is.null(ct.all)) {
        message("[pgx.makeAutoContrastsStratified] WARNING : no valid contrasts")
        return(NULL)
    }

    ## sample-wise contrasts
    rownames(ct.all) <- ct.all[,"sample"]
    ct.all <- as.matrix(ct.all[,-1,drop=FALSE])  ## drop sample column
    ct.all <- ct.all[match(rownames(df),rownames(ct.all)),]
    rownames(ct.all) <- rownames(df)
    ct.all[is.na(ct.all)] <- 0
    dim(ct.all)
    ct.all
}


##mingrp=3;slen=20;ref=NULL;fix.degenerate=FALSE;skip.hidden=TRUE
pgx.makeAutoContrasts <- function(df, mingrp=3, slen=20, ref=NULL, 
                                 fix.degenerate=FALSE, skip.hidden=TRUE )
{
    ## "Automagiccally" parse dataframe and create contrasts using
    ## default assumptions.
    ##
    ##
    shortestunique <- function(xx,slen=3) {
        dup <- sapply(1:max(nchar(xx)),
                      function(i) any(duplicated(substring(xx,1,i))))
        if(!any(!dup)) return(xx)
        k <- min(which(!dup))
        substring(xx,1,max(k,slen))
    }
    
    autoContrast1 <- function(x, ref1, slen, mingrp) {
        ## Automatically create contrast. If 2 levels, create A-vs-B,
        ## otherwise create A-vs-others.
        ##
        if(is.null(ref1)) ref1 <- NA
        x <- as.character(x)
        x <- iconv(x, "latin1", "ASCII", sub="")
        nx <- table(x)
        too.small <- names(which(nx < mingrp))
        if(length(too.small)) x[which(x %in% too.small)] <- NA
        nx <- table(x)
        if(length(nx)<2) return(NULL)
        x <- factor(x)
        if(!is.na(ref1)) x <- relevel(x, ref=ref1)
        
        xlevels <- gsub("[^[:alnum:]+-]","",levels(x))        
        levels(x) <- shortestunique(xlevels,slen=slen)
        xref <- gsub("[^[:alnum:]]","",levels(x)[1])
        nn <- length(nx)
        nn

        if(length(levels(x)) != nn) {
            ## something got wrong...
            return(NULL)
        } else if(nn<2) {
            return(NULL)
        } else if(nn == 2) {
            ct <- model.matrix(~x)[,2,drop=FALSE]
            colnames(ct) <- paste0(levels(x)[2],"_vs_",levels(x)[1])
        } else if(nn >= 3) {
            if(is.na(ref1)) {
                ct <- model.matrix(~0 + x)
                colnames(ct) <- paste0(levels(x),"_vs_other")
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
        ## ct[is.na(ct)] <- 0 ## no!!
        ct
    }
    
    ## repeat ref if too short
    if(!is.null(ref) && length(ref)<ncol(df)) ref <- Matrix::head(rep(ref,99),ncol(df))

    ## filter out 'internal/hidden' and 'group' parameters
    ##not.used <- grepl("^[.]|group",colnames(df))
    not.used <- grepl("^[.]",colnames(df))    
    if(skip.hidden && sum(not.used)>0 && sum(!not.used)>0 ) {
        df <- df[,!not.used,drop=FALSE]
    }
    
    ## first all to characters
    df.rownames <- rownames(df)
    df <- data.frame(df, check.names=FALSE)
    df <- apply(df, 2, as.character)
    
    ## trim leading/end parts that are equal
    df <- apply(df, 2, trimsame)
    
    ## try detect (fluffy) comment fields (and remove)
    countSpaces <- function(s) { sapply(gregexpr(" ", s), function(p) { sum(p>=0) } ) }    
    justComment <- function(x) {
        x <- iconv(x, "latin1", "ASCII", sub="")
        (nchar(x) > 50 || countSpaces(x)>=4)
    }
    is.comment <- sapply(df[1,], justComment)
    is.comment
    sel <- which(!is.comment)
    sel
    if(length(sel)==0) {
        cat("WARNING:: could not auto-find variables...\n")
        return(NULL)
    }
    df <- df[,sel,drop=FALSE]
    if(!is.null(ref)) ref <- ref[sel]
    df[df==""] <- NA
    df[df==" "] <- NA
    dim(df)   

    ## ----------- use type.convert to infer parameters
    df <- type.convert(data.frame(df,check.names=FALSE))

    ## ----------- convert numeric variables into bins
    ii <- which(sapply(df,class) %in% c("integer","numeric"))
    if(length(ii)) {
        for(i in ii) {
            x <- as.numeric(as.character(df[,i]))
            x <- c("low","high")[1 + 1*(x > median(x,na.rm=TRUE))]
            df[,i] <- factor(x, levels=c("low","high"))
        }        
    }
    
    ## ----------- try to detect time series (detect factors by time)
    if(0) {
        has.time <- any(grepl("hr|hour|time|day",colnames(df),ignore.case=TRUE))
        has.time
        if(has.time && NCOL(df)>1) {
            time.col <- grep("hr|hour|time|day",colnames(df),ignore.case=TRUE)[1]
            df1 <- df[,-time.col,drop=FALSE]
            time <- df[,time.col]
            vars.bytime <- names(which(apply(df1, 2, function(x) sum(table(x,time)==0)==0)))
            vars.bytime
            if(length(vars.bytime)>0) {
                cat("detected time-series variables:",vars.bytime,"\n")
                dt <- c()
                i=1
                for(i in 1:length(vars.bytime)) {
                    dt <- cbind(dt, paste(df1[,vars.bytime[i]], time, sep="_"))
                }
                colnames(dt) <- paste0(vars.bytime,"_",colnames(df)[time.col])
                Matrix::head(dt)
                df <- cbind(df, dt)
            }
        }
        Matrix::head(df)
    }
    
    ## emergency bail out...
    if(ncol(df)==0) {
        return(NULL)
    }

    ## For each phenotype parameter we 'automagically' try to create a
    ## contrast
    K <- NULL
    i=1
    for(i in 1:ncol(df)) {
        ref1 <- NA
        if(!is.null(ref)) ref1 <- ref[i]
        x <- df[,i]
        too.small <- (x %in% names(which(table(x)<mingrp)))
        x[too.small] <- NA
        x <- iconv(x, "latin1", "ASCII", sub="")  ## ???
        if(!(ref1 %in% x)) ref1 <- NA
        ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$|scrambl|none|dmso|vehicle|low|null|zero|^not"
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
    
    if(is.null(K)) {
        warning("[pgx.makeAutoContrasts] non valid contrasts")
        return(NULL)
    }
    
    rownames(K) <- df.rownames
    Matrix::head(K)
    dim(K)

    ## Now try to infer the underlying "conditions"
    K1 <- contrastAsLabels(K-0.5)
    kcode <- apply(K1,1,paste,collapse="_")
    xc <- factor(kcode, levels=unique(kcode))  ## experimental condition
    if(ncol(K1)>10) levels(xc) <- paste0("condition",1:length(levels(xc))) ## too long...
    jj <- which(!duplicated(kcode))
    length(jj)    
    K2 <- K[jj,colnames(K1),drop=FALSE]
    rownames(K2) <- xc[jj]
    dim(K2)
    Matrix::head(K2)
    is.degenerate = (length(jj) > 0.9*nrow(K1) || mean(table(xc)==1)>0.5 )
    is.degenerate
    
    ## THIS IS EXPERIMENTAL: remove 
    if(fix.degenerate && is.degenerate) {
        cat("WARNING:: contrast matrix looks degenerate. trying to remove some contrasts...\n")
        is.degenerate = TRUE
        iter=0
        while(is.degenerate && iter<100) {
            ## Now try to infer the underlying "conditions"
            kcode <- apply(K1,1,paste,collapse="_")
            xc <- factor(kcode, levels=unique(kcode))  ## experimental condition
            if(ncol(K1)>10) levels(xc) <- paste0("condition",1:length(levels(xc))) ## too long...
            jj <- which(!duplicated(kcode))
            length(jj)
            ## SPECIAL CASE!!! if comparisons are degenerate (no valid
            ## condition groups). LIMMA does not like that. Then delete
            ## phenotype with most levels one by one
            is.degenerate = (length(jj)==nrow(K1) || mean(table(xc)==1)>0.5 )
            is.degenerate        
            if(is.degenerate) {
                ptype <- sub("[:].*","",colnames(K1))
                del.ptype <- names(which.max(table(ptype)))
                del.ptype
                del <- which(ptype %in% del.ptype)
                K1 <- K1[,-del,drop=FALSE]
            }
            iter = iter+1
        }
        iter
        length(jj)    
        K2 <- K[jj,colnames(K1),drop=FALSE]
        rownames(K2) <- xc[jj]
        Matrix::head(K2)
    } else if(!fix.degenerate && is.degenerate) {
        cat("WARNING:: contrast matrix looks degenerate. going for NULL design...\n")
        ## Go for zero design (no-replicates)
        K2 <- NULL
    }
    
    ## Translate coding 0/NA/1 to -1/0/+1 coding of contrast
    K[K==0] <- -1
    K[is.na(K)] <- 0
    if(!is.null(K2)) {
        K2[K2==0] <- -1
        K2[is.na(K2)] <- 0
    }
    rownames(K) <- df.rownames
    
    list(group = xc, contr.matrix = K2, exp.matrix=K)
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

contrastAsLabels <- function(contr.matrix, as.factor=FALSE) {
    contrastAsLabels.col <- function(contr, contr.name) {
        grp1 <- gsub(".*[:]|_vs_.*","",contr.name)
        grp0 <- gsub(".*_vs_|@.*","",contr.name)
        x <- rep(NA,length(contr))
        x[which(contr<0)] <- grp0
        x[which(contr>0)] <- grp1
        if(as.factor) x <- factor(x, levels=c(grp0,grp1))
        x
    }
    K <- data.frame(contr.matrix[,0])
    i=1
    for(i in 1:ncol(contr.matrix)) {
        contr <- contr.matrix[,i]
        contr.name <- colnames(contr.matrix)[i]
        k1 <- contrastAsLabels.col(contr, contr.name)
        K <- cbind(K, k1)
    }
    ##colnames(K) <- sub("[:].*",colnames(contr.matrix))
    colnames(K) <- colnames(contr.matrix)
    rownames(K) <- rownames(contr.matrix)
    K    
}

makeContrastsFromLabelMatrix <- function(lab.matrix) {
    ct.names <- colnames(lab.matrix)
    main.grp <- sapply(strsplit(ct.names, split="_vs_"),"[",1)
    ctrl.grp <- sapply(strsplit(ct.names, split="_vs_"),"[",2)
    main.grp <- sub(".*:","",main.grp)    
    ctrl.grp <- sub("@.*","",ctrl.grp)    
    main.grp
    ctrl.grp
    contr.mat <- matrix(0, nrow(lab.matrix), ncol(lab.matrix))
    rownames(contr.mat) <- rownames(lab.matrix)
    colnames(contr.mat) <- colnames(lab.matrix)
    for(i in 1:ncol(lab.matrix)) {
        lab1 <- trimws(lab.matrix[,i])
        j1 <- which( lab1 == main.grp[i] )
        j0 <- which( lab1 == ctrl.grp[i] )        
        contr.mat[j1,i] <- +1 / length(j1)       
        contr.mat[j0,i] <- -1 / length(j0)        
    }
    contr.mat
}

    
##=====================================================================================
##=========================== END OF FILE =============================================
##=====================================================================================
