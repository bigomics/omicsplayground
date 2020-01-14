
##-----------------------------------------------------------------------------
## Contrast creation functions
##-----------------------------------------------------------------------------

##mingrp=3;slen=8;ref=NULL
pgx.makeAutoContrast <- function(df, mingrp=3, slen=8, ref=NULL)
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
        if(is.null(ref1)) ref1 <- NA
        x <- as.character(x)
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
        ct[is.na(ct)] <- 0
        ct
    }

    if(!is.null(ref) && length(ref)!=ncol(df)) ref <- head(rep(ref,99),ncol(df))
    
    K <- c()
    i=1
    for(i in 1:ncol(df)) {
        ref1 <- NA
        if(!is.null(ref)) ref1 <- ref[i]
        x <- df[,i]
        too.small <- (x %in% names(which(table(x)<mingrp)))
        x[too.small] <- NA
        if(!(ref1 %in% x)) ref1 <- NA
        ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$|scrambl|none|vehicle"
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
        warning("[pgx.makeAutoContrast] non valid contrasts")
        return(NULL)
    }
    
    rownames(K) <- rownames(df)
    head(K)

    K1 <- contrastAsLabels(K-0.5)

    ##
    is.degenerate = TRUE
    iter=0
    while(is.degenerate && iter<100) {
        ## Now try to infer the underlying "conditions"
        kcode <- apply(K1,1,paste,collapse="_")
        xc <- factor(kcode, levels=unique(kcode))  ## experimental condition
        ##levels(xc) <- paste0("condition",1:length(levels(xc)))
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
    head(K2)
    
    ## Translate coding 0/NA/1 to -1/0/+1 coding of contrast
    K[K==0] <- -1
    K[is.na(K)] <- 0
    K2[K2==0] <- -1
    K2[is.na(K2)] <- 0
        
    list(group = xc, contr.matrix = K2, exp.matrix=K)
}

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

contrastAsLabels <- function(contr.matrix) {
    contrastAsLabels.col <- function(contr, contr.name) {
        grp1 <- gsub(".*[:]|_vs_.*","",contr.name)
        grp0 <- gsub(".*_vs_","",contr.name)
        x <- rep(NA,length(contr))
        x[which(contr<0)] <- grp0
        x[which(contr>0)] <- grp1
        x
    }
    K <- c()
    i=1
    for(i in 1:ncol(contr.matrix)) {
        contr <- contr.matrix[,i]
        contr.name <- colnames(contr.matrix)[i]
        k1 <- contrastAsLabels.col( contr, contr.name)
        K <- cbind(K, k1)
    }
    ##colnames(K) <- sub("[:].*",colnames(contr.matrix))
    colnames(K) <- colnames(contr.matrix)
    rownames(K) <- rownames(contr.matrix)
    K    
}

##Y=ngs$samples;na.rm=TRUE
makeDirectContrasts2 <- function(Y, ref, na.rm=TRUE)
{
    detectGroups <- function(contr.matrix) {
        group <- apply(contr.matrix,1,paste,collapse="_")
        table(group)
        n.group <- length(unique(group))
        group <- factor(group)
        if(ncol(contr.matrix)>10) {
            levels(group) <- paste0("group",1:n.group)
        }
        group
    }

    contr.matrix <- makeDirectContrasts(Y=Y, ref=ref, na.rm=na.rm, warn=FALSE) 
    contr.matrix <- sign(contr.matrix)    
    contr.matrix0 <- contr.matrix
    no.vs <- grep("_vs_|_VS_",colnames(contr.matrix0),invert=TRUE)
    no.vs
    if(length(no.vs)>0) {
        colnames(contr.matrix0)[no.vs] <- paste0(colnames(contr.matrix0)[no.vs],":Y_vs_N")
    }    
    if(all(grepl("_vs_|_VS_",colnames(contr.matrix0)))) {
        contr.matrix0 <- contrastAsLabels(contr.matrix0)
    }
    group <- detectGroups(contr.matrix0)
    table(group)
    if(length(levels(group)) > 0.5*nrow(contr.matrix)) {
        cat("WARNING:: contrast matrix looks degenerate. consider removing a contrast.\n")
    }

    group.contr.matrix <- contr.matrix[which(!duplicated(group)),]
    rownames(group.contr.matrix) <- group[which(!duplicated(group))]
    
    list(contr.matrix=group.contr.matrix, group=group)
}

makeDirectContrasts <- function(Y, ref, na.rm=TRUE, warn=TRUE) {
    if(warn) warning("makeDirectContrasts is deprectated. please use makeDirectContrasts2()")
    contr.matrix <- c()
    i=2
    all <- c("all","other","others","rest")
    for(i in 1:ncol(Y)) {
        m1 <- NULL
        if(!is.na(ref[i]) && !ref[i] %in% all ) {
            x <- as.character(Y[,i])
            x[is.na(x)|x=="NA"] <- "_"
            m1 <- model.matrix( ~ 0 + x)
            colnames(m1) <- sub("^x","",colnames(m1))
            m1 <- m1 - m1[,ref[i]]  ## +1/-1 encoding
            m1 <- m1[,which(colnames(m1)!=ref[i]),drop=FALSE]  ## remove refvsref...
            m1 <- m1[,!colnames(m1) %in% c("NA","_"),drop=FALSE]
            colnames(m1) <- paste0(colnames(m1),"_vs_",ref[i])
        } else if(!is.na(ref[i]) && ref[i] %in% all ) {
            x <- as.character(Y[,i])
            x[is.na(x)|x=="NA"] <- "_"            
            m1 <- model.matrix( ~ 0 + x)
            colnames(m1) <- sub("^x","",colnames(m1))
            ##m1 <- m1 - m1[,ref[i]]  ## +1/-1 encoding
            m1 <- t(t(m1==1) / colSums(m1==1) - t(m1==0) / colSums(m1==0))
            ##m1 <- m1[,which(colnames(m1)!=ref[i]),drop=FALSE]
            m1 <- m1[,!colnames(m1) %in% c("NA","_"),drop=FALSE]            
            colnames(m1) <- paste0(colnames(m1),"_vs_others")
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
    sign(contr.matrix)
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
    colnames(m1) <- paste0(colnames(m1),"_vs_others")
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



##=====================================================================================
##=========================== END OF FILE =============================================
##=====================================================================================
