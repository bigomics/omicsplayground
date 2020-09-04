##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

imputeMedian <- function(X) {
    mx <- apply(X,1,median,na.rm=TRUE)
    mx[is.na(mx)] <- median(mx,na.rm=TRUE)
    impX <- X
    impX[is.na(impX)] <- 0
    impX <- impX + is.na(X) * mx
    return(impX)
}

##X=ngs$count;y=ngs$samples$group
averageByGroup <- function(X, y) {
    t(apply(X, 1, function(x) tapply(x, y, mean)))
}

gmean <- function(x) {
    ## geometric mean
    exp(mean(log(x + 1e-40)))
}

mat2hugo <- function(x) {
    x <- x[order(-apply(x,1,sd,na.rm=TRUE)),]
    hx <- alias2hugo(rownames(x))
    jj <- which(!duplicated(hx))
    x <- x[jj,]
    rownames(x) <- hx[jj]
    return(x)
}

##s=symbol    
alias2hugo.SAVE <- function(s) {
    require(org.Hs.eg.db,quietly=TRUE)
    ##eg <- sapply(lapply(s, get, env=org.Hs.egALIAS2EG),"[",1)
    s.na = which(!is.na(s) & s!="" & s!=" ")
    s1 <- s[s.na]
    eg <- sapply(mget(s1, env=org.Hs.egALIAS2EG, ifnotfound=NA),"[",1)
    eg[is.na(eg)] <- "unknown"
    symb <- sapply(mget(eg, env=org.Hs.egSYMBOL, ifnotfound=NA),"[",1)
    jj <- which(is.na(symb))
    if(length(jj)) symb[jj] <- s1[jj]
    symb
    symb0 <- rep(NA,length(s))
    symb0[s.na] <- symb
    return(symb0)
}

gx.hist <- function(gx, main="",ylim=NULL) {
    h0 <- hist(as.vector(gx), breaks=120, main=main,
               col="grey",freq=FALSE, ylim=ylim, xlab="signal")
    i = 1
    for(i in 1:ncol(gx)) {
        h1 <- hist(gx[,i], breaks=h0$breaks,plot=FALSE)
        lines( h0$mids, h1$density, col=i+1 )
    }
}

## from http://menugget.blogspot.ch/2011/09/converting-values-to-color-levels.html
#this function converts a vector of values("z") to a vector of color
#levels. One must define the number of colors. The limits of the color
#scale("zlim") or the break points for the color changes("breaks") can 
#also be defined. when breaks and zlim are defined, breaks overrides zlim.
val2col<-function(z, zlim, col = heat.colors(12), breaks){
    if(!missing(breaks)){
        if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
    }
    if(missing(breaks) & !missing(zlim)){
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
    }
    if(missing(breaks) & missing(zlim)){
        zlim <- range(z, na.rm=TRUE)
        zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
        zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
        breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
    }
    CUT <- cut(z, breaks=breaks)
    colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
    return(colorlevels)
}

##HUGO.SYMBOLS <- unique(unlist(as.list(org.Hs.egSYMBOL)))

symbol2hugo <- function(genes, remove.non.hugo=TRUE, silent=FALSE,
                        take.only.first=FALSE, split.char=";", unknown="unknown_gene")
{
    ##remove.non.hugo=TRUE;silent=FALSE;take.only.first=FALSE;split.char=";";unknown="unknown_gene"
    require(org.Hs.eg.db)
    HUGO.SYMBOLS <- unique(unlist(as.list(org.Hs.egSYMBOL)))
    ss <- as.character(genes)
    ss <- gsub("Sep 0","SEPT",ss)  # typical XLS error
    ss[is.na(ss)|ss==""] <- unknown
    ii <- which( !(ss %in% HUGO.SYMBOLS) & ss!=unknown )
    length(ii)
    if(length(ii)==0) {
        return(genes)
    }
    if(!silent) cat("trying to convert",length(ii),"aliases to HUGO\n")
    ss0 <- sapply(ss[ii],strsplit,split=split.char)
    ee0 <- lapply(ss0, function(s) unlist(mget(s,envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
    ee0 <- lapply(ee0, function(e) {e[is.na(e)|e==""|is.nan(e)] <- unknown; e})
    gg  <- lapply(ee0, function(e) unlist(mget(e, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
    if(remove.non.hugo)
        gg <- lapply(gg, intersect, HUGO.SYMBOLS)
    gg.len <- lapply(gg,length)
    sum(gg.len>1)
    if(sum(gg.len>1) && !silent ) {
        cat("warning:",sum(gg.len>1),"entrezID have multiple symbols\n")
    }
    if(sum(gg.len>1) && take.only.first) {
        gg <- sapply(gg,"[",1)
    } else {
        gg <- sapply(gg,paste,collapse=split.char)
    }
    if(!silent)
        cat("updating",length(gg),"deprecated symbols\n")
    gg[is.na(gg)|gg==""] <- unknown
    ss[ii] <- gg
    ss
}


gx.collapse2symbol <- function(X, symbol) {
    j1 <- order(-apply(X,1,sd))
    X <- X[j1,]
    symbol <- symbol[j1]
    j2 <- which(!duplicated(symbol) & !(symbol %in% c(""," ","NA",NA)) )
    X <- X[j2,]
    rownames(X) <- symbol[j2]
    return(X)
}
