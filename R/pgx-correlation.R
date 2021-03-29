##df=ngs$samples
pgx.testTraitRelationship <- function(me, df, plot=TRUE, cex=1)
{

    require(corrplot)
    require(WGCNA)

    df <- type.convert(df)
    cl <- sapply(df,class)
    cvar <- which(cl %in% c("numeric","integer"))
    dvar <- which(cl %in% c("factor","character"))    
    dc <- df[,cvar,drop=FALSE]
    dd <- df[,dvar,drop=FALSE]

    dim(dc)
    dim(dd)
    
    ## contious vs continous -> correlation
    rho.P <- NULL
    if(ncol(dc)) {
        rho.P <- matrix(NA,ncol(me),ncol(dc))
        i=1;j=2
        rho <- cor(me, dc, use="pairwise")
        rho.P <- cor.pvalue(P, nrow(me))
        dim(rho.P)
    }

    ## continous vs discrete -> ANOVA
    anova.P <- NULL
    if(ncol(dd)) {
        anova.P <- matrix(NA,ncol(me),ncol(dd))
        colnames(anova.P) <- colnames(dd)
        rownames(anova.P) <- colnames(me)
        for(i in 1:ncol(dd)) {
            y <- dd[,i]
            res <- gx.limmaF( t(me), y, fdr=1, lfc=0)
            anova.P[,i] <- res[colnames(me),"P.Value"]
        }
        anova.P
    }
    
    P <- cbind(rho.P, anova.P)
    P <- P[,colnames(df)]
        
    if(plot==TRUE) {

        require(corrplot)    
        sigP <- -log10(P+1e-8)
        sigP <- (1 - P)**1
        ##par(oma=c(0,0,0,1))
        corrplot( sigP, is.corr=FALSE, ##type="upper",
                 mar = c(0,0,0,2),
                 ##p.mat = Q, sig.level = 0.05, ##insig = "blank",
                 tl.cex = cex, tl.col="black", tl.offset = 1,
                 cl.align.text = "l", cl.offset = 0.25, cl.cex = 0.7, 
                 pch.col = "grey50")


    }
    return(P)
}


##df=ngs$samples
pgx.testPhenoCorrelation <- function(df, plot=TRUE, cex=1)
{

    require(corrplot)
    
    cl <- sapply(df,class)
    nlev = apply(df,2,function(x) length(unique(x[!is.na(x)])))
    cvar <- which(cl %in% c("numeric","integer")  & nlev>=2 )
    dvar <- which(cl %in% c("factor","character") & nlev>=2 )    
    dc <- df[,cvar,drop=FALSE]
    dd <- df[,dvar,drop=FALSE]

    dim(dc)
    dim(dd)
    
    ## discrete vs discreate -> Fisher test
    fisher.P <- NULL
    if(ncol(dd)) {
        fisher.P <- matrix(NA,ncol(dd),ncol(dd))
        i=1;j=2
        for(i in 1:(ncol(dd)-1)) {
            kk <- which( !is.na(dd[,i]) & !is.na(dd[,j]) )
            if(length(unique(dd[kk,i])) < 2 || length(unique(dd[kk,j])) < 2) next            
            for(j in (i+1):ncol(dd)) {
                tb <- table(dd[,i], dd[,j])
                fisher.P[i,j] <- fisher.test(tb, simulate.p.value=TRUE)$p.value
            }
        }
        rownames(fisher.P) <- colnames(dd)
        colnames(fisher.P) <- colnames(dd)
    }
        
    ## discrete vs continuous -> ANOVA or Kruskal-Wallace
    kruskal.P <- NULL
    if(ncol(dc)>0) {
        kruskal.P <- matrix(NA,ncol(dd),ncol(dc))
        i=1;j=2
        for(i in 1:ncol(dd)) {
            kk <- which( !is.na(dc[,j]) & !is.na(dd[,i]) )
            if(length(unique(dd[kk,i])) < 2) next
            for(j in 1:ncol(dc)) {
                kruskal.P[i,j] <- kruskal.test(dc[kk,j], dd[kk,i])$p.value
            }
        }
        rownames(kruskal.P) <- colnames(dd)
        colnames(kruskal.P) <- colnames(dc)
    }
    
    ## continuous vs continuous -> correlation test
    cor.P <- NULL
    if(ncol(dc)>1) {
        cor.P <- matrix(NA,ncol(dc),ncol(dc))
        i=1;j=2
        for(i in 1:(ncol(dc)-1)) {
            for(j in (i+1):ncol(dc)) {
                cor.P[i,j] <- cor.test(dc[,i], dc[,j])$p.value
            }
        }
        rownames(cor.P) <- colnames(dc)
        colnames(cor.P) <- colnames(dc)
    }
    
    P <- matrix(NA,ncol(df),ncol(df))
    rownames(P) <- colnames(P) <- colnames(df)

    if(!is.null(fisher.P)) {
        ii <- match(rownames(fisher.P),rownames(P))
        jj <- match(colnames(fisher.P),colnames(P))
        P[ii,jj] <- fisher.P
    }

    if(!is.null(kruskal.P)) {
        ii <- match(rownames(kruskal.P),rownames(P))
        jj <- match(colnames(kruskal.P),colnames(P))
        P[ii,jj] <- kruskal.P
    }

    if(!is.null(cor.P)) {
        ii <- match(rownames(cor.P),rownames(P))
        jj <- match(colnames(cor.P),colnames(P))
        P[ii,jj] <- cor.P
    }
    
    ij <- which(!is.na(P),arr.ind=TRUE)
    qv <- p.adjust(P[ij], method="BH")
    Q <- P
    Q[ij] <- qv
    
    P[is.na(P)] <- 0
    P <- (P + t(P))/2
    Q[is.na(Q)] <- 0
    Q <- (Q + t(Q))/2
    
    if(plot==TRUE) {

        require(corrplot)    
        logP <- -log10(P+1e-8)
        logQ <- -log10(Q+1e-8)
        diag(logQ) <- 0
        ##par(oma=c(0,0,0,1))
        corrplot( logQ, is.corr=FALSE, type="upper",
                 mar = c(0,0,0,2),
                 p.mat = Q, sig.level = 0.05, ##insig = "blank",
                 tl.cex = cex, tl.col="black", tl.offset = 1,
                 cl.align.text = "l", cl.offset = 0.25, cl.cex = 0.7, 
                 pch.col = "grey50",
                 order="hclust")

    }

    return(list(P=P, Q=Q))
}
