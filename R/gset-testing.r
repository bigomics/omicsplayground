
gset.culling <- function( gsets, rho.max=0.50, nmax=400 ) {
    if(nmax<0) nmax=length(gsets)
    if(length(gsets)==1) {
        return(names(gsets))  ## really? if no name?
    }
    if(length(gsets)>nmax) {
        cat("reducing to",nmax,"initial gene sets\n")              
        gsets <- gsets[1:nmax]
    }
    gg <- unique(unlist(gsets))
    G  <- matrix(0, nrow=length(gg), ncol=length(gsets))
    rownames(G) <- gg
    colnames(G) <- names(gsets)
    for(i in 1:ncol(G)) {
        jj <- which( gg %in% gsets[[i]] )
        G[jj,i] <- 1
    }
    cat("culling using rho.max=",rho.max,"\n")      
    G <- scale( G )
    G <- G / sqrt(nrow(G)-1) 
    cf <-  t(G) %*% G
    ##cat("dim(cf)=",dim(cf),"\n")      
    cf[which(cf > rho.max)] <- 1
    cf.wmax <- apply(cf,1,which.max)
    jj <- colnames(G)[which(cf.wmax==1:nrow(cf))]
    length(jj)
    if(length(jj)<=2) jj <- names(gsets)[1:2]
    cat("retaining",length(jj),"gene sets\n")      
    return(jj)
}

## clustering of results of functional analysis
## dup.remove=TRUE;n=60; nb.prop=TRUE; method="ward"
gset.cluster <- function(z, dup.remove=TRUE,n=60, nb.prop=TRUE, method="ward") {
    z <- z[!is.na(z$genes),]
    z <- z[order(z$p.value),]
    if(dup.remove) {
        z <- z[!duplicated(z$genes),]
    }
    z <- z[1:min(nrow(z),n),]
    z.gg <- lapply(z$genes,function(x) strsplit(as.character(x),split="\\|")[[1]])
    n0 <- sort(unique(unlist(z.gg)))
    Z <- matrix(0, nrow=nrow(z), ncol=length(n0))
    colnames(Z) <- n0
    rownames(Z) <- rownames(z)
    for(i in 1:nrow(Z)) {
        Z[i,] <- 1*(colnames(Z) %in% z.gg[[i]])
    }
    if(nb.prop==TRUE) {
        jj <- intersect(colnames(Z),colnames(GNET$shortest.path))
        B <- 1*(GNET$shortest.path[jj,jj]==1)
        ##    B <- 1*(GNET$shortest.path[jj,jj]>0 & GNET$shortest.path[jj,jj]<=1)    
        ##    B <- (1.0/GNET$shortest.path[jj,jj])**4
        B[which(is.na(B))] <- 0
        Zjj <- t(B %*% t(Z[,jj]))
        Zjj <- Zjj / (1e-12 +apply(Zjj,1,max))
        diag(Zjj) <- 0
        Z[,jj] <- Z[,jj] + 0.66*(Zjj>0)
    }
    Z <- Z + matrix(rnorm(ncol(Z)*nrow(Z))*1e-8,nrow=nrow(Z),ncol=ncol(Z))
    ##  Z <- Z / (1e-12 +apply(Z,1,max))
    Z <- pmax(pmin(Z,1),0)
    ##  heatmap(Z, scale="none",mar=c(10,15))
    ##  gx.heatmap2(Z, mar=c(15,20), scale="none", col=grey.colors(64), method="ward" )
    heatmap(Z, mar=c(15,20), scale="none", col=heat.colors(64))
    
    ## NMF
    if(0) {
        ## unstable!!
        require(NMF)
        res <- nmf(t(Z), 3, method="brunet")
        a1 <- t(NMF::basis(res))
        a2 <- NMF::consensus(res)
        a3 <- NMF::coefficients(res)
        h1 <- hclust( dist(t(a1)) )
        h2 <- hclust( dist(t(a2)) )
        h3 <- hclust( dist(t(a3)) )        
        heatmap(a1,mar=c(20,20))
        heatmap(a2,mar=c(20,20))
        heatmap(a3,mar=c(20,20))
        
        ## filtered Z
        Zf <- t(a3) %*% a1
        
        ## plot
        png(w=1600,h=1600,file="gsnmf1.png")
        heatmap(Z, Rowv=as.dendrogram(h3),
                Colv=as.dendrogram(h1),
                scale="none", mar=c(15,20) )
        dev.off()
        png(w=1600,h=1600,file="gsnmf2.png")
        heatmap(Zf, Rowv=as.dendrogram(h3),
                Colv=as.dendrogram(h1),
                scale="none", mar=c(15,20) )
        dev.off()
        png(w=1600,h=1600,file="gsnmf3.png")
        heatmap(Zf, scale="none", mar=c(15,20) )
        dev.off()
    }
    
    return(Z)
}

##========================================================================
##======================== GSA ===========================================
##========================================================================

gx.dfa <- function(gep, pheno, fdr=0.01) { 
    z0 <- gx.limma(gep, pheno, fdr=fdr, compute.means=FALSE)
    gg <- rownames(z0)
    z1 <- gx.geneset.analysis(gg, fdr=fdr )
    z1
}

gx.gsa <- function(gep, pheno, genesets=GENESETS, fdr=0.05, nperm=100,
                   resp.type="Two class unpaired" ) {
    ##fdr=0.05;nperm=100;resp.type="Two class unpaired"
    require(GSA)
    jj <- which(!is.na(pheno))
    gep <- gep[,jj]
    pheno <- pheno[jj]
    if(resp.type=="Two class unpaired") {
        pheno <- as.integer(as.factor(pheno)) ## convert to 1/2
    }
    GSA.obj <- GSA(gep, pheno, genenames=rownames(gep),genesets=genesets,
                   resp.type=resp.type, nperm=nperm)
    gsa <- GSA.listsets( GSA.obj, geneset.names=names(genesets), FDRcut=fdr )
    ss <- rbind(gsa$positive, gsa$negative)
    if(nrow(ss)==0) {
        cat("no significant genesets:: try larger FDR\n")
        return(NULL)
    }
    
    gs.name   <- names(genesets)[ as.integer(ss[,"Gene_set"])]
    gs.len    <- unlist(sapply(genesets[gs.name],length))
    gs.genes  <- genesets[gs.name]
    gs.score  <- as.numeric(ss[,"Score"])
    
    ## quickly limma to determine significant (top)genes
    if(resp.type=="Two class unpaired") {
        mm <- gx.limma( gep, pheno, fdr )
        sum(mm$adj.P.Val<fdr)
        gg.top <- rownames(mm)
    } else {
        pv <- apply(gep[,],1,function(x) cor.test(x,pheno)$p.value)
        qv <- p.adjust(pv, method="fdr")
        gg.top <- names(qv)[which(qv<fdr)]
    }
    gs.overlap <- sapply(gs.genes, function(x) length(intersect(x,rownames(gep))))
    gs.overlap.top <- sapply(gs.genes, function(x) length(intersect(x,gg.top)))  
    gs.genes2 <- sapply(gs.genes,function(x) paste(sort(intersect(gg.top,x)),collapse="|"))
    gs.overlap.top <- unlist(gs.overlap.top)
    G <- data.frame(
        ## gene.set=gs.name,
        p.value=as.numeric(as.character(ss[,"p-value"])),
        q.value=as.numeric(as.character(ss[,"FDR"])),
        score=as.numeric(as.character(ss[,"Score"])),
        set.size=gs.len,
        n.sig=gs.overlap.top,
        sig.ratio=(gs.overlap.top/gs.len),
        sig.genes=gs.genes2 )
    G <- G[gs.overlap>0 & G$n.sig>0 & G$sig.ratio>0.0,]
    G <- G[order(-abs(G$score)),]
    ##  rownames(G) <- NULL
    G
}

##========================================================================
##===================== LASSO-based =======================================
##========================================================================

## abs.rho=FALSE;sort.by="fc";pheno.method="project";correct.sign=TRUE;p.sig=0.01;calc.pv=c("rho","tt","rt","ks","perm.samples","meta");min.results=3 
run.lasso <- function( F, y, trim.coef=0.10, alpha=1.0 )
{
    require(glmnet)
    ## elastic net regression
    cat("running elastic net with alpha=",alpha,"\n")
    fit <- cv.glmnet( F, y, alpha=alpha)
    if(0){
        par(mfrow=c(2,2))
        plot(fit)
        plot(fit$glmnet.fit, "lambda", label=TRUE)
    }
    cx <- coef(fit, s=fit$lambda.min)[,1]
    cx <- cx[which(names(cx)!="(Intercept)" & cx!=0)]
    cx <- cx[order(-abs(cx))]
    cx[1:10]
    if(trim.coef>0) {
        ##    jj <- which( (cumsum(abs(cx)) / sum(abs(cx))) < 0.95 )
        jj <- which( abs(cx) / max(abs(cx)) > trim.coef )
        cat("trimming with c=",trim.coef,"to",length(jj),"coefficients\n")    
        cx <- cx[jj]
    }
    cat("number of nonzero coeff:",sum( abs(cx)>0),"\n")
    cx
}

##alpha=1.0
gset.lasso <- function( genes, genesets, alpha=1.0, trim.coef=0.10,
                       min.genes=8, g.scale=FALSE )
{
    ng <- sapply(genesets, function(x) length(intersect(genes,x)))
    genesets <- genesets[which(ng>=min.genes)]
    length(genesets)
    gg <- unique(c(genes,unlist(genesets)))
    G <- matrix(0, nrow=length(gg), ncol=length(genesets))
    cat("building regression matrix:",nrow(G),"x",ncol(G),"\n")
    rownames(G) <- gg
    for(j in 1:ncol(G)) {
        ii <- which( gg %in% genesets[[j]] )
        G[ii,j] <- 1
    }
    y <- 1*(gg %in% genes)
    names(y) <- gg
    colnames(G) <- names(genesets)
    if(g.scale) G <- as.matrix(scale(G,center=FALSE))  
    z0 <- run.lasso(G, y, trim.coef=trim.coef, alpha=alpha)
    ggx <- lapply( genesets[names(z0)], intersect, genes )
    ggx <- sapply(ggx, paste, collapse="|")
    res <- data.frame( coef=z0, genes=ggx)
    return(res)
}

##trim.coef=1.0;alpha=1.0;limma.fdr=0.20;min.fc=0.20;min.sd=0.20;min.genes=10
gx.lasso <- function(gep, pheno, genesets, trim.coef=0.10, alpha=1.0,
                     limma.fdr=0.20, min.sd=0.20, min.genes=10)
{
    ## run LIMMA to annotate with significant genes
    if(min.sd>0) {
        cat("filtering on min.sd=",min.sd,"\n")
        gep.sd <- apply(gep,1,sd)
        sum(gep.sd>min.sd)
        gep <- gep[which(gep.sd>min.sd),]
    }

    cat("running LIMMA for gene filtering\n")
    mx <- gx.limma(gep, pheno, fdr=limma.fdr)
    sig.genes <- rownames(mx)
    cat("number of significant genes:",length(sig.genes)," (fdr=",limma.fdr,")\n")
    
    ## minmal overlap with sig.genes
    ng <- sapply(genesets, function(x) length(intersect(sig.genes,x)))
    genesets <- genesets[which(ng>=min.genes)]
    length(genesets)

    ## build regression matrix
    gg <- unique(unlist(genesets))
    gg <- intersect(gg, rownames(gep))
    rho <- cor( t(gep[gg,]), pheno )[,1]
    names(rho) <- gg
    G <- matrix(0, nrow=length(gg), ncol=length(genesets))
    rownames(G) <- gg
    colnames(G) <- names(genesets)
    cat("building regression matrix:",nrow(G),"x",ncol(G),"\n")
    for(j in 1:ncol(G)) {
        ii <- which( gg %in% genesets[[j]] )
        G[ii,j] <- 1
    }
    ## G <- as.matrix(scale(G,center=FALSE))
    z0 <- run.lasso(G, rho, trim.coef=trim.coef, alpha=alpha)

    ggx <- lapply( genesets[names(z0)], intersect, sig.genes )
    ggx <- sapply(ggx, paste, collapse="|")
    data.frame( coef=z0, sig.genes=ggx)
}


##========================================================================
##======================= end of file ====================================
##========================================================================

