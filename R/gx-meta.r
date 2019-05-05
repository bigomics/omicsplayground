
GX.TEST_METHODS=c("ttest","rank.ttest","limma","limma.trend","pearson","rank.pearson")
methods=GX.TEST_METHODS

gx.runMethods.TOBEFINISHED <- function(xx, yy, ref, fdr=0.2, lfc=0.5,
                                       methods=GX.TEST_METHODS)
{
    cat("gmt.size=",length(gmt),"\n")
    yy
    if(is.null(ref))
        ref=levels(yy)[1]

    fy = as.integer(relevel(factor(yy),ref=ref)) -1
    fy

    ## calculate significant genes with LIMMA (we need all gene for GSEA-PR)
    ##lfc=0;ref=NULL
    all.tests = list()
    ## Standard Fisher exact test
    if("limma" %in% methods) {
        limma = gx.limma( xx, yy, fdr=1.0, lfc=lfc, ref=ref, trend=FALSE )  ## trend true for NGS
        limma = limma[rownames(xx),]
        all.tests[["limma"]] <- limma
    }
    if("limma.trend" %in% methods) {
        limma = gx.limma( xx, yy, fdr=1.0, lfc=lfc, ref=ref, trend=TRUE )  ## trend true for NGS
        limma = limma[rownames(xx),]
        all.tests[["limma.trend"]] <- limma
    }
    if("ttest" %in% methods) {
        require(matrixTests)
        suppressWarnings( resTT <- matrixTests::row_t_welch( xx[,which(fy==1)], xx[,which(fy==0)]) )
        resTT = resTT[rownames(xx),]
        resTT$qvalue = p.adjust( resTT$pvalue, method="fdr")
        all.tests[["ttest"]] <- resTT
    }
    if("rank.ttest" %in% methods) {
        require(matrixTests)
        xx1 = apply(xx, 2, rank)
        ##xx = scale(xx)
        suppressWarnings( resTT <- matrixTests::row_t_welch( xx1[,which(fy==1)], xx1[,which(fy==0)]) )
        resTT = resTT[rownames(xx),]
        resTT$qvalue = p.adjust( resTT$pvalue, method="fdr")
        all.tests[["rank.ttest"]] <- resTT
    }
    if("pearson" %in% methods) {
        require(matrixTests)
        suppressWarnings( resTT <- matrixTests::row_cor_pearson( xx[,which(fy==1)], xx[,which(fy==0)]) )
        resTT = resTT[rownames(xx),]
        resTT$qvalue = p.adjust( resTT$pvalue, method="fdr")
        all.tests[["pearson"]] <- resTT
    }
    if("rank.pearson" %in% methods) {
        require(matrixTests)
        xx = apply(xx, 2, rank)
        suppressWarnings( resTT <- matrixTests::row_cor_pearson( xx[,which(fy==1)], xx[,which(fy==0)]) )
        resTT = resTT[rownames(xx),]
        resTT$qvalue = p.adjust( resTT$pvalue, method="fdr")
        all.tests[["rank.pearson"]] <- resTT
    }
    names(all.tests)
    return(all.tests)
}

gx.metaTestMethods.TOBEFINISHED <- function(all.tests)
{

    ##--------------------------------------------------------------
    ## count significant terms
    ##--------------------------------------------------------------
    ##fdr = 0.25
    i=1
    sig.up = list()
    sig.dn = list()
    N1 = c()
    for(i in 1:length(all.tests)) {
        R = all.tests[[i]]
        rr = tolower(colnames(R))
        pv = R[,grep("nom p-val|^p.value|^p.value|^np$|pval",rr)[1]]
        qv = R[,grep("fdr|adj.p.val|q.value|adjusted.p|qval",rr)[1]]
        fx = R[,grep("^nes$|logfc|sign|statistic",rr)]
        sig.up1 = rownames(R)[which(pv <= 0.05 & qv <= fdr & fx >= 0)]
        sig.dn1 = rownames(R)[which(pv <= 0.05 & qv <= fdr & fx < 0)]
        n.up = length(sig.up1)
        n.down = length(sig.dn1)
        n.notsig = sum(qv > fdr, na.rm=TRUE)
        N1 = rbind(N1, c( n.down, n.notsig, n.up))
        sig.up[[i]] = sig.up1
        sig.dn[[i]] = sig.dn1
    }
    colnames(N1) = c("n.down","n.zero","n.up")
    rownames(N1) = names(all.tests)
    names(sig.up) = names(all.tests)
    names(sig.dn) = names(all.tests)
    N1

    ##--------------------------------------------------
    ## Add overlap counts
    ##--------------------------------------------------
    nmax = 4
    nmax = length(all.tests)
    nmax
    common.up = names(which(table(unlist(sig.up))==nmax))
    common.dn = names(which(table(unlist(sig.dn))==nmax))
    common.sig <- unique(c(common.dn,common.up))
    length(common.sig)
    n1 = length(common.dn)
    n2 = length(common.up)
    n0 = NA
    N1 <- rbind(N1, "**common**" = c(n1, n0, n2))
    N1
    ##N1 = cbind(comparison=rownames(N1), N1)
    ##rownames(N1) = NULL
    overlap.up = sapply(sig.up, function(x)
        sapply(sig.up, function(y) length(intersect(x,y))))
    overlap.dn = sapply(sig.dn, function(x)
        sapply(sig.dn, function(y) length(intersect(x,y))))
    overlap.up
    overlap.dn

    ##--------------------------------------------------
    ## meta analysis, aggregate p-values
    ##--------------------------------------------------
    gsets = sort(unique(unlist(lapply(all.tests,rownames))))
    length(gsets)

    i=1
    pv = c()
    qv = c()
    fx = c()
    j=1
    for(j in 1:length(all.tests)) {
        x = all.tests[[j]]
        pv.col = grep("p.value|^p$|p-val|pval|nom.p.val|nom p-val|p.val",tolower(colnames(x)))[1]
        qv.col = grep("q.value|^q$|q-val|pval|fdr.q.val|fdr q-val|q.val|adj.p",tolower(colnames(x)))[1]
        fx.col = grep("sign|nes|logfc|fc",tolower(colnames(x)))[1]
        pv = cbind(pv, x[match(gsets,rownames(x)),pv.col] )
        qv = cbind(qv, x[match(gsets,rownames(x)),qv.col] )
        fx = cbind(fx, (x[match(gsets,rownames(x)),fx.col] ))
        colnames(pv)[ncol(pv)] = names(all.tests)[j]
        colnames(qv)[ncol(qv)] = names(all.tests)[j]
        colnames(fx)[ncol(fx)] = names(all.tests)[j]
    }
    head(pv)
    rownames(pv) = gsets
    rownames(qv) = gsets
    rownames(fx) = gsets

    pv[is.na(pv)] = 0.9999
    fx[is.na(fx)] = 0
    require(metap)
    ##p.meta   = apply(pv, 1, function(p) metap::allmetap(p, method="sumlog")$p[[1]])
    p.meta = apply(pv, 1, function(p) metap::allmetap(p, method="sumz")$p[[1]])  ## stouffer
    ##p.meta = apply(pv, 1, function(p) metap::allmetap(p, method="maximump")$p[[1]])
    q.meta = p.adjust(p.meta, method="fdr")
    meta.fx = rowMeans(apply(fx, 2, scale, center=FALSE, na.rm=TRUE))
    meta = data.frame(fx=meta.fx, p=p.meta, q=q.meta)
    ##res.meta = data.frame(meta=I(meta), fx=I(fx), p=I(pv), q=I(qv))
    res.meta = data.frame(meta=meta, fx=fx, p=pv, q=qv)
    rownames(res.meta) = gsets

    res = list( outputs=all.tests, meta=res.meta,
               sig.counts=N1, sig.up=sig.up, sig.dn=sig.dn,
               common.up=common.up, common.dn=common.dn,
               overlap.up=overlap.up, overlap.dn=overlap.dn)

    return(res)
}


##======================================================================
##======================================================================
##======================================================================
