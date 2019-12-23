##=============================================================================
##==========    Platform global and user settings =============================
##=============================================================================

USER.GENESETTEST.METHODS <- NULL
USER.GENETEST.METHODS <- NULL

##=============================================================================
##==========    Platform helper functions =====================================
##=============================================================================
##s=rep("abc",100)

##check.names=FALSE;row.names=1;stringsAsFactors=FALSE
fread.csv <- function(file, check.names=FALSE, row.names=1,
                      stringsAsFactors=FALSE)
{
    require(data.table)
    df <- fread(file=file, check.names=check.names)
    x <- data.frame(df[,2:ncol(df)], stringsAsFactors=stringsAsFactors,
                    check.names=check.names)
    rownames(x) <- df[[row.names]]
    if(all(sapply(x,class)=="numeric")) x <- as.matrix(x)
    if(all(sapply(x,class)=="character")) x <- as.matrix(x)
    if(all(sapply(x,class)=="integer")) x <- as.matrix(x)
    return(x)
}

tagDuplicates <- function(s) {
    ## Tag duplicate with blanks
    ##
    jj <- which(duplicated(s))
    t <- s[jj][1]
    for(t in unique(s[jj])) {
        ii <- which(s==t)
        spaces <- paste(rep(" ",length(ii)),collapse="")
        blanks <- substring(spaces,0,0:(length(ii)-1))
        ##s[ii] <- paste(s[ii],1:length(ii),sep=".")
        s[ii] <- paste0(s[ii],blanks)
    }
    s <- gsub("[.]1$","",s)
    s
}

wrapHyperLink <- function(s, gs) {
    
    ## GEO/GSE accession
    gs = as.character(gs)
    s1 = s = as.character(s)
    jj <- grep("GSE[0-9]",gs)
    if(length(jj)) {
        acc = sub("[-_ ].*","",gsub("^.*GSE","GSE",gs[jj]))
        url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",acc)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }

    ## KEGG
    jj <- grep("HSA-[0-9][0-9]|hsa[0-9][0-9]",gs)
    if(length(jj)) {
        id = sub("^.*HSA-|.*hsa","hsa",gs[jj])
        id = sub("[-_ ].*","",id)
        url = paste0("https://www.genome.jp/kegg-bin/show_pathway?map=",id,"&show_description=show")
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }

    ## Wikipathways
    jj <- grep("_WP[0-9][0-9]",gs)
    if(length(jj)) {
        id = sub("^.*_WP","WP",gs[jj])
        url = paste0("https://www.wikipathways.org/index.php/Pathway:",id)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }

    ## MSigDB
    jj <- grep("^H:|^C[1-8]:|HALLMARK",gs)
    if(length(jj)) {
        gs1 = sub("^.*:","",gs[jj])
        url = paste0("http://software.broadinstitute.org/gsea/msigdb/cards/",gs1)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }
    
    ## GO reference
    jj <- grep("\\(GO_.*\\)$",gs)
    if(length(jj)) {
        id = gsub("^.*\\(GO_|\\)$","",gs[jj])
        url = paste0("http://amigo.geneontology.org/amigo/term/GO:",id)
        s1[jj] <- paste0("<a href='",url,"' target='_blank'>",s[jj],"</a>")
    }
    
    return(s1)
}

reverse.AvsB.1 <- function(comp) {
    paste(rev(strsplit(comp, split="_vs_|_VS_")[[1]]),collapse="_vs_")
}
reverse.AvsB <- function(comp) {
    sapply(comp,reverse.AvsB.1)
}

is.POSvsNEG <- function(ngs) {
    ## Determines from contrast matrix if notation is 'A_vs_B' or
    ## 'B_vs_A'.
    ##
    cntrmat <- ngs$model.parameters$contr.matrix
    ##ct0 <- cntrmat[,comp]        
    grp1 <- sapply(strsplit(colnames(cntrmat),split="_vs_"),"[",1)
    grp2 <- sapply(strsplit(colnames(cntrmat),split="_vs_"),"[",2)
    grp1x <- intersect(grp1,rownames(cntrmat))
    grp2x <- intersect(grp2,rownames(cntrmat))        
    grp1.sign <- mean(cntrmat[intersect(grp1,rownames(cntrmat)),which(grp1 %in% grp1x)])
    grp2.sign <- mean(cntrmat[intersect(grp2,rownames(cntrmat)),which(grp2 %in% grp2x)])
    grp1.sign
    grp2.sign
    is.PosvsNeg1 <- (grp1.sign > 0 && grp2.sign < 0)
    ##is.NegvsPos1 <- (grp2.sign > 0 && grp1.sign < 0)
    is.PosvsNeg1
    
    grp1.neg2 <- mean(grepl("neg|untr|ref|wt|ctr|control",tolower(grp1)))
    grp2.neg2 <- mean(grepl("neg|untr|ref|wt|ctr|control",tolower(grp2)))
    grp1.neg2
    grp2.neg2
    is.PosvsNeg2 <- ( grp2.neg2 > grp1.neg2)
    is.PosvsNeg2
    
    ok <- setdiff(c(is.PosvsNeg1,is.PosvsNeg2),NA)[1]
    if(is.na(ok) || length(ok)==0) ok <- TRUE  ## DEFAULT if not known !!!!!
    ok
}

is.categorical <- function(x, max.ncat=20, min.ncat=2) {
    is.factor <- any(class(x) %in% c("factor","character"))
    is.factor
    n.unique <- length(unique(setdiff(x,NA)))
    n.notna  <- length(x[!is.na(x)])
    is.id    <- (n.unique > 0.8*n.notna)
    is.id
    is.factor2 <- (is.factor & !is.id & n.unique>=min.ncat & n.unique<= max.ncat)
    is.factor2
}

pgx.getCategoricalPhenotypes <-function(df, max.ncat=20, min.ncat=2) {
    ##
    ##
    ##

    is.bad = 0
    is.bad1 <- grepl("^sample$|id$|replicate|patient|donor|individ|ratio|year|month|day|age|efs|dfs|surv|follow",tolower(colnames(df)))
    ## is.factor <- sapply(sapply(data.frame(df), class), function(s) any(s %in% c("factor","character")))
    is.bad2 <- apply(df,2,function(x) any(grepl("^sample|patient|replicate|donor|individ",x,ignore.case=TRUE)))    
    is.bad <- (is.bad1 | is.bad2)
    
    is.factor <- apply(df, 2, is.categorical)
    is.factor
    n.unique <- apply(df,2,function(x) length(unique(setdiff(x,c(NA,"NA","")))))
    n.notna  <- apply(df,2,function(x) length(x[!is.na(x)]))
    is.id    <- (n.unique > 0.8*n.notna)
    is.id
    is.factor2 <- (!is.bad & is.factor & !is.id & n.unique>=min.ncat & n.unique<= max.ncat)
    is.factor2
    colnames(df)[which(is.factor2)]
}

pgx.getMetaFoldChangeMatrix <- function(ngs, what="meta")
{
    fc0 = NULL
    qv0 = NULL
    ##ngs <- inputData()
    sel = names(ngs$gset.meta$meta)
    methods = colnames(unclass(ngs$gx.meta$meta[[1]]$fc))
    if(what %in% methods) {
        fc0 = sapply(ngs$gx.meta$meta[sel], function(x) unclass(x$fc)[,what])
        qv0 = sapply(ngs$gx.meta$meta[sel], function(x) unclass(x$q)[,what])
        rownames(fc0)=rownames(qv0)=rownames(ngs$gx.meta$meta[[1]])
    } else if(what=="meta") {
        fc0 = sapply(ngs$gx.meta$meta[sel], function(x) x$meta.fx)
        qv0 = sapply(ngs$gx.meta$meta[sel], function(x) x$meta.q)
        rownames(fc0)=rownames(qv0)=rownames(ngs$gx.meta$meta[[1]])
    } else {
        cat("WARNING:: pgx.getMetaFoldChangeMatrix: unknown method")
        return(NULL)
    }
    res = list(fc=fc0, qv=qv0)
    return(res)
}

##gene="CD4"
pgx.getGeneCorrelation <- function(gene, xref) {

    rho.genes <- unlist(lapply(xref, rownames))
    rho.genes <- sort(unique(rho.genes))
    length(rho.genes)

    R <- NULL
    ## correlation using external datasets
    k=1
    for(k in 1:length(xref)) {
        has.gene <- (toupper(gene) %in% toupper(rownames(xref[[k]])))
        if( has.gene ) {
            ## take geometric mean with correlation in TISSUE
            xx <- log2(1 + xref[[k]])
            jj <- match(toupper(gene), toupper(rownames(xx)))
            tx = xx[jj,]
            if(class(xx)=="dgCMatrix") {
                suppressWarnings( rho1 <- corSparse(xx, cbind(tx))[,1] )
            } else {
                suppressWarnings( rho1 <- cor(t(xx), tx)[,1])
            }
            if(sum(is.na(rho1))>0) rho1[which(is.na(rho1))] = 0

            rho1 = rho1[match(rho.genes,names(rho1))]
            R <- cbind(R, rho1)
            if(NCOL(R)==1) R <- matrix(R, ncol=1)
            rownames(R) <- rho.genes
            colnames(R)[ncol(R)] <- names(xref)[k]
        }
    }

    ##if(is.null(R)) return(NULL)
    dim(R)
    if(!is.null(R) && NCOL(R)>0) {
        R[is.na(R)] <- 0
        R[is.nan(R)] <- 0
        if(NCOL(R)==1) R <- matrix(R, ncol=1)
        rownames(R) <- rho.genes
        R <- R[which(rowSums(R!=0,na.rm=TRUE)>0),,drop=FALSE]
        ## geneset rho has no sign (=cosine correlation) so we use the sign of others
        if(ncol(R)>1 && "gene sets" %in% colnames(R)) {
            k <- which(colnames(R)=="gene sets")
            R[,"gene sets"] <- R[,"gene sets"] * sign(rowMeans(R[,-k,drop=FALSE],na.rm=TRUE))
        }
        head(R)
    }
    return(R)
}

ngs.save <- function(ngs, file, update.date=TRUE, light=TRUE, system=FALSE) {

    if(update.date||is.null(ngs$date)) ngs$date <- Sys.Date()

    if(light) {
        ## ------- make a light version
        ngs$gx.meta$outputs <- NULL
        ngs$gset.meta$outputs <- NULL
        ngs$model.parameters$efit <- NULL
        ngs$gmt.all <- NULL
        ngs$families <- NULL
        ngs$collections <- NULL
        ## ngs$counts <- NULL
        ngs$gset.meta$matrices <- NULL
    }
    if(system==FALSE) {
        ## remove system (is big...)
        ngs$omicsnet <- NULL
        ngs$omicsnet.reduced <- NULL
    }
    sort(sapply(ngs, object.size)) / 1e9
    sum(sapply(ngs, object.size)) / 1e9        
    
    cat(">>> saving PGX file to",file,"\n")
    save(ngs, file=file)
}

##comp=1;level="geneset";probe=rownames(ngs$gsetX)[1]

getLevels <- function(Y) {
    yy = Y[,grep("group|title|name|sample",colnames(Y),invert=TRUE),drop=FALSE]
    yy = yy[,which(apply(yy,2,function(x) length(unique(x)))<20),drop=FALSE]
    ##levels = setdiff(unique(as.vector(apply(yy,2,as.character))),c("",NA))
    levels = lapply(1:ncol(yy), function(i) unique(paste0(colnames(yy)[i],"=",yy[,i])))
    levels = sort(unlist(levels))
    return(levels)
}



selectSamplesFromSelectedLevels <- function(Y, levels)
{
    if( is.null(levels) || levels=="") {
        return(rownames(Y))
    }
    pheno <- sapply(strsplit(levels,split="="),"[[",1)
    ptype <- sapply(strsplit(levels,split="="),"[[",2)
    sel = rep(1, nrow(Y))
    i=1
    for(ph in unique(pheno)) {
        k <- which(pheno==ph)
        sel <- sel * (Y[,ph] %in% ptype[k])
    }
    rownames(Y)[which(sel==1)]
}

##fields=c("symbol","name","alias","map_location","summary")
getMyGeneInfo <- function(eg, fields=c("symbol","name","alias","map_location","summary"))
{
    require(mygene)
    ##res = getGene(eg, fields="all")[[1]]
    ##fields = c("symbol","name","alias","map_location","summary")
    info = lapply(fields, function(f) getGene(eg, fields=f)[[1]] )
    names(info) <- fields
    info = lapply(info, function(x) ifelse(length(x)==3,x[[3]],"(not available)") )
    info = sapply(info, paste, collapse=",")
    if(0) {
        rifs = getGene(eg, fields="generif")[[1]]
        collapse.rif <- function(r) paste0(r$text," (PMID=",r$pubmed,")")
        rifs = rifs[which(sapply(rifs,length)>2)]
        xrif = sapply(rifs[[3]], function(rr) collapse.rif(rr))
    }
    return(info)
}

## much faster and off-line
getHSGeneInfo <- function(eg, as.link=TRUE) {
    require(org.Hs.eg.db)
    require(KEGG.db)
    require(GO.db)

    env.list <- c("symbol"=org.Hs.egSYMBOL,
                  "name"=org.Hs.egGENENAME,
                  "map_location"=org.Hs.egMAP,
                  "OMIM" = org.Hs.egOMIM,
                  "KEGG"=org.Hs.egPATH,
                  ##"PMID"=org.Hs.egPMID,
                  "GO"=org.Hs.egGO)

    info <- lapply(env.list, function(env) mget(eg, env=env, ifnotfound=NA)[[1]])
    names(info) <- names(env.list)
    gene.symbol <- toupper(mget(as.character(eg), env=org.Hs.egSYMBOL))[1]
    info[["symbol"]] <- gene.symbol
    
    ## create link to GeneCards
    if(as.link) {
        genecards.link = "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=GENE' target='_blank'>GENE</a>"
        info[["symbol"]] <- gsub("GENE", info[["symbol"]], genecards.link)
    }

    ## create link to OMIM
    if(as.link) {
        omim.link = "<a href='https://www.omim.org/entry/OMIM' target='_blank'>OMIM</a>"
        info[["OMIM"]] <- sapply(info[["OMIM"]], function(x) gsub("OMIM",x,omim.link))
    }

    ## create link to KEGG
    kegg.link = "<a href='https://www.genome.jp/kegg-bin/show_pathway?map=hsaKEGGID&show_description=show' target='_blank'>KEGGNAME (KEGGID)</a>"
    for(i in 1:length(info[["KEGG"]])) {
        kegg.id = info[["KEGG"]][[i]]
        kegg.id = setdiff(kegg.id,NA)
        if(length(kegg.id)>0) {
            kegg.name = mget(kegg.id, env=KEGGPATHID2NAME, ifnotfound=NA)[[1]]
            if(!is.na(kegg.name) && as.link) {
                info[["KEGG"]][[i]] <- gsub("KEGGNAME",kegg.name,gsub("KEGGID",kegg.id,kegg.link))
            } else {
                info[["KEGG"]][[i]] <- kegg.name
            }
        }
    }

    ## create link to GO
    if(!is.na(info[["GO"]][1])) {
        go.evidence = c("EXP","IDA","IPI","IMP","IGI","IEP")
        amigo.link = "<a href='http://amigo.geneontology.org/amigo/term/GOID' target='_blank'>GOTERM (GOID)</a>"
        sel <- which( sapply(info[["GO"]],"[[",2) %in% go.evidence  &
                      sapply(info[["GO"]],"[[",3) %in% c("BP")  )
        sel
        info[["GO"]] <- info[["GO"]][sel]

        ## sometimes GO.db is broken...
        suppressWarnings( try.out <- try(Term(mget("GO:0000001", env=GOTERM, ifnotfound=NA)[[1]])))
        go.ok <- (class(try.out) !="try-error")
        if(go.ok && length(sel)>0) {
            i=1
            for(i in 1:length(info[["GO"]])) {
                go_id = info[["GO"]][[i]][[1]]
                go_term = Term(mget(go_id, env=GOTERM, ifnotfound=NA)[[1]])
                if(as.link) {
                    info[["GO"]][[i]] = gsub("GOTERM",go_term,gsub("GOID",go_id,amigo.link))
                } else {
                    info[["GO"]][[i]] = go_term
                }
            }
        } else {
            info[["GO"]] <- NULL
        }
    }

    return(info)
}

##skipifexists=0;perplexity=30;sv.rank=-1;prefix="C";kclust=1;ntop=1000;fromX=FALSE;prior.counts=1;mean.center=TRUE;method="tsne";determine.clusters=1;dims=2;find.clusters=TRUE
pgx.clusterSamples <- function(ngs, skipifexists=FALSE, perplexity=NULL,
                               ntop=1000, sv.rank=-1, prefix="C", 
                               fromX=FALSE, is.logx=FALSE,
                               kclust=1, prior.counts=NULL, 
                               dims=c(2,3), find.clusters=TRUE,
                               row.center=TRUE, row.scale=FALSE,
                               method=c("tsne","umap","pca") )
{


    sX <- ngs$counts
    if(fromX) sX <- 2**ngs$X

    res <- pgx.clusterSamplesFromMatrix(
        sX, perplexity=perplexity, is.logx=FALSE,
        ntop=ntop, sv.rank=sv.rank, prefix=prefix,         
        kclust=kclust, prior.counts=prior.counts, 
        dims=dims, find.clusters=find.clusters,
        row.center=row.center, row.scale=row.scale,
        method=method)

    if(!is.null(res$pos2d)) ngs$tsne2d <- res$pos2d
    if(!is.null(res$pos3d)) ngs$tsne3d <- res$pos3d
    if(!is.null(res$idx)) ngs$samples$cluster <- res$idx
    
    return(ngs)
}

pgx.clusterSamplesFromMatrix <- function(counts, perplexity=NULL,
                                         ntop=1000, sv.rank=-1, prefix="C", 
                                         fromX=FALSE, is.logx=FALSE,
                                         kclust=1, prior.counts=NULL, 
                                         dims=c(2,3), find.clusters=TRUE,
                                         row.center=TRUE, row.scale=FALSE,
                                         method=c("tsne","umap","pca") )
{
    require(Rtsne)
    require(irlba)
    ##set.seed(0)
    method <- method[1]
    
    sX <- counts
    if(is.logx) sX <- 2**sX
    
    if(is.null(prior.counts)) {
        qq <- quantile(as.vector(sX[sX>0]),probs=0.50) ## at 50%
        qq
        prior.counts <- qq
    }
    sX <- log2(prior.counts + sX)
    sX <- limma::normalizeQuantiles(sX)  ## in linear space
    if(row.center) sX <- sX - rowMeans(sX,na.rm=TRUE)
    sX = head( sX[order(-apply(sX,1,sd)),], ntop)
    ## sX = t(scale(t(sX),scale=TRUE))  ## really? or just centering?
    ##sX = t(scale(t(sX),scale=FALSE))  ## really? or just centering?
    if(row.scale) sX <- (sX / apply(sX,1,sd,na.rm=TRUE))
    
    dim(sX)
    ## some randomization is sometimes necessary if the data is 'too
    ## clean' and clusters become lines..
    ##sX = sX + 0.001*matrix(rnorm(length(sX)),nrow(sX),ncol(sX))

    ## ------------ find t-SNE clusters
    if(is.null(perplexity)) {
        perplexity = max(1,min(30, round((ncol(sX)-1)/4)))
        perplexity
    }

    ##sv.rank=20
    if(sv.rank > 0) {
        cat("performing tSNE on reduced PCA k=",sv.rank,"\n")
        svd <- irlba(sX, nv=sv.rank)
        sv <- svd$v %*% diag(svd$d[1:ncol(svd$v)])
        rownames(sv) <- colnames(sX)
        sX <- t(sv)
    }

    pos2=pos3=NULL
    if(method=="umap") {
        require(umap)
        if(2 %in% dims) {
            pos2 = umap(
                t(sX),
                n_neighbours = perplexity,
                n_components = 2,
                metric = "euclidean"
            )$layout
            colnames(pos2) <- c("umap_1","umap_2")
        }
        if(3 %in% dims) {
            pos3 = umap(
                t(sX),
                n_neighbours = perplexity,
                n_components = 3,
                metric = "euclidean"
            )$layout
            colnames(pos3) <- c("umap_1","umap_2","umap_3")
        }
    } else if(method=="tsne") {
        if(2 %in% dims) {
            pos2 = Rtsne( t(sX), dim=2, perplexity=perplexity,
                         check_duplicates=FALSE, num_threads=99)$Y
            colnames(pos2) <- c("tnse_1","tnse_2")
        }
        if(3 %in% dims) {
            pos3 = Rtsne( t(sX), dim=3, perplexity=perplexity,
                         check_duplicates=FALSE, num_threads=99)$Y
            colnames(pos3) <- c("tnse_1","tnse_2","tnse_3")
        }
    } else if(method=="pca") {
        sv.rank
        svd <- irlba(sX, nv=3)
        if(2 %in% dims) {
            pos2 = svd$v[,1:2]
            colnames(pos2) <- c("pca_1","pca_2")
        }
        if(3 %in% dims) {
            pos3 = svd$v[,1:3]
            colnames(pos3) <- c("pca_1","pca_2","pca_3")
        }
    }
    if(!is.null(pos2)) {
        rownames(pos2) = colnames(sX)
    }
    if(!is.null(pos3)) {
        rownames(pos3) = colnames(sX)
    }
    
    ## ------------ find t-SNE clusters from graph
    idx = NULL
    if(find.clusters) {
        cat("Finding clusters...\n")
        require(igraph)
        if(!is.null(pos2)) pos <- pos2
        if(!is.null(pos3)) pos <- pos3
        dist = as.dist(dist(scale(pos))) ## use 3D distance??
        
        ##dist = dist + 0.1*mean(dist)
        gr = graph_from_adjacency_matrix(1.0/dist, diag=FALSE, mode="undirected")
        ## should we iteratively cluster???
        hc <- hclustGraph(gr)  ##
        dim(hc)
        idx <- hc[,min(kclust,ncol(hc))]
        ##idx <- cluster_louvain(gr)$membership
        table(idx)
        
        ## ------------ zap small clusters to "0"
        sort(table(idx))
        min.size <- pmax(3, 0.01*length(idx))
        min.size
        small.clusters <- names(which(table(idx) < min.size))
        idx[ which(idx %in% small.clusters)] <- "0"
        sort(table(idx))
        
        ## rename levels with largest cluster first
        idx <- factor(idx, levels=names(sort(-table(idx))))
        levels(idx) <- paste0(prefix,1:length(levels(idx)))
        table(idx)
        cat("Found",length(unique(idx)),"clusters...\n")
    }

    res <- list(pos2d=pos2, pos3d=pos3, idx=idx)
    return(res)
}

##levels="gene";contrast="Bmem_activation";layout=NULL;gene="IRF4";layout="layout_with_fr";hilight=NULL
pgx.getGeneFamilies <- function(genes, FILES="../files", min.size=10, max.size=500)
{
    ##-----------------------------------------------------------------------------
    ## Gene families
    ##-----------------------------------------------------------------------------
    ##xgene = genes[,"gene_name"]
    ##xtitle = genes[,"gene_title"]
    families <- list()
    families[["<all>"]] = genes  ## X is sorted

    gmt.kea  <- read.gmt(file.path(FILES,"gmt/kinase_substrates_kea.gmt"))
    gmt.chea <- read.gmt(file.path(FILES,"gmt/tf_targets_chea.gmt"))
    families[["Kinases (KEA)"]] = names(gmt.kea)
    families[["Transcription factors (ChEA)"]] = names(gmt.chea)

    ## Read standard HGNC gene families (www.genefamilies.org)
    ##gmt.hgnc <- read.gmt(file.path(FILES,"gmt/hgnc_genefamilies.gmt"))
    gmt.hgnc <- read.gmt(file.path(FILES,"gmt/hgnc_genefamilies_EDITED.gmt"))

    gmt.hgnc.size <- sapply(gmt.hgnc,length)
    gmt.hgnc <- gmt.hgnc[ which( gmt.hgnc.size >=50 & gmt.hgnc.size <= 1000)]
    length(gmt.hgnc)
    names(gmt.hgnc) <- paste0(names(gmt.hgnc)," (HGNC)")

    families <- c(families, gmt.hgnc)

    ##xgene = as.character(genes$gene_name)
    ##gtype = as.character(genes$gene_type)
    ##gtype = gsub("TR_.*gene","TR_gene",gtype)
    ##gtype = gsub(".*pseudogene","pseudogene",gtype)
    ##sort(table(gtype))
    ##families[["<protein_coding>"]] = genes[which(gtype=="protein_coding")]
    ##families[["Antisense genes"]] = genes[which(gtype=="antisense")]
    ##families[["Pseudogenes"]] = genes[which(gtype=="pseudogene")]
    ##families[["TR genes"]] = genes[which(gtype=="TR_gene")]
    ##families[["lincRNA"]] = genes[which(gtype=="lincRNA")]
    ##families[["snoRNA"]] = genes[which(gtype=="snoRNA")]
    ##families[["Arginine related"]] = genes[grep("arginine",tolower(xtitle))]
    ##families[["CD family"]] = genes[grep("^CD[1-9]",genes)]
    families[["Interleukins (IL)"]] = genes[grep("^IL[1-9]",genes)]
    families[["Chemokines"]] = genes[grep("CCL|CCR|CXCR|CXCL|XCL|CX3",genes)]
    families[["Ribosomal proteins"]] = genes[grep("^RPS|^RPL",genes)]
    families[["Ribosomal (mitochondrial)"]] = genes[grep("^MRPL|^MRPS",genes)]
    ## families[["TR family"]] = genes[grep("^TR[ABDG][VJC]",genes)]
    families[["G-protein family"]] = genes[grep("^GN[ABG]|^GPBAR|^GPER|^FZD",genes)]
    families[["Heatshock proteins"]] = genes[grep("^HSP",genes)]
    families[["Integrin family"]] = genes[grep("^ITG",genes)]
    families[["MAPK family"]] = genes[grep("^MAP[1-9]|^MAPK",genes)]
    ##families[["Olfactory receptors"]] = genes[grep("^OR[1-9]",genes)]
    families[["Myosins"]] = genes[grep("^MYO[1-9]",genes)]
    families[["Protein phosphatase (PPP)"]] = genes[grep("^PPP",genes)]
    families[["PTP family"]] = genes[grep("^PTP",genes)]
    families[["Small nucleolar RNA"]] = genes[grep("^SNOR",genes)]
    families[["Toll-like receptors"]] = genes[grep("^TLR",genes)]
    ##families[["Zinc-fingers C2H2-type"]] = genes[grep("^ZNF",genes)]
    families[["EIF factors"]] = genes[grep("^EIF",genes)]
    families[["GPR proteins"]] = genes[grep("^GPR",genes)]
    families[["CDCC proteins"]] = genes[grep("^CDCC",genes)]
    families[["KIAA proteins"]] = genes[grep("^KIAA",genes)]
    families[["TRIM proteins"]] = genes[grep("^TRIM",genes)]
    families[["TNF proteins"]] = genes[grep("^TNF",genes)]
    families[["SLC proteins"]] = genes[grep("^SLC",genes)]
    families[["ZBTB proteins"]] = genes[grep("^ZBTB",genes)]
    families[["TMEM family"]] = genes[grep("^TMEM",genes)]
    families[["STAT family"]] = genes[grep("^STAT",genes)]
    families[["DNA/RNA polymerases"]] = genes[grep("^POL",genes)]
    families[["Proteasome"]] = genes[grep("^PSM",genes)]
    families[["IFN/IFIT family"]] = genes[grep("^IFN|^IFIT",genes)]
    families[["Nuclear receptors"]] = genes[grep("^NR[0-9]|^RXR|^ESR|^PGR$|^AR$|^HNF4|^ROR|^PPAR|^THR|^VDR", genes)]
    families[["Cytochrome family"]] = genes[grep("^CYP|^CYB|^CYC|^COX|^COA",genes)]
    families[["Micro RNA"]] = genes[grep("^MIR",genes)]

    ## add pathways?
    ##kk = grep("^BIOCARTA_",names(ngs$gmt.all))
    ##kk = head(kk,5) ## just try
    ##pathways = ngs$gmt.all[kk]
    ##families = c(families, pathways)

    ## convert to mouse???
    is.mouse = (mean(grepl("[a-z]",sub(".*:","",genes))) > 0.8)
    is.mouse
    if(is.mouse) {
        library(org.Mm.eg.db)
        mouse.genes = as.character(unlist(as.list(org.Mm.egSYMBOL)))
        names(mouse.genes) = toupper(mouse.genes)
        families <- mclapply(families, function(s) setdiff(as.character(mouse.genes[toupper(s)]),NA),
                             mc.cores=16)
    }

    ## sort all
    families = lapply(families, function(f) intersect(f,genes))

    ## ----------- filter on size
    nsize <- sapply(families, length)
    sel <- which(nsize >= min.size & nsize < max.size)
    sel <- unique(c(which(names(families)=="<all>"),sel))
    families <- families[sel]

    return(families)
}

pgx.getGeneSetCollections <- function(gsets, min.size=10, max.size=500)
{
    ##-----------------------------------------------------------------------------
    ## Gene set collections
    ##-----------------------------------------------------------------------------
    ##gsets = rownames(ngs$gsetX)

    kegg0 <- sort(gsets[grep("KEGG|.*hsa[0-9]{5}$",gsets)])
    kegg0 <- kegg0[!duplicated(getKeggID(kegg0))]

    collections = list(
        "Hallmark collection" = gsets[grep("HALLMARK",gsets)],
        "KEGG pathways" = kegg0,
        "KEGG metabolic pathways" = kegg0[grep("^00|^01",getKeggID(kegg0))],
        ##"REACTOME pathways" = gsets[grep("^REACTOME",gsets)],
        ##"Gene Ontology" = gsets[grep("^GO[_:]",gsets)],
        ##"GO biological processes" = gsets[grep("^GOBP",gsets)],
        ##"GO cellular compartment" = gsets[grep("^GOCC",gsets)],
        ##"GO molecular function" = gsets[grep("^GOMF",gsets)],
        ##"MSigDB C2" = gsets[grep("^C2",gsets)],
        ##"MSigDB C3" = gsets[grep("^C3",gsets)],
        ##"MSigDB C6" = gsets[grep("^C6",gsets)],
        ##"MSigDB C7" = gsets[grep("^C7",gsets)],
        ##"Kinase substrates (KEA)" = gsets[grep("^KEA",gsets)],
        ##"TF targets (ChEA)" = gsets[grep("^CHEA",gsets)],
        ## "Drug signatures" = gsets[grep("DRUG",gsets,ignore.case=TRUE)],
        ##"Drug signatures (up)" = gsets[grep("^DRUG_.*_UP$",gsets)],
        ##"Drug signatures (down)" = gsets[grep("^DRUG_.*_DN$",gsets)],
        ##"PID pathways" = gsets[grep("^PID_",gsets)],
        ##"PANTHER pathways" = gsets[grep("PANTHER",gsets)],
        ##"BIOCARTA pathways" = gsets[grep("^BIOCARTA_",gsets)],
        "Pathway related" = gsets[grep("pathway",gsets,ignore.case=TRUE)],
        ##"Kinase related" = gsets[grep("kinase",gsets,ignore.case=TRUE)],
        "Metabolism related" = gsets[grep("metaboli",gsets,ignore.case=TRUE)],
        "Signalling related" = gsets[grep("signal",gsets,ignore.case=TRUE)],
        "T-cell related" = gsets[grep("tcell|t-cell|t[ ]cell",gsets,ignore.case=TRUE)],
        "B-cell related" = gsets[grep("bcell]b-cell|b[ ]cell",gsets,ignore.case=TRUE)],
        "Response related" = gsets[grep("response",gsets,ignore.case=TRUE)],
        "Cancer related" = gsets[grep("cancer",gsets,ignore.case=TRUE)],
        "Immune related" = gsets[grep("immune",gsets,ignore.case=TRUE)],
        "Cell differentiation" = gsets[grep("differentiation",gsets,ignore.case=TRUE)],
        "Checkpoint related" = gsets[grep("checkpoint",gsets,ignore.case=TRUE)],
        "IL gene sets" = gsets[grep("IL[1-9]{1,2}",gsets,ignore.case=TRUE)]
    )

    collections[["<all>"]] = gsets  ## X is sorted
    collections = collections[which(sapply(collections,length)>=10)]
    collections = collections[order(names(collections))]

    ## ----------- add main collections from gene set prefixes
    gsets.db = sub(":.*","", gsets)
    gsets.groups = tapply(gsets, gsets.db, list)
    collections <- c(collections, gsets.groups)

    ## ----------- filter on size
    nsize <- sapply(collections, length)
    sel <- which( nsize >= min.size & nsize < max.size )
    collections <- collections[sel]
    return(collections)
}


##-----------------------------------------------------------------------------
## Generic module functions
##-----------------------------------------------------------------------------

filterFamily <- function(genes, family, ngs) {
    ##ngs <- isolate(inputData())
    gg = ngs$families[[10]]
    gg = ngs$families[[family]]
    ## check probe name, short probe name or gene name for match
    p0 = (toupper(sub(".*:","",rownames(genes))) %in% toupper(gg))
    p1 = (toupper(rownames(genes))  %in% toupper(gg))
    p2 = (toupper(as.character(genes$gene_name))  %in% toupper(gg))
    jj = which(p0 | p1 | p2 )
    rownames(genes)[jj]
}

filterProbes <- function(genes, gg) {
    ##genes = ngs$families[[10]]
    ##genes = ngs$families[[family]]
    ## check probe name, short probe name or gene name for match
    p0 = (toupper(sub(".*:","",rownames(genes))) %in% toupper(gg))
    p1 = (toupper(rownames(genes))  %in% toupper(gg))
    p2 = (toupper(as.character(genes$gene_name))  %in% toupper(gg))
    jj = which(p0 | p1 | p2)
    if(length(jj)==0) return(NULL)
    return( rownames(genes)[jj] )
}

computeFeatureScore <- function(X, Y, features)
{
    ##features=X=NULL
    ##Y = ngs$Y[,grep("group|sample|patient",colnames(ngs$Y),invert=TRUE)]
    sdx = apply(X,1,sd)
    names(sdx) = rownames(X)
    S = matrix(NA, nrow=length(features), ncol=ncol(Y))
    rownames(S) = names(features)
    colnames(S) = colnames(Y)
    for(k in 1:ncol(Y)) {
        grp = Y[colnames(X),k]
        grp = as.character(grp)
        score = rep(NA, length(features))
        names(score) = names(features)
        i=1
        for(i in 1:length(features)) {
            pp = features[[i]]
            ## if(input$cl_level=="gene") pp = filterProbes(ngs$GENES, features[[i]])
            pp = head(pp[order(-sdx[pp])],100)
            mx = t(apply(X[pp,], 1, function(x) tapply(x,grp,mean)))
            ##D = as.matrix(dist(t(mx)))
            D = 1 - cor(mx, use="pairwise")
            diag(D) = NA
            score[i] = mean(D,na.rm=TRUE)
        }
        S[,k] = score
    }
    if(is.null(S)) return(NULL)
    return(S)
}


##-----------------------------------------------------------------------------
## KEGG
##-----------------------------------------------------------------------------

getKeggID <- function(gsets)
{
    ## Guess KEGG id from gene set name
    ##
    require(KEGGgraph)
    require(KEGG.db)
    require(pathview)
    kegg.names <- unlist(as.list(KEGG.db::KEGGPATHID2NAME))
    kegg.ids <- names(kegg.names)
    kegg.namesUPPERCASE <- toupper(gsub("[- ]","_",gsub("[)(/,.']","",kegg.names)))
    kegg.namesUPPERCASE <- gsub("__|___","_",kegg.namesUPPERCASE)

    k1 <- grep("hsa[0-9]{5}$",gsets,value=TRUE,ignore.case=TRUE)
    names(k1) <- sub(".*_hsa0","0",k1)
    k2 <- grep("kegg",gsets,value=TRUE,ignore.case=TRUE)
    k2x <- sub(".*KEGG_","",k2)
    names(k2) <- kegg.ids[match(k2x,kegg.namesUPPERCASE)]

    kegg.gsets <- c(k1, k2)
    kegg.ids <- names(kegg.gsets)[match(gsets, kegg.gsets)]
    kegg.ids

    return(kegg.ids)
}

##-----------------------------------------------------------------------------
## Generic helper functions
##-----------------------------------------------------------------------------

##mat=ngs$X;group=ngs$samples$group;FUN=mean
averageByGroup <- function(mat, group, FUN=mean) {
    ##sum(duplicated(ngs$genes$gene_name))
    ##out <- t(apply(mat, 1, function(x) tapply(x, group, FUN)))
    out <- do.call(cbind, tapply(1:ncol(mat), group,
                                 function(i) rowMeans(mat[,i,drop=FALSE])))
    return(out)
}

makeAcronym <- function(x) {
    xp <- strsplit(x,split="[_ -]")
    sapply(xp, function(s) {
        if(length(s)==1) return(substring(s,1,2))
        toupper(paste(substring(s,1,1),collapse=""))
    })
}


relevelFactorFirst <- function(f) {
    factor(f, levels=f[!duplicated(f)])
}

extremeCorrelation <- function(query_sig, ref_set, n=200) {
    gg <- intersect(rownames(ref_set), names(query_sig))
    if(n>0) {
        gg <- gg[unique(c(head(order(query_sig),n),head(order(-query_sig),n)))]
    }
    rho <- cor( ref_set[gg,], query_sig[gg], use="pairwise")
    rho <- rho[order(-rowMeans(rho**2,na.rm=TRUE)),,drop=FALSE]
    if(NCOL(rho)==1) rho <- rho[,1]
    return(rho)
}

##s=symbol
alias2hugo <- function(s) {
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

breakstringBROKEN <- function(s, n, force=FALSE) {
    if(is.null(s) || length(s)==0) return(NULL)
    if(length(s)==1) return(breakstring1(s,n=n,force=force))
    sapply(s, breakstring1,n=n,force=force)
}

##s="breakstringBROKENbreakstringBROKENbreakstringBROKENbreakstringBROKEN";n=10
breakstring <- function(s, n, nmax=999, force=FALSE, brk='\n') {
    if(is.na(s)) return(NA)
    s <- substring(as.character(s),1,nmax)
    if(nchar(s)<n) return(s)
    b = substring(s,1,n)
    n1 <- n+1
    for(i in 1:10) {
        if(n1>nchar(s)) break
        b1 <- substring(s,n1,n1+n)
        b <- paste0(b,brk,b1)
        n1 <- n1+n+1
    }
    return(b)
}

##s="breakstringBROKENbreakstringBROKENbreakstringBROKENbreakstringBROKEN";n=10
##n=20;brk="\n"
breakstring2 <- function(s, n, brk='\n', nmax=999) {
    if(is.na(s)) return(NA)
    s <- substring(as.character(s),1,nmax)
    if(is.na(s)) return(NA)
    if(nchar(s)<n) return(s)
    a = ""
    words <- paste0(strsplit(s, split=" ")[[1]]," ")
    words
    it =0
    len1 = sum(sapply(words,nchar))
    len1
    while(len1>n && it<100 && length(words)>1) {
        len = cumsum(sapply(words,nchar))
        len
        k = which(len>n)[1]
        k
        if(k==1) k=2
        a <- paste0(a,paste0(words[1:(k-1)],collapse=""),brk)
        a
        words <- words[k:length(words)]
        words
        it = it+1
        len1 = sum(sapply(words,nchar))
        len1
    }
    a <- paste0(a,paste0(words,collapse=""),brk)
    a <- sub(paste0(brk,"$"),"",gsub(paste0(" ",brk),brk,a))
    return(a)
}

shortstring <- function(s, n, dots=1 ) {
    sapply(s, shortstring0, n=n, dots=dots)
}

shortstring0 <- function(s, n, dots=1 ) {
    ##s0 = as.character(s)
    s0 <- iconv(as.character(s),to="UTF-8")
    s0 <- gsub("[&].*[;]","",s0) ## HTML special garbage...
    jj <- which(nchar(s0)>n)
    if(length(jj)==0) return(s0)
    s <- s0[jj]
    if(dots<1) {
        n1 <- ceiling(dots*n)
        s1 <- substring(s,1, n1 )
        s2 <- substring(s, nchar(s) - n + nchar(s1),nchar(s))
        aa <- paste0(s1,"...",s2)
    } else {
        aa <- paste0(substring(s,1,n),"...")
    }
    s0[jj] <- aa
    s0
}

psort <- function(x,p.col=NULL) {
    j = grep("p.value|^p$|p-val|pval",tolower(colnames(x)))[1]
    x[order(x[,j]),]
}

color_from_middle <- function (data, color1,color2) {
    ## from https://stackoverflow.com/questions/33521828/
    max_val=max(abs(data),na.rm=TRUE)
    DT::JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",max_val,color1,max_val,color1,color2,color2,max_val,max_val))
}

tidy.dataframe <- function(Y) {
    require(tidyverse)
    ##as_tibble(Y)
    Y <- Y[,which(colMeans(is.na(Y))<1)]
    Y <- apply(Y,2,function(x) sub("NA",NA,x)) ## all characters
    Y <- Y[,which(colMeans(is.na(Y))<1)]
    suppressWarnings( num.Y <- apply(Y,2,function(x) as.numeric(as.character(x))) )
    is.numeric <- ( 0.8*colMeans(is.na(num.Y)) <= colMeans(is.na(Y)) )
    nlevel <- apply(Y,2,function(x) length(unique(x)))
    is.factor  <- (!is.numeric | (is.numeric & nlevel<=3) )
    is.factor <- ( is.factor | grepl("batch|replicat|type|clust|group",colnames(Y)) )
    new.Y <- data.frame(Y)
    new.Y[,which(is.numeric)] <- num.Y[,which(is.numeric)]
    new.Y[,which(is.factor)]  <- apply(Y[,which(is.factor)],2,
                                       function(a) factor(as.character(a)))
    new.Y <- data.frame(new.Y)
    return(new.Y)
}


is.num <- function(y) {
    suppressWarnings(numy <- as.numeric(as.character(y)))
    t1 <- !all(is.na(numy)) && is.numeric(numy)
    t2 <- length(unique(y)) > 0.33*length(y)
    (t1 && t2)
}

isanumber <- function(x) {
    x <- sub("NA",NA,x)
    x[which(x=="")] <- NA
    suppressWarnings(nx <- as.numeric(x[!is.na(x)]))
    (length(nx)>0 && mean(!is.na(nx)) > 0.5)
}

expandAnnotationMatrix <- function(A) {
    expandPhenoMatrix(A)
}

expandAnnotationMatrixSAVE <- function(A) {
    ## get expanded annotation matrix
    nlevel <- apply(A,2,function(x) length(unique(x)))
    y.isnum <- apply(A,2,is.num)
    ##kk <- (y.isnum | (!y.isnum & nlevel>1 & nlevel<=5))
    ##A <- A[,which(kk),drop=FALSE]
    head(A)
    i=1
    m1 <- list()
    for(i in 1:ncol(A)) {
        if(is.num(A[,i])) {
            m0 <- matrix(rank(A[,i]), ncol=1)
            colnames(m0) <- colnames(A)[i]
        } else {
            x <- as.character(A[,i])
            x[is.na(x)] <- "_"
            m0 <- model.matrix( ~ 0 + x)
            colnames(m0) <- sub("^x","",colnames(m0))
        }
        if(NCOL(m0)>1) {
            colnames(m0) <- paste0(colnames(A)[i],"=",colnames(m0))
        }
        m1[[i]] <- m0
    }
    names(m1) <- colnames(A)
    ##m1 <- lapply(m1, function(m) {colnames(m)=sub("^x","",colnames(m));m})
    M <- do.call( cbind, m1)
    rownames(M) <- rownames(A)
    return(M)
}

expandPhenoMatrix <- function(pheno, collapse=TRUE) {
    ## get expanded annotation matrix
    ##a1 <- pheno[,grep("group|sample",colnames(pheno),invert=TRUE),drop=FALSE]
    ##m2 <- expandAnnotationMatrix(a1)
    a1 <- tidy.dataframe(pheno)
    nlevel <- apply(a1,2,function(x) length(setdiff(unique(x),NA)))
    nterms <- colSums(!is.na(a1))
    y.class <- sapply(a1,class)
    ##y.isnum <- apply(a1,2,is.num)
    y.isnum <- (y.class == "numeric")
    nlevel
    nratio <- nlevel/nterms
    nratio
    kk <- which(y.isnum | (!y.isnum & nlevel>1 & nratio < 0.66 ))
    a1 <- a1[,kk,drop=FALSE]
    a1.isnum <- y.isnum[kk]
    head(a1)
    i=1
    m1 <- list()
    for(i in 1:ncol(a1)) {
        if( a1.isnum[i] ) {
            ##m0 <- matrix(rank(a1[,i]), ncol=1)
            suppressWarnings( x <- as.numeric(a1[,i]) )
            m0 <- matrix( (x > median(x,na.rm=TRUE)), ncol=1)
            colnames(m0) <- "high"
        } else if(nlevel[i]==2) {
            x <- as.character(a1[,i])
            x1 <- tail(sort(x),1)
            m0 <- matrix(x==x1, ncol=1)
            colnames(m0) <- x1
        } else {
            x <- as.character(a1[,i])
            x[is.na(x) | x=="NA" | x==" "] <- "_"
            m0 <- model.matrix( ~ 0 + x)
            colnames(m0) <- sub("^x","",colnames(m0))
        }
        rownames(m0) <- rownames(a1)
        ## remove "_"
        if("_" %in% colnames(m0)) {
            m0 <- m0[,-which(colnames(m0)=="_")]
        }

        m1[[i]] <- m0
    }
    names(m1) <- colnames(a1)
    if(collapse) {
        for(i in 1:length(m1)) {
            colnames(m1[[i]]) <- paste0(names(m1)[i],"=",colnames(m1[[i]]))
        }
        m1 <- do.call(cbind, m1)
    }
    return(m1)
}

correctMarchSeptemberGenes <- function(gg) {
    sep.from <- c( paste0("0",1:9,"-Sep"), paste0(1:19,"-Sep"))
    sep.to <- c( paste0("SEPT",1:9), paste0("SEPT",1:19))
    mar.from <- c( paste0("0",1:9,"-Mar"), paste0(1:19,"-Mar"))
    mar.to <- c( paste0("MARCH",1:9), paste0("MARCH",1:19))

    require(plyr)
    from <- c(sep.from, mar.from)
    to <- c(sep.to, mar.to)
    gg1 <- sub("[.-]Sep$|[.-]SEP$","-Sep",gg)
    gg1 <- sub("[.-]Mar$|[.-]MAR$","-Mar",gg1)
    jj <- which( from %in% gg1)
    gg2 <- gg1
    if(length(jj)>0) {
        cat("Found ",length(jj),"Sept/Mar genes!\n")
        length(jj)
        from[jj]
        gg2 <- mapvalues(gg1, from[jj], to[jj])
    }
    return(gg2)
}

cor.pvalue <- function(x,n) pnorm(-abs(x/((1-x**2)/(n-2))**0.5))

##k=NULL;mc.cores=2
hclustGraph <- function(g, k=NULL, mc.cores=2)
{
    ## Hierarchical clustering of graph using iterative Louvain
    ## clustering on different levels. If k=NULL iterates until
    ## convergences.
    ##
    require(parallel)
    idx = rep(1, length(V(g)))
    K = c()
    maxiter=100
    if(!is.null(k)) maxiter=k
    iter=1
    ok=1
    idx.len = -1
    while( iter <= maxiter && ok ) {
        old.len = idx.len
        newidx0 = newidx = idx
        i=idx[1]
        if(mc.cores>1 && length(unique(idx))>1) {
            idx.list = tapply(1:length(idx),idx,list)
            mc.cores
            system.time( newidx0 <- mclapply(idx.list, function(ii) {
                subg = induced_subgraph(g, ii)
                subi = cluster_louvain(subg)$membership
                return(subi)
            }, mc.cores=mc.cores) )
            newidx0 = lapply(1:length(newidx0), function(i) paste0(i,"-",newidx0[[i]]))
            newidx0 = as.vector(unlist(newidx0))
            newidx = rep(NA,length(idx))
            newidx[as.vector(unlist(idx.list))] = newidx0
        } else {
            for(i in unique(idx)) {
                ii = which(idx==i)
                subg = induced_subgraph(g, ii)
                subi = cluster_louvain(subg)$membership
                newidx[ii] = paste(i,subi,sep="-")
            }
        }
        vv = names(sort(table(newidx),decreasing=TRUE))
        idx = as.integer(factor(newidx, levels=vv))
        K = cbind(K, idx)
        idx.len = length(table(idx))
        ok = (idx.len > old.len)
        iter = iter+1
    }
    if(NCOL(K)==1) K <- matrix(K, ncol=1)
    rownames(K) = V(g)$name
    if(!ok && is.null(k)) K = K[,1:(ncol(K)-1),drop=FALSE]
    dim(K)
    ##K = K[,1:(ncol(K)-1)]
    colnames(K) <- NULL
    return(K)
}

##=====================================================================================
##=========================== END OF FILE =============================================
##=====================================================================================
