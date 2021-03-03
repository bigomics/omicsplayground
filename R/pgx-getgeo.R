##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


if(0) {

    devtools::install_github("mensxmachina/BioDataome")
    
        
    library(GEOquery)
    res <- getGEO(id)
    
    a1 <- (parsePhenoFromName(geo.getTitle(geo[[1]])))
    a2 <- (parsePhenoFromName(geo.getTitle(geo[[2]])))    
    
    OPG = "/home/kwee/bigomics/omicsplayground"
    RDIR = file.path(OPG,"R")
    FILES = file.path(OPG,"lib")
    PGX.DIR = file.path(OPG,"data")
    source(file.path(RDIR,"pgx-include.R"))


    fd <- fread("~/Downloads/GSE56192_GeneLevel_Raw_data.csv")
    X <- as.matrix(fd[,3:ncol(fd)])
    rownames(X) <- fd$gene_symbol
    dim(X)
    meta1 <- pgx.getGeoMetadata("GSE56192")         
    meta2 <- pgx.getGeoMetadata.fromGSM(colnames(X))
    

    ## BiocManager::install("GEOmetadb")    
    id = "GSE21653"  ## BRCA
    id = "GSE10846"  ## DLBCL
    id = "GSE53784" ## WINJEV  
    id = "GSE56192" ## SARS

}


##-------------------------------------------------------------------------------------
## Query GEO
##-------------------------------------------------------------------------------------
id="GSE100035"
id="GSE102908"
id="GSE102908"

pgx.getGEOseries <- function(id, convert.hugo=TRUE) {
    
    ## get data/pheno matrices
    archs.h5 = file.path(FILESX,"human_matrix.h5")
    counts <- pgx.getGEOcounts(id, archs.h5=archs.h5)
    dim(counts)

    ## get sample info
    meta <- pgx.getGeoMetadata(id)     
    dim(meta)
    
    ## conform matrices
    samples <- intersect(rownames(meta),colnames(counts))
    meta <- meta[samples,,drop=FALSE]
    counts <- counts[,samples,drop=FALSE]    

    ## convert to latest official HUGO???
    if(convert.hugo) {
        sum(duplicated(rownames(counts)))
        symbol <- alias2hugo(rownames(counts))  ## auto-detect mouse/human
        sum(duplicated(symbol))
        counts1 <- tapply(1:nrow(counts), symbol, function(i) colSums(counts[i,,drop=FALSE]))
        counts <- do.call(rbind, counts1)
        remove(counts1)
    }
    dim(counts)
    genes <- ngs.getGeneAnnotation(rownames(counts))
        
    ## get categorical phenotypes
    dim(meta)
    meta1 <- apply(meta,2,trimsame)    
    sampleinfo <- pgx.discretizePhenotypeMatrix(
        meta1, min.ncat=2, max.ncat=20, remove.dup=TRUE)
    sampleinfo <- data.frame(sampleinfo, stringsAsFactors=FALSE, check.names=FALSE)
    dim(sampleinfo)
    head(sampleinfo)
    
    ## automagically create contrast matrix
    mingrp=3;slen=15;ref=NA
    ct <- pgx.makeAutoContrasts(sampleinfo, mingrp=3, slen=20, ref=NA)
    is.null(ct)
    
    if(!is.null(ct)) {
        if("group" %in% colnames(sampleinfo)) {
            colnames(sampleinfo) <- sub("group","xgroup",colnames(sampleinfo)) ## backup..
        }
        sampleinfo <- cbind(sampleinfo, group=as.character(ct$group))
    }
    head(sampleinfo)

    contrasts <- NULL
    if(!is.null(ct$contr.matrix)) {
        contrasts <- ct$contr.matrix
    } else {
        contrasts <- ct$exp.matrix
    }
    head(contrasts)
    
    geo <- list(counts=counts, genes=genes,
                samples = sampleinfo, meta=meta, 
                contrasts=contrasts)
    return(geo)
}

pgx.getGEOcounts <- function(id, archs.h5) {

    expr <- NULL
    
    if(!is.null(archs.h5) && is.null(expr)) {
        ##h5.file = "../libx/human_matrix.h5"
        expr <- pgx.getGEOcounts.archs4(id, archs.h5)
    }
    
    if(is.null(expr)) {
        expr <- pgx.getGEOcounts.recount(id)
    }

    if(is.null(expr)) {
        ## Try with GEOquery
        expr <- pgx.getGEOcounts.GEOquery(id)
    }
    
    is.null(expr)
    if(is.null(expr)) {
        cat("WARNING:: could not get GEO expression. please download manually.\n")
        return(NULL)
    }
    dim(expr)
    max(expr)
    
    ## perform linear transformation (unlog) if required
    qx <- as.numeric(quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    qx
    is.count <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
    is.count
    if(!is.count) {
        ## expr[which(expr <= 0 | is.na(expr))] <- 0 ## really???
        expr <- 2**expr
    }
    
    dim(expr)
    expr
}

pgx.getGeoMetadata <- function(id) {
    ##
    ## load series and platform data from GEO
    ##
    
    id
    ## First try without downloading the GSEMatrix
    pheno = NULL
    pheno <- pgx.getGeoMetadata.fromGSM(id) 
    is.null(pheno)
    if(is.null(pheno)) {
        ## try from Eset
        pheno <- pgx.getGeoMetadata.fromEset(id) 
    }

    has.title <- "title" %in% colnames(pheno)
    if(has.title) {
        px <- title2pheno(pheno$title, split=NULL, trim=TRUE, summarize=TRUE)
        if(!is.null(px) && NCOL(px)>0) {
            pheno <- cbind(pheno, px)
        }
    }

    colnames(pheno) <- gsub("[ ]","_",colnames(pheno)) ## no spaces???
    dim(pheno)
    head(pheno)
    pheno
}


##-------------------------------------------------------------------------------------
## Query GEO expression
##-------------------------------------------------------------------------------------
h5.file = "../libx/human_matrix.h5"
id
pgx.getGEOcounts.archs4 <- function(id, h5.file)
{
    require(rhdf5)
    h5ls(h5.file)
    sample.series <- h5read(h5.file, "meta/Sample_series_id")
    sample.series <- strsplit(as.character(sample.series), split="Xx-xX")

    gse.series <- sort(unique(unlist(sample.series)))
    id %in% gse.series

    idx <- which(sapply(sample.series, function(s) id %in% s))
    length(idx)

    if(length(idx)==0) {
        cat("WARNING: series",id,"not in ARCHS4 matrix file\n")
        return(NULL)
    }
    X <- h5read(h5.file,"data/expression", index = list(NULL, idx))

    sample.acc <- h5read(h5.file, "meta/Sample_geo_accession")
    gene_name  <- h5read(h5.file, "meta/genes")
    head(gene_name)
    colnames(X) <- sample.acc[idx]
    rownames(X) <- gene_name
    sum(duplicated(gene_name))

    ## collapse by symbol
    jj <- !is.na(gene_name) & gene_name!=""
    X <- X[jj,]
    gene_name <- gene_name[jj]
    ## sum intensities (linear)
    X2 <- tapply(1:nrow(X),gene_name,function(ii) {
        colSums(X[ii,,drop=FALSE],na.rm=TRUE)  ## not log!!
    })
    X2 <- do.call(rbind, X2)
    
    return(X2)
}


pgx.getGEOcounts.recount <- function(id)
{
    ## Vignette recount-quickstart.html
    ## Load library
    library('recount')
    
    ## Find a project of interest
    project_info <- abstract_search(id)
    project_info$project

    if(length(project_info$project)==0) {
        cat("could not find",id,"in recount database\n")
        return(NULL)
    }
    
    ## Download the gene-level RangedSummarizedExperiment data
    outdir <- file.path(tempdir(),project_info$project)
    download_study(project_info$project, outdir=outdir)

    ## Load the data
    load(file.path(outdir,'rse_gene.Rdata'))
    
    ## Scale counts by taking into account the total coverage per sample
    rse <- scale_counts(rse_gene)

    counts <- assay(rse)
    dim(counts)
    return(counts)
}

pgx.getGEOcounts.GEOquery <- function(id) {
    ## Retrieve expression matrix, phenotype and probe annotation
    ## matrices for a certain GEO id.
    ##
    require(GEOquery)

    ## load series and platform data from GEO
    id
    ##dat <- getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)
    gse <- try( getGEO(id, GSEMatrix=TRUE, getGPL=TRUE) )
    ##gse <- try( getGEO(id, GSEMatrix=TRUE, getGPL=FALSE) )
    class(gse)
    if(class(gse)=="try-error") {
        cat("ERROR: getGEO() error\n")
        return(NULL)
    }
    length(gse)
    attr(gse, "names")

    dim(exprs(gse[[1]]))
    has.expr <- sapply(gse, function(x) nrow(exprs(x))>0)
    has.expr

    supp_file <- gse[[1]]@experimentData@other$supplementary_file
    
    if(!any(has.expr)) {
        cat("WARNING: dataset has no included expression data\n")
        if(nchar(supp_file)>5) {
            cat("Supplementary file available:",supp_file,"\n")
        }
        return(NULL)
    }

    ## select which has expression
    gse <- gse[which(has.expr)]
    length(gse)
    
    ## select preferred platform is multiple exists
    expr.list <- list()
    k=1
    for(k in 1:length(gse)) {

        ## get expression
        eset <- gse[[k]]
        ex <- exprs(eset)
        dim(ex)

        if(ncol(ex)<=3) {
            ## too small dataset
            next()
        }
        
        ## perform linear transformation (unlog) if required
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
        qx
        is.count <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) ||
            (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2) ## from GEO2R script
        is.count
        if(!is.count) {
            ## ex[which(ex <= 0 | is.na(ex))] <- 0
            ex <- 2**ex
        }
        
        ## featuredata
        has.fdata <- !is.null(fData(eset)) && NCOL(fData(eset))>0
        has.fdata
        if(has.fdata) {
            fdata <- fData(eset)
        } else {
            eset@annotation
            gpl.annot <- getGEO(eset@annotation)
            fdata <- Table(gpl.annot)
            fdata <- fdata[match(rownames(ex),fdata$ID),]
            dim(fdata)
        }

        ## collapse by symbol
        colnames(fdata)
        head(fdata)
        probe.col <- grep("symbol|unigene|ensembl|refseq",colnames(fdata),ignore.case=TRUE)[1]
        probe.col
        if(length(probe.col)==0 || is.na(probe.col)) {
            cat("[pgx.getGEOcounts.GEOquery] ERROR: could not determine probe column\n")
            next()
        }

        probe.col <- probe.col[1]
        ##probes <- fdata[,probe.col]
        fsymbol <- probe2symbol(fdata[,probe.col])
        head(fsymbol)
        jj <- !is.na(fsymbol) & fsymbol!=""
        ex <- ex[jj,]
        fsymbol <- fsymbol[jj]
        ## sum intensities (linear)
        ex2 <- tapply(1:nrow(ex),fsymbol,function(ii) {
            colSums(ex[ii,,drop=FALSE],na.rm=TRUE)  ## not log!!
        })
        ex2 <- do.call(rbind, ex2)
        expr.list[[names(gse)[k]]] <- ex2
    }
    length(expr.list)

    if(length(expr.list)==0) {
        return(NULL)
    }
   
    ## collapse all expressions
    if(length(expr.list)>1) {
        probes <- sort(unique(unlist(lapply(expr.list,rownames))))
        expr.list <- lapply(expr.list, function(x) x[match(probes,rownames(x)),] )
        expr <- do.call(rbind, expr.list)
    } else {
        expr <- expr.list[[1]]
    }
    dim(expr)

    return(expr) ## return always linear intensities
}


pgx.getGEOcounts.fromSuppl <- function(id) {
    ## Retrieve expression matrix, phenotype and probe annotation
    ## matrices for a certain GEO id.
    ##
    require(GEOquery)

    ## load series and platform data from GEO
    id
    ##dat <- getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)
    gse <- try( getGEO(id, GSEMatrix=FALSE, getGPL=FALSE) )
    class(gse)
    if(class(gse)=="try-error") {
        cat("ERROR: getGEO() error\n")
        return(NULL)
    }
    length(gse)
    attr(gse, "names")

    supp_file <- gse@header$supplementary_file
    supp_file
    
    ## Fill me!!!!
    ##
    ##
    dim(expr)

    return(expr)
}


##-------------------------------------------------------------------------------------
## Query GEO metadata
##-------------------------------------------------------------------------------------

pgx.getGeoMetadata.fromEset <- function(id) {

    ## If not succesful, try with downloading the GSEMatrix
    suppressMessages(gse <- try(getGEO(id, GSEMatrix=TRUE, getGPL=FALSE)))
    class(gse)
    if(class(gse)=="try-error") {
        res <- list(error="ERROR: pgx.getGeoMetadata.fromEset() error")
        return(res)
    }
    length(gse)
    attr(gse, "names")
    class(gse[[1]])
    nsamples <- sapply(gse, function(s) nrow(pData(phenoData(s))))
    nsamples

    gse <- gse[nsamples>=3]
    length(gse)    
    eset <- gse[[1]]
    pheno.list <- lapply(gse, pgx.getGeoMetadata.fromEset1)
    pheno.list <- lapply(pheno.list, function(m) {
        rn=rownames(m)
        m1=as.matrix(apply(m,2,as.character))
        rownames(m1)=rn
        m1
    })
    
    vars <- sort(unique(unlist(lapply(pheno.list, colnames))))
    pheno.list <- lapply(pheno.list, function(m) {
        m1=m[,match(vars,colnames(m))]
        colnames(m1)=vars
        m1
    })

    gpl <- sub("_series.*","",sub(".*-GPL","GPL",names(gse)))
    pheno.list <- lapply(1:length(pheno.list), function(i) {
        pheno.list[[i]] <- cbind(GPL=gpl[i], pheno.list[[i]])
    })
    lapply(pheno.list,dim)
    pheno <- do.call(rbind, pheno.list)
    dim(pheno)
    pheno <- data.frame(pheno,stringsAsFactor=FALSE,check.names=FALSE)
    return(pheno)
}

pgx.getGeoMetadata.fromEset1 <- function(eset) {
    ##
    ## load series and platform data from GEO
    ##
        
    ## Get summary 
    summary <- experimentData(eset)@abstract

    ## pdata object
    pdata <- pData(phenoData(eset))

    gsm.title <- as.character(pdata$title)
    gsm.source <- as.character(pdata$source_name_ch1)
    gsm.samples <- as.character(pdata$geo_accession)
    head(gsm.samples)
    
    ## Base sample_info from characteristics (ch1) column
    ch1_info <- eset.getCH1(eset)
    dim(ch1_info)
    head(ch1_info)
    
    ## We can get extra information from title
    head(gsm.title)
    is.underscored <- length(gsm.title) && all(grepl("_",gsm.title))
    is.underscored
    title_info <- NULL
    if(FALSE && is.underscored) {
        title2 <- trimws(gsm.title)
        ##title2 <- trimsame(title2, split="_",ends=TRUE)
        ##title2 <- gsub("[\\(\\),-]","",title2)
        ##title2 <- gsub("[ ]","_",title2)
        title_info <- eset.parsePhenoFromTitle(title2,split="_")
        head(title_info)
    }
    
    ## All sample_info: from characterisctis_ch1 and title
    sample_info <- data.frame(GSM=gsm.samples, title=gsm.title,
                              source = gsm.source,
                              stringsAsFactors=FALSE, check.names=FALSE)
    if(!is.null(ch1_info)) sample_info <- cbind(sample_info, ch1_info)
    if(!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

    sample_info <- data.frame(sample_info, stringsAsFactors=FALSE, check.names=FALSE)
    dim(sample_info)    
    ##rownames(sample_info) <- gsm.samples
    sample_info
}

pgx.getGeoMetadata.fromGSM <- function(id) {
    ##
    ## load series and platform data from GEO
    ##
    id
    suppressMessages(gse <- try(getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)))
    class(gse)
    if(class(gse)=="try-error") {
        res <- list(error="ERROR: getGEO() error")
        return(res)
    }
    length(gse)
    attr(gse, "names")
    slotNames(gse)
    class(gse)
    length(gse@gsms)
    
    if(length(gse@gsms)==0) {
        cat("WARNING:: no GSM information in object\n")
        return(NULL)
    }
    
    ## Summary and sample names
    summary <- gse@header$summary
    gsm.samples <- gse@header$sample_id
    ##gsm.samples <- names(dat@gsms)
    head(gsm.samples)
    
    ## Get sample_info from characteristics (ch1) column
    ch1_info <- lapply(gse@gsms, function(g) g@header$characteristics_ch1)
    ##ch1_info <- lapply(ch1_info, function(x) {colnames(x)=colnames(ch1_info[[1]]);x})
    is.null(ch1_info)
    if(!is.null(ch1_info)) {        
        ch1_info <- lapply(ch1_info, function(x) sub("^Clinical info: ","",x))        
        ch1_vars <- unique(unlist(lapply(ch1_info, function(x) trimws(sub("[:=].*","",x)))))
        ch1_info <- lapply(ch1_info, function(x) {
            xvar = trimws(sub("[:=].*","",x))
            x = trimws(sub(".*[:=] ","",x))
            names(x) = xvar
            x = x[match(ch1_vars,names(x))]
            x
        })
        ch1_info <- do.call( rbind, ch1_info )
        dim(ch1_info)
        colnames(ch1_info) <- ch1_vars
        head(ch1_info)
    }
    
    ## We can get more information from title??
    gsm.title  <- sapply(gse@gsms, function(g) g@header$title)
    gsm.source <- sapply(gse@gsms, function(g) g@header$source_name_ch1)
    gsm.gpl    <- sapply(gse@gsms, function(g) g@header$platform_id)    
    head(gsm.title)
    is.underscored <- length(gsm.title) && all(grepl("_",gsm.title))
    is.underscored
    title_info <- NULL
    ## NEED RETHINK!!!!!!!!!!!!!!!!!!
    if(FALSE && is.underscored) {
        title2 <- trimws(gsm.title)
        ##title2 <- trimsame2(title2, split="_")
        ##title2 <- gsub("[\\(\\),-]","",title2)
        ##title2 <- gsub("[ ]","_",title2)
        title_info <- eset.parsePhenoFromTitle(title2,split="_")
        dim(title_info)
        head(title_info)
    }
    
    ## All sample_info: from characterisctis_ch1 and title
    sample_info <- data.frame(GPL=gsm.gpl, GSM=gsm.samples, title=gsm.title,
                              source=gsm.source, stringsAsFactors=FALSE)
    if(!is.null(ch1_info)) sample_info <- cbind(sample_info, ch1_info)
    if(!is.null(title_info)) sample_info <- cbind(sample_info, title_info)

    sample_info <- data.frame(sample_info, stringsAsFactors=FALSE, check.names=FALSE)
    dim(sample_info)    
    ##colnames(sample_info) <- gsub("[ ]","_",colnames(sample_info))
    ##rownames(sample_info) <- gsm.samples
    sample_info
}

##source(file.path(RDIR,"ngs-functions.R"))
##source(file.path(RDIR,"pgx-functions.R"))

##-------------------------------------------------------------------------------------
## HELPER functions
##-------------------------------------------------------------------------------------

eset.getPhenoData <- function(eset,field) { pData(phenoData(eset))[,field] }
eset.getTitle <- function(eset) as.character(pData(phenoData(eset))$title)
eset.getOrganism <- function(eset) unique(as.character(eset.getPhenoData(eset,"organism_ch1")))
eset.getCH1 <- function(eset) {
    pdata <- pData(phenoData(eset))
    pdata <- pdata[,grep(":ch1$",colnames(pdata)),drop=FALSE]
    clin_info <- NULL
    has.clin <- "Clinical info:ch1" %in% colnames(pdata)
    has.clin
    if(has.clin) {
        clin_ch1 <- as.character(pdata[,"Clinical info:ch1"])
        clin.terms <- sub(":.*$","",strsplit(clin_ch1[1], split=";")[[1]])
        clin_info <- t(sapply(clin_ch1, function(s) (strsplit(s,split=";")[[1]]) ))
        clin_info <- apply(clin_info,2,function(s) sub(".*[:] ","",s))
        clin_info <- data.frame(clin_info, stringsAsFactors=FALSE, check.names=FALSE)
        rownames(clin_info) <- NULL
        dim(clin_info)
        colnames(clin_info) <- clin.terms
        pdata <- cbind( pdata[,-which(colnames(pdata)=="Clinical info:ch1"),drop=FALSE],
                       clin_info )
    }
    colnames(pdata) <- sub(":ch1","",colnames(pdata))
    ##colnames(pdata) <- gsub("[ ]","_",colnames(pdata))
    pdata
}

eset.parseCharacteristicsInfo <- function(ch, split=",") {
    terms <- sub(":.*$","",trimws(strsplit(ch, split=",")[[1]]))
    value <- t(sapply(ch, function(s) (strsplit(s,split=",")[[1]]) ))
    value <- apply(value, 2, function(s) trimws(sub(".*[:]","",s)))
    value <- data.frame(value, stringsAsFactors=FALSE, check.names=FALSE)
    rownames(value) <- NULL
    dim(value)
    colnames(value) <- terms
    value
}

##title=title2
title2pheno <- function(title, split=NULL, trim=TRUE, summarize=TRUE)
{
    ##
    ##
    ##

    ## determine the split character
    if(is.null(split)) {        
        split.chars <- c(',',';','\\|','_',' ')
        ss <- c()
        for(i in 1:length(title)) {
            ns <- sapply(split.chars, function(s) sum(gregexpr(s, title[i])[[1]]>0))
            split0 <- names(ns)[which.max(ns)]
            ns1 <- ns[setdiff(names(ns)," ")]
            if(split0 == " " && any(ns1>0)) {
                split0 <- names(ns1)[which.max(ns1)]
            }
            ss[i] <- split0
        }
        split <- names(which.max(table(ss)))
    }
    split

    ## Check if all titles have equal splitted parts (e.g. nicely formatted)
    nsplit <- sapply(title, function(tt) sum(gregexpr(split, tt)[[1]]>0))
    nsplit.equal <- all(nsplit==nsplit[1])
    nsplit.equal
    if(!nsplit.equal) {
        cat("splitted title terms not equal lengths\n")
        return(NULL)
    }
    
    ## split
    tt <- as.character(sapply(as.character(title), function(s) trimws(s)))
    ff <- sapply(as.character(tt), strsplit, split=split)

    ## cleanup
    ff <- lapply(ff, trimws)  ## trim whitespace
    ff <- lapply(ff, function(s) gsub("hours$|hour$|hrs$|hr$","h",s))  ## hours
    ff <- lapply(ff, function(s) gsub("[ ][ ]*"," ",s))  ## double space
    ##ff <- lapply(ff, function(s) gsub("([0-9]*)[ _]h","\\1h",s)) ## remove any space between hours??
    
    ## make dataframe
    F1 <- do.call(rbind, ff)
    F1
    F1[is.na(F1)] <- NA
    rownames(F1) <- NULL

    ## Guess names
    getmax.term <- function(s) {
        tt <- table(unlist(strsplit(s,split="[ _]")))
        vip.tt <- grep("hour|repl|hr|time|treat|infec|pati|sampl",names(tt))
        vip.tt
        if(length(vip.tt)) tt[vip.tt] <- 1.1 * tt[vip.tt] ## boost known keywords
        names(which.max(tt))
    }
    maxterm <- apply(F1,2,function(s) getmax.term(s))
    maxterm <- paste0("_",maxterm)
    maxterm
    colnames(F1) <- maxterm

    ## trims same words/characters on both ends   
    if(trim) F1 <- apply(F1, 2, trimsame, summarize=summarize) 
    F1
}



##title=title2
eset.parsePhenoFromTitle <- function(title, split=NULL)
{
    ##
    ##
    ##
    require(Biostrings)
    require(msa)
    ##title <- as.character(pData(phenoData(eset))$title)

    if(!all(grepl(split,title))) {
        return(NULL)
    }
        
    tt <- as.character(sapply(as.character(title), function(s) trimws(s)))
    tt <- sapply(tt, function(s) gsub("[ ]*hours|[ ]*hour|[ ]*hrs|[ ]*hr","h",s))  ## hours
    tt <- sapply(tt, function(s) gsub("([0-9]*)[ _]h","\\1h",s)) ## remove any space between hours??
    tt <- trimsame(tt,split=split,ends=TRUE)
    tt <- gsub(paste0(split,split,split),split,tt)
    tt <- gsub(paste0(split,split),split,tt)
    ff <- sapply(as.character(tt), strsplit, split=split)
    nf <- max(sapply(ff,length))
    ff <- lapply(ff, function(x) head(c(x,rep(NA,nf)),nf))

    ## cleanup
    ff <- lapply(ff, trimws)  ## trim whitespace
    ff <- lapply(ff, function(s) gsub("hours$|hour$|hrs$|hr$","h",s))  ## hours
    ff <- lapply(ff, function(s) gsub("[ ][ ]*"," ",s))  ## double space
    
    F1 <- do.call(rbind, ff)
    F1
    F1[is.na(F1)] <- NA
    lapply(apply(F1,2,table),sort,decreasing=TRUE)

    AA <- setdiff(unique(GENETIC_CODE),"*")

    G <- list()
    i=1
    for(i in 1:(ncol(F1)-1)) {

        k <- min(ncol(F1),(i+1))
        a2 <- factor(as.vector(F1[,i:k]))    
        aa.dict <- levels(a2)
        names(aa.dict) <- AA[1:length(levels(a2))]
        levels(a2) <- AA[1:length(levels(a2))]        
        F2 <- matrix(a2,nrow(F1))
        F2[is.na(F2)] <- "-"

        ff <- apply(F2,1,paste,collapse="")
        names(ff) <- paste0("tt",1:nrow(F2))
        aln <- as.character(msa(ff,type="protein"))
        aln <- aln[names(ff)]
        F.aln <- do.call(rbind,sapply(aln,strsplit,split=""))
        F.aln2 <- apply(F.aln,2,function(x) aa.dict[x])
        
        if(ncol(F.aln2) > ncol(F2) && i<(ncol(F1)-1)) {
            G[[i]] <- F.aln2[,1:(ncol(F.aln)-1)]
        } else if(i==(ncol(F1)-1)) {
            G[[i]] <- F.aln2
        } else {
            G[[i]] <- F.aln2[,1]
        }
    }
    G <- do.call( cbind, G)
    G <- G[,colMeans(is.na(G))<1,drop=FALSE]
    rownames(G) <- NULL
    ##rownames(G) <- tt
    colnames(G) <- paste0("V",1:ncol(G))    
    G
}
