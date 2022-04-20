##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


if(0) {

    devtools::install_github("mensxmachina/BioDataome")
    
    res <- GEOquery::getGEO(id)    
    a1 <- (parsePhenoFromName(geo.getTitle(geo[[1]])))
    a2 <- (parsePhenoFromName(geo.getTitle(geo[[2]])))    
    
    OPG = "/home/kwee/bigomics/omicsplayground"
    RDIR = file.path(OPG,"R")
    FILES = file.path(OPG,"lib")
    PGX.DIR = file.path(OPG,"data")
    source(file.path(RDIR,"pgx-include.R"))


    fd <- data.table::fread("~/Downloads/GSE56192_GeneLevel_Raw_data.csv")
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
id="GSE141499"
##ARCHSH5 = file.path(FILESX,"human_matrix.h5")

pgx.getGEOseries <- function(id, archs.h5="human_matrix.h5", convert.hugo=TRUE)
{
    ##
    ## Highly automagic download of GEO datasets from different
    ## sources with automatic probe/gene conversion and creating
    ## autocontrasts. The GEO series is first searched in a locally
    ## stored ARCHS4 H5 file, then if it is available at the recount
    ## database, if not it is retrieved from GEO using geoquery.
    ##
    ## id:      GEO id
    ## return:  object with counts, samples, genes.
    ##
    
    
    ## get data/pheno matrices
    geo <- pgx.getGEOcounts(id, archs.h5=archs.h5)
    counts <- geo$expr
    dim(counts)

    ## get sample info
    meta <- pgx.getGeoMetadata(id)     
    names(meta)
    
    ## conform matrices
    samples <- intersect(rownames(meta),colnames(counts))
    meta    <- meta[samples,,drop=FALSE]
    counts  <- counts[,samples,drop=FALSE]    

    ## convert to latest official HUGO???
    if(convert.hugo) {
        symbol <- alias2hugo(rownames(counts))  ## auto-detect mouse/human
        rownames(counts) <- symbol
    }

    ## sum up values duplicated symbols
    ndup <- sum(duplicated(rownames(counts)))
    if(ndup > 0) {
        counts1 <- tapply(1:nrow(counts), symbol, function(i) Matrix::colSums(counts[i,,drop=FALSE]))
        counts <- do.call(rbind, counts1)
        remove(counts1)
    }

    ## get annotation
    dim(counts)
    sum(duplicated(rownames(counts)))
    genes <- ngs.getGeneAnnotation(rownames(counts))
        
    ## get categorical phenotypes
    dim(meta)
    meta1 <- apply(meta,2,trimsame)
    rownames(meta1) <- rownames(meta)
    sampleinfo <- pgx.discretizePhenotypeMatrix(
        meta1, min.ncat=2, max.ncat=20, remove.dup=TRUE)
    sampleinfo <- data.frame(sampleinfo, stringsAsFactors=FALSE, check.names=FALSE)
    dim(sampleinfo)
    Matrix::head(sampleinfo)
    
    ## automagically create contrast matrix
    contrasts <- NULL
    if(NCOL(sampleinfo)>0) {
        mingrp=3;slen=15;ref=NA
        ct <- pgx.makeAutoContrasts(sampleinfo, mingrp=3, slen=20, ref=NA)
        is.null(ct)
        if(is.null(ct)) {
            ct <- pgx.makeAutoContrasts(sampleinfo, mingrp=2, slen=20, ref=NA)            
        }
        ## if(!is.null(ct)) {
        ##     if("group" %in% colnames(sampleinfo)) {
        ##         colnames(sampleinfo) <- sub("group","xgroup",colnames(sampleinfo)) ## backup..
        ##     }
        ##     sampleinfo <- cbind(sampleinfo, group=as.character(ct$group))
        ## }
        if(!is.null(ct$exp.matrix)) {
            contrasts <- ct$exp.matrix
        } else {
            contrasts <- ct$design %*% ct$contr.matrix
        }
        Matrix::head(contrasts)
    }

    info <- pgx.getGeoExperimentInfo(id)    

    out <- list(
        counts  = counts,
        genes   = genes,
        samples = sampleinfo,
        contrasts = contrasts,
        meta    = meta, 
        info    = info,
        source  = geo$source
    )

    return(out)
}

pgx.getGEOcounts <- function(id, archs.h5) {

    expr <- NULL
    src = ""
    
    if(!is.null(archs.h5) && is.null(expr)) {
        ##h5.file = "../libx/human_matrix.h5"
        expr <- pgx.getGEOcounts.archs4(id, archs.h5)
        if(!is.null(expr)) src = "ARCHS4"
    }
    
    if(is.null(expr)) {
        expr <- pgx.getGEOcounts.recount(id)
        if(!is.null(expr)) src = "recount"
    }

    if(is.null(expr)) {
        ## Try with GEOquery
        expr <- pgx.getGEOcounts.GEOquery(id)
        if(!is.null(expr)) src = "GEO"        
    }
    
    is.null(expr)
    if(is.null(expr)) {
        cat("WARNING:: could not get GEO expression. please download manually.\n")
        return(NULL)
    }
    dim(expr)
    max(expr)

    if(0) {
        ## already done in getGEO functions...
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
    }
    
    dim(expr)
    list( expr=expr, source=src )
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

    ## Sometimes the phenotype is coded in the title string 
    has.title <- "title" %in% colnames(pheno)
    if(has.title && NCOL(pheno)==0) {
        px <- title2pheno(pheno$title, split=NULL, trim=TRUE, summarize=TRUE)
        if(!is.null(px) && NCOL(px)>0 && is.null(pheno)) {
            pheno <- px
        }
        if(!is.null(px) && NCOL(px)>0 && !is.null(pheno)) {
            pheno <- cbind(pheno, px)
        }
    }

    colnames(pheno) <- gsub("[ ]","_",colnames(pheno)) ## no spaces???
    dim(pheno)
    Matrix::head(pheno)

    ##
    ##experimentData(eset)
    

    pheno
}


##-------------------------------------------------------------------------------------
## Query GEO expression
##-------------------------------------------------------------------------------------
h5.file = "../libx/human_matrix.h5"
id
pgx.getGEOcounts.archs4 <- function(id, h5.file)
{
    
    rhdf5::h5ls(h5.file)
    sample.series <- rhdf5::h5read(h5.file, "meta/Sample_series_id")
    sample.series <- strsplit(as.character(sample.series), split="Xx-xX")

    gse.series <- sort(unique(unlist(sample.series)))
    id %in% gse.series

    idx <- which(sapply(sample.series, function(s) id %in% s))
    length(idx)

    if(length(idx)==0) {
        cat("WARNING: series",id,"not in ARCHS4 matrix file\n")
        return(NULL)
    }
    X <- rhdf5::h5read(h5.file,"data/expression", index = list(NULL, idx))

    sample.acc <- rhdf5::h5read(h5.file, "meta/Sample_geo_accession")
    gene_name  <- rhdf5::h5read(h5.file, "meta/genes")
    Matrix::head(gene_name)
    colnames(X) <- sample.acc[idx]
    rownames(X) <- gene_name
    sum(duplicated(gene_name))

    ## collapse by symbol
    jj <- !is.na(gene_name) & gene_name!=""
    X <- X[jj,]
    gene_name <- gene_name[jj]
    ## sum intensities (linear)
    X2 <- tapply(1:nrow(X),gene_name,function(ii) {
        Matrix::colSums(X[ii,,drop=FALSE],na.rm=TRUE)  ## not log!!
    })
    X2 <- do.call(rbind, X2)
    
    return(X2)
}


id="GSE10846"
pgx.getGEOcounts.recount <- function(id)
{
    ## Vignette recount-quickstart.html
    ## Load library
    
    
    ## Find a project of interest
    project_info <- recount::abstract_search(id)
    project_info$project

    if(length(project_info$project)==0) {
        cat("could not find",id,"in recount database\n")
        return(NULL)
    }
    
    ## Download the gene-level RangedSummarizedExperiment data
    outdir <- file.path(tempdir(),project_info$project)
    recount::download_study(project_info$project, outdir=outdir)

    ## Load the data
    load(file.path(outdir,'rse_gene.Rdata'))
    
    ## Scale counts by taking into account the total coverage per sample
    rse <- recount::scale_counts(rse_gene)

    counts <- MultiAssayExperiment::assay(rse)
    dim(counts)
    return(counts)
}

pgx.getGEOcounts.GEOquery <- function(id) {
    ## Retrieve expression matrix, phenotype and probe annotation
    ## matrices for a certain GEO id.
    ##
    

    ## load series and platform data from GEO
    id
    ##dat <- GEOquery::getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)
    gse <- try( GEOquery::getGEO(id, GSEMatrix=TRUE, getGPL=TRUE) )
    ##gse <- try( GEOquery::getGEO(id, GSEMatrix=TRUE, getGPL=FALSE) )
    class(gse)
    if(class(gse)=="try-error") {
        cat("ERROR: GEOquery::getGEO() error\n")
        return(NULL)
    }
    length(gse)
    attr(gse, "names")

    dim(exprs(gse[[1]]))
    has.expr <- sapply(gse, function(x) nrow(exprs(x))>0)
    has.expr
    
    if(!any(has.expr)) {
        cat("WARNING: dataset has no included expression data\n")
        supp_file <- sapply(gse, function(g) g@experimentData@other$supplementary_file)
        supp_file    
        if(any(nchar(supp_file))>5) {
            cat("Supplementary file available: ",paste(supp_file,collapse=" "),"\n")
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
            gpl.annot <- GEOquery::getGEO(eset@annotation)
            fdata <- GEOquery::Table(gpl.annot)
            fdata <- fdata[match(rownames(ex),fdata$ID),]
            dim(fdata)
        }

        ## get symbol from featuredata
        colnames(fdata)
        fsymbol <- pgx.getSymbolFromFeatureData(fdata) 
        Matrix::head(fsymbol[!is.na(fsymbol)])

        ## collapse by symbol
        jj <- which(!is.na(fsymbol) & fsymbol!="")
        ex <- ex[jj,]
        fsymbol <- fsymbol[jj]
        ## sum intensities (linear)
        ex2 <- tapply(1:nrow(ex),fsymbol,function(ii) {
            Matrix::colSums(ex[ii,,drop=FALSE],na.rm=TRUE)  ## not log!!
        })
        ex2 <- do.call(rbind, ex2)
        expr.list[[names(gse)[k]]] <- ex2
    }
    length(expr.list)

    if(length(expr.list)==0) {
        return(NULL)
    }
   
    if(length(expr.list)>1) {
        ## merge/join all expressions
        probes <- sort(unique(unlist(lapply(expr.list,rownames))))
        samples <- sort(unique(unlist(lapply(expr.list,colnames))))
        expr.list2 <- lapply(expr.list, function(x)
            x[match(probes,rownames(x)),match(samples,colnames(x))] )
        expr.list2 <- lapply(expr.list2, function(x) {x[is.na(x)]=0;x})
        expr <- Reduce('+', expr.list2)
        colnames(expr) <- samples
        rownames(expr) <- probes        
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
    

    ## load series and platform data from GEO
    id
    ##dat <- GEOquery::getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)
    gse <- try( GEOquery::getGEO(id, GSEMatrix=FALSE, getGPL=FALSE) )
    class(gse)
    if(class(gse)=="try-error") {
        cat("ERROR: GEOquery::getGEO() error\n")
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

pgx.getSymbolFromFeatureData <- function(fdata) {

    ## extract GENE symbol from featureData. The problem is that we don't
    ## know the gene column because the column names are not always
    ## consistent. Also the actual gene symbol may be part of an
    ## annotation string instead of single symbol column.
    require(org.Hs.eg.db)
    
    colnames(fdata)
    symbol <- NULL
    
    ## If there is a symbol column, than it is easy
    SYMBOL <- as.character(unlist(as.list(org.Hs.egSYMBOL)))
    symbol.col <- grep("symbol|gene|hugo",colnames(fdata),ignore.case=TRUE)
    symbol.col
    ok.symbol <- apply( fdata[,symbol.col,drop=FALSE], 2,
                       function(g) mean(toupper(g[!is.na(g)]) %in% SYMBOL))
    ok.symbol
    if(any(ok.symbol > 0.5)) {
        k <- symbol.col[which.max(ok.symbol)]
        symbol <- fdata[,k]
        return(symbol)
    }

    ## If there is an ENTREZ column, than it is easy
    ENTREZ <- biomaRt::keys(org.Hs.egSYMBOL)
    entrez.col <- grep("entrez",colnames(fdata),ignore.case=TRUE)
    entrez.col
    entrez.match <- apply( fdata[,entrez.col,drop=FALSE], 2,
                          function(g) mean(g[!is.na(g)] %in% ENTREZ))
    entrez.match
    entrez.ok <- length(entrez.col) && entrez.match > 0.5
    entrez.ok
    if(entrez.ok) {
        k <- entrez.col[which.max(entrez.match)]
        probes <- as.character(fdata[,k])
        symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, probes, 'SYMBOL', 'ENTREZID')
        return(symbol)
    }

    ## If there is an REFSEQ column
    REFSEQ <- unlist(as.list(org.Hs.egREFSEQ))
    refseq.col <- grep("refseq",colnames(fdata),ignore.case=TRUE)
    refseq.col
    refseq.match <- apply( fdata[,refseq.col,drop=FALSE], 2,
                          function(g) mean(sub("[.].*","",g[!is.na(g)]) %in% REFSEQ))
    refseq.match
    refseq.ok <- length(refseq.col) && refseq.match > 0.5
    refseq.ok
    if(refseq.ok) {
        k <- refseq.col[which.max(refseq.match)]
        probes <- sub("[.].*","",as.character(fdata[,k]))
        symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, probes, 'SYMBOL', 'REFSEQ')
        return(symbol)
    }
    
    ## Otherwise try Ensemble ID
    gene.column <- grep("gene|mrna|transcript", colnames(fdata),ignore.case=TRUE)
    gene.column
    has.ens <- apply(fdata[,gene.column,drop=FALSE], 2, function(s) mean(grepl("ENS",s)))
    has.ens
    if( any(has.ens > 0.3) ) {
        ens.col = ifelse(max(has.ens)>0, names(which.max(has.ens)), NA)        
        ens.ann = lapply(fdata[,ens.col], function(a) trimws(strsplit(a,split="//|///")[[1]]))
        ens.probes <- sapply(ens.ann, function(s) Matrix::head(grep("^ENS",s,value=TRUE),1))
        ens.probes[sapply(ens.probes,length)==0] <- NA
        ens.probes <- unlist(ens.probes)
        symbol <- probe2symbol(ens.probes)
        return(symbol)
    }

    message("WARNING:: could not parse symbol information from featureData!")
    return(NULL)
}

pgx.getGeoExperimentInfo <- function(id) {
    suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)))
    info = gse@header
    ##info <- info[c("geo_accession","title","summary","type")]
    info
}

pgx.getGeoMetadata.fromEset <- function(id) {
    
    ## If not succesful, try with downloading the GSEMatrix
    suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix=TRUE, getGPL=FALSE)))
    ## suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)))    
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
    pheno <- data.frame(pheno, stringsAsFactors=FALSE, check.names=FALSE)
    
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
    Matrix::head(gsm.samples)
    
    ## Base sample_info from characteristics (ch1) column
    ch1_info <- eset.getCH1(eset)
    dim(ch1_info)
    Matrix::head(ch1_info)
    
    ## We can get extra information from title
    Matrix::head(gsm.title)
    is.underscored <- length(gsm.title) && all(grepl("_",gsm.title))
    is.underscored
    title_info <- NULL
    if(FALSE && is.underscored) {
        title2 <- trimws(gsm.title)
        ##title2 <- trimsame(title2, split="_",ends=TRUE)
        ##title2 <- gsub("[\\(\\),-]","",title2)
        ##title2 <- gsub("[ ]","_",title2)
        title_info <- eset.parsePhenoFromTitle(title2,split="_")
        Matrix::head(title_info)
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
    suppressMessages(gse <- try(GEOquery::getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)))
    class(gse)
    if(class(gse)=="try-error") {
        res <- list(error="ERROR: GEOquery::getGEO() error")
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
    Matrix::head(gsm.samples)
    
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
        Matrix::head(ch1_info)
    }
    
    ## We can get more information from title??
    gsm.title  <- sapply(gse@gsms, function(g) g@header$title)
    gsm.source <- sapply(gse@gsms, function(g) g@header$source_name_ch1)
    gsm.gpl    <- sapply(gse@gsms, function(g) g@header$platform_id)    
    Matrix::head(gsm.title)
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
        Matrix::head(title_info)
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
    ff <- lapply(ff, function(x) Matrix::head(c(x,rep(NA,nf)),nf))

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
        aln <- as.character(msa::msa(ff,type="protein"))
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


parse_geo_series_matrix <- function(SERIES_FILE,
                                    PLATFORM_FILE,
                                    GENE_COLUMN = "GENE_SYMBOL",
                                    EXPRESSION_OUTPUT_FILE = "expression.csv",
                                    ANNOTATION_OUTPUT_FILE = "annotation.csv",                                    
                                    write.file=TRUE)
{
    ##
    ## https://gist.github.com/SimonLarsen
    ## https://gist.github.com/SimonLarsen/66f27c188039129f3510669b992b2c99

    library(data.table)
    library(foreach)

    if(0) {
        args <- commandArgs(TRUE)    
        ## Input files ####
        ## Series Matrix File
        SERIES_FILE <- args[1]
        ## Platform data table obtained from GEO.
        PLATFORM_FILE <- args[2]
        
        ## Parameters ####
        ## Column containing gene IDs in platform file #
        ## Generally "ENTREZ_GENE_ID" or "GENE"
        GENE_COLUMN <- if(!is.na(args[3])) args[3] else "ENTREZ_GENE_ID"
        
        ## Output files ####
        EXPRESSION_OUTPUT_FILE <- if(!is.na(args[4])) args[4] else "expression.csv"
        ANNOTATION_OUTPUT_FILE <- if(!is.na(args[5])) args[5] else "annotation.csv"
    }
    
    ## Read characteristics
    con <- file(SERIES_FILE, "r")
    characteristics <- c()
    while(TRUE) {
        line <- readLines(con, n=1)
        if(length(line) == 0) {
            break
        } else if(startsWith(line, "!Sample_title")) {
            titles <- unlist(strsplit(line, "\t"))[-1]
            titles <- gsub("\\\"", "", titles)
        } else if(startsWith(line, "!Sample_characteristics")) {
            characteristics <- c(characteristics, line)
        } else if(startsWith(line, "!Sample_geo_accession")) {
            accession <- unlist(strsplit(line, "\t"))[-1]
            accession <- gsub("\\\"", "", accession)
        }
    }
    close(con)
    
    ## Parse characteristics
    anno <- data.frame(lapply(characteristics, function(x) {
        values <- unlist(strsplit(x, "\t"))[-1]
        values <- gsub("\\\"", "", values)
        parts <- strsplit(values, ": ")
        
        name <- parts[[1]][[1]]
        values <- sapply(parts, function(x) x[2])
        
        out <- list()
        out[[name]] <- values
        return(out)
    }))

    anno <- data.table(sample=accession, title=titles, anno)
    
    ## Read probe-level expression data
    D <- fread(SERIES_FILE, header=TRUE, skip="\"ID_REF\"", fill=TRUE, na.strings=c("","NA","null"))
    D <- D[1:(nrow(D)-1),] # remove table end marker
    D <- as.matrix(D[,-1],rownames=D$ID_REF)
    
    
    ## Read platform data table
    ref <- read.table(PLATFORM_FILE, header=TRUE, sep="\t", quote="", comment.char="#", fill=TRUE)[,c("ID",GENE_COLUMN)]
    colnames(ref) <- c("ID_REF","gene")
    ref$gene <- as.character(ref$gene)
    ref <- subset(ref, !is.na(gene) & gene != "")
    
    gene_split <- strsplit(ref$gene, " /// ")
    ref <- data.frame(
        ID_REF = rep(ref$ID_REF, sapply(gene_split, length)),
        gene = unlist(gene_split)
    )

    ## Merge tables to map Gene genes ids
    ##m <- data.table(merge(ref, D, all=FALSE)[,-1])
    dim(m)
    gene <- ref$gene[match(rownames(D),as.character(ref$ID_REF))]
    
    ## Aggregate duplicate genes by median
    sum(duplicated(gene))
    ##m <- m[,lapply(.SD, median(na.rm=TRUE)), by=gene]
    ex <- do.call(rbind, tapply(1:nrow(D), gene, function(i) colMeans(D[i,,drop=FALSE])))
    head(ex)
        
    ## Transpose matrix, add sample column and set gene names as column names
    ##m <- data.table(samples, transpose(m[,-1]))

    ## drop empty or duplicated columns
    sel1 <- (apply(anno, 2, function(x) mean(x %in% c(NA,"NA","null",""))) < 1)
    sel2 <- !duplicated(t(anno))
    sel  <- which( sel1 & sel2 )
    ##anno1 <- anno[,..sel]
    anno1 <- anno[,sel,with=FALSE]
    dim(anno1)
        
    ## Write results to separate expression and annotation files
    dim(ex)
    dim(anno1)
    if(write.file) {
        cat("writing",EXPRESSION_OUTPUT_FILE,"\n")
        fwrite(data.table(gene=rownames(ex),ex), file=EXPRESSION_OUTPUT_FILE, sep=",")
        ##write.csv(ex, file=EXPRESSION_OUTPUT_FILE)
        cat("writing",ANNOTATION_OUTPUT_FILE,"\n")
        fwrite(anno1, file=ANNOTATION_OUTPUT_FILE, sep=",")
    }

    res <- list()
    res$values <- ex
    res$anno <- anno1
    res

}

if(0) {

    geo <- parse_geo_series_matrix(
        SERIES_FILE = "GSE44770_series_matrix.txt",
        PLATFORM_FILE = "platform.txt",
        GENE_COLUMN = "ORF",
        EXPRESSION_OUTPUT_FILE = "expression.csv",
        ANNOTATION_OUTPUT_FILE = "samples.csv",
        write.file=TRUE
    )

    names(geo)
    head(geo$values)
    min(geo$values, na.rm=TRUE)
    head(geo$anno)

    counts <- 2**(pmax(3 + 3*geo$values,0))
    hist( log2(counts[,1]), breaks=100)
    write.csv(counts, file="counts.csv")
    
}
