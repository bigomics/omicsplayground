##
## https://github.com/denalitherapeutics/archs4
##
##

##BiocManager::install("denalitherapeutics/archs4")
##BiocManager::install("GEOmetadb")
library("archs4")
library(dplyr)
library(GEOmetadb)

RDIR="../../R/"
FILES="../../lib/"
PGX.DIR="../../data/"
source("../../R/pgx-init.R")

archs4dir <- "~/.archs4data"
archs4dir <- "~/bigomics/data/archs4data"

if(0) {
    archs4_local_data_dir_create(archs4dir)
    cwd = getwd()
    cwd
    setwd(archs4dir)
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5")
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_hiseq_transcript_v2.h5")
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5")
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_hiseq_transcript_v2.h5")
    system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz")
    system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz")
    
    create_augmented_feature_info(datadir=archs4dir)
    archs4_local_data_dir_validate(datadir=archs4dir)
    setwd(cwd)
    cwd
}


library(archs4)
a4 <- Archs4Repository(datadir=archs4dir)

## Get titles from GEO experiments
if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
rs <- dbGetQuery(con,paste("select gse,title from gse where",
                           "contributor like '%Ivo%Kwee%'",sep=" "))
rs
rs <- dbGetQuery(con,'select gse,title from gse')
GSE.TITLE <- rs$title
names(GSE.TITLE) <- rs$gse

## EXAMPLES
##

## all series ID
samplesize <- table(a4$sample_table$series_id)
all_ids = names(samplesize[ samplesize>=20 & samplesize<=200 ])
length(all_ids)
table(all_ids %in% names(GSE.TITLE))

## Select studies with relevant terms
ids <- all_ids[grep("cancer|onco|immun|aging|senesc",GSE.TITLE[all_ids],ignore.case=TRUE)]
length(ids)
##ids <- head(all_ids,100)
head(ids)
GSE.TITLE[ids]
##cc <- sample_covariates(a4)$name
cc3 = c("Sample_title","Sample_source_name_ch1","Sample_characteristics_ch1")

id = "GSE114716"
id = "GSE53784"
id = sample(ids,1)
id
length(ids)
for(id in ids) {

    pgx.file <- file.path(archs4dir,paste0(id,".pgx"))
    pgx.file
    if(file.exists(pgx.file)) {
        cat("skipping already done GEO series",id,"...\n")
        next()
    }

    cat("retrieving Archs4 data for GEO series",id,":",GSE.TITLE[id],"\n")
    aa <- pgx.getArchs4Dataset(a4, id) 
    names(aa)
    if("group" %in% colnames(aa$samples)) {
        colnames(aa$samples) <- sub("group","xgroup",colnames(aa$samples))
    }

    ## get categorical phenotypes
    df <- aa$samples
    vv <- pgx.getCategoricalPhenotypes(df, max.ncat=10, min.ncat=2) 
    cat("found discrete phenotypes:",vv,"\n")
    df <- df[,vv,drop=FALSE]
    apply(df,2,table)

    ## create contrast matrix
    res <- pgx.makeAutoContrast(df, mingrp=3, slen=15, ref=NA)
    names(res)
    head(res$contr.matrix)
    aa$samples <- cbind(aa$samples, group=res$group)
    aa$contrasts <- res$contr.matrix
    head(aa$contrasts)
    
    dim(aa$genes)
    dim(aa$counts)
    dim(aa$samples)
    dim(aa$contrasts)

    if(ncol(aa$contrasts)>0) {
    
        ## Playground pre-computation
        ngs <- pgx.upload(
            aa$counts, aa$samples, aa$contrasts,
            max.genes=5000,
            gx.methods = c("trend.limma","edger.qlf","edger.lrt"),
            gset.methods = c("fisher","gsva","fgsea"),
            ##extra.methods = c("meta.go","deconv","infer","drugs"),
            extra.methods = c("meta.go","infer","drugs"),
            lib.dir = "../../lib",
            progress=NULL)
        
        names(ngs)
        
        ## save
        ngs$name <- id
        ngs$description <- paste0("GEO Series ",id,". ",aa$title,".")
        ngs$datatype <- "RNA-seq"
        ngs$organism <- ngs.detectOrganism(ngs)
        ngs$date <- Sys.Date()

        cat("object size: ",format(object.size(ngs), units="MB"),"\n")
        ngs.save(ngs, file=pgx.file)

    }    
}

mingrp=3;slen=8;ref=NULL
pgx.makeAutoContrast <- function(df, mingrp=3, slen=8, ref=NULL) {

    shortestunique <- function(xx,slen=3) {
        k <- min(which(!sapply(1:max(nchar(xx)),function(i) any(duplicated(substring(xx,1,i))))))
        substring(xx,1,max(k,slen))
    }
    
    autoContrast1 <- function(x, ref, slen) {
        if(is.null(ref)) ref <- NA
        x <- as.character(x)
        nx <- table(x)
        too.small <- names(which(nx <= mingrp))
        if(length(too.small)) x[which(x %in% too.small)] <- NA
        nx <- table(x)
        x <- factor(x)
        if(!is.na(ref)) x <- relevel(x, ref=ref)
        xlevels <- gsub("[^[:alnum:]]","",levels(x))
        levels(x) <- shortestunique(xlevels,slen=slen)
        xref <- gsub("[^[:alnum:]]","",levels(x)[1])
        nn <- length(nx)
        nn
        if(nn<2) {
            return(NULL)
        } else if(nn == 2) {
            ct <- model.matrix(~x)[,2]
            ct <- matrix(ct, ncol=1)
            colnames(ct) <- paste0(levels(x)[2],"_vs_",levels(x)[1])
        } else if(nn >= 3) {
            if(is.na(ref)) {
                ct <- model.matrix(~0 + x)
                colnames(ct) <- paste0(levels(x),"_vs_rest")
            } else {
                ct <- model.matrix(~0 + x)
                colnames(ct) <- paste0(levels(x),"_vs_",xref)
                i=1
                for(i in 1:ncol(ct)) {
                    j <- which(!(x %in% levels(x)[c(1,i)]))
                    ct[j,i] <- NA
                }
                ct <- ct[,2:ncol(ct)] ## remove REFvsREF                
            }
        }
        ct
    }

    if(!is.null(ref) && length(ref)!=ncol(df)) ref <- head(rep(ref,99),ncol(df))
    
    K <- c()
    for(i in 1:ncol(df)) {
        refi <- NA
        if(!is.null(ref)) refi <- ref[i]
        x <- df[,i]
        if(!(refi %in% x)) refi <- NA
        ref.pattern <- "wt|contr|ctr|untreat|normal|^no$|neg|ref"
        detect.ref <- any(grepl(ref.pattern,x,ignore.case=TRUE))
        if(is.na(refi) & detect.ref) {
            refi <- grep(ref.pattern,x,ignore.case=TRUE,value=TRUE)[1]
            cat("automatic reference detected:",refi,"\n")
        }
        ct <- autoContrast1(x, ref=refi, slen=slen)
        colnames(ct) <- paste0(colnames(df)[i],":",colnames(ct))
        K <- cbind(K,ct)
    }
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


##sample_columns=cc3;feature_type="gene";row_id="symbol"
pgx.getArchs4Dataset <- function(a4, id)
{    
    ##ngs <- as.DGEList(a4, id, feature_type="gene", row_id="symbol")
    cc3 = c("Sample_title","Sample_source_name_ch1","Sample_characteristics_ch1")
    ngs <- my.as.DGEList(a4, id, sample_columns=cc3, feature_type="gene", row_id="symbol")
    is.null(ngs)
    if(is.null(ngs)) {
        warning("   failed to retrieve ",id,"\n")
        return(NULL)
    }

    names(ngs)
    lapply(ngs,dim)
    colnames(ngs$samples)

    cc3 <- intersect(colnames(ngs$samples), cc3)
    ngs$samples <- ngs$samples[,cc3,drop=FALSE]

    if("Sample_characteristics_ch1" %in% cc3) {
        ch1  <- ngs$samples[,"Sample_characteristics_ch1"]
        ch1x <- lapply(ch1, function(x) strsplit(x, split="Xx-xX")[[1]])
        pheno.values <- lapply(ch1x, function(x) sub(".*[:] ","",x))
        plen <- sapply(pheno.values,length)
        nf <- max(plen)
        nf
        pheno.names <- sub(":.*","",ch1x[[which.max(plen)]])
        pheno.names
        pheno.values <- t(sapply(pheno.values, function(s) head(c(s,rep(NA,nf)),nf)))        
        dim(pheno.values)
        colnames(pheno.values) <- gsub("[ /,;+&-]","_",pheno.names)
        head(pheno.values)
        title <- ngs$samples[,"Sample_title"]
        src   <- ngs$samples[,"Sample_source_name_ch1"]
        xpheno <- cbind(title=title, source=src, pheno.values)
        ngs$samples <- xpheno
    }

    rownames(ngs$samples) <- colnames(ngs$counts)
    ngs$samples <- data.frame(ngs$samples, stringsAsFactors=FALSE)
    ngs$counts <- as.matrix(ngs$counts)    
    ngs$title <- GSE.TITLE[id]

    return(ngs)
}

##features=NULL; sample_columns=c("Sample_title", "Sample_source_name_ch1"); feature_type="gene"; row_id="symbol";check_missing_samples=TRUE

my.as.DGEList <- function (a4, id, features = NULL,
                           sample_columns = c("Sample_title", "Sample_source_name_ch1"),
                           feature_type = c("gene", "transcript"), 
                           row_id = c("symbol", "ensembl"),
                           check_missing_samples = TRUE) 
{
    require(dplyr)
    ##assert_class(a4, "Archs4Repository")
    ##feature_type <- match.arg(feature_type)
    feature_type <- feature_type[1]
    ##row_id <- match.arg(row_id)
    row_id <- row_id[1]
    if (FALSE && !is.null(features)) {
        if (is.character(features)) {
            features <- feature_lookup(a4, features, type = type)
            if (any(is.na(features$a4name))) {
                warning("Removing 'not found' features from query")
                features <- filter(features, !is.na(a4name))
            }
        }
        assert_data_frame(features, min.rows = 1L)
        assert_subset(c("ensembl_id", "a4name"), colnames(features))
        features <- distinct(features, ensembl_id, .keep_all = TRUE)
    }

    si <- sample_info(a4, id, columns = sample_columns, check_missing_samples = check_missing_samples)
    nna <- sum(!is.na(si$sample_id))
    nna
    if(nna < 3) {
        warning("Could not retrieve enough samples or series not found")
        return(NULL)
    }
    
    si <- as.data.frame(si, stringsAsFactors = FALSE)
    si <- distinct(si, sample_id, .keep_all = TRUE)
    rownames(si) <- si$sample_id
    not.found <- filter(si, is.na(organism))
    not.found    
    if (nrow(not.found)) {
        ##stop("The following samples could not be found: ", paste(sprintf("%s::%s", 
        ##    not.found$series_id, not.found$sample_id), collapse = "; "))
        warning("The following samples could not be found: ",
                paste(sprintf("%s::%s",not.found$series_id, not.found$sample_id),
                      collapse = "; "),"\n")
    }
    
    org <- unique(si$organism)
    org
    if (length(org) != 1L) {
        stop("You are querying across species")
    }

    finfo <- feature_info(a4, feature_type, org, augmented = FALSE)
    ##finfo <- feature_info(a4, feature_type, org, augmented = TRUE)
    finfo <- as.data.frame(finfo, stringsAsFactors = FALSE)
    dim(finfo)
    head(finfo)
    
    if (feature_type == "gene") {
        ##cat("feature_type = gene\n")
        if (row_id == "ensembl") {
            finfo <- filter(finfo, !is.na(ensembl_id))
            rownames(finfo) <- finfo[["ensembl_id"]]
        }
        else {
            rownames(finfo) <- finfo[["a4name"]]
        }
    } else {
        ##cat("feature_type = transcript\n")
        dup.ensid <- duplicated(finfo[["ensembl_id"]])
        if (any(dup.ensid)) {
            warning("Duplicated ensembl identifiers when version is removed, ", 
                    "rownames maintain their versioned id")
            rownames(finfo) <- finfo[["ensembl_id_full"]]
        }
        else {
            rownames(finfo) <- finfo[["ensembl_id"]]
        }
    }

    if (!is.null(features)) {
        finfo <- subset(finfo, a4name %in% features[["a4name"]])
    }
    h5col <- paste0("sample_h5idx_", feature_type)
    isna <- is.na(si[[h5col]])
    isna
    if(any(isna)) {
        msg <- paste0("The following samples do not have ", feature_type, 
                      "-level quantitation and will not be included in the expression ", 
                      "container:\n  ", paste(sprintf("%s::%s", si$series_id[isna], 
                                                      si$sample_id[isna]), collapse = "; "))
        warning(msg, immediate. = TRUE)
        si <- si[!isna, , drop = FALSE]
    }
    
    if (nrow(si) == 0L) {
        ##stop("No samples left to assemble expression data")
        warning("No samples left to assemble expression data")
        return(NULL)
    }

    libinfo <- sample_table(a4) %>%
        select(series_id, sample_id, 
               ## lib.size = libsize, norm.factors = normfactor,
               a4libsize) %>% 
        distinct(sample_id, .keep_all = TRUE)
    libinfo <- select(si, series_id, sample_id) %>% left_join(libinfo, 
        by = c("series_id", "sample_id"))
    na.overflow <- is.na(libinfo$lib.size) | is.na(libinfo$norm.factors)
    if (any(na.overflow)) {
        warning("Removing ", sum(na.overflow), " samples due to libsize NA overflow issues")
        libinfo <- filter(libinfo, !na.overflow)
        si <- subset(si, sample_id %in% libinfo$sample_id)
    }
    
    rownames(si) <- si[["sample_id"]]
    counts <- local({
        key <- paste(org, feature_type, sep = "_")
        h5.fn <- file_path(a4, key)
        index <- list(finfo$h5idx, si[[h5col]])
        cnts <- rhdf5::h5read(h5.fn, "data/expression", index = index)
        colnames(cnts) <- rownames(si)
        rownames(cnts) <- rownames(finfo)
        cnts
    })
    out  <- suppressWarnings(edgeR::DGEList(counts, genes = finfo, samples = si))

    xref <- match(colnames(out), libinfo$sample_id)
    if (any(is.na(xref))) {
        stop("Problem matching sample_id to libinfo data.frame")
    }
    if (!all(colnames(out) == libinfo$sample_id[xref])) {
        stop("Mismatch in outgoing DGEList to libinfo data.frame")
    }
    ##out$samples$lib.size <- libinfo$lib.size[xref]
    ##out$samples$norm.factors <- libinfo$norm.factors[xref]
    out$samples$lib.size <- colSums(out$counts)
    out$samples$norm.factors <- rep(1,ncol(out$counts))

    out
}
