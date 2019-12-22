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
source("../../R/pgx-archs4.R")
source("../../R/ngs-functions.R")
source("../../R/pgx-contrasts.R")

archs4dir <- "~/.archs4data"
##archs4dir <- "~/bigomics/data/archs4data"
archs4dir <- "/data/Projects/Data/archs4data"

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
metadb <- file.path(archs4dir,'GEOmetadb.sqlite')
if(!file.exists(metadb)) getSQLiteFile(destdir=archs4dir)
con <- dbConnect(SQLite(),metadb)
rs <- dbGetQuery(con,paste("select gse,title from gse where",
                           "contributor like '%Ivo%Kwee%'",sep=" "))
rs

rs <- dbGetQuery(con,'select gse,title from gse')
GSE.TITLE <- rs$title
names(GSE.TITLE) <- rs$gse

if(0) {
    D <- read.csv("../../lib/L1000_repurposing_drugs.txt",
                  sep="\t",comment.char="#")
    drugs = as.character(D$pert_iname)
    drugs1 = grep("[[:punct:]]",drugs,value=TRUE,invert=TRUE)
    about.drug <- mclapply(GSE.TITLE[], function(a) any(sapply(drugs1,function(d)grepl(d,a))))
    about.drug <- unlist(about.drug)
    names(about.drug) <- GSE.TITLE
    table(about.drug)
}

##================================================================================
##================================= MAIN =========================================
##================================================================================

## all series ID
samplesize <- table(a4$sample_table$series_id)
all_ids = names(samplesize[ samplesize>=20 & samplesize<=200 ])
all_ids = names(samplesize[ samplesize>=10 & samplesize<=500 ])
length(all_ids)
table(all_ids %in% names(GSE.TITLE))

## Select studies with relevant terms
ids <- all_ids[grep("lymphom|leukaem|hemato",GSE.TITLE[all_ids],ignore.case=TRUE)]
ids <- all_ids[grep("prostate.cancer",GSE.TITLE[all_ids],ignore.case=TRUE)]
ids <- all_ids[grep("breast.cancer",GSE.TITLE[all_ids],ignore.case=TRUE)]
ids <- all_ids[grep("cancer|onco|tumor|tumour",GSE.TITLE[all_ids],ignore.case=TRUE)]
length(ids)

## Select studies with relevant terms about AGING
ids1 <- all_ids[grep("[ ]aging|^aging|senesc",GSE.TITLE[all_ids],ignore.case=TRUE)]
length(ids1)

## Select studies with relevant terms about Immuno AND oncology terms
ids <- all_ids[grep("cancer|onco|tumor|tumour",GSE.TITLE[all_ids],ignore.case=TRUE)]
ids2 <- ids[grep("immun",GSE.TITLE[ids],ignore.case=TRUE)]
length(ids2)

ids <- c(ids1, ids2)
##ids <- head(ids,20)
length(ids)
head(ids)
GSE.TITLE[ids]
##cc <- sample_covariates(a4)$name


id = "GSE53784"
id = ids[1]
id
length(ids)
res <- mclapply( ids[1:4], function(id) prepArchs4Dataset(id) )
res



prepArchs4Dataset <- function(id) {
    
    pgx.file <- file.path(archs4dir,paste0(id,".pgx"))
    pgx.file
    if(file.exists(pgx.file)) {
        cat("skipping already done GEO series",id,"...\n")
        return("already done")
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
    vv
    if(length(vv)==0) {
        cat("skipping. no valid phenotypes in GEO series",id,"...\n")
        return("skipped. no valid phenotypes")
    }
    cat("found discrete phenotypes:",vv,"\n")
    df <- df[,vv,drop=FALSE]
    apply(df,2,table)
    
    ## create contrast matrix
    mingrp=3;slen=15;ref=NA
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

    has.contrast <- (!is.null(aa$contrasts) && ncol(aa$contrasts)>0)
    has.contrast    
    if(!has.contrast) {
        return("skipped. no valid contrasts")
    } else {
    
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
        ngs$description <- paste0("GEO Series ",id,". ",GSE.TITLE[id],".")
        ngs$datatype <- "RNA-seq"
        ngs$organism <- ngs.detectOrganism(ngs)
        ngs$date <- Sys.Date()
        ngs$link <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",id)
        
        cat("object size: ",format(object.size(ngs), units="MB"),"\n")
        ngs.save(ngs, file=pgx.file)
        return("OK")
    }    
}
