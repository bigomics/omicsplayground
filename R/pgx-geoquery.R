
if(0) {
    ## BiocManager::install("GEOmetadb")    
    library(GEOmetadb)
    if(!file.exists('GEOmetadb.sqlite')) getSQLiteFile()
    
    con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
    ##dbDisconnect(con)
    geo_tables <- dbListTables(con)
    geo_tables
    
    rs <- dbGetQuery(con,paste("select gse,title from gse where",
                               "contributor like '%Ivo%Kwee%'",sep=" "))
    rs
    
    platforms = c(
        "GPL570", ## Affymetrix HG-U133plus2
        "GPL6244", ## Affymetrix HuGene 1.0st
        "GPL10558", ## Illumina HT-12 v4.0
        "GPL11154", ## Illumina HiSeq 2000
        "GPL16791", ## Illumina HiSeq 2500
        "GPL20401" ## Illumina HiSeq 4000
    )
    
    ## Select GSE datasets
    library(dplyr)
    db = src_sqlite('GEOmetadb.sqlite')
    gse = data.frame(tbl(db,'gse'))
    gpl = data.frame(tbl(db,'gse_gpl'))
    gsm = data.frame(tbl(db,'gse_gsm'))
    gsm0 = data.frame(tbl(db,'gsm'))
    nsamples = table(gsm$gse)
    nsamples <- nsamples[match(gpl$gse,names(nsamples))]
    
    s1 <- ( gpl$gpl %in% platforms )
    s2 <- ( nsamples >= 20 & nsamples <= 200 )
    sel.gse <- gpl$gse[ which(s1 & s2) ]
    length(sel.gse)    
    gse1 <- gse[ which(gse$gse %in% sel.gse),]
    head(gse1)[,1:3]

    library(GEOquery)
    id = gse1$gse[1]
    id = "GSE21653"  ## BRCA
    id = "GSE138126" ## 
    id = "GSE10846"  ## DLBCL
    id = "GSE108467"
    id = "GSE4824"
    ##id = "GDS507"
    id = "GSE53784" ## WINJEV
    
    ## examples
    id = "GSE10846"
    geo <- NULL
    geo <- pgx.GEOquery(id)
    is.null(geo$error)
    if(!is.null(geo$error)) {
        write(geo$error, file=file.path(SAVEDIR,paste0("geo-",id,".error")))
    } else {
        lapply(geo, dim)
        SAVEDIR = "."
        saveRDS(geo, file.path(SAVEDIR,paste0("geo-",id,".rds")))
    }
    geo <- pgx.GEOquery("ABC")

    
}

##source(file.path(RDIR,"ngs-functions.R"))
##source(file.path(RDIR,"pgx-functions.R"))

##-------------------------------------------------------------------------------------
## Query GEO
##-------------------------------------------------------------------------------------

GPL.PLATFORMS = c(
    "GPL570", ## Affymetrix HG-U133plus2
    "GPL6244", ## Affymetrix HuGene 1.0st
    "GPL10558", ## Illumina HT-12 v4.0
    "GPL11154", ## Illumina HiSeq 2000
    "GPL16791", ## Illumina HiSeq 2500
    "GPL20401" ## Illumina HiSeq 4000
)
platforms = GPL.PLATFORMS

pgx.GEOquery <- function(id, platforms=GPL.PLATFORMS) {
    ## Retrieve expression matrix, phenotype and probe annotation
    ## matrices for a certain GEO id.
    ##
    require(GEOquery)
    
    ## load series and platform data from GEO
    id
    ##geo <- getGEO(id, GSEMatrix=FALSE, getGPL=FALSE)
    geo <- try( getGEO(id, GSEMatrix=TRUE, getGPL=FALSE) )
    class(geo)
    if(class(geo)=="try-error") {
        res <- list(error="ERROR: getGEO() error")
        return(res)
    }
    length(geo)
    attr(geo, "names")
    
    idx = 1
    if (length(geo) > 1) {
        geo.gpl <- gsub(".*GPL","GPL",attr(geo, "names"))
        geo.gpl <- gsub("_series.*","",geo.gpl)
        idx <- which(geo.gpl %in% platforms)[1]
    }
    eset <- geo[[idx]]
    
    ## perform log2 transform if required
    ex <- exprs(eset)
    dim(ex)
    if(nrow(ex)==0) {
        res <- list(error="ERROR: No counts data in GEO object")
        return(res)
    }

    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    qx
    not.log <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
    not.log
    if (not.log) {
        ex[which(ex <= 0 | is.na(ex))] <- 0
        exprs(eset) <- log2(1 + ex)
    }
    
    ## get phenotype data
    pdata = pData(eset)
    colnames(pdata)
    head(pdata$title)  ## sample names
    pheno = pdata[,grep(":ch1$",colnames(pdata)),drop=FALSE]
    colnames(pheno) <- sub(":ch1$","",colnames(pheno))
    colnames(pheno)
    dim(pheno)
    if(ncol(pheno)==0) {
        res <- list(error="ERROR: No CH1 phenotype data")
        return(res)
    }

    ## unpack clinical info
    if("Clinical info" %in% colnames(pheno)) {
        clin.info <- pheno[,"Clinical info"]
        clin.terms <- sub(":.*$","",strsplit(clin.info[1], split=";")[[1]])
        pheno.clin <- t(sapply(clin.info, function(s) (strsplit(s,split=";")[[1]]) ))
        pheno.clin <- apply(pheno.clin,2,function(s) sub(".*: ","",s))
        colnames(pheno.clin) <- clin.terms
        rownames(pheno.clin) <- NULL
        pheno$"Clinical info" <- NULL
        pheno <- data.frame(pheno, pheno.clin, stringsAsFactors=FALSE, check.names=FALSE)
    }
    pheno <- data.frame(pheno, stringsAsFactors=FALSE, check.names=FALSE)
    colnames(pheno) <- gsub("[ ]","_",colnames(pheno))
    dim(pheno)
    
    ## get platform annotation
    gpl.id <- eset@annotation
    gpl.id
    gpl <- getGEO(gpl.id, destdir="/tmp")
    gpl.title <- Meta(gpl)$title
    gpl.title
    colnames(Table(gpl))
    gpl.table <- Table(gpl)
    if(any(grepl("gene.symbol|gene.title",colnames(gpl.table),ignore.case=TRUE))) {
        gpl.table <- gpl.table[,grep("^id$|gene.symbol|gene.title|gene.name",colnames(gpl.table),ignore.case=TRUE)]
        gpl.table <- gpl.table[match(rownames(eset),gpl.table$ID),,drop=FALSE]
    } else if ("mrna_assignment" %in% colnames(gpl.table)) {
        
    }
    head(gpl.table)
    rownames(gpl.table) <- gpl.table[,"ID"]
    
    ## make contrasts
    sel <- pgx.getCategoricalPhenotypes(pheno, max.ncat=5, min.ncat=2) 
    sel
    ##pheno = pheno[,sel]
    if(length(sel)==0) {
        res <- list(error="ERROR: No valid phenotypes variables")
        return(res)
    }    
    Y = pheno[,sel,drop=FALSE]
    Y = apply(Y,2,function(s) gsub("[ ]","_",s))
    ref = rep(NA,length(sel))
    ref = apply(Y,2,function(x) sort(setdiff(unique(x),c(NA,"NA")))[1] )
    nlev = apply(Y,2,function(x) length(setdiff(unique(x),c(NA,"NA"))))
    ref[nlev>=4] <- "all"
    ##ct0 <- expandPhenoMatrix(Y)
    ct1 <- sign(makeDirectContrasts(Y, ref=ref))
    ct.min <- apply(ct1,2,function(x) min(table(x[ x!=0 & !is.na(x)])))
    ct1 <- ct1[,ct.min>=3,drop=FALSE]
    dim(ct1)
    if(ncol(ct1)==0) {
        res <- list(error="ERROR: Empty contrast matrix")
        return(res)
    }    
    
    ## prepare object
    counts  = 2**exprs(eset)
    samples = pheno
    genes   = gpl.table

    ## make sure matrices are aligned
    ct1 <- ct1[match(colnames(counts),rownames(ct1)),,drop=FALSE]
    genes <- genes[match(rownames(counts),rownames(genes)),,drop=FALSE]
    
    res <- list(
        counts = counts,
        samples = pheno,
        genes = gpl.table,
        contrasts = ct1,
        info = c(geo=id, platform.id=gpl.id, platform=gpl.title),
        error = NULL
    )

    return(res)
}
