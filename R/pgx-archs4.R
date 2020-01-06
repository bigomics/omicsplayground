##
## https://github.com/denalitherapeutics/archs4
##
##

ARCHS4.DIR = "~/.archs4data"

if(0) {

    ##BiocManager::install("denalitherapeutics/archs4")
    ##BiocManager::install("GEOmetadb")
    library("archs4")
    library(dplyr)
    library(GEOmetadb)
    
    ARCHS4.DIR <- "~/.archs4data"
    archs4_local_data_dir_create(ARCHS4.DIR)
    cwd = getwd()
    cwd
    setwd(ARCHS4.DIR)
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5")
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/human_hiseq_transcript_v2.h5")
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_matrix.h5")
    system("wget https://s3.amazonaws.com/mssm-seq-matrix/mouse_hiseq_transcript_v2.h5")
    system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz")
    system("wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.gtf.gz")
    
    create_augmented_feature_info(datadir=ARCHS4.DIR)
    archs4_local_data_dir_validate(datadir=ARCHS4.DIR)
    setwd(cwd)
    cwd
}

parseGeoAnnot <- function(annot) {
    parse.annot <- function(a) {
        aa <- strsplit(a, split="Xx-xX")[[1]]
        aa <- lapply(aa, function(x) strsplit(x,split=": ")[[1]])
        array( sapply(aa, "[",2), dimnames=list(sapply(aa,"[",1)))
    }        
    lapply(annot[], parse.annot)
}

keyword="cancer"
geo.selectSeries <- function(keyword, variables=NULL)
{
    metadb <- file.path(ARCHS4.DIR,'GEOmetadb.sqlite')
    require(GEOmetadb)        
    ## Get titles from GEO experiments
    ##if(!file.exists(metadb)) getSQLiteFile(destdir=ARCHS4.DIR)
    con <- dbConnect(SQLite(),metadb)
    dbListTables(con)
    dbListFields(con,'gse')        
    series.str <- paste0("('",paste(unique(series),collapse="','"),"')")
    ##rs <- dbGetQuery(con,paste("SELECT gse,title FROM gse WHERE gse IN",series.str))
    ww <- strsplit(keyword,split="\\|")[[1]]
    gse <- c()
    for(w in ww) {
        keyword.str <- paste0("'",tolower(w),"')>0")
        rs1 <- dbGetQuery(con,paste("SELECT gse,title FROM gse WHERE INSTR(LOWER(title),",
                                    keyword.str))
        rs2 <- dbGetQuery(con,paste("SELECT gse,summary FROM gse WHERE INSTR(LOWER(summary),",
                                    keyword.str))
        gse <- unique(c(gse,unlist(rs1$gse),unlist(rs2$gse)))
    }
    length(gse)
    
    if(!is.null(variables)) {
        matrix_file = file.path(ARCHS4.DIR, "human_matrix.h5")
        ##extracted_expression_file = "/tmp/example_expression_matrix.tsv"
        file.exists(matrix_file)
        ## Retrieve information from compressed data
        annot0  = h5read(matrix_file, "meta/Sample_characteristics_ch1")
        series0 = h5read(matrix_file, "meta/Sample_series_id")        
        sel <- grep(variables, annot0)
        head(annot0[sel])
        sel <- sel[which(!duplicated(series0[sel]))]
        gse3   <- series0[sel]
        annot1 <- parseGeoAnnot(annot0[sel])
        annot1 <- sapply(annot1, function(s) paste(names(s),collapse="|"))
        names(annot1) <- gse3
        length(gse3)
        gse <- intersect(gse, gse3)
    }

    series.str <- paste0("('",paste(unique(gse),collapse="','"),"')")    
    rt <- dbGetQuery(con,paste("SELECT gse,title,summary FROM gse WHERE gse IN ",series.str))
    rt$variables <- annot1[rt$gse]
    return(rt)
}

res <- geo.selectSeries("immune","surv")
res <- geo.selectSeries("cancer|lymphom","surv")

gse.series = c("GSE76514","GSE81475")
genes = c('CD101','CD109','CD14','CD151','CD160','CD163','CD163L1','CD164','CD164L2','CD177')
studies = c("prad_tcga","ov_tcga")

pgx.getArchs4expression <- function(gse.series, genes=NULL, batch.correct=TRUE)
{
    ## Adapted from: https://amp.pharm.mssm.edu/archs4/help.html
    ##
    ##
    ## R script to download selected samples
    ## Copy code and run on a local machine to initiate download
    ## Check for dependencies and install if missing

    ## packages <- c("rhdf5", "preprocessCore", "sva")
    ## if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
    ##     print("Install required packages")
    ##     source("https://bioconductor.org/biocLite.R")
    ##     biocLite("rhdf5")
    ##     biocLite("preprocessCore")
    ##     biocLite("sva")
    ## }
    library("rhdf5")
    library("preprocessCore")
    library("sva")
    
    matrix_file = file.path(ARCHS4.DIR, "tcga_matrix.h5")
    matrix_file = file.path(ARCHS4.DIR, "mouse_matrix.h5")
    matrix_file = file.path(ARCHS4.DIR, "human_matrix.h5")
    ##extracted_expression_file = "/tmp/example_expression_matrix.tsv"
    file.exists(matrix_file)
    h5ls(matrix_file)[,1:2]
    
    ## Check if gene expression file was already downloaded, if not in current directory download file form repository
    if(!file.exists(matrix_file)){
        ##print("Downloading compressed gene expression matrix.")
        ##url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"
        ##download.file(url, matrix_file, quiet = FALSE)
        stop("pgx.getArchs4Data:: missing file = ",matrix_file)
    } 
    
    ## Retrieve information from compressed data
    samples0 = h5read(matrix_file, "meta/Sample_geo_accession")
    tissue0  = h5read(matrix_file, "meta/Sample_source_name_ch1")
    genes0   = h5read(matrix_file, "meta/genes")
    series0  = h5read(matrix_file, "meta/Sample_series_id")
    sampletitle0 = h5read(matrix_file, "meta/Sample_title")
    
    if(is.null(genes)) {
        genes <- genes0
    } else {
        genes <- intersect(genes, genes0)
    }
    
    ## Identify columns to be extracted
    ##sample_locations = which(samples %in% samp)
    sample_index = which(series0 %in% gse.series)
    length(sample_index)
    series = series0[sample_index]
    samples = samples0[sample_index]
    gene_index <- match(genes, genes0)
    tissue <- tissue0[sample_index]
    sample_title <- sampletitle0[sample_index]

    ##annot = h5read(matrix_file, "meta/Sample_characteristics_ch1")
    annot = h5read(matrix_file, "meta/Sample_characteristics_ch1",
                   index = list(sample_index))
    length(annot)
    annot <- parseGeoAnnot(annot)    
    names(annot) <- samples
    
    gse.info <- NULL
    metadb <- file.path(ARCHS4.DIR,'GEOmetadb.sqlite')
    if(file.exists(metadb)) {
        require(GEOmetadb)        
        ## Get titles from GEO experiments
        ##if(!file.exists(metadb)) getSQLiteFile(destdir=ARCHS4.DIR)
        con <- dbConnect(SQLite(),metadb)
        dbListTables(con)
        dbListFields(con,'gse')        
        series.str <- paste0("('",paste(unique(series),collapse="','"),"')")
        ##rs <- dbGetQuery(con,paste("SELECT gse,title FROM gse WHERE gse IN",series.str))
        gse.info <- dbGetQuery(con,paste("SELECT gse,title,pubmed_id,summary FROM gse WHERE gse IN",series.str))
        ## rs <- dbGetQuery(con,'SELECT ID,gse,title FROM gse')
    }

    ## extract gene expression from compressed data
    expression = h5read(
        matrix_file, "data/expression",
        index = list(gene_index, sample_index)
    )
    H5close()
    dim(expression)
    rownames(expression) = genes
    ##colnames(expression) = samples
    colnames(expression) = paste0("[",series,"] ",samples)
    
    ## normalize samples and correct for differences in gene count distribution
    expression = log2(expression+1)
    ##expression = normalize.quantiles(expression)
    expression = limma::normalizeQuantiles(expression)    
    
    ## correct batch effects in gene expression
    batch.correct = TRUE
    if(batch.correct && length(gse.series)>1 ) {
        batchid = match(series, unique(series))
        table(batchid)
        suppressMessages(
            expression <- ComBat(dat=expression, batch=batchid, par.prior=TRUE, prior.plots=FALSE)
        )
    }    
    
    res <- list()
    res$expression <- expression
    res$characteristics <- annot
    res$tissue <- tissue
    ##res$series <- series
    res$series.info <- gse.info
    ## res$samples <- samples
    
    return(res)
}

id = "GSE29943"
pgx.getArchs4Dataset <- function(id, datadir = ARCHS4.DIR, a4 = NULL)
{    
    library(archs4)
    if(is.null(a4)) a4 <- Archs4Repository(datadir=ARCHS4.DIR)  ## takes long time...
    
    ##feature_type="gene";row_id="symbol"
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
    cc3
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
        if(length(pheno.names)==1) {
            pheno.values <- t(pheno.values)
        }
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
    ngs$title <- paste("GEO Series",id)

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

    ## libinfo <- sample_table(a4) %>%
    ##     select(series_id, sample_id, 
    ##            ## lib.size = libsize, norm.factors = normfactor,
    ##            a4libsize) %>% 
    ##     distinct(sample_id, .keep_all = TRUE)
    ## libinfo <- select(si, series_id, sample_id) %>% left_join(libinfo, 
    ##     by = c("series_id", "sample_id"))

    ## na.overflow <- is.na(libinfo$lib.size) | is.na(libinfo$norm.factors)
    ## if (any(na.overflow)) {
    ##     warning("Removing ", sum(na.overflow), " samples due to libsize NA overflow issues")
    ##     libinfo <- filter(libinfo, !na.overflow)
    ##     si <- subset(si, sample_id %in% libinfo$sample_id)
    ## }
    
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

    ## xref <- match(colnames(out), libinfo$sample_id)
    ## if (any(is.na(xref))) {
    ##     stop("Problem matching sample_id to libinfo data.frame")
    ## }
    ## if (!all(colnames(out) == libinfo$sample_id[xref])) {
    ##     stop("Mismatch in outgoing DGEList to libinfo data.frame")
    ## }
    ##out$samples$lib.size <- libinfo$lib.size[xref]
    ##out$samples$norm.factors <- libinfo$norm.factors[xref]
    out$samples$lib.size <- colSums(out$counts)
    out$samples$norm.factors <- rep(1,ncol(out$counts))  ## HACK ????

    out
}
