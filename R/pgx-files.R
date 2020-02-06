
##h5.file="test.h5";chunk=100
pgx.saveMatrixH5 <- function(X, h5.file, chunk=NULL, del=TRUE )
{   
    require(rhdf5)
    if(del) unlink(h5.file)
    
    if(is.null(chunk)) {
        if(del) h5createFile(h5.file)    
        ## h5createGroup("myhdf5file.h5","foo")
        ## A = matrix(1:10,nr=5,nc=2)
        ## h5write(A, "myhdf5file.h5","foo/A")    
        h5createGroup(h5.file,"data")    
        h5write( X, h5.file, "data/matrix")
    } else {
        if(del) h5createFile(h5.file)    
        if(del) h5createGroup(h5.file,"data")
        h5createDataset(
            h5.file, "data/matrix",
            c(nrow(X),ncol(X)),
            ##storage.mode = "integer",
            chunk = chunk,
            level = 7
        )
        h5write(
            X,
            file = h5.file,
            name = "data/matrix",
            index = list(1:nrow(X),1:ncol(X))
        )
    }

    h5write( rownames(X), h5.file, "data/rownames")
    h5write( colnames(X), h5.file, "data/colnames")    
    
    h5closeAll()
}


file = "./OPTIONS_"
pgx.readOptions <- function(file = "./OPTIONS") {
    if(!file.exists(file)) return(NULL)
    opt <- read.table(file, sep="=", row.names=1)
    opt <- gsub("^[ ]*|[ ]*$","",apply(opt,1,c))
    opt <- sapply(opt,list)
    opt <- sapply(opt,strsplit,split=";")
    names(opt) <- gsub("^[ ]*|[ ]*$","",names(opt))
    opt
}

##pgx.dir=PGX.DIR
pgx.initDatasetFolder <- function(pgx.dir, verbose=TRUE, force=FALSE)
{
    i=1
    for(i in 1:length(pgx.dir)) {
        ## pgx.initDatasetFolder(pgx.dir[i], verbose=verbose, force=force)
        if(!dir.exists(pgx.dir[i])) next()

        ## skip if no pgx files
        npgx <- length(dir(pgx.dir[i],"pgx$"))
        if(npgx==0) {
            cat("no pgx files in",pgx.dir[i],"\n")
            next()
        }

        ## skip if too many...
        if(npgx > 100) {
            cat("too many files in",pgx.dir[i],". Please init manually\n")
            next()
        }
        
        info <- pgx.initDatasetFolder1 (
            pgx.dir[i],
            allfc.file = "datasets-allFC.csv",
            info.file = "datasets-info.csv",
            force = force,
            verbose = verbose)    

    }
}

if(0) {
    ##h5ls(h5.file)
    pgx.createSignatureDatabaseH5( pgx.files, h5.file, update.only=FALSE)
    h5ls(h5.file)
    pgx.addEnrichmentSignaturesH5(h5.file, X=NULL, mc.cores=8,
                                  lib.dir=FILES, methods=c("gsea"))
    h5ls(h5.file)
}


verbose=TRUE;file="datasets-allFC.csv"
pgx.readDatasetProfiles <- function(pgx.dir, file="datasets-allFC.csv",
                                    verbose=TRUE)
{
    F <- NULL
    pgx.dir <- unique(pgx.dir)
    i=1
    for(i in 1:length(pgx.dir)) {
        if(!dir.exists(pgx.dir[i])) next()
        f1 <- pgx.readDatasetProfiles1(pgx.dir[i], file=file,
                                       verbose=verbose)
        if(is.null(F)) {
            F <- f1
        } else {
            gg <- sort(unique(rownames(F),rownames(f1)))
            gg <- setdiff(gg, c(NA,""))
            F <- cbind( F[match(gg,rownames(F)),,drop=FALSE],
                       f1[match(gg,rownames(f1)),,drop=FALSE] )
            rownames(F) <- gg
        }
    }
    dim(F)
    return(F)
}

##pgx.dir=PGX.DIR;file="datasets-allFC.csv"
pgx.readDatasetProfiles1 <- function(pgx.dir, file="datasets-allFC.csv",
                                     verbose=TRUE)
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.readDatasetProfiles1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    fn <- file.path(pgx.dir,file)
    fn
    if(!file.exists(fn)) {
        stop("FATAL : could not find profiles matrix. please create first with initDatasetFolder().\n")
        ## pgx.updateDatasetProfiles(pgx.dir, file=file)
        return()
    } else {
        if(verbose) cat("Found existing dataset profiles matrix\n")
    }
    require(data.table)
    allFC <- fread.csv(file=file.path(pgx.dir, file),
                           row.names=1, check.names=FALSE)
    allFC <- as.matrix(allFC)
    if(verbose) cat("dataset profiles matrix : dim=",dim(allFC),"\n")
    dim(allFC)    
    return(allFC)
}


file="datasets-info.csv"
pgx.scanInfoFile <- function(pgx.dir, file="datasets-info.csv", verbose=TRUE)
{
    pgxinfo <- NULL
    i=1
    for(i in 1:length(pgx.dir)) {
        
        pgxinfo.file <- file.path(pgx.dir[i], file)
        if(!file.exists(pgxinfo.file)) next()  ## really?? no updating??
        
        info = fread.csv(pgxinfo.file, stringsAsFactors=FALSE, row.names=1)
        dim(info)       
        info$path <- pgx.dir[i]
        if(is.null(pgxinfo)) {
            pgxinfo <- info
        } else {
            jj <- match(colnames(pgxinfo),colnames(info))
            pgxinfo <- rbind(pgxinfo, info[,jj])
        }
    }
    dim(pgxinfo)
    return(pgxinfo)
}


##pgx.dir=PGX.DIR[1];allfc.file="datasets-allFC.csv";verbose=1;info.file="datasets-info.csv";force=1
pgx.initDatasetFolder1 <- function( pgx.dir,
                                   allfc.file = "datasets-allFC.csv",
                                   info.file = "datasets-info.csv",
                                   force=FALSE, verbose=TRUE)
{

    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.updateDatasetsMetaFiles1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    
    ## all public datasets
    pgx.dir <- pgx.dir[1]  ## only one folder!!!
    pgx.files <- dir(pgx.dir, pattern="[.]pgx$")
    pgx.files

    ##----------------------------------------------------------------------
    ## If an allFC file exists
    ##----------------------------------------------------------------------

    allfc.file1 <- file.path(pgx.dir, allfc.file)
    has.fc <- file.exists(allfc.file1)
    if(verbose && has.fc) cat("file",allfc.file1,"exists: YES\n")
    if(verbose && !has.fc) cat("file",allfc.file1,"exists: NO\n")

    info.file1 <- file.path(pgx.dir, info.file)
    has.info <- file.exists(info.file1)
    if(verbose && has.info) cat("file",info.file1,"exists: YES\n")
    if(verbose && !has.info) cat("file",info.file1,"exists: NO\n")
    
    ##----------------------------------------------------------------------
    ## If an allFC file exits, check if it is done for all PGX files
    ##----------------------------------------------------------------------

    pgxinfo <- NULL
    pgx.missing0 <- pgx.files
    pgx.missing1 <- pgx.files
    force

    allFC <-NULL
    if(!force && file.exists(allfc.file1)) {
        if(verbose) cat("checking which pgx files already done in allFC...\n")
        allFC <- read.csv(allfc.file1,row.names=1,check.names=FALSE,nrow=5)  ## just HEADER!!!
        dim(allFC)
        files.done <- gsub("\\[|\\].*","",colnames(allFC))
        files.done <- unique(paste0(files.done,".pgx"))
        pgx.missing0 <- setdiff(pgx.missing0, files.done)
    }
    dim(allFC)
    
    if(!force && has.info) {
        if(verbose) cat("checking which pgx files already in PGX info...\n")        
        pgxinfo = fread.csv(info.file1, stringsAsFactors=FALSE, row.names=1)
        dim(pgxinfo)       
        jj <- which(!(sub(".pgx$","",pgx.missing1) %in% sub(".pgx$","",pgxinfo$dataset)))
        pgx.missing1 = pgx.missing1[jj]
    }
   
    ##----------------------------------------------------------------------
    ## Check if it is done for all PGX files
    ##----------------------------------------------------------------------
    
    ## files to be done either for allFC or missing in INFO
    pgx.missing <- unique(c(pgx.missing0, pgx.missing1))
    length(pgx.missing)
    if(verbose) cat("FORCE =",force,"\n")
    
    if(length(pgx.missing)==0) {
        if(verbose) cat("no update required. use FORCE=1 for forced update.\n")
        return(NULL)
    }
    if(verbose) cat("scanning",length(pgx.missing),"PGX files in folder",pgx.dir,"\n")
        
    ##----------------------------------------------------------------------
    ## Reread allFC file. Before we only read the header.
    ##----------------------------------------------------------------------

    allFC <-NULL
    if(!force && file.exists(allfc.file1)) {
        allFC <- read.csv(allfc.file1,row.names=1,check.names=FALSE)  
    }
    dim(allFC)

    ##----------------------------------------------------------------------
    ## For all new PGX files, load the PGX file and get the meta FC
    ## matrix.
    ##----------------------------------------------------------------------

    info.cols <- NULL
    missing.FC <- list()
    pgx = pgx.missing[1]
    
    for(pgx in pgx.missing) {

        pgx
        if(verbose) cat(".")        
        try.error <- try( load( file.path(pgx.dir,pgx),verbose=0) )
        if(class(try.error)=="try-error") {
            cat(paste("error in loading PGX file:",pgx,". skipping\n"))
            next()
        }        

        ##---------------------------------------------
        ## extract the meta FC matrix
        ##---------------------------------------------
        rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        rownames(meta$fc) <- toupper(rownames(meta$fc))
        missing.FC[[pgx]] <- meta$fc

        ##---------------------------------------------
        ## compile the info for update
        ##---------------------------------------------
        
        cnd = colnames(ngs$samples)
        cnd = cnd[grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",cnd,invert=TRUE)]
        is.mouse = (mean(grepl("[a-z]",ngs$genes$gene_name))>0.8)
        organism = c("human","mouse")[1 + is.mouse]
        date = ifelse(is.null(ngs$date), "", as.character(ngs$date))

        this.info <- c(
            dataset = pgx,
            datatype = ifelse(is.null(ngs$datatype),"", ngs$datatype),
            description = ifelse(is.null(ngs$description),"", ngs$description),
            organism = organism,
            nsamples = nrow(ngs$samples),
            ngenes = nrow(ngs$X),
            nsets = nrow(ngs$gsetX),
            conditions = paste(cnd,collapse=" "),
            date = date
        )
        
        info.cols <- unique(c(info.cols, names(this.info)))
        if(!is.null(pgxinfo) && NCOL(pgxinfo)>0 && nrow(pgxinfo)>0 ) {
            this.info = this.info[match(info.cols,names(this.info))]
            names(this.info) = info.cols
            pgxinfo = pgxinfo[,match(info.cols,colnames(pgxinfo)),drop=FALSE]
            colnames(pgxinfo) = info.cols
            pgxinfo <- rbind( pgxinfo, this.info)
        } else {
            pgxinfo <- rbind(NULL,this.info)
        }
        ## pgxinfo <- rbind( pgxinfo, this.info)
    }
    if(verbose) cat("\n")
    rownames(pgxinfo) <- NULL    
    pgxinfo <- data.frame(pgxinfo)    

    ##----------------------------------------------------------------------
    ## Update the INFO meta file
    ##----------------------------------------------------------------------    
    if(force) {
        ## remove unneccessary entries if forced.
        sel <- which(pgxinfo$dataset %in% pgx.files)
        sel <- which(sub(".pgx$","",pgxinfo$dataset) %in% sub(".pgx$","",pgx.files))
        pgxinfo <- pgxinfo[sel,,drop=FALSE]
    }
    if(verbose) cat("writing new PGX.INFO file to",info.file1,"...\n")
    write.csv(pgxinfo, file = info.file1)
    Sys.chmod(info.file1, "0666")

    ##----------------------------------------------------------------------
    ## Update the ALL.FC meta file
    ##----------------------------------------------------------------------
    
    ## find most common genes
    all.gg <- toupper(as.character(unlist(sapply(missing.FC, rownames))))
    gg.tbl <- table(all.gg)
    table(gg.tbl)
        
    ## Conform the multiple metaFC matrices
    gg <- names(gg.tbl)
    length(gg)
    missing.FC <- lapply(missing.FC, function(x) {
        x <- x[match(gg,toupper(rownames(x))),,drop=FALSE]
        rownames(x) <- gg
        return(x)
    })
    
    ## append file name in front of contrast names
    id <- paste0("[",sub("[.]pgx","",names(missing.FC)),"]")
    id
    for(i in 1:length(missing.FC)) {
        colnames(missing.FC[[i]]) <- paste0(id[i]," ",colnames(missing.FC[[i]]))
    }
    allFC.new <- do.call(cbind, missing.FC)
    allFC.new <- as.matrix(allFC.new)

    if(is.null(allFC)) {
        allFC <- allFC.new
    } else {
        gg <- sort(unique(c(rownames(allFC), rownames(allFC.new))))
        j1 <- match(gg, rownames(allFC))
        j2 <- match(gg, rownames(allFC.new))
        allFC <- allFC[j1,,drop=FALSE]
        allFC.new <- allFC.new[j2,,drop=FALSE]
        allFC <- cbind(allFC, allFC.new)
        rownames(allFC) <- gg
    }
    dim(allFC)
    
    ## restrict to 8000 genes
    allfc.sd <- apply(allFC, 1, sd, na.rm=TRUE)
    allfc.nna <- rowMeans(!is.na(allFC))
    jj <- head( order(-allfc.sd * allfc.nna), 8000)
    allFC <- allFC[jj,,drop=FALSE]
    dim(allFC)

    ## check for duplicates
    allFC <- allFC[,!duplicated(colnames(allFC)),drop=FALSE]
    allFC <- allFC[,order(colnames(allFC)),drop=FALSE]
    dim(allFC)
    
    if(verbose) cat("writing all fold-changes to",allfc.file1,"...\n")
    write.csv(allFC, file=allfc.file1)
    Sys.chmod(allfc.file1, "0666")

    ##load(file="../files/allFoldChanges.rda", verbose=1)
    return(pgxinfo)
}




pgx.updateDatasetProfiles.NOTUSED <- function(pgx.dir, file="datasets-allFC.csv",
                                          force=FALSE, verbose=TRUE)
{
    i=1
    for(i in 1:length(pgx.dir)) {
        if(!dir.exists(pgx.dir)) next()
        df <- pgx.updateDatasetProfiles1(pgx.dir[i], file=file,
                                         force=force, verbose=verbose)
    }    
}

pgx.updateInfoFile.NOTUSED <- function(pgx.dir, file="datasets-info.csv", 
                                   force=FALSE, verbose=TRUE )
{
    pgxinfo <- NULL
    i=1
    for(i in 1:length(pgx.dir)) {
        if(!dir.exists(pgx.dir[i])) next()
        info <- pgx.updateInfoFile1(pgx.dir[i], file=file, force=force, verbose=verbose)
        dim(info)
        info$path <- pgx.dir[i]
        if(is.null(pgxinfo)) {
            pgxinfo <- info
        } else {
            jj <- match(colnames(pgxinfo),colnames(info))
            pgxinfo <- rbind(pgxinfo, info[,jj])
        }
    }
    dim(pgxinfo)
    return(pgxinfo)
}


##pgx.dir=PGX.DIR;file="datasets-allFC.csv";verbose=1
pgx.updateDatasetProfiles1.NOTUSED <- function(pgx.dir, file="datasets-allFC.csv",
                                       force=FALSE, verbose=TRUE)
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.updateDatasetProfiles1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    
    ## all public datasets
    pgx.files <- dir(pgx.dir, pattern="[.]pgx$")
    ##pub.id <- sub("-.*","",pgx.files)
    ##pgx.files <- pgx.files[!duplicated(pub.id)]
    pgx.files

    allfc.file <- file.path(pgx.dir,file)
    has.fc <- file.exists(allfc.file)
    if(verbose && has.fc) cat("checking if allFC file",allfc.file,"exists: YES\n")
    if(verbose && !has.fc) cat("checking if allFC file",allfc.file,"exists: NO\n")
    
    ##----------------------------------------------------------------------
    ## If an allFC file exitss, check if it is done for all PGX files
    ##----------------------------------------------------------------------

    allFC <-NULL
    if(!force && file.exists(allfc.file)) {
        if(verbose) cat("checking which pgx files already done...\n")
        allFC <- read.csv(allfc.file,row.names=1,check.names=FALSE,nrow=5)
        dim(allFC)
        files.done <- gsub("\\[|\\].*","",colnames(allFC))
        files.done <- unique(paste0(files.done,".pgx"))
        pgx.files <- setdiff(pgx.files, files.done)
    }
    length(pgx.files)
    
    if(length(pgx.files)==0) {
        if(verbose) cat("all done. no update required.\n")
        return(NULL)
    }

    if(!force && file.exists(allfc.file)) {
        if(verbose) cat("checking which pgx files already done...\n")
        allFC <- fread.csv(allfc.file,row.names=1,check.names=FALSE)
        dim(allFC)
    }
    
    if(verbose) cat("scanning",length(pgx.files),"PGX files in folder\n")
    
    ##----------------------------------------------------------------------
    ## For all new PGX files, load the PGX file and get the meta FC
    ## matrix.
    ##----------------------------------------------------------------------

    FC <- list()
    pgx=pgx.files[2]
    for(pgx in pgx.files) {
        if(verbose) cat(".")        
        try.error <- try( load(file.path(pgx.dir,pgx),verbose=0) )
        if(class(try.error)=="try-error") {
            cat(paste("error in loading",pgx,"!"))
            next()
        }        
        rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
        ##meta.fx <- sapply(ngs$gx.meta$meta,function(x) x$meta.fx)
        ##rownames(meta.fx) <- toupper(rownames(ngs$gx.meta$meta[[1]]))
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        rownames(meta$fc) <- toupper(rownames(meta$fc))
        FC[[pgx]] <- meta$fc
    }
    if(verbose) cat("\n")
    
    ## find common genes
    all.gg <- toupper(as.character(unlist(sapply(FC, rownames))))
    gg.tbl <- table(all.gg)
    table(gg.tbl)
        
    ## take 8000 most frequent genes
    gg <- head(names(sort(-gg.tbl)),8000)
    length(gg)
    FC <- lapply(FC, function(x) {
        x <- x[match(gg,toupper(rownames(x))),,drop=FALSE]
        rownames(x) <- gg
        return(x)
    })
    
    ## append file name in front of contrast names
    FC.id <- paste0("[",sub("[.]pgx","",names(FC)),"]")
    for(i in 1:length(FC)) {
        colnames(FC[[i]]) <- paste0(FC.id[i]," ",colnames(FC[[i]]))
    }
    allFC.new <- do.call(cbind, FC)
    allFC.new <- as.matrix(allFC.new)
    if(!is.null(allFC)) {
        jj <- match(rownames(allFC), rownames(allFC.new))
        allFC.new <- allFC.new[jj,]
        rownames(allFC.new) <- rownames(allFC)
    }
    allFC <- cbind(allFC, allFC.new)
    
    dim(allFC)
    if(verbose) cat("writing all fold-changes to",allfc.file,"...\n")
    write.csv(allFC, file=allfc.file)
    Sys.chmod(allfc.file, "0666")

    ##load(file="../files/allFoldChanges.rda", verbose=1)
    ##return(allFC)
}

##pgx.dir=PGX.DIR;file="datasets-info.csv";force=FALSE;verbose=1
pgx.updateInfoFile1.NOTUSED <- function(pgx.dir, file="datasets-info.csv", 
                                    force=FALSE, verbose=TRUE )
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.updateInfoFile1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }

    require(shiny)
    if(verbose) cat(">>> updating data sets info file in:",pgx.dir,"\n")

    pgx.dir <- sub("/$","",pgx.dir)
    pgx.files  <- dir(pgx.dir, pattern="[.]pgx$", full.names=FALSE)

    pgxinfo.file <- file.path(pgx.dir, file)
    pgxinfo <- c()
    has.info <- file.exists(pgxinfo.file)
    if(verbose && has.info) cat("checking if PGX-info file",pgxinfo.file,"exists: YES\n")
    if(verbose && !has.info) cat("checking if PGX-info file",pgxinfo.file,"exists: NO\n")
    
    if(!force && has.info) {
        if(verbose) cat("File exists. appending to existing info file\n")

        pgxinfo = fread.csv(pgxinfo.file, stringsAsFactors=FALSE, row.names=1)
        dim(pgxinfo)       
        jj <- which(!(sub(".pgx$","",pgx.files) %in% sub(".pgx$","",pgxinfo$dataset)))
        jj
        if(length(jj)==0) {
            if(verbose) cat("nothing to scan. all up to date. \n")
            return(pgxinfo)
        }
        pgx.files = pgx.files[jj]
    }
   
    ##pgx.files =head(pgx.files,3)
    f = pgx.files[1]
    cols <- NULL
    i=1
    if(verbose) cat("scanning",length(pgx.files),"PGX files")
    for(i in 1:length(pgx.files)) {
        ##if(verbose) cat("scanning info from",pgx.files[i],"\n")
        if(verbose) cat(".")        

        try.error <- try( load(file.path(pgx.dir,pgx.files[i]), verbose=0 ))
        if(class(try.error)=="try-error") {
            warning(paste("error in loading",pgx.files[i],"!"))
            next()
        }        
        
        cnd = colnames(ngs$samples)
        cnd = cnd[grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",cnd,invert=TRUE)]
        is.mouse = (mean(grepl("[a-z]",ngs$genes$gene_name))>0.8)
        organism = c("human","mouse")[1 + is.mouse]
        date = ifelse(is.null(ngs$date), "", as.character(ngs$date))

        this.info <- c(
            dataset = pgx.files[i],
            datatype = ifelse(is.null(ngs$datatype),"", ngs$datatype),
            description = ifelse(is.null(ngs$description),"", ngs$description),
            organism = organism,
            nsamples = nrow(ngs$samples),
            ngenes = nrow(ngs$X),
            nsets = nrow(ngs$gsetX),
            conditions = paste(cnd,collapse=" "),
            date = date
        )

        cols <- unique(c(cols, names(this.info)))
        if(!is.null(pgxinfo) && NCOL(pgxinfo)>0 ) {
            this.info = this.info[match(cols,names(this.info))]
            names(this.info) = cols
            pgxinfo = pgxinfo[,match(cols,colnames(pgxinfo)),drop=FALSE]
            colnames(pgxinfo) = cols
        }
        ##cat("i=",i,": ",length(this.info),"\n")
        pgxinfo <- rbind( pgxinfo, this.info)
    }
    if(verbose) cat("\n")
    rownames(pgxinfo) <- NULL
    
    pgxinfo <- data.frame(pgxinfo)    
    if(verbose) cat("writing pgx info to",pgxinfo.file,"...\n")
    write.csv(pgxinfo, file=pgxinfo.file)
    Sys.chmod(pgxinfo.file, "0666")
    return(pgxinfo)
}

