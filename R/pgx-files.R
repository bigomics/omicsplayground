
pgx.initDatasetFolder <- function(pgx.dir, verbose=TRUE, force=FALSE)
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.initDatasetFolder] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }

    if(verbose) cat("init of data folder (info and allFC) ...\n")
    r1 <- pgx.updateDatasetProfiles(pgx.dir, file="datasets-allFC.csv",
                                  force=force, verbose=verbose)
    r2 <- pgx.updateInfoFile(pgx.dir, file="datasets-info.csv",
                             force=force, verbose=verbose) 
}

##pgx.dir=PGX.DIR;file="datasets-allFC.csv"
pgx.readDatasetProfiles <- function(pgx.dir, file="datasets-allFC.csv",
                                    verbose=TRUE)
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.initDatasetFolder] FATAL ERROR : folder",pgx.dir,"does not exist"))
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
    if(verbose) cat("Reading dataset profiles matrix...\n")
    require(data.table)
    allFC <- fread.csv(file=file.path(pgx.dir, file),
                           row.names=1, check.names=FALSE)
    allFC <- as.matrix(allFC)
    if(verbose) cat("dataset profiles matrix : dim=",dim(allFC),"\n")
    dim(allFC)    
    return(allFC)
}

##pgx.dir=PGX.DIR;file="datasets-allFC.csv";verbose=1
pgx.updateDatasetProfiles <- function(pgx.dir, file="datasets-allFC.csv",
                                      force=FALSE, verbose=TRUE)
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.initDatasetFolder] FATAL ERROR : folder",pgx.dir,"does not exist"))
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
    
    FC <- list()
    pgx=pgx.files[2]
    for(pgx in pgx.files) {
        if(verbose) cat(".")        
        try.error <- try( load(file.path(pgx.dir,pgx),verbose=0) )
        if(class(try.error)=="try-error") {
            warning(paste("error in loading",pgx,"!"))
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
pgx.updateInfoFile <- function(pgx.dir, file="datasets-info.csv", 
                               force=FALSE, verbose=TRUE )
{
    if(!dir.exists(pgx.dir)) {
        stop(paste("[pgx.initDatasetFolder] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }

    require(shiny)
    if(verbose) cat(">>> updating data sets info file... (pgx.updateInfoFile) \n")

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

