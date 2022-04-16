##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##access.dirs=FILESX

##filter.get="omicsplayground";from=NULL;to=NULL;unique=TRUE
pgx.parseAccessLogs <- function(logs.dir, from=NULL, to=NULL,
                                filter.get=NULL, unique = TRUE)
{
    
    ##logs.dir <- c(FILESX, file.path(FILESX,"apache2"),"/var/www/html/logs","/var/log/apache2")
    logs.dir <- logs.dir[dir.exists(logs.dir)]
    if(length(logs.dir)==0) return(NULL)
    logs.dir
    access.files <- lapply(logs.dir, dir, pattern="access.log",
                           full.names=TRUE, recursive=TRUE )
    access.files
    access.files <- unlist(access.files)
    ##access.files <- grep("shinyproxy",access.files,value=TRUE,invert=TRUE)
    access.files
    
    access.logs <- lapply(access.files, function(f)
        suppressMessages(suppressWarnings(try(read.table(f)))))
    ##access.logs <- lapply(access.files, function(f)
    ##    suppressMessages(suppressWarnings(try(data.table::fread(f,sep=" ")))))
    access.logs <- access.logs[sapply(access.logs,class)!="try-error"]
    access.logs <- access.logs[sapply(access.logs,nrow)>0]
    length(access.logs)
    if(length(access.logs)==0) return(NULL)
    
    i=3
    for(i in 1:length(access.logs)) {
        df <- data.frame(access.logs[[i]])
        cols <- c("V1","V4","V6")
        if(mean(grepl(":80$",df[,"V1"])) > 0.9) {
            cols <- c("V2","V5","V7")
        }
        access.logs[[i]] <- df[,cols]
        colnames(access.logs[[i]]) <- c("ip","date","get")
    }    
    remove(df)
    
    ## Filter access log
    acc <- do.call(rbind, access.logs)
    dim(acc)

    if(!is.null(filter.get)) {
        sel <- grep(filter.get,acc[,"get"])
        acc <- acc[sel,]
    }
    dim(acc)
    Matrix::head(acc)
    
    ## Extract visiting period
    ## if the operating system is not windows set the timezone to LC_TIME
    if(Sys.info()["sysname"] != "Windows") {
        Sys.setlocale("LC_TIME","en_US.UTF-8")
    }   
    ##Sys.setlocale("LC_TIME","C") ## just to make sure
    acc$date <- gsub("[:].*|\\[","",as.character(acc[,"date"]))
    acc$date <- as.Date(acc$date, format = "%d/%b/%Y")
    acc <- acc[order(acc$date),]

    from.date <- Matrix::head(acc$date,1)
    to.date <- Matrix::tail(acc$date,1)
    from.to <- paste(from.date,"-",to.date)
    from.to
    
    ## Extract IP
    acc.ip <- as.character(acc[,"ip"])
    ##loc <- ip_api(unique(acc.ip))
    unique.ip <- unique(acc.ip)
    ## unique.hostname <- ip_to_hostname(unique.ip)  ## takes loooonnnggg... time!!
    ## names(unique.hostname)  <- unique.ip
    
    ## create lookup-table for IP to country
    
    file <- system.file("extdata","GeoLite2-Country.mmdb", package = "rgeolocate")
    loc  <- rgeolocate::maxmind(unique.ip, file, c("country_code","country_name"))
    loc$ip <- unique.ip
    ##file <- file.path(lib.dir,"GeoLite2-City.mmdb")
    ##loc <- rgeolocate::maxmind(ip, file, c("country_code", "country_name", "city_name"))
    loc$country_name[which(loc$ip %in% c("127.0.0.1"))] <- "<local.ip>"
    ##loc$country_code[which(loc$ip %in% c("127.0.0.1"))] <- "<local.ip>"
    loc$country_name[is.na(loc$country_name)] <- "(unknown)"
    
    country_codes <- unique(loc$country_code)
    names(country_codes) <- loc[match(country_codes,loc$country_code),"country_name"]
    country_codes["(unknown)"] = "(unknown)"

    ## now map IP to country_code
    acc$country_code <- loc$country_code[match(acc.ip,loc$ip)]
    acc$country_name <- loc$country_name[match(acc.ip,loc$ip)]
    acc$country_code[is.na(acc$country_code)] <- "(unknown)"
    acc$country_name[is.na(acc$country_name)] <- "(unknown)"
    Matrix::tail(sort(table(acc$country_code)),40)

    if(0) {
        getDodgy <- function(acc0,n=100) {
            ii <- grep("omicsplayground",acc0[,"get"],invert=TRUE)
            ii <- ii[which(nchar(as.character(acc0[ii,"get"])) > 20)]
            Matrix::head(acc0[ii,],n=n)
        }
        Matrix::tail(sort(table(acc$country_code)),40)
        ii <- which(acc.cc=="US")
        getDodgy(acc[ii,],100)
        ii <- which(acc.cc=="CN")
        getDodgy(acc[ii,],100)        
    }
    
    ## cumulative table: counting total number of hits since start
    acc$days <- acc$date - from.date + 1
    ndays <- max(acc$days)    
    ncountries <- length(unique(acc$country_code))
    M <- matrix(0, nrow=ncountries, ncol=ndays)
    rownames(M) <- sort(unique(acc$country_code)) 
    d = 1
    for(d in 1:ndays) {
        jj <- which(acc$days <= d)
        if(unique) jj <- jj[which(!duplicated(acc$ip[jj]))]  ## unique per day
        tt <- table(as.character(acc$country_code[jj]))
        tt <- tt[match(rownames(M),names(tt))]
        tt[is.na(tt)] <- 0
        M[,d] <- tt
    }
    M[is.na(M)] <- 0
    colnames(M) <- as.character(from.date + 1:ncol(M) -1 )
    
    ## weekly table: counting hits per week since start
    week.since <- seq(from.date, max(acc$date), 7)
    W <- matrix(0, nrow=ncountries, ncol=length(week.since)-1)
    rownames(W) <- sort(unique(acc$country_code))
    colnames(W) <- paste0("week.",1:ncol(W))
    i = 1
    acc0 <- acc
    ## acc0 <- acc[!duplicated(acc$ip),] ## unique, new visitors
    for(i in 1:(ncol(W)-1)) {
        jj <- which(acc0$date >= week.since[i] & acc0$date < week.since[i+1] )
        if(unique) jj <- jj[which(!duplicated(acc0$ip[jj]))]
        tt <- table(as.character(acc0$country_code[jj]))
        tt <- tt[match(rownames(W),names(tt))]
        tt[is.na(tt)] <- 0
        W[,i] <- tt
    }
    
    if(0) {
        barplot(Matrix::colSums(W,na.rm=TRUE),las=3, main="unique.IP")
    }
    
    ## final table
    if(unique) {
        tt <- table(as.character(loc$country_name))  ## unique IP
    } else {
        tt <- table(as.character(acc$country_name))  ## full access log
    }
    sum(tt)
    cc <- country_codes[names(tt)]
    names(cc) <- names(tt)
    df <- data.frame( country_name = names(tt),
                     country_code = cc,
                     count = (as.integer(tt)))
    df <- df[order(-df$count),]
    sum(df$count,na.rm=TRUE)
    M <- M[match(df$country_code, rownames(M)),]
    W <- W[match(df$country_code, rownames(W)),]
    
    res <- list(visitors=df, period=from.to, from=from.date, to=to.date,
                table=M, weekly=W, access=acc)
    return(res)
}

h5exists <- function(h5.file, obj) {
        
    xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
    obj %in% gsub("^/|^//","",xobjs)
}

##h5.file="test.h5";chunk=100
pgx.saveMatrixH5 <- function(X, h5.file, chunk=NULL, del=TRUE )
{   
    
    if(del) unlink(h5.file)
    
    if(is.null(chunk)) {
        if(del) rhdf5::h5createFile(h5.file)    
        ## rhdf5::h5createGroup("myhdf5file.h5","foo")
        ## A = matrix(1:10,nr=5,nc=2)
        ## rhdf5::h5write(A, "myhdf5file.h5","foo/A")    
        rhdf5::h5createGroup(h5.file,"data")    
        rhdf5::h5write( X, h5.file, "data/matrix")
    } else {
        if(del) rhdf5::h5createFile(h5.file)    
        if(del) rhdf5::h5createGroup(h5.file,"data")
        if(h5exists(h5.file,"data/matrix")) rhdf5::h5delete(h5.file, "data/matrix")        
        rhdf5::h5createDataset(
            h5.file, "data/matrix",
            c(nrow(X),ncol(X)),
            ##storage.mode = "integer",
            chunk = chunk,
            level = 7
        )
        rhdf5::h5write(
            X,
            file = h5.file,
            name = "data/matrix",
            index = list(1:nrow(X),1:ncol(X))
        )
    }

    rhdf5::h5write( rownames(X), h5.file, "data/rownames")
    rhdf5::h5write( colnames(X), h5.file, "data/colnames")    
    
    rhdf5::h5closeAll()
}

pgx.readOptions <- function(file = "./OPTIONS") {
    if(!file.exists(file)) return(NULL)
    opt <- read.table(file, sep="=", row.names=1)
    opt <- gsub("^[ ]*|[ ]*$","",apply(opt,1,c)) ## strip leading/post spaces
    opt <- sapply(opt,list)
    opt <- sapply(opt,strsplit,split="[;,]")
    is.bool <- sapply(opt, function(x) all(tolower(x) %in% c("true","false")))
    is.bool
    opt[is.bool] <- sapply(opt[is.bool], function(x) tolower(x) %in% c("true"))
    names(opt) <- trimws(names(opt))
    opt
}

##pgx.dir=PGX.DIR;verbose=TRUE;force=FALSE
pgx.initDatasetFolder <- function(pgx.dir, verbose=TRUE, force=FALSE)
{
    ## Initialized information file for multiple folders.
    ##
    ##
    i=1
    for(i in 1:length(pgx.dir)) {
        ## pgx.initDatasetFolder(pgx.dir[i], verbose=verbose, force=force)
        if(!dir.exists(pgx.dir[i])) next()

        ## skip if no pgx files
        npgx <- length(dir(pgx.dir[i],"pgx$"))
        if(npgx==0) {
            message("[initDatasetFolder] no pgx files in",pgx.dir[i])
            next()
        }

        ## skip if too many...
        if(npgx > 100 && force==FALSE) {
            message("[initDatasetFolder] too many files in",pgx.dir[i],". Please init manually")
            next()
        }
        
        info <- pgx.initDatasetFolder1 (
            pgx.dir[i],
            allfc.file = "datasets-allFC.csv",
            info.file = "datasets-info.csv",
            force = force,
            verbose = verbose
        )    

    }

}

if(0) {
    ##h5ls(h5.file)
    pgx.createSignatureDatabaseH5( pgx.files, h5.file, update.only=FALSE)
    rhdf5::h5ls(h5.file)
    pgx.addEnrichmentSignaturesH5(h5.file, X=NULL, mc.cores=8,
                                  lib.dir=FILES, methods=c("gsea"))
    rhdf5::h5ls(h5.file)
}


##verbose=TRUE;file="datasets-allFC.csv"
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
        stop(paste("[readDatasetProfiles1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    fn <- file.path(pgx.dir,file)
    fn
    if(!file.exists(fn)) {
        stop("FATAL : could not find profiles matrix. please create first with initDatasetFolder().\n")
        ## pgx.updateDatasetProfiles(pgx.dir, file=file)
        return()
    } else {
        if(verbose) message("[readDatasetProfiles1] Found existing dataset profiles matrix")
    }
    
    allFC <- fread.csv(file=file.path(pgx.dir, file),
                           row.names=1, check.names=FALSE)
    allFC <- as.matrix(allFC)
    if(verbose) message("[readDatasetProfiles1] dataset profiles matrix : dim=",dim(allFC))
    dim(allFC)    
    return(allFC)
}


file="datasets-info.csv";force=FALSE
pgx.scanInfoFile <- function(pgx.dir, file="datasets-info.csv", force=FALSE, verbose=TRUE)
{
    pgxinfo <- NULL
    i=1
    for(i in 1:length(pgx.dir)) {        
        pgx.initDatasetFolder1(pgx.dir[i], force=force, verbose=TRUE)        
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

##pgx.dir1=PGX.DIR[1];allfc.file="datasets-allFC.csv";verbose=1;info.file="datasets-info.csv";force=0
pgx.initDatasetFolder1 <- function( pgx.dir1,
                                   allfc.file = "datasets-allFC.csv",
                                   info.file = "datasets-info.csv",
                                   force=FALSE, verbose=TRUE)
{
    ##
    ## Initialize file information file for SINGLE folder
    ##
    ##

    if(!dir.exists(pgx.dir1)) {
        stop(paste("[initDatasetFolder1] FATAL ERROR : folder",pgx.dir,"does not exist"))
    }
    
    ## all public datasets
    pgx.dir1 <- pgx.dir1[1]  ## only one folder!!!
    pgx.files <- dir(pgx.dir1, pattern="[.]pgx$")
    pgx.files

    ##----------------------------------------------------------------------
    ## If an allFC file exists
    ##----------------------------------------------------------------------

    allfc.file1 <- file.path(pgx.dir1, allfc.file)
    has.fc <- file.exists(allfc.file1)
    if(verbose && has.fc) message("[initDatasetFolder1] file ",allfc.file1," exists: YES")
    if(verbose && !has.fc) message("[initDatasetFolder1] file ",allfc.file1," exists: NO")

    info.file1 <- file.path(pgx.dir1, info.file)
    has.info <- file.exists(info.file1)
    if(verbose && has.info) message("[initDatasetFolder1] file ",info.file1," exists: YES")
    if(verbose && !has.info) message("[initDatasetFolder1] file ",info.file1," exists: NO")
    
    ##----------------------------------------------------------------------
    ## If an allFC file exits, check if it is done for all PGX files
    ##----------------------------------------------------------------------

    pgxinfo <- NULL
    pgx.missing0 <- pgx.files
    pgx.missing1 <- pgx.files
    force

    allFC <-NULL
    if(!force && file.exists(allfc.file1)) {
        if(verbose) message("[initDatasetFolder1] checking which pgx files already done in allFC...")
        allFC <- read.csv(allfc.file1,row.names=1,check.names=FALSE,nrows=5)  ## just HEADER!!!
        dim(allFC)
        files.done <- gsub("\\[|\\].*","",colnames(allFC))
        files.done <- unique(paste0(files.done,".pgx"))
        pgx.missing0 <- setdiff(pgx.missing0, files.done)
    }
    dim(allFC)
    
    if(!force && has.info) {
        if(verbose) message("[initDatasetFolder1] checking which pgx files already in PGX info...")        
        pgxinfo = fread.csv(info.file1, stringsAsFactors=FALSE, row.names=1)
        dim(pgxinfo)       
        jj <- which(!(sub(".pgx$","",pgx.missing1) %in% sub(".pgx$","",pgxinfo$dataset)))
        jj
        pgx.missing1 = pgx.missing1[jj]
    }
   
    ##----------------------------------------------------------------------
    ## Check if it is done for all PGX files
    ##----------------------------------------------------------------------
    
    ## files to be done either for allFC or missing in INFO
    pgx.missing <- unique(c(pgx.missing0, pgx.missing1))
    length(pgx.missing)
    if(verbose) message("[initDatasetFolder1] FORCE = ",force)
    
    if(length(pgx.missing)==0) {
        if(verbose) message("[initDatasetFolder1] no update required. use FORCE=1 for forced update.")
        return(NULL)
    }
    if(verbose) message("[initDatasetFolder1] scanning ",length(pgx.missing)," PGX files in folder ",pgx.dir1)
        
    ##----------------------------------------------------------------------
    ## Reread allFC file. Before we only read the header.
    ##----------------------------------------------------------------------
    allFC <-NULL
    if(!force && file.exists(allfc.file1)) {
        ##allFC <- read.csv(allfc.file1,row.names=1,check.names=FALSE)
        allFC <- fread.csv(allfc.file1,row.names=1,check.names=FALSE)  
    }
    dim(allFC)

    ##----------------------------------------------------------------------
    ## For all new PGX files, load the PGX file and get the meta FC
    ## matrix.
    ##----------------------------------------------------------------------

    info.cols <- NULL
    missing.FC <- list()
    message("[initDatasetFolder1] missing pgx = ",pgx.missing)
    pgxfile = pgx.missing[1]
    
    ngs <- NULL
    for(pgxfile in pgx.missing) {

        pgxfile
        try.error <- try( load(file.path(pgx.dir1,pgxfile),verbose=0) )
        if(class(try.error)=="try-error") {
            message(paste("[initDatasetFolder1] ERROR in loading PGX file:",pgxfile,". skipping\n"))
            next()
        }                

        message("[initDatasetFolder1] pgxfile = ",pgxfile)
        message("[initDatasetFolder1] names(ngs) = ", paste(names(ngs),collapse=' '))

        if(!pgx.checkObject(ngs)) {
            message(paste("[initDatasetFolder1] INVALID PGX object",pgxfile,". Skipping"))
            next()            
        }

        ## check if name exists
        ##if(is.null(ngs$name)) ngs$name <- sub(".pgx$","",pgxfile)
        ngs$name <- sub(".pgx$","",pgxfile)  ## force filename as name
        
        ##---------------------------------------------
        ## extract the meta FC matrix
        ##---------------------------------------------
        ## rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        rownames(meta$fc) <- toupper(rownames(meta$fc))
        missing.FC[[pgxfile]] <- meta$fc

        ##---------------------------------------------
        ## compile the info for update
        ##---------------------------------------------
        pgxinfo <- pgx.updateInfoPGX(pgxinfo, ngs)
        Matrix::tail(pgxinfo)
        ## pgxinfo <- rbind( pgxinfo, this.info)

    }

    ngs <- NULL
    rownames(pgxinfo) <- NULL    
    pgxinfo <- data.frame(pgxinfo)    

    if(length(missing.FC)==0) {
        ## no valid new files
        return(pgxinfo)
    }

    
    ##----------------------------------------------------------------------
    ## Update the INFO meta file
    ##----------------------------------------------------------------------    
    if(force) {
        ## remove unneccessary entries if forced.
        sel <- which(pgxinfo$dataset %in% pgx.files)
        sel <- which(sub(".pgx$","",pgxinfo$dataset) %in% sub(".pgx$","",pgx.files))
        pgxinfo <- pgxinfo[sel,,drop=FALSE]
    }
    if(verbose) message("[initDatasetFolder1] writing updated PGX.INFO file to ",info.file1,"...")
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
    jj <- Matrix::head( order(-allfc.sd * allfc.nna), 8000)
    allFC <- allFC[jj,,drop=FALSE]
    dim(allFC)
    
    ## check for duplicates
    allFC <- allFC[,!duplicated(colnames(allFC)),drop=FALSE]
    allFC <- allFC[,order(colnames(allFC)),drop=FALSE]
    dim(allFC)
    
    if(verbose) message("[initDatasetFolder1] writing updated all fold-changes to",allfc.file1,"...")
        write.csv(allFC, file=allfc.file1)
    Sys.chmod(allfc.file1, "0666")
    
    ##load(file="../files/allFoldChanges.rda", verbose=1)
    return(pgxinfo)
}


pgx.updateInfoPGX <- function(pgxinfo, ngs, remove.old=TRUE)
{

    cond = grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
                colnames(ngs$samples),invert=TRUE,value=TRUE)
    cond = grep("title|source|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
                colnames(ngs$samples),invert=TRUE,value=TRUE)
    cond

    is.mouse = (mean(grepl("[a-z]",ngs$genes$gene_name))>0.8)
    organism = c("human","mouse")[1 + is.mouse]
    if("organism" %in% names(ngs)) organism <- ngs$organism
    
    this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    date = ifelse(is.null(ngs$date), this.date, as.character(ngs$date))
    dataset.name <- ngs$name
    ## dataset.name <- ifelse(is.null(ngs$name), pgxfile, ngs$name)
        
    this.info <- c(
        dataset = dataset.name,
        ## author = "", ## add author? maintainer? owner??
        collection = ngs$collection,
        datatype = ifelse(is.null(ngs$datatype),"", ngs$datatype),
        description = ifelse(is.null(ngs$description),"", ngs$description),
        organism = organism,
        nsamples = nrow(ngs$samples),
        ngenes = nrow(ngs$X),
        nsets = nrow(ngs$gsetX),
        conditions = paste(cond,collapse=" "),
        date = as.character(date),
        path = NULL
    )

    ## force to be character...
    !is.null(pgxinfo) && NCOL(pgxinfo)>0 && nrow(pgxinfo)>0
    
    if(!is.null(pgxinfo) && NCOL(pgxinfo)>0 && nrow(pgxinfo)>0 )
    {
        if("date" %in% colnames(pgxinfo)) {
            pgxinfo$date <- as.character(pgxinfo$date)
        }
        which.factor <- which(sapply(pgxinfo,is.factor))
        which.factor
        for(i in which.factor) {
            pgxinfo[,i] <- as.character(pgxinfo[,i])
        }
        
        ## remove existing entries??
        if(remove.old && nrow(pgxinfo)>0 ) {
            d1 <- sub("[.]pgx$","",pgxinfo$dataset)
            d2 <- sub("[.]pgx$","",this.info["dataset"])
            if( !is.null(d2) && !is.na(d2) && 
                d2 %in% d1 && d2!="") {
                sel <- which(d1!=d2)
                pgxinfo <- pgxinfo[sel,,drop=FALSE]
            }
        }
        
        ## merge with same columns
        info.cols <- colnames(pgxinfo)
        info.cols <- unique(c(info.cols, names(this.info)))
        this.info = this.info[match(info.cols,names(this.info))]
        names(this.info) = info.cols
        pgxinfo1 <- pgxinfo
        for(f in setdiff(info.cols,colnames(pgxinfo1))) {
            pgxinfo1[[f]] <- NA
        }
        match(info.cols,colnames(pgxinfo1))
        pgxinfo1 = pgxinfo1[,match(info.cols,colnames(pgxinfo1)),drop=FALSE]
        colnames(pgxinfo1) = info.cols        
        pgxinfo <- rbind( pgxinfo1, this.info)
    } else {
        pgxinfo <- data.frame(rbind(this.info))
    }

    pgxinfo
}