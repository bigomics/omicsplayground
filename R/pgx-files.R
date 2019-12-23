##pgx.dir=PGX.DIR;file="datasets-allFC.csv"

pgx.forceInitDatasetFolder <- function(datadir) {
    cat("force init of data folder (info and allFC) ...\n")
    pgx.scanDatasetProfiles(datadir, file="datasets-allFC.csv")
    pgx.updateInfoFile(datadir, file="datasets-info.csv", force=TRUE) 
}

pgx.readDatasetProfiles <- function(pgx.dir, file="datasets-allFC.csv") {
    fn <- file.path(pgx.dir,file)
    fn
    if(!file.exists(fn)) {
        pgx.scanDatasetProfiles(pgx.dir, file=file)
    }
    allFC <- read.csv(file=file.path(pgx.dir, file),
                      row.names=1, check.names=FALSE)
    allFC <- as.matrix(allFC)
    dim(allFC)    
    return(allFC)
}

pgx.scanDatasetProfiles <- function(pgx.dir, file="datasets-allFC.csv") {
    
    ## all public datasets
    pgx.files <- dir(pgx.dir, pattern="[.]pgx$")
    pub.id <- sub("-.*","",pgx.files)
    pgx.files <- pgx.files[!duplicated(pub.id)]
    pgx.files

    allfc.file <- file.path(pgx.dir,file)
    
    aa <- c()
    if(file.exists(allfc.file)) {
        aa <- read.csv(allfc.file,row.names=1,check.names=FALSE)
        dim(aa)
    }
    
    FC <- list()
    pgx=pgx.files[2]
    for(pgx in pgx.files) {
        load(file.path(pgx.dir,pgx),verbose=0)
        rownames(ngs$X) <- toupper(sub(".*:","",rownames(ngs$X)))
        ##meta.fx <- sapply(ngs$gx.meta$meta,function(x) x$meta.fx)
        ##rownames(meta.fx) <- toupper(rownames(ngs$gx.meta$meta[[1]]))
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        rownames(meta$fc) <- toupper(rownames(meta$fc))
        FC[[pgx]] <- meta$fc
    }
    
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

    ## append file/experiment name in front of contrast names
    FC.id <- paste0("[",sub("-.*","",names(FC)),"]")
    for(i in 1:length(FC)) {
        colnames(FC[[i]]) <- paste0(FC.id[i]," ",colnames(FC[[i]]))
    }

    allFC <- do.call(cbind, FC)
    allFC <- as.matrix(allFC)
    dim(allFC)
    write.csv(allFC, file=allfc.file)
    ##load(file="../files/allFoldChanges.rda", verbose=1)
    ##return(allFC)
}


##pgx.dir=PGX.DIR;file="datasets-info.csv"
pgx.updateInfoFile <- function(pgx.dir, file="datasets-info.csv", force=FALSE) {

    ## what files to use, heavy or light objects
    pgxinfo <- c()
    pgxinfo.file <- file.path(pgx.dir, file)

    if(!force && file.exists(pgxinfo.file)) {
        ## info file exists, check and update        
        info = read.csv(pgxinfo.file, stringsAsFactors=FALSE, row.names=1)
        pgx.files = dir(pgx.dir, pattern=".pgx$")
        all(pgx.files %in% info$dataset)
        newfiles <- setdiff(pgx.files, info$dataset)
        cat("[pgx.updateInfo] nr of new files=",length(newfiles),"\n")
        
        if(length(newfiles)>0) {
            ## scan files not yet in info file, then update info file
            info <- pgx.scanInfo(pgx.dir=pgx.dir, inc.progress=FALSE, pgx=info, verbose=FALSE)
            rownames(info) <- NULL
            Sys.chmod(pgxinfo.file, mode="0666")
            write.csv(info, file=pgxinfo.file)
        }
        ##info = info[match(pgx.files,info$dataset),] 
        rownames(info) <- NULL
        pgxinfo <- info  
    } else {
        ## no info file, scan and create from scratch
        info <- pgx.scanInfo(pgx.dir=pgx.dir, inc.progress=FALSE, verbose=FALSE)
        Sys.chmod(pgxinfo.file, mode="0666")
        write.csv(info, file=pgxinfo.file)
        pgxinfo <- info
    }
    return(pgxinfo)
}

pgx.scanInfo <- function(pgx.dir, inc.progress=FALSE,
                         pgx=NULL, verbose=TRUE )
{
    require(shiny)
    ##cat(">>> scanning available data sets...\n")
    pgx.dir <- sub("/$","",pgx.dir)
    pgx.files  <- dir(pgx.dir, pattern="[.]pgx$", full.names=TRUE)
    pgx.files0 <- dir(pgx.dir, pattern="[.]pgx$", full.names=FALSE)
    if(!is.null(pgx)) {
        jj <- which(!sub(".pgx$","",pgx.files0) %in% sub(".pgx$","",pgx$dataset))
        jj
        if(length(jj)==0) {
            if(verbose) cat("[pgx.scanInfo] all up to date\n")
            return(pgx)
        }
        pgx.files = pgx.files[jj]
        pgx.files0 = pgx.files0[jj]
    }

    ##pgx.files =head(pgx.files,3)
    f = pgx.files[1]
    pgx.info <- c()
    cols <- NULL
    i=1
    for(i in 1:length(pgx.files)) {
        if(verbose) cat("scanning info from",pgx.files[i],"\n")
        load( pgx.files[i] )
        cnd = colnames(ngs$samples)
        cnd = cnd[grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",cnd,invert=TRUE)]
        is.mouse = (mean(grepl("[a-z]",ngs$genes$gene_name))>0.8)
        organism = c("human","mouse")[1 + is.mouse]
        date = ifelse(is.null(ngs$date), "", as.character(ngs$date))

        this.info <- c(
            dataset = pgx.files0[i],
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
        if(!is.null(pgx.info)) {
            this.info = this.info[match(cols,names(this.info))]
            names(this.info) = cols
            pgx.info = pgx.info[,match(cols,colnames(pgx.info)),drop=FALSE]
            colnames(pgx.info) = cols
        }
        ##cat("i=",i,": ",length(this.info),"\n")
        pgx.info <- rbind( pgx.info, this.info)
        if(inc.progress) incProgress( 1/length(pgx.files) )
    }
    rownames(pgx.info) <- NULL
    
    if(is.null(pgx)) {
        if(verbose) cat(">>> found",nrow(pgx.info),"data sets\n")
        pgx.info <- data.frame(pgx.info)
        rownames(pgx.info) <- NULL
    } else {
        if(verbose) cat(">>> updated",nrow(pgx.info),"data sets\n")
        pgx.info <- data.frame(pgx.info)
        pgx.info <- pgx.info[,match(colnames(pgx),colnames(pgx.info))]
        pgx.info = rbind(pgx, pgx.info)
        rownames(pgx.info) <- NULL
    }
    ##write.csv(pgx.info, file="../pgx/pgx-info.csv")
    return(pgx.info)
}

