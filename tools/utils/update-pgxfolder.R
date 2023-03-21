RDIR = "../R"
FILES = "../lib"
FILESX = "../libx"
PGX.DIR = "../data"
FILES
source(file.path(RDIR,"pgx-include.R"))

pgx.files <- dir(".", pattern=".pgx$")
##pgx.files <- grep("X.pgx$",pgx.files,invert=TRUE,value=TRUE)
##pgx.files <- dir("../data.BAK/", pattern=".pgx",full.names=TRUE)

##pgx.files <- grep('155249',pgx.files,value=TRUE)
##pgx.files <- grep('birabresib',pgx.files,value=TRUE)
##pgx.files <- grep('Nnm.pgx',pgx.files,value=TRUE)
grep('ibet',pgx.files)

pgx.files
pgx.file = pgx.files[1]
pgx.file = "geiger2016-arginine.pgx"
pgx.file = "GSE22886-immune.pgx"
pgx.file

for(pgx.file in pgx.files) {
    
    cat("*********** updating",pgx.file,"**********\n")    
    load(pgx.file, verbose=1)
    object.size(ngs)/1e6

    names(ngs)
    names(ngs$connectivity)
    names(ngs$drugs)

    has.r <- ("clust" %in% names(ngs$drugs[[1]]))
    has.r
    if(has.r) next()
    
    if(0 && "cluster.genes" %in% names(ngs)) {
        cat("already done. skipping...\n")
        next()
    }

    if(0) {
        message("clustering genes...")
        ngs = pgx.clusterGenes(ngs, methods='umap', dims=c(2,3), level='gene')
        ngs = pgx.clusterGenes(ngs, methods='umap', dims=c(2,3), level='geneset')
        names(ngs$cluster.genes$pos)
    }

    if(0) {
        ##extra <- c("meta.go","deconv","infer","drugs","wordcloud")
        extra <- c("drugs","connectivity")
        extra <- c("drugs-combo")
        extra <- c("connectivity")
        ##ngs$connectivity <- NULL
        ##ngs$drugs <- NULL
        sigdb = "../libx/sigdb-gtex.h5"
        sigdb = c("../libx/sigdb-lincs.h5","../libx/sigdb-creeds.h5","../libx/sigdb-drugsx.h5")
        all.db <- dir("../libx","sigdb.*h5$")
        db1 <- setdiff(all.db, names(ngs$connectivity))
        db1
        ## if(length(db1)==0) next()    
        sigdb = paste0("../libx/",db1)
        sigdb = NULL
        ngs <- compute.extra(ngs, extra, lib.dir=FILES, sigdb=sigdb) 
        ngs$connectivity[["sigdb-lincs.h5"]] <- NULL  ## old, now split into cp and gt
        names(ngs$connectivity)
    }
    
    if(1) {
        extra <- c("drugs")    
        ngs <- compute.extra(ngs, extra, lib.dir=FILES, libx.dir=FILESX)     
        names(ngs$drugs)
        names(ngs$drugs[[1]])
        lapply(ngs$drugs, names)
    }
    
    ##------------------ save new object -------------------
    names(ngs)
    pgx.file    
    ngs.save(ngs, file=pgx.file)    
    
}










