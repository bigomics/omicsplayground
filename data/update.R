library(knitr)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(matrixTests)
library(kableExtra)
library(knitr)
library(rhdf5)

FILES = "../lib"
FILESX = "../libx"
RDIR  = "../R"
PGX.DIR = "../data"
OUT.DIR = "."
source("../R/pgx-include.R")
##source("options.R")
FILES

pgx.files <- dir(".", pattern=".pgx")
##pgx.files <- grep("X.pgx$",pgx.files,invert=TRUE,value=TRUE)
##pgx.files <- dir("../data.BAK/", pattern=".pgx",full.names=TRUE)
pgx.files
pgx.file = pgx.files[7]
pgx.file

for(pgx.file in pgx.files) {
    
    cat("*********** updating",pgx.file,"**********\n")    
    load(pgx.file, verbose=1)
    object.size(ngs)/1e6
    
    if(0 && "connectivity" %in% names(ngs)) {
        cat("already done. skipping...\n")
        next()
    }
    
    ##extra <- c("meta.go","deconv","infer","drugs","wordcloud")
    extra <- c("connectivity")
    ngs$connectivity <- NULL
    ngs <- compute.extra(ngs, extra, lib.dir=c(FILES,FILESX) )     
    names(ngs$connectivity)
    sort(sapply(ngs$connectivity,object.size))
    
    ##------------------ save new object -------------------
    names(ngs)
    pgx.file0 <- sub(".*[/]","",pgx.file)
    pgx.file0 <- file.path(OUT.DIR,pgx.file0)
    pgx.file0
    ngs.save(ngs, file=pgx.file0)    
}










