##
##
##
## Mar2019/IK: adding new scale/normalize
##
##
##

library(knitr)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(matrixTests)
library(kableExtra)
library(knitr)

FILES = "../lib/"
RDIR  = "../R/"
source("../R/pgx-include.R")
##source("options.R")
FILES
MAX.GENES
MAX.GENES=5000

pgx.files <- dir(".", pattern=".pgx")
pgx.files <- dir("../data.BAK/", pattern=".pgx",full.names=TRUE)
pgx.files
pgx.file = pgx.files[1]
pgx.file

for(pgx.file in pgx.files) {
    
    load(pgx.file, verbose=1)
    object.size(ngs)/1e6
    
    ##extra <- c("meta.go","deconv","infer","drugs","wordcloud")
    ##extra <- c("wordcloud")
    ##ngs <- compute.extra(ngs, extra, lib.dir=FILES)     

    i=1
    for(i in 1:length(ngs$gx.meta$meta)) {
        colnames(ngs$gx.meta$meta[[i]])
        rownames(ngs$gx.meta$meta[[i]]$p) <- NULL
        rownames(ngs$gx.meta$meta[[i]]$q) <- NULL
        rownames(ngs$gx.meta$meta[[i]]$fc) <- NULL
    }
    
    for(i in 1:length(ngs$gset.meta$meta)) {
        colnames(ngs$gset.meta$meta[[i]])
        rownames(ngs$gset.meta$meta[[i]]$p) <- NULL
        rownames(ngs$gset.meta$meta[[i]]$q) <- NULL
        rownames(ngs$gset.meta$meta[[i]]$fc) <- NULL        
    }
    
    object.size(ngs)/1e6

    names(ngs)
    pgx.file0 <- sub(".*[/]","",pgx.file)
    ngs.save(ngs, file=pgx.file0)
}








