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

pgx.files <- dir(".", pattern=".pgx")
##pgx.files <- dir("../data.BAK/", pattern=".pgx",full.names=TRUE)
pgx.files
pgx.file = pgx.files[3]
pgx.file

for(pgx.file in pgx.files) {
    
    load(pgx.file, verbose=1)
    object.size(ngs)/1e6
    
    ##extra <- c("meta.go","deconv","infer","drugs","wordcloud")
    extra <- c("connectivity")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES)     
    names(ngs$connectivity)
    
    ##------------------ save new object -------------------
    names(ngs)
    pgx.file0 <- sub(".*[/]","",pgx.file)
    pgx.file0
    ngs.save(ngs, file=pgx.file0)
    
}








