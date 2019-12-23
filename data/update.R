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

source("../R/gset-gsea.r")
source("../R/gset-meta.r")
source("../R/pgx-functions.R")
source("../R/pgx-deconv.R")
source("../R/pgx-proteomics.R")

source("../R/pgx-graph.R")
source("../R/pgx-drugs.R")
source("../R/pgx-wordcloud.R")
source("../R/compute2-genes.R")
source("../R/compute2-genesets.R")
source("../R/compute2-extra.R")    

source("options.R")
FILES
MAX.GENES
MAX.GENES=5000

pgx.files <- dir(".", pattern=".pgx")
pgx.file = pgx.files[3]
pgx.file

for(pgx.file in pgx.files) {
    
    load(pgx.file, verbose=1)

    extra <- c("meta.go","deconv","infer","drugs","wordcloud")
    extra <- c("wordcloud")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES)     
    names(ngs)
    ngs.save(ngs, file=pgx.file)
}








