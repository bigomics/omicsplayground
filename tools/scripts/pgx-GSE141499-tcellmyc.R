##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##

OPG="~/Playground/omicsplayground/"
RDIR = file.path(OPG,"R")
FILES = file.path(OPG,"lib")
FILESX = file.path(OPG,"libx")
PGX.DIR = file.path(OPG,"data")
source(file.path(RDIR,"pgx-include.R"))
source(file.path(RDIR,"pgx-getgeo.R"))

##------------------------------------------------------------
## Get GEO data
##------------------------------------------------------------

## load series and platform data from GEO
archs.h5 = file.path(FILESX,"human_matrix.h5")
geo <- pgx.getGEOseries(id="GSE141499", archs.h5=archs.h5)
names(geo)
head(geo$samples)
head(geo$meta)

##------------------------------------------------------------
## Set PGX information
##------------------------------------------------------------

pgx = pgx.createPGX(
    counts = geo$counts,
    samples = geo$samples,
    contrasts = geo$contrasts
)

gx.methods    = c("trend.limma","edger.qlf","deseq2.wald")
gset.methods  = c("fisher","gsva","fgsea")
extra.methods = c("deconv","wordcloud","connectivity")
extra.methods  = c("meta.go","infer","drugs","deconv","wordcloud","connectivity")

pgx <- pgx.computePGX(
    pgx,
    max.genes = 20000,
    max.genesets = 5000, 
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    use.design = TRUE,        ## no.design+prune are combined 
    prune.samples = FALSE,  ##
    do.cluster = TRUE,                
    progress = NULL,
    lib.dir = FILES 
)

extra.methods <- c("connectivity")
pgx <- compute.extra(pgx, extra.methods, lib.dir=FILES) 

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------

rda.file="../data/GSE141499-tcellmyc.pgx"
ngs=pgx   ## still old...
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$organism = "mouse"
ngs$description = "GSE141499. Homeostasis and transitional activation of regulatory T cells require c-Myc (Zheng et al, Sci Adv 2020)."

names(ngs)
ngs$timings

rda.file
ngs.save(ngs, file=rda.file)












