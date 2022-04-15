##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##
##
##

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")

##source("options.R")
MAX.GENES = 8000
MAX.GENESETS = 8000
MAX.GENES = 999999
MAX.GENESETS = 999999
BATCH.CORRECT=FALSE
BATCH.CORRECT=TRUE
SUBSAMPLE=TRUE

## run all available methods 
GENE.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","fgsea")

rda.file="../data/GSE98588-dlbcl.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE98588 DLBCL (Chapuy et al, Nat Med. 2018)"


##------------------------------------------------------------
## Read data
##------------------------------------------------------------
if(0) {
    BiocManager::install("hgu133plus2.db")
    BiocManager::install("GEOquery")
}
library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE98588", GSEMatrix =TRUE, AnnotGPL=TRUE)
length(geo)
attr(geo, "names")

pdata = pData(geo[[1]])
head(pdata)

sampleTable <- data.frame(pdata[,c("title","geo_accession","source_name_ch1")])

## merge data sets
X = exprs(geo[[1]])
max(X)  ## check max for checking logarithm
X = X[order(-apply(X,1,sd)),]
dim(X)

## convert affymetrix ID to GENE symbol
orf = fData(geo[[1]])$ORF
symbol = probe2symbol(orf)
sel <- which(!duplicated(symbol) & !is.na(symbol) & symbol!="")
X = X[sel,]
rownames(X) <- symbol[sel]
dim(X)

sum(is.na(X))
sum(rowSums(X)==0)
dim(X)

## conform tables
table(rownames(sampleTable) == colnames(X))

##-------------------------------------------------------------------
## create ngs object
##-------------------------------------------------------------------

samples = sampleTable
counts = 2**X

## Create contrasts 
df <- 
contr.matrix <- sampleTable[,0]


ngs <- pgx.createPGX(
    counts = counts,
    samples = samples,
    contrasts = contr.matrix,
    is.logx = FALSE,
    do.cluster = TRUE,
    do.clustergenes = TRUE,
    batch.correct = FALSE,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = FALSE,
    only.hugo = FALSE,
    convert.hugo = FALSE,
    only.proteincoding = FALSE
)
names(ngs)

gx.methods    = c("trend.limma")
gset.methods  = c("fisher")
extra.methods  = c("meta.go","infer")
extra.methods  = c("drugs")

gx.methods    = c("ttest","ttest.welch","trend.limma","edger.qlf","deseq2.wald")
gset.methods  = c("fisher","gsva","fgsea","spearman","camera","fry")
extra.methods  = c("meta.go","infer","drugs","deconv","wordcloud","connectivity")

ngs <- pgx.computePGX(
    ngs,
    max.genes = 20000,
    max.genesets = 5000, 
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    use.design = TRUE,      ## no.design+prune are combined 
    prune.samples = FALSE,  ##
    do.cluster = TRUE,                
    progress = NULL,
    lib.dir = FILES 
)

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------

rda.file="../data/geiger2016-arginine-test.pgx"
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "LC/MS proteomics"
ngs$description = "Proteome profiles of activated  vs resting human naive T cells at different times (Geiger et al., Cell 2016)."
ngs$organism = 'human'

rda.file
ngs.save(ngs, file=rda.file)


##===================================================================
##========================= END OF FILE =============================
##===================================================================








