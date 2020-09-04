##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
## Author: BigOmics Analytics (IK)
## Date:   2020
## 

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")
FILES

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/geiger2016-arginineX.pgx"
rda.file="../data/geiger2016-arginine.pgx"
rda.file

## load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "LC/MS proteomics"
ngs$description = "Proteome profiles of activated  vs resting human naive T cells at different times (Geiger et al., Cell 2016)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
##devtools::install_github("bartongroup/Proteus", build_opts= c("--no-resave-data", "--no-manual"), build_vignettes=FALSE)    
library(org.Hs.eg.db)
library(proteus)    
metadataFile <- "../ext/arginine/meta.txt"
proteinGroupsFile <- "../ext/arginine/proteinGroups.txt"
## Read the proteinGroups file
prot <- proteus.readProteinGroups(
    file=proteinGroupsFile, meta=metadataFile,
    collapse.gene=1, unit="intensity", is.log2=TRUE, use.LFQ=FALSE)

##-------------------------------------------------------------------
## scale/normalize counts
##-------------------------------------------------------------------    
## impute missing values
norm.counts <- prot.imputeMissing(
    prot$tab, groups=prot$metadata$condition,
    method="group.median", zero.na=TRUE)

## scale/normalize
norm.counts <- prot.normalizeCounts(
    norm.counts, scaling=1e6,
    qnormalize=TRUE, prior.count=0, plot=0)
##hist(log2(1e-8+norm.counts), breaks=100)

##-------------------------------------------------------------------
## create ngs object
##-------------------------------------------------------------------
ngs$samples = prot$metadata
colnames(ngs$samples) <- sub("condition","group",colnames(ngs$samples))
##ngs$counts = data$counts
ngs$counts = norm.counts
colnames(ngs$counts)==ngs$samples$sample
short.name <- sub(".*_tcell_","",colnames(ngs$counts))
rownames(ngs$samples)=colnames(ngs$counts)=short.name

## relevel factors??
ngs$samples$group <- relevelFactorFirst(ngs$samples$group)
ngs$samples$activated <- relevelFactorFirst(ngs$samples$activated)
ngs$samples$time <- relevelFactorFirst(ngs$samples$time)   

require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
ngs$genes = data.frame( gene_name = prot$gene,
                       gene_alias = prot$gene.names,
                       gene_title = GENE.TITLE[prot$gene] )

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
sum(duplicated(ngs$genes$gene_name))
x1 = apply( ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
ngs$genes = ngs$genes[match(rownames(x1), ngs$genes$gene_name),]
ngs$counts = x1
rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
remove(x1)

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=3)
head(ngs$samples)
table(ngs$samples$cluster)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
group.levels <- levels(ngs$samples$group)
group.levels
## 10 contrasts in total
contr.matrix <- makeContrasts(
    act12h_vs_notact = act12h - notact,
    act24h_vs_notact = act24h - notact,
    act48h_vs_notact = act48h - notact,
    act72h_vs_notact = act72h - notact,
    act96h_vs_notact = act96h - notact,
    act_vs_notact = (act96h + act72h + act48h + act24h + act12h)/5 - notact,
    act48h_vs_act12h = act48h - act12h,
    act72h_vs_act12h = act72h - act12h,
    act72h_vs_act48h = act72h - act48h,
    act96h_vs_act48h = act96h - act48h,
    act96h_vs_act72h = act96h - act72h,
    levels = group.levels)
contr.matrix
##contr.matrix = contr.matrix[,1:3]

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
ngs$timings <- c()

GENE.METHODS=c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","camera","fgsea")

GENE.METHODS=c("ttest","ttest.welch", ## "ttest.rank",
                   "voom.limma","trend.limma","notrend.limma",
                   "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                    "camera", "fry","fgsea") ## no GSEA, too slow...

MAX.GENES = 8000
MAX.GENESETS = 8000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features=MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENES,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("drugs-combo")
extra <- c("connectivity")
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

if(0) {
    sigdb = NULL
    sigdb = c("../libx/sigdb-lincs.h5","../libx/sigdb-virome.h5")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES, sigdb=sigdb) 
    names(ngs$connectivity)
}

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)

##===================================================================
##========================= END OF FILE =============================
##===================================================================
