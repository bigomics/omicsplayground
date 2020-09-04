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
FILES

DOWNSAMPLE=75
COMPARE.CLUSTERS=FALSE
##COMPARE.CLUSTERS=TRUE
FILTER.GENES=TRUE

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE92332-scintestine.pgx"
##if(COMPARE.CLUSTERS) rda.file <- sub(".pgx$",paste0("-vsCLUST.pgx"),rda.file)
##if(DOWNSAMPLE>0) rda.file <- sub(".pgx$",paste0("-s",DOWNSAMPLE,".pgx"),rda.file)
##if(SMALL>0) rda.file <- sub(".pgx$",paste0("-",EXT,".pgx"),rda.file)
rda.file


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "scRNA-seq"
ngs$organism = "mouse"
ngs$description = "GSE92332 data set. A single-cell survey of the small intestinal epithelium (Haber et al., Cell 2017)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(data.table)

##--------------------------------------------------------------
## Read SC counts
##--------------------------------------------------------------
## load series and platform data from GEO
gse.file = "/tmp/GSE92332_atlas_UMIcounts.txt.gz"   ## 10x data???
if(!file.exists(gse.file)) {
    system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_atlas_UMIcounts.txt.gz -P /tmp")
}
gse = fread(gse.file,nrow=-1000)
dim(gse)
head(gse)[,1:10]
head(colnames(gse),10)
counts = as.matrix(gse[,2:ncol(gse)])
rownames(counts) = as.character(gse[[1]])
head(counts)[,1:10]
summary(colSums(counts))

##-------------------------------------------------------------------
## Dowsample samples
##-------------------------------------------------------------------
DOWNSAMPLE
if(DOWNSAMPLE>0) {
    ## sample each category
    tissue <- sub(".*_","",colnames(counts))
    table(tissue)
    jj <- order(-colSums(counts))
    jj <- tapply(jj, tissue[jj], function(x) head(x,DOWNSAMPLE))
    ##jj <- tapply(jj, tissue[jj], function(x) head(x,20))
    jj <- unlist(jj)
    table(tissue[jj])
    counts <- counts[,jj]
    dim(counts)
    summary(colSums(counts))
}

##--------------------------------------------------------------
## gene annotation
##--------------------------------------------------------------
if(!require("org.Mm.eg.db")) {
    BiocManager::install("org.Mm.eg.db", update=FALSE,ask=FALSE)
}
require(org.Mm.eg.db)
gene <- as.character(gse[[1]])
genes <- select(org.Mm.eg.db, gene, c("ENTREZID","GENENAME","CHR"), "SYMBOL")
colnames(genes) <- c("gene_name","entrez_id","gene_title","chr")
head(genes)
dim(genes)
sum(duplicated(genes$gene_name))
head(genes)
genes <- genes[match(gene, genes$gene_name),]
dim(genes)
head(genes)

## add gene type information
##require("EnsDb.Hsapiens.v86")
if(!require("EnsDb.Mmusculus.v79"))  {
    BiocManager::install("EnsDb.Mmusculus.v79", update=FALSE,ask=FALSE)
}
require("EnsDb.Mmusculus.v79")
daf <- transcripts(EnsDb.Mmusculus.v79,
                   columns = c("gene_name", "gene_biotype"),
                   return.type="DataFrame")
head(daf)
genes$gene_biotype <- daf$gene_biotype[match(genes$gene_name,daf$gene_name)]
table(genes$gene_biotype)
head(genes)

##--------------------------------------------------------------
## Prepare sample table
##--------------------------------------------------------------

## geo <- getGEO("GSE92332", GSEMatrix =TRUE, AnnotGPL=TRUE)
## str(geo[[1]])
## pdata = pData(geo[[1]])
## head(pdata)
## table(pdata$source_name_ch1)
## table(pdata$"tissue:ch1")
## table(pdata$"treatment:ch1")
tissue <- sub(".*_","",colnames(counts))
batch  <- sub("_.*","",colnames(counts))
sampleTable <- data.frame(batch=batch, tissue=tissue)
rownames(sampleTable) <- colnames(counts)
head(sampleTable)

##-------------------------------------------------------------------
## Now create an PGX object
##-------------------------------------------------------------------
ngs$counts <- round(counts)
ngs$samples <- sampleTable
ngs$genes = genes
##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
##ngs$samples$batch <- NULL
##ngs$samples$batch <- as.integer(lib.size2)

## tagged rownames???
row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])
rownames(ngs$genes) = rownames(ngs$counts) = row.id
names(ngs)

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
sum(duplicated(ngs$genes$gene_name))
dim(ngs$counts)
ngs <- ngs.collapseByGene(ngs)
head(rownames(ngs$counts))
head(rownames(ngs$genes))

##-------------------------------------------------------------------
## gene filtering
##-------------------------------------------------------------------
sort(table(ngs$genes$gene_biotype))

FILTER.GENES
if(FILTER.GENES) {
    is.sex <- (ngs$genes$chr %in% c("chrX","chrY","X","Y",23,24))
    prot.coding <- (ngs$genes$gene_biotype == "protein_coding" &
                    !is.na(ngs$genes$gene_biotype) )
    imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
    keep <- (prot.coding | imm.gene)
    table(ngs$genes$chr[which(keep)])
    sort(table(ngs$genes$gene_biotype[which(keep)]))
    ##keep <- ( keep & !is.sex & prot.coding)
    ##keep <- (keep & !is.sex & (prot.coding | imm.gene))
    ##keep <- (keep & rowMeans(ngs$counts >= 3) >= 0.001)
    table(keep)
    ngs$counts <- ngs$counts[keep,]
    ngs$genes  <- ngs$genes[keep,]
}

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE,
                          perplexity=30, kclust=1)
table(ngs$samples$cluster)

##-------------------------------------------------------------------
## take top varying
##-------------------------------------------------------------------
MAX.GENES
dim(ngs$counts)
if(TRUE && MAX.GENES>0) {
    cat("shrinking data matrices: n=",MAX.GENES,"\n")
    logcpm = edgeR::cpm(ngs$counts, log=TRUE)
    jj <- head( order(-apply(logcpm,1,sd)), MAX.GENES )  ## how many genes?
    head(jj)
    ##bX <- bX[jj,]
    ngs$counts <- ngs$counts[jj,]
    ngs$genes  <- ngs$genes[jj,]
}
dim(ngs$counts)
dim(ngs$genes)
ngs$timings <- c()


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

ngs$samples$group <- as.character(ngs$samples$tissue)
table(ngs$samples$cluster, ngs$samples$group)

## ----------------- test genes ------------------------------------------
## COMPARE.CLUSTERS
## clusters <- ngs$samples$tissue
## if(COMPARE.CLUSTERS) {
##     clusters <- ngs$samples$cluster
## }
## ## make model matrix for group vs. rest
## table(clusters)
## contr.matrix <- makeClusterContrasts(clusters, full=FALSE)
## contr.matrix
## ##contr.matrix = contr.matrix[,1:3]

## combine TISSUE and CLUSTER contrast in single matrix
tissue.contr  <- makeClusterContrasts(ngs$samples$tissue, full=FALSE, by.sample=TRUE)
cluster.contr <- makeClusterContrasts(ngs$samples$cluster, full=FALSE, by.sample=TRUE)
contr.matrix <- cbind(tissue.contr, cluster.contr)
rownames(contr.matrix) <- rownames(ngs$samples)
contr.matrix <- normalizeContrasts(contr.matrix)
head(contr.matrix)


##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

ngs$timings <- c()

GENE.METHODS=c("ttest.welch","trend.limma","edger.qlf")
GENESET.METHODS = c("fisher","gsva","camera","fgsea")

MAX.GENES = 8000
MAX.GENESETS = 8000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features=MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("drugs-combo")
extra <- c("connectivity")
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

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


