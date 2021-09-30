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
FILES
MAX.GENES = 8000

COMPARE="group"
COMPARE="clusters"
COMPARE="pheno"

rda.file="../data/GSE110496-dengue.pgx"
##rda.file = sub(".pgx$",paste0("-vs",COMPARE,".pgx"),rda.file)
rda.file
##load(file=rda.file, verbose=1)

ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "scRNA-seq"
ngs$description ="GSE110496 dengue/zika virus scRNA-seq data set."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

library(Biobase)
library(GEOquery)
library(data.table)

##--------------------------------------------------------------
## Download from GEO
##--------------------------------------------------------------
geo <- getGEO("GSE110496", GSEMatrix=TRUE, AnnotGPL=TRUE)
length(geo)
attr(geo, "names")
xp <- exprs(geo[[1]])
dim(xp)
pdata = pData(geo[[1]])
dim(pdata)
head(pdata)

##--------------------------------------------------------------
## Download actual expression data from GEO FTP
##--------------------------------------------------------------
url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110496/suppl/GSE110496_RAW.tar"
system("(mkdir -p /tmp/GSE110496)")
system(paste0("(cd /tmp/GSE110496 && rm -f GSE110496_RAW.tar && wget ",url,")"))
system("(cd /tmp/GSE110496 && tar xvf GSE110496_RAW.tar)")
curwd <- getwd()
curwd
    
## concatenate file into one count matrix
library(parallel)
files <- dir("/tmp/GSE110496", pattern=".tsv.gz",full.names=TRUE)
head(files)
xx <- mclapply(files[], function(f) read.table(f,header=TRUE)[,2])
counts <- do.call(cbind, xx)
rownames(counts) <- read.table(files[1],header=TRUE)[,1]
colnames(counts) <- sub("_counts.*","",files[1:ncol(counts)])
colnames(counts) <- sub(".*/","",colnames(counts))

dim(counts)
head(counts)[,1:5]

## convert ENSEMBLE to symbol
require("EnsDb.Hsapiens.v86")
daf <- transcripts(EnsDb.Hsapiens.v86,
                   columns = c("gene_name", "gene_biotype", "gene_id"),
                   return.type="DataFrame")
head(summary(as.integer(table(as.character(daf$gene_id)))))
daf$symbol <- alias2hugo(daf$gene_name)
head(daf)
        
## do we have immune cell coding genes?
imm.genes <- grep("^IGH|^IGJ|^IGK|^IGL|^TRA[VJCD]|^TRB[VJCD]|^TRD[VJCD]|^TRG[VJCD]",
                  daf$symbol,value=TRUE)
imm.genes

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
head(rownames(counts))
sum(duplicated(rownames(counts)))
##x1 = apply(counts, 2, function(x) tapply(x, rownames(counts), sum))
gene <- daf$symbol[ match(rownames(counts), daf$gene_id) ]
x1 = tapply(1:nrow(counts), gene, function(i) colSums(counts[i,,drop=FALSE]))
x1 <- do.call(rbind, x1)
x1[1:3,1:3]
counts = x1
remove(x1)
dim(counts)

##--------------------------------------------------------------
## gene annotation
##--------------------------------------------------------------
require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)

biotype = daf$gene_biotype[ match(rownames(counts), daf$symbol) ]
gene_id = daf$gene_id[ match(rownames(counts), daf$symbol) ]
genes = data.frame( gene_id = gene_id,
                   gene_name = rownames(counts),
                   gene_title = GENE.TITLE[rownames(counts)],
                   gene_biotype = biotype )
rownames(genes) = rownames(counts)
head(genes)

##--------------------------------------------------------------
## Prepare sample table
##--------------------------------------------------------------

## get main sample annotation from rows 1-3
pdata = pData(geo[[1]])
dim(pdata)
head(pdata)
pdata <- pdata[,grep(":ch1$",colnames(pdata))]
colnames(pdata)
colnames(pdata) <- sub(":ch1$","",colnames(pdata))
colnames(pdata) <- sub("\\[h\\]$","",colnames(pdata))

sampleTable = pdata[,grep("molecule|time|virus|moi",colnames(pdata))]
rownames(sampleTable) = colnames(counts)
colnames(sampleTable) = sub("_molecules","",colnames(sampleTable))
apply(sampleTable, 2, table)

##-------------------------------------------------------------------
## Now create an PGX object
##-------------------------------------------------------------------
ngs$counts  <- round(counts)
ngs$samples <- sampleTable
ngs$genes = genes
##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
ngs$samples$treatment <- c("ctrl","infected")[1 + 1*(ngs$samples$moi>0)]
ngs$samples$batch <- NULL
##ngs$samples$batch <- as.integer(lib.size2)
ngs$samples$.counts <- colSums(counts)

## tagged rownames???
##row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])
row.id = ngs$genes[,"gene_name"]
rownames(ngs$genes) = rownames(ngs$counts) = row.id
names(ngs)

if(0) {
    jj <- sample(1:ncol(counts),400)
    ngs$counts  <- ngs$counts[,jj]
    ngs$samples <- ngs$samples[jj,]
}

##-------------------------------------------------------------------
## remove outliers (low/high counts)
##-------------------------------------------------------------------
hist(log10(1+colSums(ngs$counts)), breaks=100)
logcounts <- log10(1+colSums(ngs$counts))
zcounts <- (logcounts - mean(logcounts)) / sd(logcounts)
keep <- (abs(zcounts) <= 3 & ngs$samples$moi %in% c(0,1))
##  keep <- (abs(zcounts) <= 3)
colSums(counts)[which(!keep)]
ngs$counts <- ngs$counts[,keep]
ngs$samples <- ngs$samples[keep,]        
dim(counts)

##-------------------------------------------------------------------
## gene filtering
##-------------------------------------------------------------------
##keep <- rep(TRUE,nrow(ngs$counts))
##keep <- filterByExpr(ngs)  ## default edgeR filter
keep <- (rowMeans( edgeR::cpm(ngs$counts) > 1) >= 0.20)
##keep <- (rowMeans( ngs$counts >= 3) >= 0.01)
table(keep)
ngs$counts <- ngs$counts[keep,]
ngs$genes  <- ngs$genes[keep,]
dim(ngs$genes)

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, prefix="C")
head(ngs$samples)
    
##-------------------------------------------------------------------
## take top varying
##-------------------------------------------------------------------
MAX.GENES
if(FALSE && MAX.GENES>0) {
    cat("shrinking data matrices: n=",MAX.GENES,"\n")
    logcpm = edgeR::cpm(ngs$counts, log=TRUE)
    jj <- head( order(-apply(logcpm,1,sd)), MAX.GENES )  ## how many genes?
    head(jj)
    ngs$counts <- ngs$counts[jj,]
    ngs$genes  <- ngs$genes[jj,]
}


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

dim(ngs$counts)
ngs$timings <- c()

ngs$samples$group <- paste0(ngs$samples$virus,"_",ngs$samples$moi)
ngs$samples$group <- paste0(ngs$samples$virus,"_",ngs$samples$treatment,
                            "_",ngs$samples$time)   
table(ngs$samples$group)

COMPARE="groups"
##COMPARE="pheno"
##COMPARE="clusters"
if(COMPARE=="pheno") {
    ## use phenotype directly
    
    head(ngs$samples)
    ##table(ngs$samples$cell.type)
    table(ngs$samples$virus)
    table(ngs$samples$time)
    head(ngs$samples)
    
    ##contr.matrix <- makeDirectContrasts(
    ##    ngs$samples[,c("malignant","cluster","P2RX7","PDCD1","CD274","CD8A")],
    ##    ref=c("no","cl1","neg","neg","neg","neg") )
    contr.matrix <- makeDirectContrasts(
        Y = ngs$samples[,c("virus","time","moi")],
        ref = c("dengue","4","0") )
    head(contr.matrix)
    colnames(contr.matrix) <- sub(".*:","",colnames(contr.matrix))  ## strip prefix??
    head(contr.matrix)
    ##apply(contr.matrix,2,table)
    
} else if(COMPARE=="groups") {
    
    levels = unique(as.character(ngs$samples$group))
    levels
    table(ngs$samples$group)
    ##contr.matrix <- makeFullContrasts(levels)
    
    contr.matrix <- limma::makeContrasts(
        dengue4_vs_ctrl4 = dengue_infected_4 - dengue_ctrl_4,
        dengue12_vs_ctrl12 = dengue_infected_12 - dengue_ctrl_12,
        dengue24_vs_ctrl24 = dengue_infected_24 - dengue_ctrl_24,
        dengue48_vs_ctrl48 = dengue_infected_48 - dengue_ctrl_48,
        zika4_vs_ctrl4 = zika_infected_4 - zika_ctrl_4,
        zika12_vs_ctrl12 = zika_infected_12 - zika_ctrl_12,
        zika24_vs_ctrl24 = zika_infected_24 - zika_ctrl_24,
        zika48_vs_ctrl48 = zika_infected_48 - zika_ctrl_48,
        levels = levels)
    
    dim(contr.matrix)
    contr.matrix
    
    ##contr.matrix = contr.matrix[,1:3]
} else {
    stop("COMPARE error")
}


##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
head(contr.matrix)

GENE.METHODS=c("trend.limma","edger.qlf","edger.lrt")
GENESET.METHODS=c("gsva","camera","fgsea")

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

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## Save PGX object
##-------------------------------------------------------------------
rda.file
ngs$drugs$combo <- NULL  ## save space??
ngs.save(ngs, file=rda.file)



