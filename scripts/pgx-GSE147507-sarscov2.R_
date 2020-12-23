##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE147507-sarscov2.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
ngs$description = "GSE147507. Transcriptional response of human lung epithelial cells to SARS-CoV-2 infection (Blanco-Melo et al, 2020)."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
library(Biobase)
library(GEOquery)
library(limma)
library(hgu133plus2.db)

## load series and platform data from GEO
geo <- getGEO("GSE147507", GSEMatrix=TRUE, getGPL=TRUE)
attr(geo, "names")

## get the counts table
system("wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147507/suppl/GSE147507_RawReadCounts.tsv.gz -O /tmp/GSE147507_RawReadCounts.tsv.gz")
X <- read.csv("/tmp/GSE147507_RawReadCounts.tsv.gz",
              sep="\t",row.names=1,check.names=FALSE)
X <- as.matrix(X)

if(0) {
    require(Rtsne)
    X <- X*1.0001
    res <- Rtsne( t(head(X,200)), perplexity=6)
    res <- Rtsne( t(head(X,800)), perplexity=6)        
}
    
## Get sample info
pdata = pData(geo[[1]])
head(pdata)
sampleTable <- pdata[,grep(":ch1",colnames(pdata))]
head(sampleTable)
##colnames(sampleTable) <- c("code","infected","time","replicate")
colnames(sampleTable) <- gsub("[ ]","_",sub(":ch1$","",colnames(sampleTable)))
sampleTable$time_point <- gsub(" after treatment|rs","",sampleTable$time)
sampleTable$treatment <- gsub("[ ]treatment|[ ]infected.*","",sampleTable$treatment)
sampleTable$cell_type <- NULL
##ngs$samples$induce_IFNB <- sub("YES","IFNB",ngs$samples$induce_IFNB)
head(sampleTable)

sampleTable$group <- with(sampleTable, paste(treatment,cell_line,time_point,sep="_"))
sampleTable$group <- gsub("[-]","",sampleTable$group)
sampleTable$group[c(13:14)] <- paste0(sampleTable$group[c(13:14)],"_rsv")

##colnames(X) <- rownames(sampleTable) <- tt

## match counts and sampleTable
idx <- gsub("[-]","_",sub("_S[0-9]*","",pdata$description))
colnames(X) <- gsub("[-]","_",colnames(X))
idx %in% colnames(X)
X <- X[,match(idx,colnames(X))]
rownames(sampleTable) <- colnames(X) <- sub("_index.*","",idx)

##-------------------------------------------------------------------
## gene annotation
##-------------------------------------------------------------------
require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
gene_title <- GENE.TITLE[rownames(X)]

## get chromosome locations
chrloc = sapply(as.list(org.Hs.egMAP),"[",1)
names(chrloc) = gene.symbol
chrloc <- chrloc[rownames(X)]

genes = data.frame( gene_name=rownames(X),
                   gene_title=gene_title,
                   chr=chrloc)
##genes = apply(genes,2,as.character)
rownames(genes) <- rownames(X)  ## !
head(genes)

## take out duplicated
sum(duplicated(rownames(X)))

##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
library(limma)
ngs$counts <- as.matrix(X)  ## treat as counts
ngs$samples <- data.frame(sampleTable)
ngs$genes = genes

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters 
##-------------------------------------------------------------------
dim(ngs$counts)
ngs <- pgx.clusterSamples(ngs, method="pca", perplexity=NULL,
                          skipifexists=FALSE, prefix="C")
head(ngs$samples)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
##load(file=rda.file, verbose=1)    
head(ngs$samples)
levels = unique(ngs$samples$group)
levels

contr.matrix <- makeContrasts(
    SARSCoV2NHBE_vs_Mock = SARSCoV2_NHBE_24h - Mock_NHBE_24h,
    SARSCoV2A549_vs_Mock = SARSCoV2_A549_24h - Mock_A549_24h,                
    RSVA549_vs_Mock = RSV_A549_24h - Mock_A549_24h_rsv,
    IAVA549_vs_Mock = IAV_A549_9h - Mock_A549_9h,                
    levels = levels)
contr.matrix
    
##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

rda.file
ngs$timings <- c()

GENE.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","fgsea","camera") ## no GSEA, too slow...

MAX.GENES = 20000
MAX.GENESETS = 10000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features = MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("connectivity")
extra <- c("deconv")
extra <- c("meta.go","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)












