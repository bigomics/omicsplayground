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

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/tcga-prad.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
ngs$description = "TCGA prostate cancer data set. Gene expression from patients with Gleason score. Data from cBioPortal."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

## Please download from https://amp.pharm.mssm.edu/archs4/download.html
TCGA_MATRIX = "../libx/tcga_matrix.h5"
if(file.exists(TCGA_MATRIX)) {
    
    tcga <- pgx.getTCGAdataset(study="prad_tcga", genes=NULL,
                               matrix_file=TCGA_MATRIX, from.h5=TRUE)
    names(tcga)
    X <- tcga$X[[1]]
    clin <- tcga$clin[[1]]
    
} else {
    ## get data from cBioportal
    ##
    download.dir = "/tmp/prad_tcga"
    if(!dir.exists(download.dir)) {
        system(paste("mkdir -p ",download.dir))
        system(paste("wget -c http://download.cbioportal.org/prad_tcga.tar.gz -P ",download.dir))
        system(paste("(cd ",download.dir," && tar xvfz prad_tcga.tar.gz)"))
    } else {
        ## should exists..
    }
    getwd()
    dir(download.dir)
    
    clin <- read.csv(file.path(download.dir,"data_bcr_clinical_data_patient.txt"),
                     skip=4, sep="\t")
    rownames(clin) <- clin$PATIENT_ID
    
    X.data   <- read.csv(file.path(download.dir,"data_RNA_Seq_v2_expression_median.txt"),
                         sep="\t", check.names=FALSE)
    X <- X.data[,3:ncol(X.data)]
    gene <- X.data$Hugo_Symbol
    sum(duplicated(X.data$Hugo_Symbol))
    X <- apply(X, 2, function(x) tapply(x,gene,sum))
    dim(X)
    colnames(X) <- sub("-01$","",colnames(X))
    
}

sel1 <- c(
    ##"GLEASON_PATTERN_PRIMARY","GLEASON_PATTERN_SECONDARY",
    "GLEASON_SCORE",
    "TUMOR_STATUS",
    "BIOCHEMICAL_RECURRENCE_INDICATOR",
    "RADIATION_TREATMENT_ADJUVANT",
    "TREATMENT_OUTCOME_FIRST_COURSE",
    "CLIN_T_STAGE")
sel2 <- c("OS_MONTHS","OS_STATUS")
clin1 <- apply(clin[,sel1], 2, function(x) {x[grep("^\\[",x)]=NA;x})
##clin <- data.frame(cbind(clin1, clin[,sel2]))
clin <- data.frame(clin1)
clin$CLIN_T_STAGE <- sub("[abc]","",clin$CLIN_T_STAGE)  ## simplify
clin$GLEASON_SCORE <- gsub("[ ]","",paste0("G",clin$GLEASON_SCORE))  ## simplify

samples <- sort(intersect(colnames(X),rownames(clin)))
X <- X[,samples]
sampleTable <- clin[samples,]        

##-------------------------------------------------------------------
## gene annotation
##-------------------------------------------------------------------
require(org.Hs.eg.db)
GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
gene <- sub("\\[.*\\]","",rownames(X))
gene_title <- GENE.TITLE[gene]

## get chromosome locations
chrloc = as.list(org.Hs.egCHRLOC)
names(chrloc) = gene.symbol
chrloc <- chrloc[gene]
##loc <- abs(as.numeric(unlist(sapply(chrloc, "[", 1))))
chrom <- sapply(chrloc, function(s) names(s)[1])
chrom[sapply(chrom,is.null)] <- NA
chrom <- as.vector(unlist(chrom))

dtype <- gsub("\\[|\\].*","",rownames(X))
genes = data.frame(
    ##data_type = "mrna",
    gene_name=gene,
    gene_title=gene_title,
    chr=chrom)
##genes = apply(genes,2,as.character)
head(genes)
rownames(genes) = rownames(X)

##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
ngs$counts  <- as.matrix(X)  ## treat as is
ngs$samples <- sampleTable
ngs$genes   <- genes
##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
##ngs$samples$batch <- NULL  ##???
##ngs$samples$batch <- as.integer(lib.size2)

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
sum(duplicated(rownames(X)))
##ngs <- ngs.collapseByGene(ngs)
##dim(ngs$counts)

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs$tsne2d=ngs$tsne3d=NULL
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, prefix="C",
                          perplexity=30)
head(ngs$samples)

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
head(ngs$samples)
##score <- paste0("G",as.character(ngs$samples$GLEASON_SCORE))
score <- as.character(ngs$samples$GLEASON_SCORE)
##score <- as.character(ngs$samples$CLIN_T_STAGE)
score <- gsub("[ ]","",score)
score[is.na(score)] <- "NA"
table(score)
ngs$samples$group <- score

levels = setdiff(unique(ngs$samples$group),NA)
levels
colnames(ngs$samples)

## contr.matrix <- limma::makeContrasts(
##     T2_vs_T1 = T2 - T1,
##     T3_vs_T1 = T3 - T1,
##     T4_vs_T1 = T4 - T1,
##     T3_vs_T2 = T3 - T2,
##     T4_vs_T3 = T4 - T3,
##     levels = levels)

contr.matrix <- limma::makeContrasts(
    G10_vs_G9 = G10 - G9,
    G9_vs_G8 = G9 - G8,
    G8_vs_G7 = G8 - G7,
    G7_vs_G6 = G7 - G6,
    
    G10_vs_G6 = G10 - G6,
    G9_vs_G6 = G9 - G6,
    G8_vs_G6 = G8 - G6,
    
    levels = levels)

dim(contr.matrix)
head(contr.matrix)


##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

GENE.METHODS=c("ttest","ttest.welch","ttest.rank",
               "voom.limma","trend.limma","notrend.limma",
               "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                    "camera", "fry","fgsea") ## no GSEA, too slow...
GENE.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...
MAX.GENES = 20000
MAX.GENESETS = 5000
    
## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("drugs-combo")
##extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
extra <- c("meta.go","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
names(ngs$drugs)
ngs$timings    


##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)

##===================================================================
##========================= END OF FILE =============================
##===================================================================






