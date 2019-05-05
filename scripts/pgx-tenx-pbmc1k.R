##rm(list=setdiff(ls(),run.param))
library(knitr)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(matrixTests)
library(kableExtra)
library(knitr)

source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/ngs-cook.r")
source("../R/ngs-fit.r")
source("../R/gset-fisher.r")
source("../R/gset-gsea.r")
source("../R/gset-meta.r")
source("../R/pgx-functions.R")
source("../R/ngs-cook.r")

##FAST=TRUE
SMALL=8000
FAST=FALSE
EXT="8x"
if(1) {
    SMALL=8000
    FAST=TRUE
    EXT="8k"
}
rda.file="../pgx/tenx-pbmc1k.pgx"
##rda.file="./sallusto2019-th1star.pgx"
if(SMALL>0) rda.file = sub(".pgx$",paste0("-",EXT,".pgx"),rda.file)

PROCESS.DATA=1
DIFF.EXPRESSION=1
COMPUTE.EXTRA=1

SAMPLECORRECT=FALSE
EXCLUDE.SAMPLES = NULL
DOWNSAMPLE = 2000
##EXCLUDE.SAMPLES = c("20-Th1star-ut-r3")

rda.file
##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "scRNA-seq"
ngs$description = "10x example data."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## load 10x results
    tenx <- readRDS(file="../../Sallusto_10x/pmbc1k_counts.rds")
    names(tenx)
    rownames(tenx$samples) <- tenx$samples$barcode
    rownames(tenx$genes)   <- rownames(tenx$counts)
    tenx$samples$barcode <- NULL

    ##--------------------------------------------------------------
    ## gene annotation
    ##--------------------------------------------------------------
    require(org.Hs.eg.db)
    GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    CHR = sapply(as.list(org.Hs.egCHR),"[",1)  ## some have multiple chroms..
    names(CHR) = gene.symbol
    table(CHR)

    gene <- as.character(tenx$genes$gene_name)
    tenx$genes$gene_title <- GENE.TITLE[gene]
    tenx$genes$chr <- CHR[gene]
    head(tenx$genes)

    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    ngs$counts  <- round(tenx$counts)
    ngs$samples <- tenx$samples
    ngs$genes   <-  tenx$genes
    ##ngs$samples$clone = paste(sampleTable$cell.type, sampleTable$clone,sep="")
    ##ngs$samples$lib.size <- colSums(tenx$counts / 1e6)  ## summed intensity as lib.size
    ##ngs$samples$batch <- NULL
    ##ngs$samples$batch <- as.integer(lib.size2)

    ##-------------------------------------------------------------------
    ## sample QC filtering
    ##-------------------------------------------------------------------
    if(!is.null(EXCLUDE.SAMPLES)) {
        intersect(EXCLUDE.SAMPLES, colnames(ngs$counts))
        kk <- which(!(colnames(ngs$counts) %in% EXCLUDE.SAMPLES))
        ngs$counts <- ngs$counts[,kk]
        ngs$samples <- ngs$samples[colnames(ngs$counts),]
    }
    dim(ngs$counts)

    if(DOWNSAMPLE > 0 && ncol(ngs$counts)>DOWNSAMPLE ) {
        kk <- sample(ncol(ngs$counts),DOWNSAMPLE)
        ngs$counts <- ngs$counts[,kk]
        ngs$samples <- ngs$samples[colnames(ngs$counts),]
    }
    dim(ngs$counts)

    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    sum(duplicated(ngs$genes$gene_name))
    sum(is.na(ngs$counts))
    x1 = apply( ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum, na.rm=TRUE))
    ngs$counts = as.matrix(x1)
    ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    remove(x1)

    ##-------------------------------------------------------------------
    ## gene filtering
    ##-------------------------------------------------------------------
    summary(colSums(ngs$counts))
    dim(ngs$counts)
    sum(is.na(ngs$counts))
    ##keep <- rep(TRUE,nrow(ngs$counts))
    ##keep <- filterByExpr(ngs)  ## default edgeR filter
    ##cpm <- edgeR::cpm(ngs$counts)
    keep <- (rowMeans(ngs$counts >= 3) >= 0.01)
    table(keep)
    ##keep <- ( keep & !(ngs$genes$chr %in% c("chrX","chrY")) &
    ##          ngs$genes$gene_biotype == "protein_coding")
    ##sort(table( ngs$genes$gene_biotype))
    keep <- ( keep & ngs$genes$molecule_type == "Gene")
    table(keep)
    ngs$counts <- ngs$counts[keep,]
    ngs$genes  <- ngs$genes[keep,]

    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(
        ngs, skipifexists=FALSE, perplexity=30,
        sv.rank=200, gamma=1)
    head(ngs$samples)

    if(0) {
        head(ngs$samples)
        table(ngs$samples$cluster)
        klr <- rainbow(6)[factor(ngs$samples$cluster)]
        plot(ngs$tsne2d, col=klr)
    }

    ##-------------------------------------------------------------------
    ## take top varying
    ##-------------------------------------------------------------------
    SMALL
    if(SMALL>0) {
        cat("shrinking data matrices: n=",SMALL,"\n")
        logcpm = edgeR::cpm(ngs$counts, log=TRUE)
        jj <- head( order(-apply(logcpm,1,sd)), SMALL )  ## how many genes?
        head(jj)
        ##bX <- bX[jj,]
        ngs$counts <- ngs$counts[jj,]
        ngs$genes  <- ngs$genes[jj,]
    }
    rda.file
    save(ngs, file=rda.file)

}


if(DIFF.EXPRESSION) {
    load(file=rda.file, verbose=1)

    ## make model matrix for group vs. rest
    ngs$samples$group <- ngs$samples$cluster
    clusters <- ngs$samples$cluster
    table(clusters)
    contr.matrix <- makeClusterContrasts(clusters, full=TRUE)

    ## use phenotype directly
    ##contr.matrix <- makeDirectContrasts(
    ##    Y=ngs$samples[,c("treatment","cluster","cluster")],
    ##    ref=c("ut",NA,""))
    contr.matrix

    ## -------- test genes
    FAST
    source("../R/pgx-testgenes.R")
    source("../R/pgx-testgenesets.R")
    save(ngs, file=rda.file)
}

if(COMPUTE.EXTRA) {
    ## -------- extra stuff
    load(file=rda.file, verbose=1)
    FILES = "../files/"
    source("../R/pgx-extra.R")
    save(ngs, file=rda.file)
    
}

rda.file
cat(">>> output file is:",rda.file,"\n")
ngs.save(ngs, file=rda.file)





