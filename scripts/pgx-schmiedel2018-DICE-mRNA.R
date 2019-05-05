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
source("../R/pgx-graph.R")
source("../R/pgx-functions.R")

PROCESS.DATA=1
DIFF.EXPRESSION=1
COMPUTE.EXTRA=1

SMALL=8000
FAST=TRUE
EXT="8k"

if(0) {
    SMALL=8000
    FAST=FALSE
    EXT="8x"
}
rda.file="../pgx/schmiedel2018-DICE-mRNA.pgx"
if(SMALL>0) rda.file = sub(".pgx$",paste0("-",EXT,".pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*pgx/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "RNA-seq"
ngs$description ="DICE project (Schmiedel et al, Cell 2018). Transcriptomic data generated from 13 immune cell types isolated from 106 leukapheresis samples provided by 91 healthy subjects. The cell types  include three innate immune cell types (CD14high CD16— classical monocytes, CD14— CD16+ non-classical monocytes, CD56dim CD16+ NK cells), four adaptive immune cell types (naive B cells, naive CD4+ T cells, naive CD8+ T cells, and naive regulatory T cells [TREG]), six CD4+ memory or more differentiated T cell subsets (TH1, TH1/17, TH17, TH2, follicular helper T cell [TFH], and memory TREG), and two activated cell types (naive CD4+ and CD8+ T cells that were stimulated ex vivo)."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    library(Biobase)
    library(GEOquery)
    library(limma)
    library(data.table)

    ## load series and platform data from GEO
    csv.files <- dir("~/Projects/Data/DICE", pattern=".csv", full.names=TRUE)
    csv.names <- dir("~/Projects/Data/DICE", pattern=".csv", full.names=FALSE)
    csv.names

    csv.data <- lapply(csv.files, fread)
    head(csv.data[[1]])[,1:20]
    head(csv.data[[2]])[,1:20]

    geneInfo <- csv.data[[1]][,1:3]
    csv.data <- lapply(csv.data, function(x) x[,4:ncol(x)])
    head(csv.data[[1]])[,1:20]

    cells <- sub("_TPM.csv","",csv.names)
    cell.type <- c()
    for(i in 1:length(csv.data)) {
        colnames(csv.data[[i]]) <- paste0(cells[i],"_",colnames(csv.data[[i]]))
        cell.type <- c(cell.type, rep(cells[i], ncol(csv.data[[i]])))
    }
    table(cell.type)

    ##------------------------ make object --------------------------------

    ngs <- c()
    csv.data <- lapply(csv.data, as.data.frame)
    ngs$counts <- as.matrix(do.call(cbind, csv.data))
    dim(ngs$counts)

    ngs$samples = data.frame(
        sample = colnames(ngs$counts),
        cell.type = cell.type
    )
    rownames(ngs$samples) <- colnames(ngs$counts)

    ## for statistics
    ngs$samples$group <- factor(ngs$samples$cell.type)
    ##ngs$design.matrix = model.matrix( ~ ngs$samples$group )

    ##-------------------------------------------------------------------
    ## gene annotation
    ##-------------------------------------------------------------------
    colnames(geneInfo) <- c("feature","transcript_length","annotation")
    gene_name <- sub(";.*","",geneInfo$annotation)
    gene_biotype <- sub(".*;","",geneInfo$annotation)
    ngs$genes <- data.frame(
        gene_name=gene_name,
        gene_biotype=gene_biotype,
        geneInfo)
    head(ngs$genes)

    ## check duplicated genes
    sum(duplicated(ngs$genes$feature))
    rownames(ngs$genes)  <- ngs$genes$feature
    rownames(ngs$counts) <- ngs$genes$feature
    dim(ngs$counts)
    sum(duplicated(ngs$genes$gene_name))

    ## take out duplicates and convert to gene symbole
    gene <- as.character(ngs$genes$gene_name)
    jj <- order(-rowSums(ngs$counts))
    jj <- jj[which(!duplicated(gene[jj]))]
    length(jj)
    ngs$counts <- ngs$counts[jj,]
    ngs$genes  <- ngs$genes[jj,]
    rownames(ngs$counts) <- ngs$genes$gene_name
    rownames(ngs$genes ) <- ngs$genes$gene_name

    ## gene annotation
    require(org.Hs.eg.db)
    GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    head(GENE.TITLE)
    ngs$genes$gene_title <- GENE.TITLE[ngs$genes$gene_name]

    ##-------------------------------------------------------------------
    ## gene filtering
    ##-------------------------------------------------------------------
    ##keep <- rep(TRUE,nrow(ngs$counts))
    ##keep <- filterByExpr(ngs)  ## default edgeR filter
    logcounts <- edgeR::cpm(ngs$counts, log=TRUE)
    keep1 <- ( rowSums(logcounts >= 3) >=3 & apply(logcounts,1,sd)>0)
    keep2 <- grepl("protein_coding|^TR|^IG", ngs$genes$gene_biotype)
    ##keep <- (keep1 & keep2)
    keep <- (keep2)
    table(keep)

    ngs$counts <- ngs$counts[keep,]
    ngs$genes  <- ngs$genes[keep,]

    ##-------------------------------------------------------------------
    ## take out empty samples
    ##-------------------------------------------------------------------
    sdx <- apply(log(1+ngs$counts),2,sd)
    keep.samples <- (sdx > 0)
    table(keep.samples)
    ngs$counts <- ngs$counts[,keep.samples]
    ngs$samples <- ngs$samples[colnames(ngs$counts),]

    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE)
    head(ngs$samples)

    ##-------------------------------------------------------------------
    ## take top varying
    ##-------------------------------------------------------------------

    if(FALSE && SMALL>0) {
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

if(1) {

    ## check
    head(colnames(ngs$counts))
    head(rownames(ngs$samples))
    if(!all(rownames(ngs$samples) == colnames(ngs$counts))) {
        stop("sample names do not match")
    }

    head(rownames(ngs$counts))
    head(rownames(ngs$genes))
    if(!all(rownames(ngs$counts) == rownames(ngs$genes))) {
        stop("gene names do not match")
    }

    must.geneannot <- c("gene_name","gene_title")
    if(!all(must.geneannot %in% colnames(ngs$gene))) {
        stop("missing columns in ngs$genes")
    }

    must.sampleannot <- c("group")
    if(!all(must.sampleannot %in% colnames(ngs$samples))) {
        stop("missing columns in ngs$samples")
    }

}

if(DIFF.EXPRESSION) {

    load(file=rda.file, verbose=1)
    levels <- sort(unique(ngs$samples$group))
    contr.matrix <- makeFullContrasts(levels)
    dim(contr.matrix)
    colnames(contr.matrix)
    source("../R/compute-testgenes.R")
    source("../R/compute-testgenesets.R")
    source("../R/compute-extra.R")
    save(ngs, file=rda.file)
}

rda.file
cat(">>> output file is:",rda.file,"\n")
ngs.save(ngs, file=rda.file)




