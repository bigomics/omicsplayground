rm(list=setdiff(ls(),run.param))
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
source("../R/xcr-graph.r")
source("../R/pgx-functions.R")

FILES="../files/"
RDIR="../R/"

PROCESS.DATA=1
DIFF.EXPRESSION=1
TEST.GENESETS=1
EXTRA.STUFF=1
COMPARE.CLUSTERS=FALSE
COMPARE.CLUSTERS=TRUE
##DOWNSAMPLE=100

SMALL=8000
FAST=TRUE
EXT="8k"

if(0) {
    SMALL=2000
    FAST=TRUE
    EXT="2k"
}

rda.file="../pgx/GSE99795-prostateSC.pgx"
if(COMPARE.CLUSTERS) rda.file <- sub(".pgx$",paste0("-vsCLUST.pgx"),rda.file)
##if(DOWNSAMPLE>0) rda.file <- sub(".pgx$",paste0("-s",DOWNSAMPLE,".pgx"),rda.file)
if(SMALL>0) rda.file <- sub(".pgx$",paste0("-",EXT,".pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*pgx/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "scRNA-seq"
ngs$description = "GSE99795 data set. Single-cell RNA-seq reveals a subpopulation of prostate cancer cells with enhanced cell cycle-related transcription and attenuated androgen response (Horning et al., 2018)."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    ##BiocManager::install("GEOquery", version = "3.8")
    library(Biobase)
    library(GEOquery)
    library(data.table)

    ##--------------------------------------------------------------
    ## Read SC counts
    ##--------------------------------------------------------------
    ## load series and platform data from GEO
    ##system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE99nnn/GSE99795/suppl/GSE99795_raw.txt.gz; mv GSE99795* ../data/GSE")
    gse = fread("../data/GSE/GSE99795_raw.txt.gz",nrow=-1000)
    ##counts = fread("../downloads/GSE98638_HCC.TCell.S5063.count.txt.gz",nrow=-1000)
    dim(gse)
    head(gse)[,1:10]
    head(colnames(gse),10)
    counts = as.matrix(gse[,9:ncol(gse)])
    head(counts)[,1:10]
    summary(colSums(counts))

    ##--------------------------------------------------------------
    ## gene annotation
    ##--------------------------------------------------------------
    require(org.Hs.eg.db)
    refseq.id <- as.character(gse[[1]])
    genes <- select(org.Hs.eg.db, refseq.id, c("ENTREZID","SYMBOL","GENENAME","CHR"), "REFSEQ")
    colnames(genes) <- c("refseq_id","entrez_id","gene_name","gene_title","chr")
    head(genes)
    dim(genes)
    sum(duplicated(genes$gene_name))
    head(genes)
    genes <- genes[match(refseq.id, genes$refseq_id),]
    dim(genes)

    ## add gene type information
    require("EnsDb.Hsapiens.v86")
    daf <- transcripts(EnsDb.Hsapiens.v86,
                       columns = c("gene_name", "gene_biotype"),
                       return.type="DataFrame")
    head(daf)
    genes$gene_biotype <- daf$gene_biotype[match(genes$gene_name,daf$gene_name)]
    table(genes$gene_biotype)
    head(genes)

    ##--------------------------------------------------------------
    ## Prepare sample table
    ##--------------------------------------------------------------

    geo <- getGEO("GSE99795", GSEMatrix =TRUE, AnnotGPL=TRUE)
    pdata = pData(geo[[1]])
    head(pdata)
    clinvar <- pdata[,grep(":ch1$",colnames(pdata))]
    head(clinvar)

    sampleTable <- data.frame(clinvar[,0])
    ##sampleTable$title <- as.character(pdata$title)
    sampleTable$sample <- gsub(".*\\[|\\]$","",as.character(pdata$title))
    tp <- sub(" hour untreated ","hr_UT_",as.character(clinvar$"time point:ch1"))
    tp <- sub(" hour androgen-treated ","hr_AG_",tp)
    tp <- gsub("_bulk|_cell","",tp)
    table(tp)
    sampleTable$time_point <- tp
    ##sampleTable$cell_type <- as.character(clinvar$"cell type:ch1")
    sampleTable$sample_type <- sub(".*[ ]","",as.character(as.character(clinvar$"time point:ch1")))
    sampleTable$group <- paste0(sampleTable$sample_type,"_",sampleTable$time_point)

    head(sampleTable)
    count.names <- sub("\\/.*","",colnames(counts))
    all( count.names %in% sampleTable$sample)
    sampleTable <- sampleTable[match(count.names, sampleTable$sample),]
    colnames(counts) <- rownames(sampleTable) <- count.names
    head(sampleTable)
    sampleTable$sample <- NULL

    table(sampleTable$time_point)

    ##-------------------------------------------------------------------
    ## Now create an PGX object
    ##-------------------------------------------------------------------
    if(is.null(sampleTable$group)) stop("samples need group")
    table(sampleTable$group)
    ngs$counts <- round(counts)
    ngs$samples <- sampleTable
    ngs$genes = genes
    ##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
    ngs$samples$batch <- NULL
    ##ngs$samples$batch <- as.integer(lib.size2)

    ## tagged rownames???
    row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])
    rownames(ngs$genes) = rownames(ngs$counts) = row.id
    names(ngs)

    ##-------------------------------------------------------------------
    ## sample QC filtering
    ##-------------------------------------------------------------------
    ##

    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    sum(duplicated(ngs$genes$gene_name))
    x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
    dim(x1)
    ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    ngs$counts = x1
    rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    remove(x1)

    ##-------------------------------------------------------------------
    ## gene filtering
    ##-------------------------------------------------------------------
    sort(table(ngs$genes$gene_biotype))

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
                              perplexity=30, kclust=2)
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
    dim(ngs$counts)
    ngs$timings <- c()

    rda.file
    save(ngs, file=rda.file)
}


if(DIFF.EXPRESSION) {
    load(file=rda.file, verbose=1)

    table(ngs$samples$cluster, ngs$samples$group)

    ## ----------------- test genes ------------------------------------------
    COMPARE.CLUSTERS
    if(COMPARE.CLUSTERS) {
        ## make model matrix for group vs. rest
        table(ngs$samples$cluster)
        clusters <- ngs$samples$cluster
        table(clusters)
        contr.matrix <- makeClusterContrasts(clusters)
        contr.matrix

    } else {
        table(ngs$samples$group)
        levels = unique(ngs$samples$group)
        levels

        contr.matrix <- makeContrasts(
            UT_12hr_vs_UT_0hr = cell_12hr_UT - cell_0hr_UT,
            AG_12hr_vs_UT_0hr = cell_12hr_AG - cell_0hr_UT,
            AG_12hr_vs_UT_12hr = cell_12hr_AG - cell_12hr_UT,
            bulk_vs_cell = (bulk_0hr_UT + bulk_12hr_UT + bulk_12hr_AG)/3 -
                (cell_0hr_UT + cell_12hr_UT +  cell_12hr_AG) / 3,
            levels = levels)

        dim(contr.matrix)
        contr.matrix

        ##contr.matrix = contr.matrix[,1:3]
    }

    ## USER.GENETEST.METHODS=c("trend.limma","deseq2","edger.qlf")
    USER.GENESETTEST.METHODS=c("fisher","gsva","camera","fgsea")
    source("../R/compute-testgenes.R")
    source("../R/compute-testgenesets.R")
    source("../R/compute-extra.R")

}

rda.file
ngs.save(ngs, file=rda.file)



