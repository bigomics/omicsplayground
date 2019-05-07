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

source("options.R")
rda.file="../pgx/GSE88808-prostate.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*pgx/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "mRNA array"
ngs$description = "GSE88808 data set. Gleason-score matched tumor and adjacent normal samples were collected to compare gene expression differences in early-onset versus late-onset prostate cancer patients (Ding, PLOS Genet 2016)."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ################################################################
    ## Differential expression analysis with limma
    ## BiocManager::install("GEOquery", version = "3.8")
    library(Biobase)
    library(GEOquery)
    library(data.table)

    ##--------------------------------------------------------------
    ## Read SC counts
    ##--------------------------------------------------------------
    ## load series and platform data from GEO
    geo <- getGEO("GSE88808", GSEMatrix =TRUE, AnnotGPL=FALSE)
    length(geo)
    attr(geo, "names")
    counts = 2**exprs(geo[[1]])
    head(counts)[,1:10]
    summary(colSums(counts))

    ##--------------------------------------------------------------
    ## gene annotation
    ##--------------------------------------------------------------
    require(org.Hs.eg.db)
    GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    head(GENE.TITLE)
    gene_title <- GENE.TITLE[rownames(counts)]
    genes = data.frame( gene_name=rownames(counts), gene_title=gene_title)
    rownames(genes) <- rownames(counts)

    ##--------------------------------------------------------------
    ## Prepare sample table
    ##--------------------------------------------------------------
    ##geo <- getGEO("GSE99795", GSEMatrix =TRUE, AnnotGPL=TRUE)
    pdata = pData(geo[[1]])
    head(pdata)
    sampleTable <- pdata[,grep(":ch1$",colnames(pdata))]
    colnames(sampleTable) <- sub(":ch1","",colnames(sampleTable))
    head(sampleTable)

    age <- as.numeric(sampleTable$age)
    summary(age)
    sampleTable$age.group <- c("young","old")[ 1 + 1*(age >= median(age)) ]
    head(sampleTable)

    ##-------------------------------------------------------------------
    ## Now create an PGX object
    ##-------------------------------------------------------------------
    ngs$counts <- round(counts)
    ngs$samples <- sampleTable
    ngs$genes = data.frame(genes)
    rownames(ngs$genes) <- genes$gene_name

    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE,
                              perplexity=30, kclust=2)
    head(ngs$samples)

    rda.file
    save(ngs, file=rda.file)
}


if(DIFF.EXPRESSION) {
    load(file=rda.file, verbose=1)

    ##sampleTable$group <- paste0(sampleTable$Stage,"_",sampleTable$tissue)
    ngs$samples$group <- paste0(ngs$samples$age.group,"_",ngs$samples$tissue)
    table(ngs$samples$group)
    levels = unique(ngs$samples$group)
    levels

    contr.matrix <- makeContrasts(
        old_tumor_vs_old_normal = old_tumor - old_normal,
        young_tumor_vs_young_normal = young_tumor - young_normal,
        old_normal_vs_young_normal = old_normal - young_normal,
        old_tumor_vs_young_tumor = old_tumor - young_tumor,
        levels = levels)

    dim(contr.matrix)
    contr.matrix

    USER.GENETEST.METHODS=c("trend.limma","deseq2.wald","edger.qlf")
    ## USER.GENETEST.METHODS=c("trend.limma","deseq2","edger.qlf")
    USER.GENESETTEST.METHODS=c("fisher","gsva","camera","fgsea")
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")
    source("../R/compute-extra.R")

}

rda.file
ngs.save(ngs, file=rda.file)



