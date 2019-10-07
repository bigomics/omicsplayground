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
source("../R/pgx-deconv.R")
source("../R/pgx-proteomics.R")
source("../R/pgx-drugs.R")

source("../scripts/options.R")
MAX.GENES

name="uploaded"
rda.file="../pgx/upload.pgx"
rda.file

if(0) {
    counts  = read.csv("../exampledata/counts.csv", row.names=1)
    samples = read.csv("../exampledata/samples.csv", row.names=1)
    genes   = read.csv("../exampledata/genes.csv", row.names=1)
    contrasts = read.csv("../exampledata/contrasts.csv", row.names=1)
}

pgx.upload <- function(counts, samples, genes, progress=NULL) {


    library(org.Hs.eg.db)
    
    ##load(file=rda.file, verbose=1)
    ngs <- list()  ## empty object
    ngs$name = "uploaded"
    ngs$date = date()
    ngs$datatype = "uploaded datatype"
    ngs$description = "uploaded description"

    if(0) {
        colnames(counts) == rownames(samples)
        rownames(counts) == rownames(genes)
    }
    
    ##-------------------------------------------------------------------
    ## create ngs object
    ##-------------------------------------------------------------------
    ngs$samples = samples
    ngs$counts  = counts
    ngs$genes   = data.frame(genes)
    ngs$contrasts = contrasts
    
    require(org.Hs.eg.db)
    GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    
    gene = rownames(genes)
    ngs$genes$gene_title = gene_title = GENE.TITLE[gene]
    ##ngs$genes$chr  = gene_title = GENE.CHR[gene]
    ##ngs$genes$pos  = gene_title = GENE.POS[gene]
    
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
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=3)
    head(ngs$samples)
    table(ngs$samples$cluster)
    

    ##======================================================================
    ##======================================================================
    ##======================================================================


    ## make model matrix for group vs. rest
    ##contr.matrix <- makeClusterContrasts(ngs$samples$cluster)
    contr.matrix <- ngs$contrasts
    
    rda.file
    ngs$timings <- c()
    
    USER.GENETEST.METHODS=c("ttest","ttest.welch","ttest.rank",
                            "voom.limma","trend.limma","notrend.limma",
                            "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    USER.GENESETTEST.METHODS = c("fisher","gsva","ssgsea","spearman",
                                 "camera", "fry","fgsea") ## no GSEA, too slow...
    USER.GENETEST.METHODS=c("ttest.welch","ttest.rank","trend.limma")
    USER.GENESETTEST.METHODS = c("fisher","gsva","spearman")

    if(0) {
        source("../R/compute-genes.R")
        source("../R/compute-genesets.R")
        source("../R/compute-extra.R")
    } else {

        ## new callling methods
        source("../R/compute2-genes.R")
        source("../R/compute2-genesets.R")
        source("../R/compute2-extra.R")

        ## ------------------ gene level tests ---------------------
        test.methods = c("trend.limma","ttest.welch","ttest")
        test.methods = USER.GENETEST.METHODS
        ngs <- compute.testGenes(
            ngs, contr.matrix, max.features=MAX.GENES,
            test.methods=test.methods)
        head(ngs$gx.meta$meta[[1]])        
        
        ## ------------------ gene set tests -----------------------
        test.methods = c("gsva","camera","fgsea")
        test.methods = USER.GENESETTEST.METHODS
        ngs <- compute.testGenesets(
            ngs, max.features=MAX.GENES,
            test.methods=test.methods)
        head(ngs$gset.meta$meta[[1]])

        ## ------------------ extra analyses ---------------------
        extra <- c("meta.go","deconv","infer","drugs")
        extra <- c("meta.go","infer","drugs")
        ngs <- compute.extra(ngs, extra)

    }

    ngs$timings
    return(ngs)
}

rda.file
ngs.save(ngs, file=rda.file)







