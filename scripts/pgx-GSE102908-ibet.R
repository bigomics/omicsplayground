RDIR = "../R"
FILES = "../lib"
FILESX = "../libx"
PGX.DIR = "../data"
source("../R/pgx-include.R")
source("../R/pgx-getgeo.R")
##source("options.R")

rda.file="../data/GSE102908-ibet.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE102908. RNAseq data from SCCOHT1 and OVCAR8 ovarian cancer cells treated with BET inhibitors. Analysis of mRNA profile of 2 cell lines exposed to DMSO, OTX015 for 4 and 24 hours in duplicate."

PROCESS.DATA = TRUE
DIFF.EXPRESSION = TRUE

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    library(Biobase)
    library(GEOquery)
    library(limma)

    ## load series and platform data from GEO
    ngs <- pgx.getGEOseries("GSE102908")
    names(ngs)
    
    contr.matrix <- ngs$contrasts
    colnames(contr.matrix)[1] = "treatment:OTX015_vs_DMSO"
    ct1 <- makeContrasts(
        "treatment@SCCOHT1:OTX015_vs_DMSO" =
            ( OTX015100nM_SCCOHT1_NA + OTX015100nM_SCCOHT1_4h )/2 -
            DMSO001_SCCOHT1_control,
        "treatment@OVCAR8:OTX015_vs_DMSO" =
            ( OTX015100nM_OVCAR8_NA + OTX015100nM_OVCAR8_4h )/2 -
            DMSO001_OVCAR8_control,        
        levels = rownames(contr.matrix)
    )
    contr.matrix <- cbind(contr.matrix, ct1)
    colnames(contr.matrix)
        
    ## Pre-calculate t-SNE for and get clusters
    ngs <- pgx.clusterSamples(ngs, perplexity=NULL, skipifexists=FALSE, prefix="C")
    head(ngs$samples)
    
    rda.file
    ngs$timings <- c()
    
    GENETEST.METHODS=c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
    GENESET.METHODS = c("fisher","gsva","fgsea","camera") ## no GSEA, too slow...

    MAX.GENES = 20000
    MAX.GENESETS = 5000
    
    ## new callling methods
    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = MAX.GENES,
        test.methods = GENETEST.METHODS)
    
    ngs <- compute.testGenesets (
        ngs, max.features = MAX.GENESETS,
        test.methods = GENESET.METHODS,
        lib.dir=FILES)
    
    extra <- c("connectivity")
    extra <- c("drugs-combo")
    extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
    extra <- c("meta.go","deconv","infer","drugs-combo","wordcloud","connectivity")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES) 
    
    names(ngs)
    ngs$timings
    
}

## save
rda.file
ngs.save(ngs, file=rda.file)

if(0) {
    
    source("../R/pgx-include.R")
    extra <- c("connectivity")
    sigdb = NULL
    sigdb = c("../libx/sigdb-l1000.h5")
    sigdb = "../libx/sigdb-lincs-cp2.h5"
    sigdb = c("../libx/sigdb-lincs-cp.h5","../libx/sigdb-lincs-gt.h5")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES, sigdb=sigdb) 
    names(ngs$connectivity)
    
    rda.file
    ngs.save(ngs, file=rda.file)
    
    cp <- ngs$connectivity[["sigdb-lincs-cp.h5"]]
    head(cp[[6]])

}













