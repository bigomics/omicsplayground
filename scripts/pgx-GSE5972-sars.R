## TO BE FINISHED
##
##

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")

rda.file="../data-extra/GSE5972-sars.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE5972. Gene expression profiling of patients with severe acute respiratory syndrome (SARS). Interferon-mediated immunopathological events are associated with atypical innate and adaptive immune responses in patients with severe acute respiratory syndrome (Cameron, J Virol 2007)."

PROCESS.DATA = TRUE
DIFF.EXPRESSION = TRUE

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    library(Biobase)
    library(GEOquery)
    library(limma)
    ##library(hgu133plus2.db)

    ## load series and platform data from GEO
    geo <- getGEO("GSE5972", GSEMatrix=TRUE, getGPL=TRUE)
    attr(geo, "names")
    length(geo)
    
    X <- exprs(geo[[1]])
    dim(X)
    max(X)

    ## extract GENE symbol from featureData
    colnames(featureData(geo[[1]])@data)
    gene.symbol <- as.character(featureData(geo[[1]])@data$GENE_SYMBOL)
    ##gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
    gene.symbol[10000 + 1:10]    
    jj <- which( !gene.symbol %in% c(NA,"-",""))
    X <- X[jj,]
    rownames(X) <- gene.symbol[jj]

    ## take out duplicated
    sum(duplicated(rownames(X)))
    X1 <- tapply(1:nrow(X), rownames(X), function(i) colSums(2**X[i,,drop=FALSE],na.rm=TRUE))
    X <- log2(do.call(rbind, X1))
    head(X)[,1:4]
    remove(X1)
    
    ## Get sample info from title
    pdata = pData(geo[[1]])
    head(pdata)
    title <- as.character(unlist(sapply(geo, function(d) pData(d)$title)))
    gsm <- unlist(sapply(geo, function(d) pData(d)$geo_accession))
    names(title) <- gsm

    tt.pheno <- as.character(unlist(sapply(geo, function(d) pData(d)$characteristics_ch2)))
    
    ## clean up titles.... :(
    head(title)
    tt <- sub("lowMOI","_lowMOI",title)
    tt <- sub("HighMOI","_HighMOI",tt)
    tt <- sub("hrs$","hr_0",tt)
    tt <- sub("MOCK_MRC5_([1-3])","MOCK_MRC5_mock_0hr_\\1",tt)
    tt <- sub("MRC5_","",tt)
    head(tt)

    sampleTable <- do.call(rbind,strsplit(tt,split="_"))
    head(sampleTable)
    colnames(sampleTable) <- c("infected","treatment","time","replicate")
    colnames(X) <- rownames(sampleTable) <- tt
    head(sampleTable)
    
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
    rownames(genes) <- rownames(X)
    head(genes)
        
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    library(limma)
    ##X <- limma::normalizeQuantiles(X)
    ngs$counts <- X  ## treat as counts
    ngs$samples <- data.frame(sampleTable)
    ngs$genes = genes
    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, perplexity=2, skipifexists=FALSE, prefix="C")
    head(ngs$samples)

}


if(DIFF.EXPRESSION) {

    ##load(file=rda.file, verbose=1)
    
    head(ngs$samples)
    grp <- paste(ngs$samples$infected,ngs$samples$time,sep="_")
    ngs$samples$group <- grp
    levels = unique(ngs$samples$group)
    levels
    
    contr.matrix <- makeContrasts(
                
        MERS_24hr_vs_MOCK_24hr = MERS_24hr - MOCK_24hr,
        MERS_48hr_vs_MOCK_48hr = MERS_48hr - MOCK_48hr,

        SARS_24hr_vs_MOCK_24hr = SARS_24hr - MOCK_24hr,
        SARS_48hr_vs_MOCK_48hr = SARS_48hr - MOCK_48hr,

        MOCK_24hr_vs_MOCK_0hr = MOCK_24hr - MOCK_0hr,
        MOCK_48hr_vs_MOCK_0hr = MOCK_48hr - MOCK_0hr,
        
        levels = levels)
    contr.matrix
    
    rda.file
    ngs$timings <- c()
    
    GENETEST.METHODS=c("ttest","ttest.welch","ttest.rank",
                       "voom.limma","trend.limma","notrend.limma",
                       "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                        "camera", "fry","fgsea") ## no GSEA, too slow...
    GENETEST.METHODS = c("trend.limma","edger.qlf","deseq2.wald")
    GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...

    MAX.GENES = 20000
    MAX.GENESETS = 5000
    
    ## new callling methods
    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = MAX.GENES,
        test.methods = GENETEST.METHODS)
    
    ngs <- compute.testGenesets (
        ngs, max.features=MAX.GENESETS,
        test.methods = GENESET.METHODS,
        lib.dir=FILES)

    extra <- c("connectivity")
    extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES) 
    
    names(ngs)
    ngs$timings


}

## save
rda.file
ngs.save(ngs, file=rda.file)












