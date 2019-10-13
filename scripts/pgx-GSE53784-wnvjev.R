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
source("../R/ngs-functions.R")

source("options.R")

BATCH.CORRECT=1
rda.file="../data/GSE53784-wnvjev.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "RNA-seq"
ngs$description = "GSE53784 (Clarke et al., MBio 2014). Gene expression in the brain following WNV or JEV infection. WNV- or JEV-infected (N=3) vs. mock-infected (N=3) mouse brain."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    library(Biobase)
    library(GEOquery)
    library(limma)
    library(hgu133plus2.db)

    ## load series and platform data from GEO
    gset <- getGEO("GSE53784", GSEMatrix=TRUE, getGPL=TRUE)
    attr(gset, "names")
    X <- exprs(gset[[1]])
    head(X)
    
    ## extract GENE symbol from featureData
    gene.annot <- featureData(gset[[1]])@data$gene
    gene.symbol <- gsub("[ ]","",sapply(strsplit(gene.annot,split="//"),"[",2))
    gene.symbol[10000 + 1:10]    
    jj <- which( !gene.symbol %in% c(NA,"-",""))
    X <- X[jj,]
    rownames(X) <- gene.symbol[jj]
    
    ## Get sample info
    pdata = pData(gset)
    head(pdata)
    tt <- as.character(pdata$title)
    treatment <- sub("_.*","",tt)
    replicate <- sub(".*_","",tt)
    sampleTable <- data.frame(sample=tt, treatment=treatment)
    colnames(X) <- rownames(sampleTable) <- tt

    ## conform tables
    sample.names <- as.character(sampleTable$sample)
    rownames(sampleTable) = colnames(X) = sample.names
    
    ##-------------------------------------------------------------------
    ## gene annotation
    ##-------------------------------------------------------------------
    require(org.Mm.eg.db)
    GENE.TITLE = unlist(as.list(org.Mm.egGENENAME))
    gene.symbol = unlist(as.list(org.Mm.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    head(GENE.TITLE)
    gene_title <- GENE.TITLE[rownames(X)]

    ## get chromosome locations
    chrloc = as.list(org.Mm.egCHRLOC)
    names(chrloc) = gene.symbol
    chrloc <- chrloc[rownames(X)]
    loc <- sapply(chrloc, "[", 1)
    chrom <- sapply(chrloc, function(s) names(s)[1])
    loc[sapply(loc,is.null)] <- NA
    chrom[sapply(chrom,is.null)] <- NA
    chrom <- as.vector(unlist(chrom))
    loc   <- as.vector(unlist(loc))

    genes = data.frame( gene_name=rownames(X),
                       gene_title=gene_title,
                       chr=chrom, pos=loc)
    ##genes = apply(genes,2,as.character)
    head(genes)

    jj <- order(-apply(X,1,sd))
    X <- X[jj,]
    genes <- genes[jj,]
    
    jj <- which(!duplicated(genes$gene_name) & !is.na(genes$gene_name))
    X <- X[jj,]
    genes <- genes[jj,]
    rownames(X) <- rownames(genes) <- genes$gene_name
        
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    library(limma)
    X <- limma::normalizeQuantiles(X)
    ngs$counts <- 2**X  ## treat as counts
    ngs$samples <- sampleTable
    ngs$genes = genes
    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, perplexity=2, skipifexists=FALSE, prefix="C")
    head(ngs$samples)

    ##-------------------------------------------------------------------
    ## save
    ##-------------------------------------------------------------------
    
    rda.file
    save(ngs, file=rda.file)
}


if(DIFF.EXPRESSION) {

    load(file=rda.file, verbose=1)
    
    head(ngs$samples)
    ngs$samples$group <- ngs$samples$treatment
    levels = unique(ngs$samples$group)
    levels

    contr.matrix <- makeContrasts(
        JEV_vs_MOCK = JEV - MOCK,
        WNV_vs_Mock = WNV - Mock,
        levels = levels)
    contr.matrix
    
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")
    source("../R/compute-extra.R")
}

## save
rda.file
ngs.save(ngs, file=rda.file)












