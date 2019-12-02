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
source("../R/ngs-functions.R")

source("options.R")
MAX.GENES

## run all available methods 
USER.GENETEST.METHODS = "*"
USER.GENESETTEST.METHODS = c("gsva","fisher","camera","fgsea","fry","spearman")

rda.file="../data/GSE73024-debio1437.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE73024 data set (Nakanishi et al, 2015). Transcriptome analysis of response to a selective FGFR inhibitor (CH5183284/Debio 1347), a mitogen-activated protein kinase kinase (MEK) inhibitor, or a phosphoinositide 3-kinase (PI3K) inhibitor."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    if(0) {
        BiocManager::install("hgu133plus2.db")
        BiocManager::install("GEOquery")
    }
    library(Biobase)
    library(GEOquery)
    library(limma)
    library(hgu133plus2.db)
    
    ## load series and platform data from GEO
    gset <- getGEO("GSE73024", GSEMatrix =TRUE, AnnotGPL=TRUE)
    length(gset)
    attr(gset, "names")

    pdata = pData(gset[[1]])
    head(pdata)
    clinvar <- pdata[,grep(":ch1$",colnames(pdata))]
    head(clinvar)

    tt  <- pdata$title
    src <- pdata$source_name_ch1
    sampleTable <- data.frame(sample=tt, cell.line=src, clinvar)
    colnames(sampleTable) <- gsub("[.]ch1$","",colnames(sampleTable))
    colnames(sampleTable) <- gsub("compound[.]","",colnames(sampleTable))
    colnames(sampleTable) <- gsub("treatment[.]","",colnames(sampleTable))
    head(sampleTable)
    
    ## extract data
    X = exprs(gset[[1]])
    max(X)  ## check max for checking logarithm
    sdx <- apply(log2(X),1,sd)
    X <- X[order(-sdx),]
    
    ## convert affymetrix ID to GENE symbol
    affx  = sapply(as.list(hgu133plus2SYMBOL),"[[",1)
    symbol = affx[rownames(X)]
    symbol = alias2hugo(symbol)
    rownames(X) = symbol
    X = X[which(!duplicated(symbol) & !is.na(symbol) & symbol!=""),]
    dim(X)
    X = X[which(rowMeans(is.na(X))==0), ]  ## no missing values
    sum(is.na(X))
    dim(X)

    ## conform tables
    table(rownames(sampleTable) == colnames(X))

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
    chrloc = as.list(org.Hs.egCHRLOC)
    names(chrloc) = gene.symbol
    chrloc <- chrloc[rownames(X)]
    loc <- abs(sapply(chrloc, "[", 1))
    chrom <- sapply(chrloc, function(s) names(s)[1])
    chrom[sapply(chrom,is.null)] <- NA
    chrom <- as.vector(unlist(chrom))

    genes = data.frame( gene_name=rownames(X),
                       gene_title=gene_title,
                       chr=chrom, pos=loc)
    ##genes = apply(genes,2,as.character)
    head(genes)
    rownames(genes) = rownames(X)

    ##--------------------------------------------------------------------
    ## check if batch correction is needed
    ##--------------------------------------------------------------------
    BATCH.CORRECT
    if(BATCH.CORRECT) {
        require(sva)
        batch <- sampleTable$Chemotherapy
        design = model.matrix( ~ as.character(sampleTable$sample))
        ##bX = ComBat(X, batch=as.character(batch), mod=design)
        bX = removeBatchEffect(X, batch=as.character(batch),
                               ##batch2=as.character(sampleTable$gender),
                               design=design)
        X = normalizeQuantiles(bX)
    }
    
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    ngs$counts <- as.matrix(2**X)  ## treat as counts
    ngs$samples <- sampleTable
    ngs$genes = genes
    ##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
    ##ngs$samples$batch <- NULL  ##???
    ##ngs$samples$batch <- as.integer(lib.size2)

    ## tagged rownames
    ##row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])  
    row.id = ngs$genes[,"gene_name"]
    rownames(ngs$genes) = rownames(ngs$counts) = row.id
    names(ngs)

    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    sum(duplicated(ngs$genes$gene_name))

    ## x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
    ## ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    ## ngs$counts = x1
    ## dim(x1)
    ## rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    ## remove(x1)
    ngs <- ngs.collapseByGene(ngs)
    dim(ngs$counts)
        
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, prefix="C")
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
    ngs$samples$group <- as.character(ngs$samples$dlbcl.type)
    levels = unique(ngs$samples$group)
    levels

    contr.matrix <- makeContrasts(
        ABC_vs_GCB = ABC - GCB,
        GCB_vs_ABC = GCB - ABC,
        levels = c("ABC","GCB"))
    
    contr.matrix <- makeDirectContrasts(
        Y=ngs$samples[,c("dlbcl.type","gender")],
        ref=c("GCB","male"))
    dim(contr.matrix)
    head(contr.matrix)
    colnames(contr.matrix) <- sub(".*:","",colnames(contr.matrix))  ## strip prefix 
    head(contr.matrix)
    
    ##contr.matrix = contr.matrix[,1:3]
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")
    source("../R/compute-extra.R")
}

## save
rda.file
ngs.save(ngs, file=rda.file)








