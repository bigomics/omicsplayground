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

BATCH.CORRECT=1
rda.file="../data/GSE32591-lupusnephritis.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE32591 data set (Berthier et al, 2011). transcriptome of microdissected renal biopsies from patients with lupus nephritis (LN)"

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ##   Differential expression analysis with limma
    library(Biobase)
    library(GEOquery)
    library(limma)
    library(hgu133plus2.db)

    ## load series and platform data from GEO
    gset <- getGEO("GSE32591", GSEMatrix =TRUE, getGPL=FALSE)
    attr(gset, "names")
    gset <- gset[[1]]

    pdata = pData(gset)
    head(pdata)
    clinvar <- pdata[,grep(":ch1$",colnames(pdata))]
    colnames(clinvar) <- sub(":ch1$","",colnames(clinvar))
    head(clinvar)
    colnames(clinvar) <- c("disease","patient","tissue")
    clinvar$disease <- sub("LN patient","LN",clinvar$disease)
    clinvar$tissue  <- sub(" from kidney biopsy","",clinvar$tissue)
    clinvar$tissue  <- sub("Tubulointerstitium","Tubulo",clinvar$tissue)
    clinvar$patient  <- paste0("pt",clinvar$patient)  
    head(clinvar)
    sampleTable <- clinvar
    
    ## merge data sets
    X = exprs(gset)
    max(X)  ## check max for checking logarithm
    X = X[order(-apply(X,1,sd)),]
    dim(X)
    head(X)[,1:4]
    
    ## convert affymetrix ID to GENE symbol
    symbol <- unlist(as.list(org.Hs.egSYMBOL))
    symbol = alias2hugo(symbol[rownames(X)])
    rownames(X) = symbol
    X = X[which(!duplicated(symbol) & !is.na(symbol) & symbol!=""),]
    X = X[which(rowMeans(is.na(X))==0), ]  ## no missing values
    sum(is.na(X))
    dim(X)

    ## conform tables
    table(rownames(sampleTable) == colnames(X))
    
    ##-------------------------------------------------------------------
    ## sample QC filtering
    ##-------------------------------------------------------------------
    pd <- table(sampleTable$disease,sampleTable$patient)
    pd
    is.paired <- sampleTable$patient %in% names(which(colSums(pd==0)==0))
    sampleTable <- sampleTable[is.paired,]
    X  <- X[,rownames(sampleTable)]
    dim(X)

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
        batch <- sampleTable$patient
        design = model.matrix( ~ as.character(sampleTable$disease))
        ##bX = ComBat(X, batch=as.character(batch), mod=design)
        bX = removeBatchEffect(X, batch=as.character(batch),
                               ##batch2=as.character(sampleTable$gender),
                               design=design)
        X = normalizeQuantiles(bX)
    }
    
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    ngs$counts <- (2**X)  ## treat as counts
    ngs$samples <- sampleTable
    ngs$genes = genes

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
    ngs <- pgx.clusterSamples(ngs, perplexity=5, skipifexists=FALSE, prefix="C")
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
    ngs$samples$group <- paste0(ngs$samples$tissue,"_",ngs$samples$disease)
    levels = unique(ngs$samples$group)
    levels

    contr.matrix <- makeContrasts(
        Tubulo_LN_vs_control = Tubulo_LN - Tubulo_control,
        Glomeruli_LN_vs_control = Glomeruli_LN - Glomeruli_control,
        LN_vs_control = (Tubulo_LN + Glomeruli_LN) -
            (Tubulo_control + Glomeruli_control),
        levels = levels)
    ##contr.matrix <- makeDirectContrasts(
    ##    Y=ngs$samples[,c("dlbcl.type","gender","cluster")],
    ##    ref=c("GCB","male",NA))
    contr.matrix

    USER.GENETEST.METHODS=c("trend.limma","deseq2.wald","edger.qlf")
    USER.GENESETTEST.METHODS=c("gsva","fisher","camera","fgsea")
    ##USER.GENETEST.METHODS="*"
    ##USER.GENESETTEST.METHODS=c("gsva","fisher","camera","fgsea","fry","spearman")

    ##contr.matrix = contr.matrix[,1:3]
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")
    source("../R/compute-extra.R")
}

## save
rda.file
ngs.save(ngs, file=rda.file)












