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
source("../R/ngs-functions.R")
source("../R/gset-fisher.r")
source("../R/gset-gsea.r")
source("../R/gset-meta.r")
source("../R/pgx-drugs.R")
source("../R/pgx-graph.R")
source("../R/pgx-functions.R")

source("options.R")
MAX.GENES
MAX.GENES = 2000
BATCH.CORRECT=TRUE

## run all available methods 
USER.GENETEST.METHODS = c("trend.limma","edger.qlf","deseq2.wald")
USER.GENESETTEST.METHODS = c("gsva","fisher","camera","fgsea","fry","spearman")

rda.file="../pgx/tcga-brca_pub-GX.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "RNA-seq"
ngs$datatype = "multi-omics"
ngs$description = "BRCA-TCGA public data set (from cBioPortal)."


## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ## Old OMX data
    ##

    ##load("~/OMX/shiny-omx/omxdata/tcga/brca_tcga_pub-omx.rda",verbose=1)
    load("~/Projects/Data/tcga-omx/brca_tcga_pub-omx.rda",verbose=1)
    names(omx)
    names(omx$level[[1]]$mat)
    lapply(omx$level[[1]]$mat,dim)
    mat <- omx$level[[1]]$mat[1]
    ##mat <- omx$level[[1]]$mat[1:5]
    names(mat)

    ## impute missing values
    sapply(mat,function(x) sum(is.na(x)))
    mat <- lapply(mat, imputeMedian)
    sapply(mat,function(x) sum(is.na(x)))
    
    if(length(mat)==1) {
        X <- 2**mat[["gx"]]
        samples <- colnames(mat[["gx"]])
    } else {
        ## decide necessary transformation
        ##mat[["gx"]] <- 2**mat[["gx"]]  ## assume counts
        ##mat[["px"]] <- 2**mat[["px"]]  ## assume counts
        ##mat[["cn"]] <- (mat[["cn"]])  ## ?
        ##mat[["me"]] <- (mat[["me"]])  ## ?
        ##mat[["mt"]] <- (mat[["mt"]])  ## ?
                
        samples <- Reduce(intersect, sapply(mat,colnames))
        ##samples <- Reduce(union, sapply(mat[1:5],colnames))
        samples
        names(mat)
        lapply(mat,dim)
        for(i in 1:length(mat)) {
            mat[[i]] <- mat[[i]][,match(samples,colnames(mat[[i]]))]
            colnames(mat[[i]]) <- samples
            prefix <- paste0("[",names(mat)[i],"]")
            rownames(mat[[i]]) <- paste0(prefix,rownames(mat[[i]]))
        }
        head(rownames(mat[[2]]))        
        X <- do.call(rbind, mat)
        dim(X)
    }
    
    samples <- sort(intersect(samples,rownames(omx$pheno)))
    colnames(omx$pheno) <- sub("HER2_FINAL_STATUS","HER2_STATUS",colnames(omx$pheno))
    sampleTable <- omx$pheno[samples,c("PR_STATUS","ER_STATUS","HER2_STATUS")]
    
    pam50 <- omx$pheno[samples,grep("PAM50",colnames(omx$pheno))]
    pam50 <- sub("PAM50_SUBTYPE:","",colnames(pam50)[max.col(pam50)])
    pam50 <- gsub("[-]|enriched|like","",pam50)
    sampleTable <- cbind(sampleTable, PAM50=pam50)
    head(sampleTable)
    X <- X[,samples]
    
    ##-------------------------------------------------------------------
    ## gene annotation
    ##-------------------------------------------------------------------
    require(org.Hs.eg.db)
    GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    head(GENE.TITLE)
    gene <- sub("\\[.*\\]","",rownames(X))
    gene_title <- GENE.TITLE[gene]

    ## get chromosome locations
    chrloc = as.list(org.Hs.egCHRLOC)
    names(chrloc) = gene.symbol
    chrloc <- chrloc[gene]
    ##loc <- abs(as.numeric(unlist(sapply(chrloc, "[", 1))))
    chrom <- sapply(chrloc, function(s) names(s)[1])
    chrom[sapply(chrom,is.null)] <- NA
    chrom <- as.vector(unlist(chrom))
    
    dtype <- gsub("\\[|\\].*","",rownames(X))
    table(dtype)
    genes = data.frame(
        data_type = dtype,
        gene_name=gene,
        gene_title=gene_title,
        chr=chrom)
    ##genes = apply(genes,2,as.character)
    head(genes)
    rownames(genes) = rownames(X)
    
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    ngs$counts  <- as.matrix(X)  ## treat as is
    ngs$samples <- sampleTable
    ngs$genes   <- genes
    ##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
    ##ngs$samples$batch <- NULL  ##???
    ##ngs$samples$batch <- as.integer(lib.size2)

    ## tagged rownames
    ##row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])  
    ##row.id = ngs$genes[,"gene_name"]
    ##rownames(ngs$genes) = rownames(ngs$counts) = row.id
    ##names(ngs)

    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    sum(duplicated(rownames(X)))
    ## x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
    ## ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
    ## ngs$counts = x1
    ## dim(x1)
    ## rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
    ## remove(x1)
    ##ngs <- ngs.collapseByGene(ngs)
    ##dim(ngs$counts)
        
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs$X <- as.matrix(X) ## cluster will skip log
    ngs$X <- NULL
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, prefix="C",
                              perplexity=5)
    head(ngs$samples)
    
    ##-------------------------------------------------------------------
    ## save
    ##-------------------------------------------------------------------    
    rda.file
    save(ngs, file=rda.file)
}


if(DIFF.EXPRESSION) {
    rda.file
    load(file=rda.file, verbose=1)
    
    head(ngs$samples)
    ngs$samples$group <- as.character(ngs$samples$HER2_STATUS)
    ngs$samples$group <- as.character(ngs$samples$PAM50)
    table(ngs$samples$group)
    levels = unique(ngs$samples$group)
    levels
    colnames(ngs$samples)
    
    contr.matrix <- makeContrasts(
        LuminalA_vs_normal = LuminalA - Normal,
        LuminalB_vs_normal = LuminalB - Normal,
        LuminalB_vs_LuminalA = LuminalB - LuminalA,
        Basal_vs_normal = Basal - Normal,
        HER2_vs_normal = HER2 - Normal,
        levels=levels)
    dim(contr.matrix)
    head(contr.matrix)

    ## some methods (yet) cannot handle direct contrasts!!!
    ## contr.matrix <- makeDirectContrasts(
    ##     ngs$samples[,c("PR_STATUS","ER_STATUS","HER2_STATUS")],
    ##     ref = c("0","0","0"))
    ## dim(contr.matrix)
    ## head(contr.matrix)
    ##contr.matrix = contr.matrix[,1:3]
    
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")    
    source("../R/compute-extra.R")
    
}

## save
rda.file
ngs.save(ngs, file=rda.file)


if(0) {
    load(file=rda.file, verbose=1)
    table(ngs$genes[rownames(ngs$X),]$data_type)
    jj <- match(rownames(ngs$X),rownames(ngs$genes))
    table(ngs$genes[jj,]$data_type)

    
}





