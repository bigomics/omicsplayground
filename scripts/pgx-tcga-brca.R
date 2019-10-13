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
MAX.GENES = 4000
BATCH.CORRECT=TRUE

## run all available methods 
USER.GENETEST.METHODS = c("trend.limma","edger.qlf","deseq2.wald")
USER.GENESETTEST.METHODS = c("gsva","fisher","camera","fgsea","fry","spearman")
USER.GENESETTEST.METHODS = c("gsva","camera","fgsea")

rda.file="../data/tcga-brca_pub-gx.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "RNA-seq"
##ngs$datatype = "multi-omics"
ngs$description = "TCGA breast cancer data set. Gene expression from 526 patients annotated with PAM50 classification. Data from cBioPortal."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ## get data

    ## Uncomment to download data from cBioPortal
    ## system("mkdir -p ../downloads/brca_tcga_pub")
    ## system("wget http://download.cbioportal.org/brca_tcga_pub.tar.gz -P ../downloads/brca_tcga_pub/")
    ## system("(cd ../downloads/brca_tcga_pub/ && tar xvfz brca_tcga_pub.tar.gz)")
    download.dir = "../downloads/brca_tcga_pub"
    dir(download.dir)
    
    clin1 <- read.csv(file.path(download.dir,"data_clinical_sample.txt"),
                      skip=4, sep="\t", stringsAsFactors=FALSE)
    clin2 <- read.csv(file.path(download.dir,"data_clinical_patient.txt"), 
                      skip=4, sep="\t", stringsAsFactors=FALSE)
    dim(clin1)
    dim(clin2)
    colnames(clin1)
    colnames(clin2)
    pat.id1 <- sub("-01$","",clin1$PATIENT_ID)
    clin2 <- clin2[match(pat.id1, clin2$PATIENT_ID),]
    table( pat.id1 == clin2$PATIENT_ID)
    clin <- cbind(clin1, clin2[,-1])
    rownames(clin) <- clin1$PATIENT_ID
    sel <- c(
        ##"GLEASON_PATTERN_PRIMARY","GLEASON_PATTERN_SECONDARY",
        "PAM50_SUBTYPE","TUMOR_STAGE","NODES", 
        "ER_STATUS","PR_STATUS", "HER2_STATUS",
        "METASTASIS","OS_STATUS", "OS_MONTHS"
    )
    clin <- clin[,sel]
    clin$PAM50_SUBTYPE <- gsub("[- ]|like|enriched","",clin$PAM50_SUBTYPE)
    table(clin$PAM50_SUBTYPE)
    clin$ER_STATUS[!grepl("Pos|Neg",clin$ER_STATUS)] <- NA
    clin$PR_STATUS[!grepl("Pos|Neg",clin$PR_STATUS)] <- NA
    clin$HER2_STATUS[!grepl("Pos|Neg",clin$HER2_STATUS)] <- NA        
    
    X.data   <- read.csv(file.path(download.dir,"data_mRNA_median_Zscores.txt"),
                         sep="\t", check.names=FALSE)
    head(colnames(X.data))
    X <- as.matrix(X.data[,3:ncol(X.data)])
    gene <- X.data$Hugo_Symbol
    rownames(X) <- gene
    sum(duplicated(X.data$Hugo_Symbol))
    ##X <- apply(X, 2, function(x) tapply(x,gene,sum))
    X <- X[rowMeans(is.na(X))<0.5, ]
    dim(X)
    X <- imputeMedian(X)
    sum(is.na(X))
    
    ##colnames(X) <- sub("-01$","",colnames(X))        
    samples <- sort(intersect(colnames(X),rownames(clin)))
    samples <- samples[!is.na(clin[samples,"PAM50_SUBTYPE"])]
    length(samples)
    X <- X[,samples]
    sampleTable <- clin[samples,]        
    
    ## scale back to counts...
    counts <- 10*edgeR::cpm(2**X, log=FALSE)  ## 10 million
    summary(colSums(counts))
    
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
    
    ##dtype <- gsub("\\[|\\].*","",rownames(X))
    ##table(dtype)
    genes = data.frame(
        ##data_type = dtype,
        gene_name=gene,
        gene_title=gene_title,
        chr=chrom)
    ##genes = apply(genes,2,as.character)
    rownames(genes) = rownames(X)
    head(genes)
    
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    ngs$counts  <- as.matrix(counts)  ## treat as is
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
                              perplexity=30)
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
    ##ngs$samples$group <- as.character(ngs$samples$HER2_STATUS)
    pam50 <- as.character(ngs$samples$PAM50_SUBTYPE)
    ngs$samples$group <- pam50
    table(ngs$samples$group)
    levels = setdiff(unique(ngs$samples$group),c("",NA))
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






