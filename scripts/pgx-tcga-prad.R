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
USER.GENESETTEST.METHODS = c("gsva","camera","fgsea")

rda.file="../pgx/tcga-prad-gx.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "RNA-seq"
ngs$description = "TCGA prostate cancer data set (from cBioPortal)."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ## get data
    if(0) {
        system("mkdir -p ../downloads/prad/")
        system("wget http://download.cbioportal.org/prad_tcga.tar.gz -P ../downloads/prad")
        system("(cd ../downloads/prad/ && tar xvfz prad_tcga.tar.gz)")
    }
    download.dir = "../downloads/prad"
    ##download.dir = "../../pub/cbio/prad_tcga"
    dir(download.dir)
    
    clin <- read.csv(file.path(download.dir,"data_bcr_clinical_data_patient.txt"),
                     skip=4, sep="\t")
    rownames(clin) <- clin$PATIENT_ID
    sel1 <- c(
        ##"GLEASON_PATTERN_PRIMARY","GLEASON_PATTERN_SECONDARY",
        "GLEASON_SCORE",
        "TUMOR_STATUS",
        "BIOCHEMICAL_RECURRENCE_INDICATOR",
        "RADIATION_TREATMENT_ADJUVANT",
        "TREATMENT_OUTCOME_FIRST_COURSE",
        "CLIN_T_STAGE")
    sel2 <- c("OS_MONTHS","OS_STATUS")
    clin1 <- apply(clin[,sel1], 2, function(x) {x[grep("^\\[",x)]=NA;x})
    ##clin <- data.frame(cbind(clin1, clin[,sel2]))
    clin <- data.frame(clin1)
    clin$CLIN_T_STAGE <- sub("[abc]","",clin$CLIN_T_STAGE)  ## simplify
    clin$GLEASON_SCORE <- gsub("[ ]","",paste0("G",clin$GLEASON_SCORE))  ## simplify
    
    X.data   <- read.csv(file.path(download.dir,"data_RNA_Seq_v2_expression_median.txt"),
                         sep="\t", check.names=FALSE)
    X <- X.data[,3:ncol(X.data)]
    gene <- X.data$Hugo_Symbol
    sum(duplicated(X.data$Hugo_Symbol))
    X <- apply(X, 2, function(x) tapply(x,gene,sum))
    dim(X)
    colnames(X) <- sub("-01$","",colnames(X))
    
    samples <- sort(intersect(colnames(X),rownames(clin)))
    X <- X[,samples]
    sampleTable <- clin[samples,]        
            
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
        ##data_type = "mrna",
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
    ngs$tsne2d=ngs$tsne3d=NULL
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
    ##score <- paste0("G",as.character(ngs$samples$GLEASON_SCORE))
    score <- as.character(ngs$samples$GLEASON_SCORE)
    ##score <- as.character(ngs$samples$CLIN_T_STAGE)
    score <- gsub("[ ]","",score)
    score[is.na(score)] <- "NA"
    table(score)
    ngs$samples$group <- score

    levels = setdiff(unique(ngs$samples$group),NA)
    levels
    colnames(ngs$samples)
    
    ## contr.matrix <- makeContrasts(
    ##     T2_vs_T1 = T2 - T1,
    ##     T3_vs_T1 = T3 - T1,
    ##     T4_vs_T1 = T4 - T1,
    ##     T3_vs_T2 = T3 - T2,
    ##     T4_vs_T3 = T4 - T3,
    ##     levels = levels)

    contr.matrix <- makeContrasts(
        G10_vs_G9 = G10 - G9,
        G9_vs_G8 = G9 - G8,
        G8_vs_G7 = G8 - G7,
        G7_vs_G6 = G7 - G6,
        G10_vs_G6 = G10 - G6,
        G9_vs_G6 = G9 - G6,
        G8_vs_G6 = G8 - G6,
        levels = levels)
    
    dim(contr.matrix)
    head(contr.matrix)
    
    ## some methods (yet) cannot handle direct contrasts!!!
    ## contr.matrix <- makeDirectContrasts(
    ##     ngs$samples[,c("PR_STATUS","ER_STATUS","HER2_STATUS")],
    ##     ref = c("0","0","0"))
    ## dim(contr.matrix)
    ## head(contr.matrix)
    ##contr.matrix = contr.matrix[,1:3]
    
    ##USER.GENETEST.METHODS = c("trend.limma","ttest.welch")
    ##USER.GENESETTEST.METHODS = c("gsva","camera")
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





