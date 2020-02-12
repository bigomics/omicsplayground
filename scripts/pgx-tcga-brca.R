##rm(list=setdiff(ls(),run.param))

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")

MAX.GENES
MAX.GENES = 8000
MAX.GENESETS = 8000
BATCH.CORRECT=TRUE

## run all available methods 
GENETEST.METHODS = c("trend.limma","edger.qlf","deseq2.wald")
GENESETTEST.METHODS = c("gsva","fisher","camera","fgsea","fry","spearman")
GENESETTEST.METHODS = c("gsva","camera","fgsea")

rda.file="../data/tcga-brca_pub-gx.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "RNA-seq"
##ngs$datatype = "multi-omics"
ngs$description = "TCGA breast cancer data set. Gene expression from 526 patients annotated with PAM50 classification. Data from cBioPortal."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## ##############################################################
    ## get data: download data from cBioPortal
    ## 
    download.dir = "/tmp/brca_tcga_pub"
    if(!dir.exists(download.dir)) {
        system(paste("mkdir -p ",download.dir))
        system(paste("wget -c http://download.cbioportal.org/brca_tcga_pub.tar.gz -P ",download.dir))
        system(paste("(cd ",download.dir," && tar xvfz brca_tcga_pub.tar.gz)"))
    } else {
        ## should exists..
    }
    getwd()
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

    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    sum(duplicated(rownames(X)))
    ##ngs <- ngs.collapseByGene(ngs)
    ##dim(ngs$counts)
        
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
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

    GENETEST.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
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

    extra <- c("drugs-combo")
    extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
    ngs <- compute.extra(ngs, extra, lib.dir=FILES) 
    
    names(ngs)
    names(ngs$drugs)
    ngs$timings    

}

## save
rda.file
ngs.save(ngs, file=rda.file)






