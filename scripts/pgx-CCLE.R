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
source("../scripts/options.R")


USER.GENETEST.METHODS=c("ttest","ttest.welch","trend.limma")
USER.GENESETTEST.METHODS = c("fisher","gsva","spearman","camera","fgsea") ## no GSEA

BATCH.CORRECT=1
rda.file="../pgx/CCLE-drugsensitivity.pgx"
##if(BATCH.CORRECT) rda.file = sub(".pgx$",paste0("-BC.pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$datatype = "RNA-seq"
ngs$description = "CCLE lung cancer subset."


## READ/PARSE DATA
if(PROCESS.DATA) {
    
    library(Biobase)
    library(PharmacoGx)
    data("CCLEsmall")
    
    availablePSets()
    ##CCLE <- downloadPSet("CCLE") ## 245.5 MB!!    
    pset = CCLEsmall
    ##pset = CCLE
    
    genes <- fNames(pset, "rna")
    slotNames(pset)
    names(pset@molecularProfiles)
    ic50 <- summarizeSensitivityProfiles(
        pSet=pset, sensitivity.measure='ic50_published',
        summary.stat="median", verbose=FALSE)
    auc <- summarizeSensitivityProfiles(
        pSet=pset, sensitivity.measure='auc_published',
        summary.stat="median", verbose=FALSE)

    cnv <- summarizeMolecularProfiles(
        pset, cellNames(pset),
        mDataType="cnv", verbose=FALSE)
    mut <- summarizeMolecularProfiles(
        pset, cellNames(pset),
        mDataType="mutation", summary.stat="and",
        verbose=FALSE)
    rna <- summarizeMolecularProfiles(
        pset, cellNames(pset),
        mDataType="rnaseq", verbose=FALSE)
    table( colnames(ic50)==colnames(cnv))
    cnv <- as.matrix(cnv)
    mut <- 1*(as.matrix(mut)=="1")
    rna <- as.matrix(rna)
    mut[is.na(mut)] <- 0
    mut <- mut[rowSums(mut)>=3,]
    ##mut <- mut + 1e-3*matrix(rnorm(length(mut)),nrow(mut),ncol(mut))
    dim(mut)
    dim(rna)

    head(rownames(cnv))
    head(rownames(mut))
    head(rownames(rna))

    ## translate ENSEMBL to symbol
    if( grepl("ENSG",rownames(rna)[1])) {
        library(org.Hs.eg.db)
        symbols <- mapIds(org.Hs.eg.db, keys=rownames(rna),
                          keytype="ENSEMBL", column="SYMBOL")
        jj <- which(!is.na(symbols))
        rna <- rna[jj,]
        rownames(rna) <- symbols[jj]
        rna <- rna[order(-apply(rna,1,sd)),]
        rna <- rna[!duplicated(rownames(rna)),]
    }
    
    rownames(ic50)
    dim(ic50)
    ##drug = "Crizotinib"
    
    has.all <- ( colMeans(is.na(cnv)) < 0.05 &
                 colMeans(is.na(mut)) < 0.05 &
                 colMeans(is.na(rna)) < 0.05 &
                 colMeans(is.na(ic50)) < 0.5 )
    table(has.all)
    sel <- which(has.all)
    data <- list( MUT=mut[,sel], CNV=cnv[,sel], RNA=rna[,sel])

    ## drugs sensitivity
    jj <- 1:nrow(ic50)
    ##jj <- head(order(rowMeans(is.na(ic50))),10) ## only 10?
    Y <- ic50[jj,sel]
    Y <- apply(Y,2,function(y) factor(c("NOT","RESP")[1 + 1*(y<median(y))]))
    rownames(Y) <- rownames(ic50)
    lapply(data, dim)
    
    if(0) {
        library(ggplot2)
        par(mfrow=c(2,2))
        boxplot( data$CNV[1,] ~ Y[1,])
        points( as.integer(Y), data$CNV["ERBB2",], pch=20, cex=4)
        boxplot( data$RNA["ERBB2",] ~ Y)
        points( as.integer(Y), data$RNA["ERBB2",], pch=20, cex=4)
    }

    medianImpute <- function(X) {
        mx <- apply(X,1,median,na.rm=TRUE)
        X0=X
        X0[is.na(X0)]=0
        X0 + is.na(X)*mx
    }
    data <- lapply(data, medianImpute)
    lapply(data, function(x) sum(is.na(x)))

    ##-------------------------------------------------------------------
    ## Concatenate data types
    ##-------------------------------------------------------------------
    for(i in 1:length(data)) {
        rownames(data[[i]]) <- paste0(names(data)[i],":",rownames(data[[i]]))
    }
    X <- do.call(rbind, data)
    dim(X)
        
    ##-------------------------------------------------------------------
    ## gene annotation
    ##-------------------------------------------------------------------
    require(org.Hs.eg.db)
    GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    head(GENE.TITLE)

    gene_name <- sub(".*:","",rownames(X))
    gene_title <- GENE.TITLE[gene_name]

    ## get chromosome locations
    chrloc = as.list(org.Hs.egCHRLOC)
    names(chrloc) = gene.symbol
    chrloc <- chrloc[gene_name]
    loc <- sapply(chrloc, "[", 1)
    chrom <- sapply(chrloc, function(s) names(s)[1])
    loc[sapply(loc,is.null)] <- NA
    chrom[sapply(chrom,is.null)] <- NA
    chrom <- as.vector(unlist(chrom))
    loc   <- as.vector(unlist(loc))

    data_type <- sub(":.*","",rownames(X))
    genes = data.frame( gene_name=gene_name,
                       data_type = data_type,
                       gene_title=gene_title,
                       chr=chrom, pos=loc)
    ##genes = apply(genes,2,as.character)
    rownames(genes) <- rownames(X)
    head(genes)
    
    ##-------------------------------------------------------------------
    ## Now create an DGEList object  (see tximport Vignette)
    ##-------------------------------------------------------------------
    library(limma)
    ##X <- limma::normalizeQuantiles(X)
    ngs$counts <- 2**X  ## treat as counts
    ngs$samples <- data.frame(t(Y))
    ngs$genes = genes
    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, perplexity=20, skipifexists=FALSE, prefix="C")
    head(ngs$samples)

    if(0) {


        omx.type <- sub(":.*","",rownames(X))
        unique(omx.type)
        S <- list()
        for(tp in unique(omx.type)) {
            x1 <- X[which(omx.type == tmp),]
            sx <- computeGeneSetExpression(x1, method=NULL)
            S[[tp]] <- sx
        }



    }


    
    ##-------------------------------------------------------------------
    ## save
    ##-------------------------------------------------------------------
    
    rda.file
    save(ngs, file=rda.file)
}


if(DIFF.EXPRESSION) {

    load(file=rda.file, verbose=1)
    
    head(ngs$samples)
    ngs$samples$group <- NULL

    contr.matrix <- makeDirectContrasts(
        ngs$samples, ref=c(rep("NOT",24),"all"))        
    dim(contr.matrix)
    
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")
    source("../R/compute-extra.R")

}

## save
rda.file
ngs.save(ngs, file=rda.file)













