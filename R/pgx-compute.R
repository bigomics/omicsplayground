if(0) {

    ##source("../scripts/options.R")
    ## MAX.GENES
    DIR = "../data/exampledata/"
    DIR = "~/bigomics/projects/vogel2019-tcell/data"
    counts  = as.matrix(read.csv(file.path(DIR,"counts.csv"), row.names=1))
    samples = read.csv(file.path(DIR,"samples.csv"), row.names=1, stringsAsFactors=FALSE)
    ##genes   = read.csv(file.path(DIR,"genes.csv"), row.names=1, stringsAsFactors=FALSE)
    contrasts = as.matrix(read.csv(file.path(DIR,"contrasts.csv"), row.names=1))
    progress = NULL

    counts = aa$counts
    samples = aa$samples
    contrasts = aa$contrasts

}

if(0) {
    max.genes=max.genesets=25000;lib.dir=FILES;progress=NULL;only.hugo=1;extra.methods=c("meta.go","deconv","infer","drugs","wordcloud")
    gx.methods=c("ttest.welch","trend.limma");gset.methods=c("fisher","gsva");
}

pgx.createPGX <- function(counts, samples, contrasts, ## genes, 
                          only.hugo=TRUE, only.proteincoding=TRUE)
{
    
    ##if(!is.null(progress)) progress$inc(0.01, detail = "creating PGX object")
    
    if(!"group" %in% colnames(samples)) {
        stop("samples information must have 'group' column\n")
        return(NULL)
    }

    ##-------------------------------------------------------------------
    ## convert to gene symbol
    ##-------------------------------------------------------------------
    symbol <- probe2symbol(rownames(counts), type=NULL)
    ##symbol <- alias2symbol(symbol)  ## to latest HUGO
    jj <- which(!is.na(symbol))
    counts <- as.matrix(counts[jj,])
    rownames(counts) <- symbol[jj]

    ##-------------------------------------------------------------------
    ## check counts
    ##-------------------------------------------------------------------
    is.log <- (min(counts,na.rm=TRUE) < 0 || max(counts,na.rm=TRUE) < 100)
    is.log
    if(is.log) {
        counts <- 2**counts  ## undo logarithm
        cat("[pgx.createPGX] undo logarithm\n")
    }

    mean.counts <- mean(colSums(counts,na.rm=TRUE))
    is.toobig <- log10(mean.counts) > 10
    is.toobig
    counts_multiplier = 1
    if(is.toobig) {
        ## scale to about 10 million reads
        ##progress$inc(0.01, detail = "scaling down counts")
        unit <- 10**(round(log10(mean.counts)) - 7)  
        counts <- counts / unit
        counts_multiplier = unit
    }
    counts_multiplier
    cat("[pgx.createPGX] count_multiplier= ",counts_multiplier,"\n")
    
    ##-------------------------------------------------------------------
    ## create ngs object
    ##-------------------------------------------------------------------
    ##load(file=rda.file, verbose=1)
    ngs <- list()  ## empty object
    ngs$name = "data set"
    ngs$date = date()
    ngs$datatype = "unknown"
    ngs$description = "data set"

    ngs$samples = data.frame(samples, check.names=FALSE)
    ngs$counts  = as.matrix(counts)
    ##ngs$genes   = data.frame(genes)
    ngs$contrasts = as.matrix(contrasts)
    ngs$counts_multiplier = counts_multiplier
    
    ##-------------------------------------------------------------------
    ## create gene annotation if not given
    ##-------------------------------------------------------------------
    is.mouse <- (mean(grepl("[a-z]",rownames(counts))) > 0.9)
    org = ifelse(is.mouse, "mouse", "human")
    org
    cat("[pgx.createPGX] detected organism: ",org,"\n")

    gene.symbol <- NULL
    if(org == "human") {
        require(org.Hs.eg.db)
        GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
        gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
        gene.map <- sapply(as.list(org.Hs.egMAP),"[",1)
        names(GENE.TITLE) = gene.symbol
        names(gene.map) = gene.symbol

        ## get gene biotype
        require("EnsDb.Hsapiens.v86")
        daf <- transcripts(EnsDb.Hsapiens.v86,
                           columns = c("gene_name", "gene_biotype"),
                           return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
        
    }
    if(org == "mouse") {
        require(org.Mm.eg.db)
        GENE.TITLE = unlist(as.list(org.Mm.egGENENAME))
        gene.symbol = unlist(as.list(org.Mm.egSYMBOL))
        gene.map <- sapply(as.list(org.Mm.egCHR),"[",1)
        names(GENE.TITLE) = gene.symbol
        names(gene.map) = gene.symbol

        ## get gene biotype
        require("EnsDb.Mmusculus.v79")
        daf <- transcripts(EnsDb.Mmusculus.v79,
                           columns = c("gene_name", "gene_biotype"),
                           return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
    }
    
    gene <- rownames(ngs$counts)
    gene1 <- sapply(gene, function(s) strsplit(s,split="[;,]")[[1]][1])
    gene1 <- alias2hugo(gene1)  ## convert to HUGO
    ngs$genes = data.frame( gene_name = gene1,
                           gene_alias = gene,
                           chr = gene.map[gene1],
                           gene_title = GENE.TITLE[gene1],
                           gene_biotype = GENE.BIOTYPE[gene1] )
    ##rownames(ngs$genes) <- gene1
    
    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    ndup <- sum(duplicated(ngs$genes$gene_name))
    ndup
    if(ndup>0) {        
        gene = as.character(ngs$genes$gene_name)        
        x1 = tapply(1:nrow(ngs$counts), gene, function(i) colSums(ngs$counts[i,,drop=FALSE]))
        x1 <- do.call(rbind, x1)
        ngs$genes = ngs$genes[match(rownames(x1), ngs$genes$gene_name),]
        ngs$counts = x1
        rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
        remove(x1)
    }

    ##-------------------------------------------------------------------
    ## Filter genes?
    ##-------------------------------------------------------------------
    if(org == "mouse") {
        has.name <- !is.na(ngs$genes$gene_name)
        is.official = TRUE
        if(only.hugo) is.official <- (ngs$genes$gene_name %in% gene.symbol)
        rik.genes <- grepl("Rik",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        has.chrloc <- !is.na(ngs$genes$chr)
        is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
        
        keep <- (has.name & !rik.genes & is.official & has.chrloc)
        if(only.proteincoding) keep <- keep & is.protcoding
        
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]        
    }
    if(org == "human") {
        has.name <- !is.na(ngs$genes$gene_name)
        is.official = TRUE
        if(only.hugo) is.official <- (ngs$genes$gene_name %in% gene.symbol)
        ##orf.genes <- grepl("ORF",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        has.chrloc <- !is.na(ngs$genes$chr)
        is.protcoding <- ngs$genes$gene_biotype %in% c("protein_coding")
        
        keep <- (has.name & is.official & has.chrloc)
        if(only.proteincoding) keep <- keep & is.protcoding
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]        
    }
    
    cat("DBG [pgx-upload] 1: dim(ngs$counts)=",dim(ngs$counts),"\n")
    cat("DBG [pgx-upload] 1: sum.is.na(ngs$counts)=",sum(is.na(ngs$counts)),"\n")
    cat("DBG [pgx-upload] 1: dim(ngs$samples)=",dim(ngs$samples),"\n")
    cat("DBG [pgx-upload] 1: dim(ngs$genes)=",dim(ngs$genes),"\n")
    
    return(ngs)
}

if(0) {
    max.genes=19999;max.genesets=9999
    gx.methods = c("ttest.welch","trend.limma","edger.qlf")
    gset.methods = c("fisher","gsva","fgsea")
    extra.methods = c("meta.go","deconv","infer","drugs","wordcloud")
    extra.methods = c("wordcloud")    
}

pgx.computePGX <- function(ngs, 
                           max.genes = 19999, max.genesets = 9999, 
                           gx.methods = c("ttest.welch","trend.limma","edger.qlf"),
                           gset.methods = c("fisher","gsva","fgsea"),
                           extra.methods = c("meta.go","deconv","infer","drugs","wordcloud"),
                           lib.dir = "../lib", do.cluster=TRUE,
                           progress=NULL)
{

    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    if(do.cluster) {
        if(!is.null(progress)) progress$inc(0.01, detail = "clustering")
        perplexity <- max(1,min(30,round(ncol(ngs$counts)/4)))
        perplexity=NULL
        ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=perplexity)
        head(ngs$samples)
        table(ngs$samples$cluster)
    }
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    ## contrast matrix
    contr.matrix <- ngs$contrasts

    ## create groups if not exists
    if(!"group" %in% colnames(ngs$samples)) {
        ct = apply(ngs$samples,2,function(px) mean(px %in% rownames(ngs$contrasts),na.rm=TRUE))
        group.col = which(ct > 0.95)
        if(length(group.col)>0) {
            grp.var = colnames(ngs$samples)[group.col[1]]
            cat(paste0("INFO [pgx-upload] assigning '",grp.var,"' variable as 'group' column\n"))
            ## colnames(ngs$samples)[group.col[1]] <- "group"
            ngs$samples$group <- ngs$samples[,group.col[1]] 
        } else {
            stop("sample annotation file must have 'group' column\n")
        }
    }
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    ngs$timings <- c()
    
    GENETEST.METHODS = c("ttest","ttest.welch","ttest.rank",
                         "voom.limma","trend.limma","notrend.limma",
                         "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    GENESETTEST.METHODS = c("fisher","gsva","ssgsea","spearman",
                            "camera", "fry","fgsea") ## no GSEA, too slow...
    
    ## ------------------ gene level tests ---------------------
    if(!is.null(progress)) progress$inc(0.1, detail = "testing genes")
    
    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = max.genes,
        test.methods = gx.methods)
    head(ngs$gx.meta$meta[[1]])        
    
    ## ------------------ gene set tests -----------------------
    if(!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

    max.features=max.genes;test.methods=gset.methods
    ngs <- compute.testGenesets(
        ngs, max.features = max.genesets,
        test.methods = gset.methods,
        lib.dir = lib.dir )
    head(ngs$gset.meta$meta[[1]])
    
    ## ------------------ extra analyses ---------------------
    if(!is.null(progress)) progress$inc(0.3, detail = "extra modules")
    ##extra <- c("meta.go","deconv","infer","drugs")
    ##extra <- c("meta.go","infer","drugs")
    ngs <- compute.extra(ngs, extra=extra.methods, lib.dir=lib.dir)
    
    ngs$timings
    return(ngs)
}



pgx.computeObjectPGX <- function(counts, samples, contrasts, ## genes, 
                                 ##gx.methods = c("trend.limma","edger.qlf","deseq2.wald"),
                                 max.genes = 9999, max.genesets = 9999, only.hugo=TRUE,
                                 gx.methods = c("ttest.welch","trend.limma","edger.qlf"),
                                 gset.methods = c("fisher","gsva","fgsea"),
                                 extra.methods = c("meta.go","deconv","infer","drugs","wordcloud"),
                                 lib.dir = "../lib", do.cluster=TRUE,
                                 progress=NULL)
{

    library(org.Hs.eg.db)
    cat("****  DEPRECATED: please use createPGX() and computePGX() ****\n")
    
    if(!is.null(progress)) progress$inc(0.01, detail = "creating object")
    
    ##-------------------------------------------------------------------
    ## convert to gene symbol
    ##-------------------------------------------------------------------
    if(1) {
        symbol <- probe2symbol(rownames(counts), type=NULL)
        jj <- which(!is.na(symbol))
        counts <- as.matrix(counts[jj,])
        rownames(counts) <- symbol[jj]
    }

    ##-------------------------------------------------------------------
    ## create ngs object
    ##-------------------------------------------------------------------
    ##load(file=rda.file, verbose=1)
    ngs <- list()  ## empty object
    ngs$name = "data set"
    ngs$date = date()
    ngs$datatype = "unknown"
    ngs$description = "data set"

    ngs$samples = data.frame(samples, check.names=FALSE)
    ngs$counts  = as.matrix(counts)
    ##ngs$genes   = data.frame(genes)
    ngs$contrasts = contrasts
    
    ##-------------------------------------------------------------------
    ## check counts
    ##-------------------------------------------------------------------
    is.log <- (min(ngs$counts,na.rm=TRUE) < 0 || max(ngs$counts,na.rm=TRUE) < 100)
    is.log
    if(is.log) {
        ngs$counts <- 2**ngs$counts  ## undo logarithm
        cat("[pgx.computeObjectPGX] undo logarithm\n")
    }

    mean.counts <- mean(colSums(ngs$counts,na.rm=TRUE))
    is.toobig <- log10(mean.counts) > 10
    is.toobig
    ngs$counts_multiplier = 1
    if(is.toobig) {
        ## scale to about 10 million reads
        ##progress$inc(0.01, detail = "scaling down counts")
        unit <- 10**(round(log10(mean.counts)) - 7)  
        ngs$counts <- ngs$counts / unit
        ngs$counts_multiplier = unit
    }
    ngs$counts_multiplier
    cat("[pgx.computeObjectPGX] count_multiplier= ",ngs$counts_multiplier,"\n")
    
    ##-------------------------------------------------------------------
    ## create gene annotation if not given
    ##-------------------------------------------------------------------
    is.mouse <- (mean(grepl("[a-z]",rownames(counts))) > 0.9)
    org = ifelse(is.mouse, "mouse", "human")
    org
    gene.symbol <- NULL
    if(org == "human") {
        require(org.Hs.eg.db)
        GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
        gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
        gene.map <- sapply(as.list(org.Hs.egMAP),"[",1)
        names(GENE.TITLE) = gene.symbol
        names(gene.map) = gene.symbol

        ## get gene biotype
        require("EnsDb.Hsapiens.v86")
        daf <- transcripts(EnsDb.Hsapiens.v86,
                           columns = c("gene_name", "gene_biotype"),
                           return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
        
    }
    if(org == "mouse") {
        require(org.Mm.eg.db)
        GENE.TITLE = unlist(as.list(org.Mm.egGENENAME))
        gene.symbol = unlist(as.list(org.Mm.egSYMBOL))
        gene.map <- sapply(as.list(org.Mm.egCHR),"[",1)
        names(GENE.TITLE) = gene.symbol
        names(gene.map) = gene.symbol

        ## get gene biotype
        require("EnsDb.Mmusculus.v79")
        daf <- transcripts(EnsDb.Mmusculus.v79,
                           columns = c("gene_name", "gene_biotype"),
                           return.type="DataFrame")
        GENE.BIOTYPE = daf$gene_biotype
        names(GENE.BIOTYPE) = daf$gene_name
    }
    if(!is.null(progress)) {
        aa <- paste("detected organism: ",org)
        progress$inc(0.01, detail = aa)
    }
    cat("[pgx.computeObjectPGX] detected organism: ",org,"\n")
    
    gene <- rownames(ngs$counts)
    gene1 <- sapply(gene, function(s) strsplit(s,split="[;,]")[[1]][1])
    gene1 <- alias2hugo(gene1)  ## convert to HUGO
    ngs$genes = data.frame( gene_name = gene1,
                           gene_alias = gene,
                           chr = gene.map[gene1],
                           gene_title = GENE.TITLE[gene1],
                           gene_biotype = GENE.BIOTYPE[gene1] )
    ##rownames(ngs$genes) <- gene1
    
    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    ndup <- sum(duplicated(ngs$genes$gene_name))
    ndup
    if(ndup>0) {        
        gene = as.character(ngs$genes$gene_name)        
        x1 = tapply(1:nrow(ngs$counts), gene, function(i) colSums(ngs$counts[i,,drop=FALSE]))
        x1 <- do.call(rbind, x1)
        ngs$genes = ngs$genes[match(rownames(x1), ngs$genes$gene_name),]
        ngs$counts = x1
        rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
        remove(x1)
    }

    ##-------------------------------------------------------------------
    ## Filter genes?
    ##-------------------------------------------------------------------
    if(org == "mouse") {
        has.name <- !is.na(ngs$genes$gene_name)
        is.official = TRUE
        if(only.hugo) is.official <- (ngs$genes$gene_name %in% gene.symbol)
        rik.genes <- grepl("Rik",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        has.chrloc <- !is.na(ngs$genes$chr)
        
        keep <- (has.name & !rik.genes & is.official & has.chrloc)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]        
    }
    if(org == "human") {
        has.name <- !is.na(ngs$genes$gene_name)
        is.official = TRUE
        if(only.hugo) is.official <- (ngs$genes$gene_name %in% gene.symbol)
        ##orf.genes <- grepl("ORF",ngs$genes$gene_name)
        ##imm.gene <- grepl("^TR_|^IG_",ngs$genes$gene_biotype)
        has.chrloc <- !is.na(ngs$genes$chr)

        keep <- (has.name & is.official & has.chrloc)
        ngs$counts <- ngs$counts[keep,]
        ngs$genes  <- ngs$genes[keep,]        
    }
    
    cat("DBG [pgx-upload] 1: dim(ngs$counts)=",dim(ngs$counts),"\n")
    cat("DBG [pgx-upload] 1: sum.is.na(ngs$counts)=",sum(is.na(ngs$counts)),"\n")
    cat("DBG [pgx-upload] 1: dim(ngs$samples)=",dim(ngs$samples),"\n")
    cat("DBG [pgx-upload] 1: dim(ngs$genes)=",dim(ngs$genes),"\n")

    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    if(do.cluster) {
        if(!is.null(progress)) progress$inc(0.01, detail = "clustering")
        perplexity <- max(2,min(30,round(ncol(ngs$counts)/4)))
        perplexity
        ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, perplexity=perplexity)
        head(ngs$samples)
        table(ngs$samples$cluster)
    }
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    ## contrast matrix
    contr.matrix <- ngs$contrasts

    ## create groups if not exists
    if(!"group" %in% colnames(ngs$samples)) {
        ct = apply(ngs$samples,2,function(px) mean(px %in% rownames(ngs$contrasts),na.rm=TRUE))
        group.col = which(ct > 0.95)
        if(length(group.col)>0) {
            grp.var = colnames(ngs$samples)[group.col[1]]
            cat(paste0("INFO [pgx-upload] assigning '",grp.var,"' variable as 'group' column\n"))
            ## colnames(ngs$samples)[group.col[1]] <- "group"
            ngs$samples$group <- ngs$samples[,group.col[1]] 
        } else {
            stop("sample annotation file must have 'group' column\n")
        }
    }
    
    ##======================================================================
    ##======================================================================
    ##======================================================================

    ngs$timings <- c()
    
    GENETEST.METHODS = c("ttest","ttest.welch","ttest.rank",
                         "voom.limma","trend.limma","notrend.limma",
                         "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
    GENESETTEST.METHODS = c("fisher","gsva","ssgsea","spearman",
                            "camera", "fry","fgsea") ## no GSEA, too slow...
    
    ## ------------------ gene level tests ---------------------
    if(!is.null(progress)) progress$inc(0.1, detail = "testing genes")
    
    ngs <- compute.testGenes(
        ngs, contr.matrix,
        max.features = max.genes,
        test.methods = gx.methods)
    head(ngs$gx.meta$meta[[1]])        
    
    ## ------------------ gene set tests -----------------------
    if(!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

    max.features=max.genes;test.methods=gset.methods
    system.time(
        ngs <- compute.testGenesets(
            ngs, max.features = max.genesets,
            test.methods = gset.methods,
            lib.dir = lib.dir )
    )
    head(ngs$gset.meta$meta[[1]])
    
    ## ------------------ extra analyses ---------------------
    if(!is.null(progress)) progress$inc(0.3, detail = "extra modules")
    ##extra <- c("meta.go","deconv","infer","drugs")
    ##extra <- c("meta.go","infer","drugs")
    ngs <- compute.extra(ngs, extra=extra.methods, lib.dir=lib.dir)
    
    ngs$timings
    return(ngs)
}








