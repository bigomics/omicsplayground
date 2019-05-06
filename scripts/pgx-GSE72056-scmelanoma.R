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
source("../R/gx-combat.r")
source("../R/ngs-cook.r")
source("../R/ngs-fit.r")
source("../R/gset-fisher.r")
source("../R/gset-gsea.r")
source("../R/gset-meta.r")
source("../R/pgx-graph.R")
source("../R/pgx-functions.R")

source("options.R")
SMALL
COMPARE="group"
COMPARE="clusters"
COMPARE="pheno"
DOWNSAMPLE=0
DOWNSAMPLE=100

rda.file="../pgx/GSE72056-scmelanoma.pgx"
rda.file = sub(".pgx$",paste0("-vs",COMPARE,".pgx"),rda.file)
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*pgx/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "scRNA-seq"
ngs$description ="GSE72056 melanoma scRNA-seq data set. Single-cell RNA sequencing of 4645 single cells isolated from 19 patients, profiling malignant, immune, stromal, and endothelial cells. Ref: 'Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq', Tirosh et al, Science 2016."

## READ/PARSE DATA
if(PROCESS.DATA) {

    ## #############################################################
    ##   Differential expression analysis with limma
    ## BiocManager::install("GEOquery", version = "3.8")
    library(Biobase)
    library(GEOquery)
    library(data.table)

    ##--------------------------------------------------------------
    ## Read SC counts
    ##--------------------------------------------------------------
    if(0) {
        system("wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72056/suppl/GSE72056_melanoma_single_cell_revised_v2.txt.gz -P /tmp")
    }
    counts = fread("/tmp/GSE72056_melanoma_single_cell_revised_v2.txt.gz",nrow=-1000)
    dim(counts)
    head(counts)[,1:4]
    annot  <- data.frame(counts[1:3,2:ncol(counts)], check.names=FALSE)
    rownames(annot) <- as.character(counts[1:3,][[1]])

    counts <- counts[4:nrow(counts),]
    counts = counts[!is.na(counts[[1]]),]
    gene = counts[[1]]
    counts = as.matrix(counts[,2:ncol(counts)])
    rownames(counts) = gene
    head(counts)[,1:5]

    ## From processing info: "Expression levels of genes were
    ## quantified as Ei,j=log2(TPMi,j/10+1), where TPMi,j refers to
    ## transcript-per-million (TPM) for gene i in sample j, as
    ## calculated by RSEM v1.2.3 in paired-end mode."
    counts <- round(10 * (2^counts - 1))
    summary(colSums(counts))  ## should be about one million...
    dim(counts)
    
    ## correct 'excelized' genes... :(
    rownames(counts) <- sub("1-Dec","DEC1",rownames(counts))
    rownames(counts) <- alias2hugo(rownames(counts))

    ## do we have immune cell coding genes?
    imm.genes <- grep("^IGH|^IGJ|^IGK|^IGL|^TRA[VJCD]|^TRB[VJCD]|^TRD[VJCD]|^TRG[VJCD]",
                      rownames(counts),value=TRUE)
    imm.genes

    ##-------------------------------------------------------------------
    ## collapse multiple row for genes by summing up counts
    ##-------------------------------------------------------------------
    head(rownames(counts))
    sum(duplicated(rownames(counts)))
    ##x1 = apply(counts, 2, function(x) tapply(x, rownames(counts), sum))
    x1 = tapply(1:nrow(counts), rownames(counts), function(i) colSums(counts[i,,drop=FALSE]))
    x1 <- do.call(rbind, x1)
    x1[1:3,1:3]
    counts = x1
    remove(x1)

    ##--------------------------------------------------------------
    ## gene annotation
    ##--------------------------------------------------------------
    require(org.Hs.eg.db)
    GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
    gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = gene.symbol
    head(GENE.TITLE)
    genes = data.frame( gene_name=rownames(counts),
                       gene_title=GENE.TITLE[rownames(counts)] )
    rownames(genes) = rownames(counts)
    head(genes)

    ##--------------------------------------------------------------
    ## Prepare sample table
    ##--------------------------------------------------------------

    ## get main sample annotation from rows 1-3
    aa <- as.character(annot[2,colnames(counts)])
    malignant <- c("0"="unresolved","1"="no","2"="yes")[aa]
    bb <- as.character(annot[3,colnames(counts)])
    cell.type <- c("0"="unclassified","1"="Tcell","2"="Bcell","3"="Macrophage",
                   "4"="endothelial","5"="CAF","6"="NK")[bb]

    ## get sample annotation from name
    nn <- colnames(counts)
    cd45.status <- c("_","neg","pos")[1 + 1*grepl("CD45[-_]neg|45neg",nn,ignore.case=TRUE) +
                                     2*grepl("cd45[-_]pos|45pos",nn,ignore.case=TRUE) ]
    cd90.status <- c("_","pos")[1 + 1*grepl("cd90[-_]pos|90pos",nn,ignore.case=TRUE) ]
    pd1.status <- c("_","neg","pos")[1 + 1*grepl("PD1[-_]neg|PD1neg",nn,ignore.case=TRUE) +
                                     2*grepl("PD1[-_]pos|pd1pos",nn,ignore.case=TRUE) ]
    pd1L.status <- c("_","neg","pos")[1 + 1*grepl("PDL1[-_]neg|PDL1neg",nn,ignore.case=TRUE) +
                                     2*grepl("PDL1[-_]pos|pdl1pos",nn,ignore.case=TRUE) ]
    cd8.status <- c("_","pos")[1 + 1*grepl("CD8",nn,ignore.case=TRUE) ]

    xbraf <- counts["BRAF",]
    ##braf.status <- c("neg","low","high")[ 1 + 1*( xbraf > 10) + 1*( xbraf > 150)]
    braf.status <- c("neg","pos")[ 1 + 1*( xbraf > 10)]
    table(braf.status)   
    table(cell.type, braf.status)
    braf.status.BT <- braf.status
    braf.status.BT[!(cell.type %in% c("Bcell","Tcell"))] <- NA
    table(cell.type, braf.status.BT)
    
    patient <- toupper(gsub("[-_ ]","",substring(nn,1,4)))
    table(patient)

    group <- paste0(cell.type,"_",malignant)
    table(group)
    sampleTable = data.frame(
        group = as.character(group),
        patient = patient,
        cell.type = cell.type,
        malignant = malignant,
        CD45 = cd45.status,
        BRAF = braf.status,
        BRAF.BT = braf.status.BT
        ##CD90 = cd90.status,
        ##PD1 = pd1.status,
        ##PD1L = pd1L.status,
        ##CD8  = cd8.status
    )
    rownames(sampleTable) = colnames(counts)
    apply(sampleTable, 2, table)

    ##-------------------------------------------------------------------
    ## subsampling to decrease number of samples
    ##-------------------------------------------------------------------
    DOWNSAMPLE
    if(DOWNSAMPLE>0) {
        cat("downsampling samples DOWNSAMPLE=",DOWNSAMPLE,"...\n")
        dim(counts)
        table(sampleTable$group)
        summary(colSums(counts))
        sum(duplicated(colnames(counts)))

        ##counts <- counts[,sample(ncol(counts),100)]
        wt <- 0.5 + 1*(sampleTable$BRAF=="pos")
        table(sampleTable$group)
        ##sel <- tapply( 1:ncol(counts), sampleTable$group,
        ##              function(ii) head(sample(ii,prob=wt[ii]),DOWNSAMPLE))
        sel <- list()
        for(u in unique(sampleTable$group)) {
            ii <- which(sampleTable$group == u)
            if(length(ii) <= 5) next
            sel[[u]] <- head(sample(ii, prob=wt[ii]),DOWNSAMPLE)
        }
        sel <- as.vector(unlist(sel))
        length(sel)
        table(sampleTable$cell.type[sel], sampleTable$BRAF[sel])
        
        counts <- counts[,sel]
        sampleTable <- sampleTable[colnames(counts),]
        table(sampleTable$group)
        dim(counts)
    }

    table(sampleTable$cell.type, sampleTable$BRAF)
    table(sampleTable$cell.type, sampleTable$BRAF.BT)

    ##--------------------------------------------------------------
    ## Pooling??
    ##--------------------------------------------------------------
    if(0) {
        table(sampleTable$group)
        summary(colSums(counts))
        ##counts <- counts[,sample(ncol(counts),100)]
        counts <- counts[,sample(ncol(counts),1000)]
    }

    ##-------------------------------------------------------------------
    ## Now create an PGX object
    ##-------------------------------------------------------------------
    if(is.null(sampleTable$group)) stop("samples need group")
    table(sampleTable$group)
    ngs$counts <- round(counts)
    ngs$samples <- sampleTable
    ngs$genes = genes
    ##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
    ngs$samples$batch <- NULL
    ##ngs$samples$batch <- as.integer(lib.size2)

    ## tagged rownames???
    ##row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])
    row.id = ngs$genes[,"gene_name"]
    rownames(ngs$genes) = rownames(ngs$counts) = row.id
    names(ngs)

    ##-------------------------------------------------------------------
    ## gene filtering
    ##-------------------------------------------------------------------
    ##keep <- rep(TRUE,nrow(ngs$counts))
    ##keep <- filterByExpr(ngs)  ## default edgeR filter
    keep <- (rowMeans( edgeR::cpm(ngs$counts) > 1) >= 0.05)
    ##keep <- (rowMeans( ngs$counts >= 3) >= 0.01)
    table(keep)
    ngs$counts <- ngs$counts[keep,]
    ngs$genes  <- ngs$genes[keep,]
    dim(ngs$genes)
    
    ##-------------------------------------------------------------------
    ## Pre-calculate t-SNE for and get clusters early so we can use it
    ## for doing differential analysis.
    ##-------------------------------------------------------------------
    ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE, prefix="C")
    head(ngs$samples)
    
    dim(ngs$counts)
    ngs$timings <- c()
    rda.file
    save(ngs, file=rda.file)
}


if(DIFF.EXPRESSION) {

    load(file=rda.file, verbose=1)
    COMPARE

    ## ----------------- test genes ------------------------------------------
    if(COMPARE=="pheno") {
        ## use phenotype directly

        head(ngs$samples)
        ##table(ngs$samples$cell.type)
        table(ngs$samples$malignant)
        head(ngs$samples)

        ##contr.matrix <- makeDirectContrasts(
        ##    ngs$samples[,c("malignant","cluster","P2RX7","PDCD1","CD274","CD8A")],
        ##    ref=c("no","cl1","neg","neg","neg","neg") )
        contr.matrix <- makeDirectContrasts(
            Y = ngs$samples[,c("malignant","cluster","BRAF","BRAF.BT")],
            ref = c("no","all","neg","neg") )
        head(contr.matrix)
        ##apply(contr.matrix,2,table)

    } else if(COMPARE=="clusters") {
        
        ## make model matrix for group vs. rest
        clusters <- ngs$samples$cluster
        table(clusters)
        contr.matrix <- makeClusterContrasts(clusters)
        contr.matrix
        ngs$samples$orig.group <- ngs$samples$group
        ngs$samples$group <- clusters
        
    } else if(COMPARE=="groups") {

        levels = unique(as.character(ngs$samples$group))
        levels
        table(ngs$samples$group)
        contr.matrix <- makeFullContrasts(levels)
        dim(contr.matrix)
        contr.matrix

        ##contr.matrix = contr.matrix[,1:3]
    } else {
        stop("COMPARE error")
    }
    
    head(contr.matrix)
    FAST
    USER.GENETEST.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
    USER.GENETEST.METHODS=c("trend.limma","edger.qlf","edger.lrt")
    USER.GENESETTEST.METHODS=c("gsva","camera","fgsea")
    ##USER.GENETEST.METHODS="*"
    ##USER.GENESETTEST.METHODS="*"
    
    ngs$timings <- c()
    source("../R/compute-genes.R")
    source("../R/compute-genesets.R")
    source("../R/compute-extra.R")

}

rda.file
ngs.save(ngs, file=rda.file)







