##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
##source("options.R")
FILES
MAX.GENES = 8000
MAX.GENESETS = 8000


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/GSE22886-immune.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE22886 data set (Abbas et al, 2005). Twelve different types of human leukocytes from peripheral blood and bone marrow, treated to induce activation and/or differentiation, and profiled their gene expression before and after treatment. The twelve cell types are: B cells, CD14+ cells, CD4+ CD45RO+ CD45RA- T cells, CD4+ T cells, CD8+ T cells, IgG/IgA memory B cells, IgM memory B cells, Monocytes, NK cells, Neutrophils, Plasma cells from bone marrow, and Plasma cells from PBMC."

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

library(GEOquery)
if(1) {
    if(!require("hgu133plus2.db")) BiocManager::install("hgu133plus2.db")
    if(!require("hgu133a.db")) BiocManager::install("hgu133a.db")
    if(!require("hgu133b.db")) BiocManager::install("hgu133b.db")
}
library(hgu133plus2.db)
library(hgu133a.db)
library(hgu133b.db)

## load series and platform data from GEO
gset <- getGEO("GSE22886", GSEMatrix =TRUE, AnnotGPL=TRUE)
length(gset)
attr(gset, "names")

## Note: there are two batches: 133A and 133B arrays
y <- list()
X <- NULL
i=1
i=2
for(i in 1:2) {

    pdata = pData(gset[[i]])
    head(pdata)
    sample = sub(" \\[.*","",pdata[,"title"])
    replicate = sub(".*-","",sample)

    cell.type = pdata[,"cell type:ch1"]
    table(cell.type)
    cell.type = gsub("[ +-/]","",cell.type)
    tissue = pdata[,"tissue:ch1"]
    cell.type2 = sub("-[0-9]$|-[0-9][0-9]$","",sample)
    cell.type2 = sub("[0-9][0-9]hour$","",cell.type2)
    cell.type2 = sub("Bcell-naïve","Bcell-naive",cell.type2)
    cell.type2 = sub("From","",cell.type2)
    cell.type2 = gsub("[-+_ ]","",cell.type2)
    state = c("resting","activated")[1 + 1*(grepl("activ|stimula|differen",pdata$description) ) ]
    table(state)

    major.idx =
        1*grepl("Bcell|plasma",cell.type,ignore.case=TRUE) +
        2*grepl("CD8|CD4|Tcell",cell.type,ignore.case=TRUE) +
        3*grepl("NK",cell.type,ignore.case=TRUE) +
        4*grepl("mono|dendrit|macroph|CD14cell",cell.type,ignore.case=TRUE) +
        5*grepl("neutro|eosino|baso",cell.type,ignore.case=TRUE)
    major.type = c(" ","B","T","NK","Monocyte","Granulocyte")[1 + major.idx]
    table(major.type)


    y[[i]] = data.frame( sample = sample,
                        replicate = replicate,
                        cell.family = major.type,
                        cell.type = cell.type,
                        group = cell.type2,
                        ## title = pdata$title,
                        ##batch = ifelse(i==1,"133A","133B"),
                        tissue = tissue,
                        state = state )
    rownames(y[[i]]) = rownames(pdata)
}
table(y[[1]][,1] == y[[2]][,1])
sampleTable <- y[[1]]
rownames(sampleTable) = y[[1]][,"sample"]
##write.csv(sampleTable, file=sub(".pgx","-sampleTable.csv",rda.file))

## merge data sets
x1 = log2(exprs(gset[[1]]))
x2 = log2(exprs(gset[[2]]))
colnames(x1) = y[[1]]$sample
colnames(x2) = y[[2]]$sample
table(colnames(x1) == colnames(x2))
X = rbind(x1, x2)
X = X[order(-apply(X,1,sd)),]
dim(X)

## convert affymetrix ID to GENE symbol
affx  = sapply(as.list(hgu133plus2SYMBOL),"[[",1)
symbol = affx[rownames(X)]
hugo = alias2hugo(symbol)
rownames(X) = symbol
X = X[which(!duplicated(symbol) & !is.na(symbol) & symbol!=""),]
dim(X)
X = X[which(rowMeans(is.na(X))==0), ]  ## no missing values
sum(is.na(X))
dim(X)

## conform tables
table(rownames(sampleTable) == colnames(X))

## gene annotation
require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
head(GENE.TITLE)
gene_title <- GENE.TITLE[rownames(X)]
genes = data.frame( gene_name=rownames(X), gene_title=gene_title)
##genes = apply(genes,2,as.character)
head(genes)
rownames(genes) = rownames(X)

## treat RMA as "pseudo-counts"
counts = round(2**X)
sum(is.nan(as.matrix(counts)))
min(counts)
##counts[is.nan(counts)] = 0
##counts = round(counts / 1e6)  ## too big for R integers, divide by 1M

## sample annotation
##short.names <- paste(sampleTable$cell.type,sampleTable$batch,sampleTable$repl,sep="_")
short.names <- paste(colnames(counts),sampleTable$cell.type,sep="_")
short.names <- colnames(counts)
##short.names <- sub("Roche_","R",short.names)
##short.names <- sub("HUG_","H",short.names)
short.names <- gsub("[ ]","",short.names)
short.names
colnames(counts) <- short.names
rownames(sampleTable) <- short.names
head(sampleTable)
sampleTable$sample <- NULL
head(sampleTable)

##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
if(is.null(sampleTable$group)) stop("samples need group")
table(sampleTable$group)
ngs$counts <- round(counts)
ngs$samples <- sampleTable
ngs$genes = genes
##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
ngs$samples$batch <- NULL
##ngs$samples$batch <- as.integer(lib.size2)

## tagged rownames
row.id = paste0("tag",1:nrow(ngs$genes),":",ngs$genes[,"gene_name"])
rownames(ngs$genes) = rownames(ngs$counts) = row.id
names(ngs)

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
sum(duplicated(ngs$genes$gene_name))
x1 = apply(ngs$counts, 2, function(x) tapply(x, ngs$genes$gene_name, sum))
ngs$genes = ngs$genes[match(rownames(x1),ngs$genes$gene_name),]
ngs$counts = x1
rownames(ngs$genes) = rownames(ngs$counts) = rownames(x1)
remove(x1)

##-------------------------------------------------------------------
## gene filtering
##-------------------------------------------------------------------
##keep <- rep(TRUE,nrow(ngs$counts))
##keep <- filterByExpr(ngs)  ## default edgeR filter
if(0) {
    keep <- (rowSums(edgeR::cpm(ngs$counts, log=TRUE) > 1) >=3)
    table(keep)
    ngs$counts <- ngs$counts[keep,]
    ngs$genes  <- ngs$genes[keep,]
}

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters 
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE)
head(ngs$samples)


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------
levels = levels(ngs$samples$group)
contr.matrix <- limma::makeContrasts(
    
    BcellMemoryIgGIgA_vs_Bcellnaive = BcellMemoryIgGIgA - Bcellnaive,
    BcellMemoryIgM_vs_Bcellnaive = BcellMemoryIgM - Bcellnaive,
    BcellMemoryIgGIgA_vs_BcellMemoryIgM = BcellMemoryIgGIgA - BcellMemoryIgM,
    PlasmaCellBoneMarrow_vs_Bcellnaive = PlasmaCellBoneMarrow - Bcellnaive,
    PlasmaCellPBMC_vs_Bcellnaive = PlasmaCellPBMC - Bcellnaive,
    PlasmaCellBoneMarrow_vs_PlasmaCellPBMC = PlasmaCellBoneMarrow - PlasmaCellPBMC,
    
    NeutrophilResting_vs_MonocyteDay0 = NeutrophilResting - MonocyteDay0,
    MonocyteDay1_vs_MonocyteDay0 = MonocyteDay1 - MonocyteDay0,
    MonocyteDay7_vs_MonocyteDay0 = MonocyteDay7 - MonocyteDay0,
    DendriticCellLPSstimulated_vs_MonocyteDay0 = DendriticCellLPSstimulated - MonocyteDay0,
    DendriticCellControl_vs_MonocyteDay0 = DendriticCellControl - MonocyteDay0,
    DendriticCellLPSstimulated_vs_DendriticCellControl = DendriticCellLPSstimulated - DendriticCellControl,
    
    CD4TcellTh1restimulated_vs_CD4TcellN0 = CD4TcellTh1restimulated - CD4TcellN0,
    CD4TcellTh2restimulated_vs_CD4TcellN0 = CD4TcellTh2restimulated - CD4TcellN0,
    MemoryTcellROactivated_vs_CD4TcellN0 = MemoryTcellROactivated - CD4TcellN0,
    MemoryTcellROunactivated_vs_CD4TcellN0 = MemoryTcellROunactivated - CD4TcellN0,
    CD8TcellN0_vs_CD4TcellN0 = CD8TcellN0 - CD4TcellN0,
    NKcellcontrol_vs_CD4TcellN0 =  NKcellcontrol - CD4TcellN0,
    NKcellIL15stimulated_vs_NKcellcontrol = NKcellIL15stimulated - NKcellcontrol,
    NKcellIL2stimulated_vs_NKcellcontrol = NKcellIL2stimulated - NKcellcontrol,
    
    Bcellnaive_vs_CD4TcellN0 = Bcellnaive - CD4TcellN0,
    MonocyteDay0_vs_CD4TcellN0 = MonocyteDay0 - CD4TcellN0,
    MonocyteDay0_vs_Bcellnaive = MonocyteDay0 - Bcellnaive,
    
    levels = levels)
t(contr.matrix)


##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------

##contr.matrix = contr.matrix[,1:3]
##source("../R/compute-genes.R")
##source("../R/compute-genesets.R")
##source("../R/compute-extra.R")

GENE.METHODS=c("trend.limma","edger.qlf","edger.lrt")
GENESET.METHODS=c("gsva","fisher","camera","fgsea")
    
## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features = MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir = FILES)

extra <- c("connectivity")
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------

rda.file
ngs$drugs$combo <- NULL  ## save space!!
ngs.save(ngs, file=rda.file)




