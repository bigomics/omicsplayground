##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
## Author: BigOmics Analytics (IK)
## Date:   2020
## 

RDIR = "../R"
FILES = "../lib"
PGX.DIR = "../data"
source("../R/pgx-include.R")
FILES


##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------

rda.file="../data/rieckmann2017-immprot.pgx"
rda.file

##load(file=rda.file, verbose=1)
ngs <- list()  ## empty object
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$date = date()
ngs$datatype = "LC-MS proteomics"
ngs$description = "Mass-spectrometry-based proteomics of 28 primary human hematopoietic cell populations in steady and activated states (Rieckmann et al, Nat Immunol. 2017). "

##------------------------------------------------------------
## Read data
##------------------------------------------------------------

## Split data file
D = read.csv("../ext/immprot/ni.3693-S5-copynumber.csv", check.names=FALSE)
D = data.frame(D) ##
colnames(D)

## gene annotation
genes = D[,c("Gene.names","Majority.protein.IDs")]
##genes = apply(genes,2,as.character)
head(genes)
colnames(genes) = c("gene_name","protein.ids")
gg = as.character(genes$gene_name)
gg = sapply(gg, function(x) strsplit(x,split=";")[[1]][1]) ## take just first gene??
gg <- correctMarchSeptemberGenes(gg)
genes$gene_name <- gg

require(org.Hs.eg.db)
GENE.TITLE = unlist(as.list(org.Hs.egGENENAME))
gene.symbol = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = gene.symbol
genes$gene_title <- GENE.TITLE[as.character(genes$gene_name)]
genes <- genes[,c("gene_name","gene_title")]

## give rownames
rownames(D) = paste0("tag",1:nrow(D),":",genes$gene_name)  ## add one gene to tags
rownames(genes) = rownames(D)

## extract data blocks
counts = D[,grep("^CopyNumber_",colnames(D))]  ##
sum(is.nan(as.matrix(counts)))
##counts[is.nan(counts)] = 0
##counts = round(counts / 1e6)  ## too big for R integers, divide by 1M
summary(colSums(counts))
##LFQ = as.matrix(D[,grep("^LFQ",colnames(D))])
## imputed = D[,grep("^Imputed",colnames(D))]
## copynumber = D[,grep("^Copy number",colnames(D))]
## concentration = D[,grep("^Concentration",colnames(D))]
## abundance = D[,grep("^Abundance",colnames(D))]
## welch.pvalue = 10**(-D[,grep("^-Log welch p value",colnames(D))])
## welch.pvalue = exp(-D[,grep("^-Log welch p value",colnames(D))])
## welch.difference = as.matrix(D[,grep("^welch Difference",colnames(D))])

## sample annotation
sample.names = sub("CopyNumber_","",colnames(counts))
sampleTable = data.frame(t(sapply(sample.names, function(x) strsplit(x,split="_")[[1]])))
colnames(sampleTable) = c("cell.type","repl","state")

sampleTable$state = sub("steady.state","S",sampleTable$state)
sampleTable$state = sub("activated","A",sampleTable$state)

## read new names
subtype0 <- as.character(sampleTable[,1])
subtype <- subtype0
subtype = sub("Basophil","BS.ph",subtype)
subtype = sub("Eosinophil","ES.ph",subtype)
subtype = sub("Neutrophil","NT.ph",subtype)
subtype = sub("nonclassical","nc",subtype)
subtype = sub("classical","cl",subtype)
subtype = sub("intermediate","im",subtype)
subtype = sub("memory","mem",subtype)
subtype = sub("naive","nav",subtype)
subtype = sub("bright","bri",subtype)
subtype = sub("EM$","em",subtype)
subtype = sub("CM$","cm",subtype)
subtype = sub("EMRA$","emra",subtype)
subtype <- gsub("[.]","",subtype)
table(subtype)
table(subtype0)

ct.idx =
    1*grepl("Bme|Bna|plasma",subtype,ignore.case=TRUE) +
    2*grepl("T8|T4|Th|Treg",subtype,ignore.case=TRUE) +
    3*grepl("NK",subtype,ignore.case=TRUE) +
    4*grepl("mono|dendrit|macroph|DC|MO",subtype,ignore.case=TRUE) +
    5*grepl("neutro|eosino|baso|ph",subtype,ignore.case=TRUE)
cell.type = c(" ","B","T","NK","Monocyte","Granulocyte")[1 + ct.idx]
table(cell.type)
cbind(cell.type, subtype)

sampleTable$cell.type <- cell.type
sampleTable$subtype <- subtype

short.names <- apply(sampleTable[,c("subtype","repl","state")],1,paste,collapse="_")
colnames(counts) <- short.names
rownames(sampleTable) <- short.names

group <- apply(sampleTable[,c("subtype","state")],1,paste,collapse="_")
sampleTable$group = group
table(group)

##-------------------------------------------------------------------
## collapse multiple row for genes by summing up counts
##-------------------------------------------------------------------
sum(duplicated(genes$gene_name))
x1 = apply(counts, 2, function(x) tapply(x, genes$gene_name, sum))
genes = genes[match(rownames(x1),genes$gene_name),]
counts = x1
##row.id = paste0("tag",1:nrow(raw),":",raw$genes[,"gene_name"])
rownames(genes) = rownames(counts) = rownames(x1)
remove(x1)

##-------------------------------------------------------------------
## Now create an DGEList object  (see tximport Vignette)
##-------------------------------------------------------------------
library(limma)
library(edgeR)
if(is.null(sampleTable$group)) stop("samples need group")
table(sampleTable$group)

##raw <- DGEList(round(data$counts/1e6), group=NULL)  ## we like integer counts...
##raw <- DGEList(round(data$counts/1e0), group=NULL)  ## we like integer counts...
##raw <- DGEList(round(counts), group=NULL)  ## we like integer counts...
ngs$counts <- round(counts)
ngs$samples <- sampleTable
ngs$genes = genes
##raw$samples$activated <- NULL
##lib.size <- colSums(data$counts / 1e6)  ## get original summed intensity as lib.size
ngs$samples$batch <- NULL
##raw$samples$batch <- as.integer(lib.size2)

##save(raw, file="./files/NGS_rawCounts.rda")

##-------------------------------------------------------------------
## sample filtering
##-------------------------------------------------------------------

##-------------------------------------------------------------------
## gene filtering
##-------------------------------------------------------------------
##keep <- rep(TRUE,nrow(ngs$counts))
##keep <- filterByExpr(ngs)  ## default edgeR filter
keep <- (rowSums(edgeR::cpm(ngs$counts) > 1) >=3)
table(keep)
ngs$counts <- ngs$counts[keep,]
ngs$genes  <- ngs$genes[keep,]

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters early so we can use it
## for doing differential analysis.
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, skipifexists=FALSE)
head(ngs$samples)

##-------------------------------------------------------------------
## take top varying
##-------------------------------------------------------------------
if(1 && MAX.GENES>0) {
    cat("shrinking data matrices: n=",MAX.GENES,"\n")
    logcpm = edgeR::cpm(ngs$counts, log=TRUE)
    jj <- head( order(-apply(logcpm,1,sd)), MAX.GENES )  ## how many genes?
    head(jj)
    ##bX <- bX[jj,]
    ngs$counts <- ngs$counts[jj,]
    ngs$genes  <- ngs$genes[jj,]
}

##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

group.levels <- unique(ngs$samples$group)
group.levels

## 10 contrasts in total
contr.matrix <- limma::makeContrasts(
    
    Bmem_activation = Bmem_A - Bmem_S,
    Bnav_activation = Bnav_A - Bnav_S,
    
    ##Bplasma_activation = Bplasma_A - Bplasma_S,
    ##BSph_activation = BSph_A - BSph_S,
    ##ESph_activation = ESph_A - ESph_S,
    mDC_activation  = mDC_A - mDC_S,
    
    MOcl_activation = MOcl_A - MOcl_S,
    ##MOim_activation = MOim_A - MOim_S,
    ##MOnc_activation = MOnc_A - MOnc_S,
    
    mTregs_activation = mTregs_A - mTregs_S,
    NKbri_activation = NKbri_A - NKbri_S,
    NKdim_activation = NKdim_A - NKdim_S,
    
    ##NTph_activation = NTph_A - NTph_S,
    nTregs_activation = nTregs_A - nTregs_S,
    pDC_activation = pDC_A - pDC_S,
    
    T4cm_activation = T4cm_A - T4cm_S,
    T4em_activation = T4em_A - T4em_S,
    T4emra_activation = T4emra_A - T4emra_S,
    T4nav_activation = T4nav_A - T4nav_S,
    
    T8cm_activation = T8cm_A - T8cm_S,
    T8em_activation = T8em_A - T8em_S,
    T8emra_activation = T8emra_A - T8emra_S,
    T8nav_activation = T8nav_A - T8nav_S,
    
    ##Th1_activation = Th1_A - Th1_S,
    ##Th17_activation = Th17_A - Th17_S,
    ##Th2_activation = Th2_A - Th2_S,
    
    Bmem_vs_Bnav = Bmem_S - Bnav_S,
    BSph_vs_Bnav = BSph_S - Bnav_S,
    ESph_vs_Bnav = ESph_S - Bnav_S,
    NTph_vs_Bnav = NTph_S - Bnav_S,
    mDC_vs_Bnav = mDC_S - Bnav_S,
    MOcl_vs_Bnav = MOcl_S - Bnav_S,
    MOim_vs_Bnav = MOim_S - Bnav_S,
    MOnc_vs_Bnav = MOnc_S - Bnav_S,
    pDC_vs_Bnav = pDC_S - Bnav_S,
    
    mTregs_vs_T4nav = mTregs_S - T4nav_S,
    nTregs_vs_T4nav = nTregs_S - T4nav_S,
    NKbri_vs_T4nav = NKbri_S - T4nav_S,
    NKdim_vs_T4nav = NKdim_S - T4nav_S,
    
    T4cm_vs_T4nav = T4cm_S - T4nav_S,
    T4em_vs_T4nav = T4em_S - T4nav_S,
    T4emra_vs_T4nav = T4emra_S - T4nav_S,
    
    T8cm_vs_T4nav = T8cm_S - T8nav_S,
    T8em_vs_T4nav = T8em_S - T8nav_S,
    T8emra_vs_T4nav = T8emra_S - T8nav_S,
    
    Th1_vs_T4nav = Th1_S - T4nav_S,
    Th17_vs_T4nav = Th17_S - T4nav_S,
    Th2_vs_T4nav = Th2_S - T4nav_S,
    
    T4nav_vs_Bnav = T4nav_S - Bnav_S,
    T8nav_vs_Bnav = T8nav_S - Bnav_S,
    Th1_vs_Bnav = Th1_S - Bnav_S,
    Th17_vs_Bnav = Th17_S - Bnav_S,
    Th2_vs_Bnav = Th2_S - Bnav_S,
    
    levels = group.levels )
t(contr.matrix)

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
ngs$timings <- c()

GENE.METHODS=c("ttest","ttest.welch", ## "ttest.rank",
                   "voom.limma","trend.limma","notrend.limma",
                   "edger.qlf","edger.lrt","deseq2.wald","deseq2.lrt")
GENESET.METHODS = c("fisher","gsva","ssgsea","spearman",
                    "camera", "fry","fgsea") ## no GSEA, too slow...

GENE.METHODS=c("ttest.welch","trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","camera","fgsea")

MAX.GENES = 8000
MAX.GENESETS = 8000

## new callling methods
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features=MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENES,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("drugs-combo")
extra <- c("connectivity")
extra <- c("meta.go","deconv","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
ngs$timings

##-------------------------------------------------------------------
## save PGX object
##-------------------------------------------------------------------
rda.file
ngs.save(ngs, file=rda.file)

##===================================================================
##========================= END OF FILE =============================
##===================================================================





