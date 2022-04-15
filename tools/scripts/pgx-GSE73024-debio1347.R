##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

##===================================================================
##================ Script to build PGX object =======================
##===================================================================
##
##
##
##

RDIR = "../R"
FILES = "../lib"
FILESX = "../libx"
PGX.DIR = "../data"
source("../R/pgx-include.R")
source("../R/pgx-getgeo.R")
##source("options.R")

##------------------------------------------------------------
## Set data set name
##------------------------------------------------------------

rda.file="../data/GSE73024-debio1437.pgx"
rda.file
##load(file=rda.file, verbose=1)

##------------------------------------------------------------
## Read data
##------------------------------------------------------------
## load series and platform data from GEO
ngs <- pgx.getGEOseries("GSE73024")
names(ngs)
sum(duplicated(rownames(ngs$genes)))

##------------------------------------------------------------
## Set data set information
##------------------------------------------------------------
ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "GSE73024 data set (Nakanishi et al, 2015). Transcriptome analysis of response to a selective FGFR inhibitor (CH5183284/Debio1347), a mitogen-activated protein kinase kinase (MEK) inhibitor, or a phosphoinositide 3-kinase (PI3K) inhibitor."

## get cancertype
ccle <- read.csv("../lib/CCLE_rna_celline.csv", nrow=10,
                 check.names=FALSE, row.names=1)
ctype <- sub(".*,[ ]","",colnames(ccle))
names(ctype) <- toupper(gsub("[- ]","",sub(",[ ].*","",colnames(ccle))))
ctype <- c(ctype, "UMUC14"="renal pelvis carcinoma")
ctype <- c(ctype, "HSC39"="signet ring cell gastric adenocarcinoma")
ctype <- c(ctype, "NCI716"="cecum adenocarcinoma")

## Rename/simplify sample info
head(ngs$samples)
ngs$samples$agent <- gsub(" inhibitor|[_]","",ngs$samples$agent)
ngs$samples$source <- gsub("[- ]","",ngs$samples$source)
ngs$samples$tissue <- ctype[ngs$samples$source]
ngs$samples$cancertype <- makeAcronym(ngs$samples$tissue)

jj <- which(!duplicated(ngs$samples$source))
cancertype <- ngs$samples$cancertype[jj]
names(cancertype) <- ngs$samples$source[jj]

if(1) {
    ## batch correction for cells
    batch <- as.character(ngs$samples$source)
    design = model.matrix( ~ as.character(ngs$samples$agent))
    ##quantile(ngs$counts[ngs$counts>0], probs=0.01)
    bX <- log2(1 + ngs$counts)
    bX = ComBat(bX, batch=batch, mod=design)
    ##bX = removeBatchEffect(bX, batch=batch, design=design)
    bX = normalizeQuantiles(bX)
    ngs$counts <- 2**bX
}

##-------------------------------------------------------------------
## Pre-calculate t-SNE for and get clusters
##-------------------------------------------------------------------
ngs <- pgx.clusterSamples(ngs, perplexity=NULL, skipifexists=FALSE, prefix="C")


##-------------------------------------------------------------------
## Create contrasts 
##-------------------------------------------------------------------

## not very good autocontrasts...
colnames(ngs$contrasts) 
head(ngs$samples)

if(0) {
    ct1 <- makeDirectContrasts( ngs$samples[,"agent",drop=FALSE], ref="DMSO", na.rm=TRUE)
    ngs$contrasts <- ct1$contr.matrix
    ngs$samples$group <- ct1$group
} else {

    ngs$samples$group <- paste(ngs$samples$agent, ngs$samples$source, sep="_")
    levels = unique(ngs$samples$group)
    levels

    ngs$contrasts <- limma::makeContrasts(
        
        "NCIH1581:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_NCIH1581 - DMSO_NCIH1581,
        "NCIH1581:AZD4547FGFR_vs_DMSO" = AZD4547FGFR_NCIH1581 - DMSO_NCIH1581,
        "NCIH1581:CH4987655MEK_vs_DMSO" = CH4987655MEK_NCIH1581 - DMSO_NCIH1581,
        "NCIH1581:CH5132799PI3K_vs_DMSO" = CH5132799PI3K_NCIH1581 - DMSO_NCIH1581,
        
        "AN3CA:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_AN3CA - DMSO_AN3CA,
        "DMS114:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_DMS114 - DMSO_DMS114,
        "HSC39:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_HSC39 - DMSO_HSC39,
        "KATOIII:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_KATOIII - DMSO_KATOIII,
        "KMS11:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_KMS11 - DMSO_KMS11,
        "MFE280:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_MFE280 - DMSO_MFE280,
        "MFE296:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_MFE296 - DMSO_MFE296,
        "NCI716:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_NCI716 - DMSO_NCI716,
        "NCIH520:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_NCIH520 - DMSO_NCIH520,
        "NCIN87:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_NCIN87 - DMSO_NCIN87,
        "RT4:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_RT4 - DMSO_RT4,
        "SNU16:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_SNU16 - DMSO_SNU16,
        "UMUC14:CH5183284FGFR_vs_DMSO" = CH5183284FGFR_UMUC14 - DMSO_UMUC14,

        "*:CH5183284FGFR_vs_DMSO" =
            ((CH5183284FGFR_NCIH1581 - DMSO_NCIH1581) +
            (CH5183284FGFR_AN3CA - DMSO_AN3CA) +
            (CH5183284FGFR_DMS114 - DMSO_DMS114) +
            (CH5183284FGFR_HSC39 - DMSO_HSC39) +                         
            (CH5183284FGFR_KATOIII - DMSO_KATOIII) +
            (CH5183284FGFR_KMS11 - DMSO_KMS11) +
            (CH5183284FGFR_MFE280 - DMSO_MFE280) +                         
            (CH5183284FGFR_MFE296 - DMSO_MFE296) +
            (CH5183284FGFR_NCI716 - DMSO_NCI716) +
            (CH5183284FGFR_NCIH520 - DMSO_NCIH520) +                         
            (CH5183284FGFR_NCIN87 - DMSO_NCIN87) +
            (CH5183284FGFR_RT4 - DMSO_RT4) +
            (CH5183284FGFR_SNU16 - DMSO_SNU16) +
            (CH5183284FGFR_UMUC14 - DMSO_UMUC14)) / 14,
        levels = levels)
    
}

## prepend cancertype
cl <- sub(":.*","",colnames(ngs$contrasts))
cl.cancertype <- cancertype[cl]
cl.cancertype[is.na(cl.cancertype)] = "*"
colnames(ngs$contrasts) <- paste0(cl.cancertype,"_",colnames(ngs$contrasts))
colnames(ngs$contrasts) <- sub("[*]_","",colnames(ngs$contrasts))

##-------------------------------------------------------------------
## Start computations
##-------------------------------------------------------------------
rda.file
ngs$timings <- c()

GENE.METHODS=c("trend.limma","edger.qlf","deseq2.wald")
GENESET.METHODS = c("fisher","gsva","fgsea") ## no GSEA, too slow...

MAX.GENES = 20000
MAX.GENESETS = 5000

## new callling methods
contr.matrix = ngs$contrasts
ngs <- compute.testGenes(
    ngs, contr.matrix,
    max.features = MAX.GENES,
    test.methods = GENE.METHODS)

ngs <- compute.testGenesets (
    ngs, max.features=MAX.GENESETS,
    test.methods = GENESET.METHODS,
    lib.dir=FILES)

extra <- c("connectivity")
extra <- c("drugs")
extra <- c("meta.go","infer","drugs","wordcloud","connectivity")

sigdb = NULL
##sigdb = c("../libx/sigdb-lincs-cp.h5","../libx/sigdb-lincs-gt.h5")
ngs <- compute.extra(ngs, extra, lib.dir=FILES, sigdb=sigdb) 

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






