## code to prepare internal `1sysdata.rda` dataset goes here

source("../R/setdirs.R")
source(file.path(RDIR,"pgx-functions.R"))
source(file.path(RDIR,"gset-gsea.r"))
    
## All gene families in Human UPPER CASE
require(org.Hs.eg.db)
GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
GENE.SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
names(GENE.TITLE) = GENE.SYMBOL
##GSET.PREFIX.REGEX = paste(paste0("^",GSET.PREFIXES,"_"),collapse="|")
GSET.PREFIX.REGEX = "^BIOCARTA_|^C2_|^C3_|^C7_|^CHEA_|^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^KEA_|^KEGG_|^PID_|^REACTOME_|^ST_"
GENE.SUMMARY = read.csv(file.path(FILES,"gene-summary.csv"),row.names=1)
GENE.SUMMARY = array(GENE.SUMMARY[,1], dimnames=list(rownames(GENE.SUMMARY)))

## GENExGENE <- readRDS(file=file.path(FILES,"GENExGENE-cosSparseKNN500-XL.rds"))
GSETxGENE <- readRDS(file.path(FILES,"gset-sparseG-XL.rds"))
load(file.path(FILES,"gmt-all.rda"),verbose=1)
GSETS = gmt.all;remove(gmt.all)

message("[INIT] parsing gene families...")
FAMILIES <- pgx.getGeneFamilies(GENE.SYMBOL, FILES=FILES, min.size=10, max.size=9999)
fam.file <- file.path(FILES,"custom-families.gmt")
if(file.exists(fam.file)) {
    custom.gmt = read.gmt(file.path(FILES,"custom-families.gmt"),add.source=TRUE)
    names(custom.gmt)
    FAMILIES= c(FAMILIES, custom.gmt)
}
FAMILIES[["<all>"]] <- GENE.SYMBOL
f1 <- FAMILIES
names(f1) <- paste0("FAMILY:",names(f1))
names(f1) <- sub("FAMILY:<all>","<all>",names(f1))
GSETS <- c(GSETS,f1)

## convert to integer list (more efficient)
message("[INIT] converting GSETS to list of integers...")
GSET.GENES <- sort(unique(unlist(GSETS)))  ## slow...
iGSETS <- parallel::mclapply(GSETS, function(a) match(a,GSET.GENES))  ## slow...
names(iGSETS) <- names(GSETS)
getGSETS <- function(gs) {
    lapply(iGSETS[gs],function(i) GSET.GENES[i])
}

message("[INIT] parsing collections...")
COLLECTIONS <- pgx.getGeneSetCollections(names(GSETS), min.size=10, max.size=99999)
COLLECTIONS <- COLLECTIONS[order(names(COLLECTIONS))]

##-----------------------------------------------------------------------------
## TISSUE/REFERENCE data sets
##-----------------------------------------------------------------------------
load(file.path(FILES,"sig/rna_tissue.rda"),verbose=TRUE)  ## TISSUE and TISSUE.grp

##-----------------------------------------------------------------------------
## Immune cell markers
##-----------------------------------------------------------------------------

## Really needed???
IMMPROT <- read.csv(file.path(FILES,"sig/ImmProt-signature.csv"),row.names=1)
IMMPROT_MARKERS <- rownames(read.csv(file.path(FILES,"sig/immprot-signature1000.csv"),row.names=1))
DICE_MARKERS <- rownames(read.csv(file.path(FILES,"sig/DICE-signature1000.csv"),row.names=1))
LM22 <- read.csv(file.path(FILES,"sig/LM22.txt"),sep="\t",row.names=1)
LM22_MARKERS <- rownames(LM22)

##-----------------------------------------------------------------------------
## Colors
##-----------------------------------------------------------------------------

COLORS = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
COLORS = rep(c(ggsci::pal_npg("nrc", alpha = 0.7)(10),
               ggsci::pal_aaas("default", alpha = 0.7)(10),
               ggsci::pal_d3("category10", alpha = 0.7)(10)),99)
BLUERED <- colorRampPalette(c("royalblue3","grey90","indianred3"))
PURPLEYELLOW <- colorRampPalette(c("purple","purple3","black","yellow3","yellow"))
PURPLEYELLOW <- colorRampPalette(c("purple","purple4","black","yellow4","yellow"))



##usethis::use_data(
save(
    BLUERED,
    COLLECTIONS,
    COLORS,
    DICE_MARKERS,
    FAMILIES,
    GENE.TITLE,
    GENE.SUMMARY,
    GENE.SYMBOL,
    getGSETS,
    GSET.PREFIX.REGEX,
    GSET.GENES,
    GSETxGENE,
    ## GSETS,   ## too big!
    iGSETS,
    IMMPROT,
    IMMPROT_MARKERS,
    LM22,
    LM22_MARKERS,
    PURPLEYELLOW,
    TISSUE,
    TISSUE.grp,
    ## internal = TRUE, overwrite = TRUE, compress = "gzip"  ## fastest
    file = "../lib/sysdata.rda", compress = "gzip", compression_level = 9
)


if(FALSE) {

    load("../lib/sysdata.rda",verbose=TRUE)
    load("../cache/global-init.rda",verbose=TRUE)    

}
