
require(org.Hs.eg.db)
symbol <- unlist(as.list(org.Hs.egSYMBOL))
NGENES = 1000
FILES = "../files"

source("../R/gx-util.r")

##----------------------------------------------------------------------
##------------------- CCLE CELL LINE signature -------------------------
##----------------------------------------------------------------------
require(data.table)
if(!file.exists("rna_celline.tsv")) {
    system("wget https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2770/static/E-MTAB-2770-atlasExperimentSummary.Rdata")
    system("wget https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-2770/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv")
    system("mv tpms.tsv E-MTAB-2770-tpms.tsv")
}
##load("E-MTAB-2770-atlasExperimentSummary.Rdata",verbose=1)
cell.dat = read.csv("E-MTAB-2770-tpms.tsv",sep="\t",check.names=FALSE,skip=4)
head(cell.dat)[,1:4]
CELL = cell.dat[,3:ncol(cell.dat)]
gene = cell.dat$"Gene Name"
CELL[is.na(CELL)] <- 0
CELL <- apply(CELL, 2, function(x) tapply(x,gene,sum,na.rm=TRUE))
dim(CELL)
head(CELL)[,1:4]
CELL <- CELL[rownames(CELL) %in% symbol,]
dim(CELL)
sdx <- apply(log2(100+CELL),1,sd)
topCELL <- head(CELL[order(-sdx),],NGENES)
dim(topCELL)
topCELL <- pmax(round(topCELL),0)
write.csv(topCELL, file=file.path(FILES,"CCLE_rna_celline.csv"))

## collapse to cancer type
ctype <- sub(".*[,] ","",colnames(CELL))
table(ctype)
tCELL <- t(apply(CELL+1, 1, function(x) tapply(x, ctype, gmean))-1)
sdx <- apply(log2(100+tCELL),1,sd)
topCELL <- head(tCELL[order(-sdx),],NGENES)
dim(topCELL)
topCELL <- pmax(round(topCELL),0)
write.csv(topCELL, file=file.path(FILES,"CCLE_rna_cancertype.csv"))

remove(cell.dat)
remove(CELL)
remove(tCELL)
remove(topCELL)

##----------------------------------------------------------------------
##--------Human Protein atlas TISSUE signature -------------------------
##----------------------------------------------------------------------
tissue.dat = read.csv(file.path(FILES,"rna_tissue.tsv"),sep="\t")
TISSUE = matrix(tissue.dat$Value, ncol=37, byrow=TRUE)
colnames(TISSUE) = tissue.dat$Sample[1:37]
rownames(TISSUE) = tissue.dat$Gene.name[seq(1,nrow(tissue.dat),37)]
TISSUE = TISSUE[order(-apply(TISSUE,1,sd)),]
TISSUE = TISSUE[!duplicated(rownames(TISSUE)),]
remove(tissue.dat)
dim(TISSUE)
sdx <- apply(log2(100+TISSUE),1,sd)
TISSUE <- head(TISSUE[order(-sdx),],NGENES)

hc = hclust(as.dist(1 - cor(log(10+TISSUE))))  ## should be done earlier...
##hc = hclust(dist(t(log(10+TISSUE))))
##plot(hc, hang=-1)
TISSUE = TISSUE[,hc$order]
colnames(TISSUE)
TISSUE.grp = cutree(hc,8)
remove(hc)
save(TISSUE, TISSUE.grp, file=file.path(FILES,"rna_tissue.rda"))

##----------------------------------------------------------------------
##--------Human Protein atlas CELL LINE signature ----------------------
##----------------------------------------------------------------------
require(data.table)
if(!file.exists("rna_celline.tsv")) {
    system("wget https://www.proteinatlas.org/download/rna_celline.tsv.zip")
    system("unzip rna_celline.tsv.zip")
}
cell.dat = read.csv("rna_celline.tsv",sep="\t",check.names=FALSE)
head(cell.dat)
nsamples <- length(unique(cell.dat$Sample))
nsamples
CELL = matrix(cell.dat$Value, ncol=nsamples, byrow=TRUE)
colnames(CELL) = cell.dat$Sample[1:nsamples]
gene = cell.dat$"Gene name"[seq(1,nrow(cell.dat),nsamples)]
CELL[is.na(CELL)] <- 0
CELL <- apply(CELL, 2, function(x) tapply(x,gene,sum))
dim(CELL)
head(CELL)[,1:4]
sdx <- apply(log2(100+CELL),1,sd)
CELL <- head(CELL[order(-sdx),],NGENES)
head(CELL)[,1:4]
write.csv(CELL, file=file.path(FILES,"HPA_rna_celline.csv"))
remove(cell.dat)
remove(CELL)

##----------------------------------------------------------------------
## ---------------------GTex TISSUE signature --------------------------
##----------------------------------------------------------------------
if(!file.exists("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct")) {
    system("wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz .")
    system("gunzip -f GTEx*gz")
}
gtex.data = read.csv("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct",
                     sep="\t", skip=2, header=TRUE, check.names=FALSE)
gtex <- gtex.data[,3:ncol(gtex.data)]
gtex.gene <- gtex.data[,"Description"]
gtex.tissue <- sub("Cells - |","",colnames(gtex))
gtex.tissue <- sub(" - .*","",gtex.tissue)
gtex <- t(apply(gtex+1, 1, function(x) tapply(x,gtex.tissue,gmean))-1)
gtex[is.na(gtex)] <- 0
gtex <- apply(gtex, 2, function(x) tapply(x,gtex.gene,sum))
gtex <- gtex[(rownames(gtex) %in% symbol),]
dim(gtex)
sdx <- apply(log2(100+gtex),1,sd)
gtex <- head(gtex[order(-sdx),],NGENES)
dim(gtex)
write.csv(gtex, file.path(FILES,"GTEx_rna_tissue_tpm.csv"))


##----------------------------------------------------------------------
## ------------------------ DICE signature -----------------------------
##----------------------------------------------------------------------
load("../pgx/schmiedel2018-DICE-mRNA-8k-LT.pgx",verbose=1)
grp <- ngs$samples$group
table(grp)
sum(is.na(ngs$counts))
sigx <- t(apply( ngs$counts+1, 1, function(x) tapply(x, grp, gmean))-1)
dim(sigx)
head(sigx)
sdx <- apply(log2(100+sigx),1,sd)
sigx <- head(sigx[order(-sdx),],NGENES)
write.csv(sigx, file=file.path(FILES,"DICE-signature.csv"))

##----------------------------------------------------------------------
## --------------------- ImmProt signature -----------------------------
##----------------------------------------------------------------------

load("../pgx/rieckmann2017-immprot-8k.ngs",verbose=1)
grp <- ngs$samples$group
table(grp)
sum(is.na(ngs$counts))
sigx <- t(apply( ngs$counts+1, 1, function(x) tapply(x, grp, gmean))-1)
dim(sigx)
head(sigx)
sdx <- apply(log2(100+sigx),1,sd)
sigx <- head(sigx[order(-sdx),],NGENES)
write.csv(sigx, file=file.path(FILES,"ImmProt-signature.csv"))
