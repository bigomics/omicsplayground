
RDIR = "../R"
FILES = "../lib"
FILESX = "../libx"
PGX.DIR = "../data"
source("../R/pgx-include.R")
source("../R/pgx-getgeo.R")
##source("options.R")

rda.file="../data/GSE45042-corona2.pgx"
rda.file

##load(file=rda.file, verbose=1)

## load series and platform data from GEO
ngs <- pgx.getGEOseries("GSE45042")
names(ngs)

ngs$name = gsub("^.*/|[.]pgx$","",rda.file)
ngs$datatype = "mRNA (microarray)"
ngs$description = "Cell host-response to infection with novel human coronavirus EMC predict potential antivirals and important differences with SARS-coronavirus."

grep("SSX",rownames(ngs$counts),value=TRUE)
aa <- title2pheno(ngs$samples$source, split="_")
colnames(aa) <- c("cell.line","treatment","time")
head(aa)
ngs$samples <- data.frame( source=ngs$samples$source, aa[,2:3])

## Pre-calculate t-SNE 
ngs <- pgx.clusterSamples(ngs, perplexity=NULL, skipifexists=FALSE, prefix="C")

## not very good autocontrasts...
colnames(ngs$contrasts) 

head(ngs$samples)
grp <- sub("ECL001_","",ngs$samples$source)
table(grp)
ngs$samples$group <- grp
levels = unique(ngs$samples$group)
levels
    
contr.matrix <- makeContrasts(
    EMC_0h_vs_Mock_0h = EMC_0hr - Mock_0hr,
    EMC_3h_vs_Mock_3h = EMC_3hr - Mock_3hr,
    EMC_7h_vs_Mock_7h = EMC_7hr - Mock_7hr,
    EMC_12h_vs_Mock_12h = EMC_12hr - Mock_12hr,
    EMC_18h_vs_Mock_24h = EMC_18hr - Mock_24hr,
    EMC_24h_vs_Mock_24h = EMC_24hr - Mock_24hr,
    levels = levels)

contr.matrix
##contr.matrix = contr.matrix[,1:3]

rda.file
ngs$timings <- c()

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

extra <- c("connectivity")
extra <- c("drugs")
extra <- c("meta.go","infer","drugs","wordcloud","connectivity")
ngs <- compute.extra(ngs, extra, lib.dir=FILES) 

names(ngs)
ngs$timings

## save
rda.file
ngs.save(ngs, file=rda.file)
 



