source("../R/gx-util.r")
source("../R/gset-gsea.r")
source("../R/gset-fisher.r")

load(file="../pgx/GSE10846-dlbcl-mRNA-8k.pgx",verbose=1)
##BiocManager::install("CAFE", version = "3.8")
library(DNAcopy)

head(ngs$X)
head(ngs$genes)

X <- ngs$X[,1:4]
CNA.object <- CNA( X, ngs$genes$chr, ngs$genes$pos,
                  data.type="logratio", sampleid=colnames(X))
smoothed.CNA.object <- smooth.CNA(CNA.object)

segments <- segment(smoothed.CNA.object, verbose=1)
plot(segments, plot.type="w")

