
require(devtools)
##BiocManager::install("LRBase.Hsa.eg.db")
##install_github("rikenbit/scTensor")

require("LRBase.Hsa.eg.db")
require("scTensor")

columns(LRBase.Hsa.eg.db)
keytypes(LRBase.Hsa.eg.db)
key_HSA <- keys(LRBase.Hsa.eg.db, keytype="GENEID_L")
head(select(LRBase.Hsa.eg.db, keys=key_HSA[1:2],
            columns=c("GENEID_L", "GENEID_R"), keytype="GENEID_L"))


data(GermMale)
data(labelGermMale)
data(tsneGermMale)

sce <- SingleCellExperiment(assays=list(counts = GermMale))
reducedDims(sce) <- SimpleList(TSNE=tsneGermMale$Y)
plot(reducedDims(sce)[[1]], col=labelGermMale, pch=16, cex=2,
  xlab="Dim1", ylab="Dim2", main="Germline, Male, GSE86146")
legend("topleft", legend=c(paste0("FGC_", 1:3), paste0("Soma_", 1:4)),
  col=c("#9E0142", "#D53E4F", "#F46D43", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"),
  pch=16)


library("LRBase.Hsa.eg.db")
library("MeSH.Hsa.eg.db")
## source("scTensor/R/scTensor-internal.R")
data(GermMale)
data(labelGermMale)
data(tsneGermMale)
sce <- SingleCellExperiment(assays=list(counts = GermMale))
reducedDims(sce) <- SimpleList(TSNE=tsneGermMale$Y)
cellCellSetting(sce, LRBase.Hsa.eg.db, labelGermMale, names(labelGermMale))
rks <- cellCellRanks(sce)
cellCellDecomp(sce, ranks=rks$selected)
cellCellReport(sce, reducedDimNames="TSNE",
    title="Cell-cell interaction within Germline_Male, GSE86146",
    author="Koki Tsuyuzaki", thr=30, html.open=TRUE)




