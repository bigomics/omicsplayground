source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/gset-gsea.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")

##BiocManager::install("meshes")
library(meshes)
library(DOSE)

load("../pgx/geiger2018-arginine2-8k-LT.pgx")

data(geneList, package="DOSE")
de <- names(geneList)[1:100]
## minGSSize = 200 for only speed up the compilation of the vignette
x <- enrichMeSH(de, MeSHDb = "MeSH.Hsa.eg.db", database='gendoo', category = 'C', minGSSize=200)
head(x)


