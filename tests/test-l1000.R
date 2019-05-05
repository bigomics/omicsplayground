source("../R/gx-heatmap.r")
source("../R/gx-limma.r")
source("../R/gx-util.r")
source("../R/gset-gsea.r")
source("../R/xcr-graph.r")
source("../R/pgx-graph.R")

library(Rtsne.multicore)
library(qlcMatrix)

lincs.dn <- read.gmt("../../Data/gmt/mol_lincs_dn.gmt")
lincs.up <- read.gmt("../../Data/gmt/mol_lincs_up.gmt")

genes <- sort(unique(c(unlist(lincs.dn),unlist(lincs.up))))
length(genes)

require(Matrix)
idx.up <- lapply(lincs.up, function(g) match(g, genes))
idx.dn <- lapply(lincs.dn, function(g) match(g, genes))
row.up <- lapply(1:length(lincs.up), function(i) rep(i, length(lincs.up[[i]])))
row.dn <- lapply(1:length(lincs.dn), function(i) rep(i, length(lincs.dn[[i]])))

i <- c( unlist(idx.up), unlist(idx.dn))
j <- c( unlist(row.up), unlist(row.dn))
x <- c( rep(+1, length(unlist(row.up))), rep(-1, length(unlist(row.dn))))
G <- sparseMatrix( i=i, j=j, x=x, dims=c(length(genes),length(lincs.up)))
rownames(G) <- genes
colnames(G) <- names(lincs.up)
saveRDS(G, file="../files/lincs1000-sparsematrix.rds")

