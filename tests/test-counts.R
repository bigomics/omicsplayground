

load(file="../pgx/geiger2018-arginine-8k.pgx",verbose=1)
ngs1 <- ngs
load(file="../pgx/geiger2018-arginine2-4k.pgx",verbose=1)
ngs2 <- ngs

source("../shiny-dev/pgx-init.R", local=TRUE)

genes <- FAMILIES[["Ribosomal proteins"]]
genes <- FAMILIES[["Kinases (KEA)"]]

dim(ngs1$genes)
dim(ngs1$counts)

g1 <- which(ngs1$genes$gene_name %in% genes)
colSums(ngs1$counts[g1,])
colSums(ngs1$counts[g1,]) / colSums(ngs1$counts[,])

g2 <- which(ngs2$genes$gene_name %in% genes)
colSums(ngs2$counts[g2,])
colSums(ngs2$counts[g2,]) / colSums(ngs2$counts[,])
