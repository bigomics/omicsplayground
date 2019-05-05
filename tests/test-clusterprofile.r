library(AnnotationHub)
library(clusterProfiler)
library(org.Hs.eg.db)
data(geneList)
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
        toType = c("ENSEMBL", "SYMBOL"),
        OrgDb = org.Hs.eg.db)

symbol <- gene.df[,"SYMBOL"]

load("../pgx/GSE10846-dlbcl-mRNA-8k-LT.pgx",verbose=1)
fc <- ngs$gx.meta$meta[[1]]$meta.fx
names(fc) <- rownames(ngs$gx.meta$meta[[1]])
fc <- sort(fc,decreasing=TRUE)
gene <- head(names(sort(-abs(fc))),500)

ego <- enrichGO(gene  = gene,
                universe = names(fc),
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.25,
                qvalueCutoff  = 0.25)
head(ego@result)
summary(ego@result$qvalue)

dotplot(ego, showCategory=30)
emapplot(ego, vertex.label.cex=2.5)

cnetplot(ego, foldChange=fc)
plotGOgraph(ego, firstSigNodes = 20)

gsecc <- gseGO(geneList=fc, ont="BP", OrgDb=org.Hs.eg.db,
               pvalueCutoff = 0.25,
               keyType = "SYMBOL", verbose=FALSE)
head(gsecc@result)
gseaplot(gsecc, geneSetID="GO:0048646")

plotGOgraph(gsecc, firstSigNodes = 5)

