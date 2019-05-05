require(KEGGgraph)
require(KEGG.db)
require(pathview)
kegg.names <- unlist(as.list(KEGG.db::KEGGPATHID2NAME))
kegg.ids <- names(kegg.names)


##load data
data(gse16873.d)
data(demo.paths)

##KEGG view: gene data only
i <- 1
pv.out <- pathview(
    gene.data = gse16873.d[, 1],
    pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
    kegg.native = TRUE)
str(pv.out)
head(pv.out$plot.data.gene)
##result PNG file in current directory



data(gse16873.d);fc=gse16873.d[,1]
id="00073"
id="04110" ## CELL CYCLE

pathview( gene.data = fc, pathway.id=id, ##gene.idtype="SYMBOL",
         species = "hsa", out.suffix="pathview",
         limit = list(gene=3, cpd=1),
         ##kegg.dir="../files/kegg-xml",
         kegg.native=TRUE, same.layer=FALSE )


library(org.Hs.eg.db)
kegg <- org.Hs.egPATH2EG
mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])

symbol <- unlist(as.list(org.Hs.egSYMBOL))

kegg3 <- lapply(kegg2, function(ee) as.vector(symbol[ee]))

