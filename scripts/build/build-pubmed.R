library(org.Hs.eg.db)

if(0) {
    pmid   <- as.list(org.Hs.egPMID2EG)
    symbol <- as.list(org.Hs.egSYMBOL)
    eg <- names(symbol)
    symbol <- sapply(symbol,"[",1)
    names(symbol) <- as.character(eg)
    ngene <- sapply(pmid,length)
    pmid <- pmid[which(ngene <=10)]
    
    ## collapse duplicates
    pmid.gg <- sapply(pmid,paste,collapse=",")
    idx <- tapply(names(pmid), pmid.gg, function(x) x)
    pmid <- pmid[sapply(idx,"[",1)]
    idx <- lapply(idx, function(x) paste0("PMID:",x))
    names(pmid) <- sapply(idx,paste,collapse=",")
    
    ## build PMID2SYMBOL matrix
    require(parallel)
    require(Matrix)
    require(qlcMatrix)
    idx0 <- mclapply(1:length(pmid), function(i) cbind(i, which(eg %in% pmid[[i]])))
    idx <- do.call(rbind,idx0)
    dim(idx)
    P <- sparseMatrix( i=idx[,1], j=idx[,2], x=rep(1,nrow(idx)),
                      dims=c(length(idx0), length(eg)) )
    rownames(P) <- names(pmid)
    colnames(P) <- symbol
    P <- P[,which(Matrix::colSums(P)>0)]
    dim(P)
    saveRDS(P, file="PMID2SYMBOL_sparsematrix.rds")
}


##----------------------------------------------------------------------
##-------------------- create exact GRAPH --------------------------------
##----------------------------------------------------------------------
library(igraph)
library(scran)
require(Rtsne.multicore)
P <- readRDS(file="PMID2SYMBOL_sparsematrix.rds")
dim(P)
P <- P[which(Matrix::rowSums(P) <= 10),]
P <- P[which(Matrix::rowSums(P) >= 2),]
P <- P[,which(Matrix::colSums(P)>0)]
##P <- head(P,4000)
##P <- P[sample(nrow(P),5000),]
P[1:10,1:10]
dim(P)

## creat graph from overlap
M <- P[,] %*% t(P[,])
dim(M)
diag(M) <- 0
gr <- graph_from_adjacency_matrix(M, mode="undirected",
                                  diag=FALSE, weighted=TRUE)
gr <- subgraph.edges(gr, which(E(gr)$weight>0))
##gr

require(visNetwork)
ee <- get.edges(gr, E(gr))
dim(ee)
##p.genes <- apply(P[V(gr)$name,],1,function(x) colnames(P)[which(x!=0)])
p.genes <- mclapply( V(gr)$name, function(v) colnames(P)[which(P[v,]!=0)] )
##shared.genes <- apply(ee[,],1,function(e) names(which(Matrix::colSums(P[e,])==2)))
shared.genes <- mapply( intersect, p.genes[ee[,1]], p.genes[ee[,2]] )
nshared <- sapply(shared.genes, length)
shared.genes2 <- mclapply(shared.genes,paste,collapse=",")
shared.genes2 <- unlist(shared.genes2)
head(shared.genes2)
E(gr)$genes  <- shared.genes
E(gr)$weight <- nshared

table(nshared>0)
gr <- subgraph.edges(gr, which(nshared>0))
gr
P <- P[V(gr)$name,]

plink <- mclapply(V(gr)$name,function(x) gsub("PMID:","",strsplit(x,split=",")[[1]]))
nref <- sapply(plink,length)
##vgenes <- apply(P[V(gr)$name,],1,function(x) paste(names(which(x!=0)),collapse=","))
vgenes <- mclapply(1:nrow(P),function(i) names(which(P[i,]!=0)))
vgenes2 <- unlist(sapply(vgenes,paste,collapse=","))
V(gr)$size <- nref
V(gr)$genes <- vgenes
V(gr)$pmid  <- plink

saveRDS(gr, file="PMID2SYMBOL_xgraph.rds")

##----------------------------------------------------------------------
##-------------------- extract and show GRAPH --------------------------
##----------------------------------------------------------------------

gr <- readRDS(file="PMID2SYMBOL_xgraph.rds")

pmid.extractGene <- function(gr, gene, nmin=3) {
    jj <- c()
    for(g in gene) jj <- c(jj, grep(g, V(gr)$genes))
    jj <- unique(jj)
    gr1 <- induced_subgraph(gr, jj)
    nshared <- sapply(E(gr1)$genes,length)
    ee <- which(!( nshared==1 & sapply(E(gr1)$genes,"[",1)==gene))
    gr1 <- subgraph.edges(gr1, ee, delete.vertices=TRUE)
    cmp <- components(gr1)
    jj <- which(cmp$membership %in% which(cmp$csize >= nmin))
    gr1 <- induced_subgraph(gr1, jj)
    gr1
}

pubmedlink <- function(s) {
    paste0("<a href='https://www.ncbi.nlm.nih.gov/pubmed/",s,
           "' target='_blank'>PMID:",s,"</a>") }

gene <- "CYP1A1"
gene <- "P2RX7"
gene <- "RORC"
gene <- "PDCD1"
gene <- c("PDCD1","RORC","P2RX7")
gene <- "P53"
gr1 <- pmid.extractGene(gr, gene, nmin=3)
gr1
E(gr1)$title  <- sapply(E(gr1)$genes,paste,collapse=",")

vlink <- lapply(V(gr1)$pmid, function(p) paste(sapply(p,pubmedlink),collapse=" "))
vgenes <- sapply(V(gr1)$genes, paste, collapse=",")
V(gr1)$title  <- paste0(vgenes,"<br>",vlink)
head(V(gr1)$title)
nref <- V(gr1)$size
V(gr1)$size <- 10 * log(1+nref)
V(gr1)$name  <- vgenes

pos1 <- layout_with_graphopt(gr1)
##pos1 <- layout_with_fr(gr1)

data <- toVisNetworkData(gr1)
visNetwork(nodes = data$nodes, edges = data$edges,
           height=1200, width=1600) %>%
    visIgraphLayout("layout.norm", layoutMatrix=pos1,
                    physics=FALSE) 


