

source("../R/pgx-pubmed.R")
NCORE <- detectCores(all.tests = FALSE, logical = TRUE)/2
NCORE
BLUERED <- colorRampPalette(
    rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0",
          "#92C5DE", "#4393C3", "#2166AC", "#053061")))

P <- pmid.buildMatrix()
saveRDS(P, file="PMID2SYMBOL_sparsematrix.rds")

gr <- pmid.buildGraph(P)
saveRDS(gr, file="../files/PMID2SYMBOL_xgraph.rds")

if(0) {    
    gr <- readRDS(file="PMID2SYMBOL_xgraph.rds")
    ##gr <- readRDS(file="../files/PMID2SYMBOL_xgraph.rds")
    length(V(gr))
    length(E(gr))
    
    gene <- "CYP1A1"
    gene <- "P2RX7"
    gene <- "RORC"
    gene <- "PDCD1"
    gene <- c("PDCD1","RORC","P2RX7")
    gene <- "TP53"
    gene <- "NAT2"
    gr1 <- pmid.extractGene(gr, gene, nmin=3)
    gr1
    
    library(RISmed)
    ids <- unique(unlist(V(gr1)$pmid))
    ids
    cat("retrieving PMID titles...\n")
    system.time( tt <- ArticleTitle(EUtilsGet(ids)) )
    names(tt) <- ids

    E(gr1)$title  <- sapply(E(gr1)$genes,paste,collapse=",")
    ##vlink <- lapply(V(gr1)$pmid, function(p) paste(sapply(p,pubmedlink),collapse=" "))
    vlink <- lapply(V(gr1)$pmid, function(p) paste(paste0(tt[p]," ",sapply(p,pubmedlink),"<br><br>"),
                                                   collapse=" "))
    vgenes <- sapply(V(gr1)$genes, paste, collapse="+")
    V(gr1)$title  <- paste0(vgenes,"<p>",vlink)
    head(V(gr1)$title)
    nref <- sapply(V(gr1)$pmid,length)
    V(gr1)$size <- 10 * log(1+nref)**0.5
    V(gr1)$name  <- vgenes

    ##gr1 <- mst(gr1, weights=1/E(gr1)$weight)
    ##pos1 <- layout_with_graphopt(gr1)
    pos1 <- layout_with_fr(gr1, weights=E(gr1)$weight**2)
    ##pos1 <- layout_with_fr(gr1)
    
    klust <- cluster_louvain(gr1)$membership
    V(gr1)$color <- rep(rainbow(6),99)[klust]
        
    library(visNetwork)
    data <- toVisNetworkData(gr1)
    visNetwork(nodes = data$nodes, edges = data$edges,
               height=1200, width=1800) %>%
        visIgraphLayout("layout.norm", layoutMatrix=pos1,
                        physics=FALSE) 


    
    
}

