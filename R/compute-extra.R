##
##
## This should contain computation not strictly needed for the BASIC
## functionality of the Playground
##
##
##

SAVED.PARAM <- c("ngs",ls())

##rm(list=ls())
library(knitr)
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(matrixTests)
library(kableExtra)
library(knitr)

source(file.path(RDIR,"gx-limma.r"))
source(file.path(RDIR,"gx-util.r"))
source(file.path(RDIR,"gset-fisher.r"))
source(file.path(RDIR,"gset-gsea.r"))
source(file.path(RDIR,"gset-meta.r"))
source(file.path(RDIR,"pgx-graph.R"))
source(file.path(RDIR,"pgx-functions.R"))
source(file.path(RDIR,"pgx-deconv.R"))
source(file.path(RDIR,"pgx-drugs.R"))

##FILES = "../files"
##EXTRA.STUFF=1
##rda.file="../files/geiger2018-liver-fltBC.pgx"
##rda.file="../files/rieckmann2017-immprot.pgx"

##--------------------  cleanup deprecated objects --------------------
ngs$tsne2d.genes <- NULL
ngs$tsne3d.genes <- NULL
ngs$tsne2d.gset <- NULL
ngs$tsne3d.gset <- NULL
ngs$gset_tsne_graph <- NULL
ngs$genes_tsne_graph <- NULL
ngs$families <- NULL
ngs$collections <- NULL

## -------------------GO union core graph ----------------------------------------
if(1) {
    ##source(file.path(RDIR,"pgx-graph.R", local=TRUE)
    cat(">>> Computing GO core graph...\n")
    ngs$meta.go <- pgx.computeCoreGOgraph(ngs, fdr=0.05)
}

## ------------------ Omics graphs --------------------------------
if(1){
    cat(">>> computing OmicsGraphs for",rda.file,"\n")

    source(file.path(RDIR,"xcr-graph.r"))
    source(file.path(RDIR,"pgx-graph.R"))

    ## gr1$layout <- gr1$layout[V(gr1)$name,]  ## uncomment to keep entire layout
    ngs$omicsnet <- pgx.createOmicsGraph(ngs)
    ngs$pathscores <- pgx.computePathscores(ngs$omicsnet, strict.pos=FALSE)

    ## compute reduced graph
    ngs$omicsnet.reduced <- pgx.reduceOmicsGraph(ngs)
    ngs$pathscores.reduced <- pgx.computePathscores(ngs$omicsnet.reduced, strict.pos=FALSE)
    ##save(ngs, file=rda.file)

}

## -------------- deconvolution analysis --------------------------------
if(1) {
    
    cat(">>> computing deconvolution for",rda.file,"\n")
    source(file.path(RDIR,"pgx-deconv.R"))
    ##load(file=rda.file,verbose=1)
    
    ## list of reference matrices
    refmat <- list()
    require(FARDEEP)
    readSIG <- function(f) read.csv(file.path(FILES,f), row.names=1, check.names=FALSE)   
    refmat[["Immune cell (LM22)"]] <- LM22
    refmat[["Immune cell (ImmProt)"]] <- readSIG("immprot-signature1000.csv")
    refmat[["Immune cell (DICE)"]] <- readSIG("DICE-signature1000.csv")
    refmat[["Immune cell (ImmunoStates)"]] <- readSIG("ImmunoStates_matrix.csv")
    refmat[["Tissue (HPA)"]] <- readSIG("rna_tissue_matrix.csv")
    refmat[["Tissue (GTEx)"]] <- readSIG("GTEx_rna_tissue_tpm.csv")
    refmat[["Cell line (HPA)"]] <- readSIG("HPA_rna_celline.csv")
    refmat[["Cell line (CCLE)"]] <- readSIG("CCLE_rna_celline.csv")
    refmat[["Cancer type (CCLE)"]] <- readSIG("CCLE_rna_cancertype.csv")
    
    ## list of methods to compute
    methods = DECONV.METHODS
    methods = c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor","CIBERSORT","EPIC","FARDEEP")
    ##methods <- c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor")
    methods <- c("DCQ","DeconRNAseq","NNLM","cor")    
    ## methods <- c("NNLM","cor")
    
    counts <- ngs$counts
    rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
    res <- pgx.multiDeconvolution( counts, refmat=refmat, method=methods)
    ngs$deconv <- res$results
    rownames(res$timings) <- paste0("[deconvolution]",rownames(res$timings))
    res$timings
    ngs$timings <- rbind(ngs$timings, res$timings)

    remove(refmat)
}

## -------------- infer sample characteristics --------------------------------
if(1){
    cat(">>> adding characteristics for",rda.file,"\n")

    pp <- rownames(ngs$counts)
    is.mouse = (mean(grepl("[a-z]",sub(".*:","",pp))) > 0.8)
    is.mouse
    if(!is.mouse) {

        if(1) {
            cat("estimating cell cycle...\n")
            ngs$samples$cell.cycle <- NULL
            ngs$samples$.cell.cycle <- NULL
            counts <- ngs$counts
            rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
            ngs$samples$.cell_cycle <- pgx.inferCellCyclePhase(counts)
            table(ngs$samples$.cell_cycle)
        }
        if(!(".gender" %in% colnames(ngs$samples) )) {
            cat("estimating gender...\n")
            ngs$samples$.gender <- NULL
            X <- log2(1+ngs$counts)
            gene_name <- ngs$genes[rownames(ngs$counts),"gene_name"]
            ngs$samples$.gender <- pgx.inferGender( X, gene_name )
            table(ngs$samples$.gender)
        } else {
            cat("gender already estimated. skipping...\n")
        }
        head(ngs$samples)

    }

}

## -------------- drug enrichment
if(1) {
    
    source(file.path(RDIR,"pgx-drugs.R"))
    ##source(file.path(RDIR,"pgx-graph.R", local=TRUE)
    cat(">>> Computing drug enrichment...\n")

    X <- readRDS(file=file.path(FILES,"l1000_es_5685drugs.rds"))
    ## X <- readRDS(file=file.path(FILES,"l1000_es_8221drugs.rds"))
    x.drugs <- gsub("_.*$","",colnames(X))
    length(table(x.drugs))
    dim(X)
    
    NPRUNE=-1
    NPRUNE=250
    res.mono <- pgx.computeDrugEnrichment(
        ngs, X, x.drugs, methods=c("GSEA","cor"),
        nprune=NPRUNE, contrast=NULL )

    res.combo <- pgx.computeComboEnrichment(
        ngs, X, x.drugs, res.mono=res.mono,
        ntop=15, nsample=80, nprune=NPRUNE)
    names(res.combo)

    dim(res.mono[["GSEA"]]$X)

    ngs$drugs <- NULL
    ngs$drugs[["mono"]] <- res.mono[["GSEA"]]
    ngs$drugs[["combo"]] <- res.combo
    names(ngs$drugs)

    ## SHOULD MAYBE BE DONE IN PREPROCESSING....
    annot0 <- read.csv(file.path(FILES,"L1000_repurposing_drugs.txt"),
                  sep="\t", comment.char="#")   
    rownames(annot0) <- annot0$pert_iname
    ##annot0$pert_iname <- NULL
    ngs$drugs$annot <- annot0
    
    remove(X)
    remove(x.drugs)

}



##-------------------- cleanup -------------------------------------
rm( list=setdiff( ls(), c("ngs",SAVED.PARAM)))

