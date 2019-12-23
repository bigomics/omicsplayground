##
##
## This should contain computation not strictly needed for the BASIC
## functionality of the Playground
##
##
##

##extra <- c("meta.go","deconv","infer","drugs")
compute.extra <- function(ngs, extra, lib.dir) {

    ## detect if it is single or multi-omics
    single.omics <- !any(grepl("\\[",rownames(ngs$counts)))
    single.omics
    if(single.omics) {
        cat(">>> computing extra for SINGLE-OMICS\n")
        rna.counts <- ngs$counts
    } else {
        cat(">>> computing extra for MULTI-OMICS\n")
        data.type <- gsub("\\[|\\].*","",rownames(ngs$counts))
        jj <- which(data.type %in% c("gx","mrna"))
        length(jj)
        if(length(jj)==0) {
            stop("FATAL. could not find gx/mrna values.")
        }
        rna.counts <- ngs$counts[jj,]
        ##rownames(rna.counts) <- gsub(".*:|.*\\]","",rownames(rna.counts))
        is.logged <- ( min(rna.counts, na.rm=TRUE) < 0 ||
                       max(rna.counts, na.rm=TRUE) < 50 )
        if(is.logged) {
            cat("expression data seems log. undoing logarithm\n")
            rna.counts <- 2**rna.counts
        }
    }

    if("meta.go" %in% extra) {
        cat(">>> Computing GO core graph...\n")
        ngs$meta.go <- pgx.computeCoreGOgraph(ngs, fdr=0.05)
    }

    if("deconv" %in% extra) {
        cat(">>> computing deconvolution\n")
        ngs <- compute.deconvolution(
            ngs, lib.dir=lib.dir, rna.counts=rna.counts,
            full=FALSE) 
    }

    if("infer" %in% extra) {
        cat(">>> inferring extra phenotypes...\n")
        ngs <- compute.cellcycle.gender(ngs, rna.counts=rna.counts)
    }

    if("drugs" %in% extra) {
        cat(">>> Computing drug enrichment...\n")
        ngs <- compute.drugEnrichment(ngs, lib.dir=lib.dir) 
    }
    
    if("graph" %in% extra) {
        cat(">>> computing OmicsGraphs...\n")
        ngs <- compute.omicsGraphs(ngs) 
    }
    
    if("wordcloud" %in% extra) {
        cat(">>> computing WordCloud statistics...\n")
        res <- pgx.calculateWordFreq(ngs, progress=NULL, pg.unit=1)        
        ngs$wordcloud <- res
        remove(res)
    }
    
    return(ngs)
}


## -------------- deconvolution analysis --------------------------------
compute.deconvolution <- function(ngs, lib.dir, rna.counts, full=FALSE) {
    
    ## list of reference matrices
    refmat <- list()
    readSIG <- function(f) read.csv(file.path(lib.dir,f), row.names=1, check.names=FALSE)
    LM22 <- read.csv(file.path(lib.dir,"LM22.txt"),sep="\t",row.names=1)
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
    ##methods = DECONV.METHODS
    methods = c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor","CIBERSORT","EPIC","FARDEEP")
    ##methods <- c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor")
    methods <- c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor")
    ##methods <- c("DCQ","I-NNLS","NNLM","cor")
    ## methods <- c("NNLM","cor")
    ##if(ncol(ngs$counts)>100) methods <- setdiff(methods,"CIBERSORT")  ## too slow...

    ## list of reference matrices
    refmat <- list()
    readSIG <- function(f) read.csv(file.path(lib.dir,f), row.names=1, check.names=FALSE)
    LM22 <- read.csv(file.path(lib.dir,"LM22.txt"),sep="\t",row.names=1)
    refmat[["Immune cell (LM22)"]] <- LM22
    refmat[["Immune cell (ImmProt)"]] <- readSIG("immprot-signature1000.csv")
    refmat[["Immune cell (DICE)"]] <- readSIG("DICE-signature1000.csv")
    refmat[["Immune cell (ImmunoStates)"]] <- readSIG("ImmunoStates_matrix.csv")
    refmat[["Tissue (HPA)"]]       <- readSIG("rna_tissue_matrix.csv")
    refmat[["Tissue (GTEx)"]]      <- readSIG("GTEx_rna_tissue_tpm.csv")
    refmat[["Cell line (HPA)"]]    <- readSIG("HPA_rna_celline.csv")
    refmat[["Cell line (CCLE)"]] <- readSIG("CCLE_rna_celline.csv")
    refmat[["Cancer type (CCLE)"]] <- readSIG("CCLE_rna_cancertype.csv")

    ## list of methods to compute
    ##methods = DECONV.METHODS
    methods = c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor","CIBERSORT","EPIC")
    ## methods <- c("NNLM","cor")

    if(full==FALSE) {
        ## Fast methods, subset of references
        sel = c("Immune cell (LM22)","Immune cell (ImmunoStates)",
                "Tissue (GTEx)","Cell line (HPA)","Cancer type (CCLE)")
        refmat <- refmat[intersect(sel,names(refmat))]
        methods <- c("DCQ","DeconRNAseq","I-NNLS","NNLM","cor")        
    }
    
    ##counts <- ngs$counts
    counts <- rna.counts
    rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
    res <- pgx.multiDeconvolution(counts, refmat=refmat, method=methods)
    ngs$deconv <- res$results
    rownames(res$timings) <- paste0("[deconvolution]",rownames(res$timings))
    res$timings
    ngs$timings <- rbind(ngs$timings, res$timings)

    remove(refmat)
    remove(res)

    return(ngs)
}

## -------------- infer sample characteristics --------------------------------
compute.cellcycle.gender <- function(ngs, rna.counts) {
    pp <- rownames(rna.counts)
    is.mouse = (mean(grepl("[a-z]",gsub(".*:|.*\\]","",pp))) > 0.8)
    is.mouse
    if(!is.mouse) {
        if(1) {
            cat("estimating cell cycle (using Seurat)...\n")
            ngs$samples$cell.cycle <- NULL
            ngs$samples$.cell.cycle <- NULL
            ##counts <- ngs$counts
            counts <- rna.counts
            rownames(counts) <- toupper(ngs$genes[rownames(counts),"gene_name"])
            res <- try( pgx.inferCellCyclePhase(counts) )  ## can give bins error
            if(class(res)!="try-error") {
                ngs$samples$.cell_cycle <- res
                table(ngs$samples$.cell_cycle)
            }
        }
        if(!(".gender" %in% colnames(ngs$samples) )) {
            cat("estimating gender...\n")
            ngs$samples$.gender <- NULL
            X <- log2(1+rna.counts)
            gene_name <- ngs$genes[rownames(X),"gene_name"]
            ngs$samples$.gender <- pgx.inferGender( X, gene_name )
            table(ngs$samples$.gender)
        } else {
            cat("gender already estimated. skipping...\n")
        }
        head(ngs$samples)
    }
    return(ngs)
}

compute.drugEnrichment <- function(ngs, lib.dir) {
    ## -------------- drug enrichment
    ##source(file.path(RDIR,"pgx-drugs.R"))
    ##source(file.path(RDIR,"pgx-graph.R", local=TRUE)
    X <- readRDS(file=file.path(lib.dir,"l1000_es.rds"))
    x.drugs <- gsub("_.*$","",colnames(X))
    length(table(x.drugs))
    dim(X)

    NPRUNE=-1
    NPRUNE=250
    res.mono <- pgx.computeDrugEnrichment(
        ngs, X, x.drugs, methods=c("GSEA","cor"),
        nprune=NPRUNE, contrast=NULL )

    if(is.null(res.mono)) {
        cat("[compute.drugEnrichment] WARNING:: pgx.computeDrugEnrichment failed!\n")
        return(ngs)
    }
    
    res.combo <- pgx.computeComboEnrichment(
        ngs, X, x.drugs, res.mono=res.mono,
        contrasts = NULL,
        ntop=15, nsample=80, nprune=NPRUNE)
    names(res.combo)

    dim(res.mono[["GSEA"]]$X)

    ngs$drugs <- NULL
    ngs$drugs[["mono"]] <- res.mono[["GSEA"]]
    ngs$drugs[["combo"]] <- res.combo
    names(ngs$drugs)

    ## attach annotation
    annot0 <- read.csv(file.path(lib.dir,"L1000_repurposing_drugs.txt"),
                       sep="\t", comment.char="#")
    rownames(annot0) <- annot0$pert_iname
    ##annot0$pert_iname <- NULL
    ngs$drugs$annot <- annot0

    remove(X)
    remove(x.drugs)
    return(ngs)
}

## ------------------ Omics graphs --------------------------------
compute.omicsGraphs <- function(ngs) {
    ## gr1$layout <- gr1$layout[V(gr1)$name,]  ## uncomment to keep entire layout
    ngs$omicsnet <- pgx.createOmicsGraph(ngs)
    ngs$pathscores <- pgx.computePathscores(ngs$omicsnet, strict.pos=FALSE)

    ## compute reduced graph
    ngs$omicsnet.reduced <- pgx.reduceOmicsGraph(ngs)
    ngs$pathscores.reduced <- pgx.computePathscores(ngs$omicsnet.reduced, strict.pos=FALSE)
    ##save(ngs, file=rda.file)
    return(ngs)
}

