##
##
## This should contain computation not strictly needed for the BASIC
## functionality of the Playground
##
##
##

##extra <- c("meta.go","deconv","infer","drugs")

compute.extra <- function(ngs, extra, lib.dir, sigdb=NULL) {
    
    if(is.null(sigdb)) {
        sigdb = c(
            file.path(lib.dir,"../data/datasets-allFC.csv"),
            file.path(lib.dir,"sigdb-archs4.h5"),
            file.path(lib.dir,"sigdb-creeds.h5")
        )
    }
    
    timings <- c()
    
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
        tt <- system.time({
            ngs$meta.go <- pgx.computeCoreGOgraph(ngs, fdr=0.05)
        })
        timings <- rbind(timings, c("meta.go", tt))
    }

    if("deconv" %in% extra) {
        cat(">>> computing deconvolution\n")
        tt <- system.time({
            ngs <- compute.deconvolution(
                ngs, lib.dir=lib.dir, rna.counts=rna.counts,
                full=FALSE) 
        })
        timings <- rbind(timings, c("deconv", tt))
    }

    if("infer" %in% extra) {
        cat(">>> inferring extra phenotypes...\n")
        tt <- system.time({
            ngs <- compute.cellcycle.gender(ngs, rna.counts=rna.counts)
        })
        timings <- rbind(timings, c("infer", tt))
    }

    if("drugs" %in% extra) {
        cat(">>> Computing drug enrichment (single)...\n")
        ngs$drugs <- NULL  ## reset??

        tt <- system.time({
            ngs <- compute.drugActivityEnrichment(ngs, lib.dir=lib.dir, combo=FALSE) 
        })
        timings <- rbind(timings, c("drugs", tt))
        
        tt <- system.time({
            ngs <- compute.drugSensitivityEnrichment(
                ngs, lib.dir=lib.dir, ref.db = c("CTRPv2","GDSC"), combo=FALSE) 
        })
        timings <- rbind(timings, c("drugs-sx", tt))

    }

    if("drugs-combo" %in% extra) {
        cat(">>> Computing drug enrichment (combo)...\n")
        tt <- system.time({
            ngs <- compute.drugActivityEnrichment(
                ngs, lib.dir=lib.dir, combo=TRUE) 
        })
        timings <- rbind(timings, c("drugs-combo", tt))

        tt <- system.time({
            ngs <- compute.drugSensitivityEnrichment(
                ngs, lib.dir=lib.dir, ref.db = c("CTRPv2","GDSC"), combo=TRUE) 
        })
        timings <- rbind(timings, c("drugs-sx-combo", tt))
    }
    
    if("graph" %in% extra) {
        cat(">>> computing OmicsGraphs...\n")
        tt <- system.time({
            ngs <- compute.omicsGraphs(ngs) 
        })
        timings <- rbind(timings, c("graph", tt))
    }
    
    if("wordcloud" %in% extra) {
        cat(">>> computing WordCloud statistics...\n")
        tt <- system.time({
            res <- pgx.calculateWordFreq(ngs, progress=NULL, pg.unit=1)        
        })
        timings <- rbind(timings, c("wordcloud", tt))
        ngs$wordcloud <- res
        remove(res)
    }

    if("connectivity" %in% extra) {
        cat(">>> computing connectivity scores...\n")
        ngs$connectivity <- NULL  ## clean up
        ## sigdb.list = c(
        ##     file.path(PGX.DIR,"datasets-allFC.csv"),
        ##     file.path(FILES,"sigdb-archs4.h5")
        ## )
        
        for(db in sigdb) {
            if(file.exists(db)) {
                ntop = 9999
                cat("computing scores for sigDB",db,"\n")                
                tt <- system.time({
                    scores <- pgx.computeConnectivityScores(
                        ngs, db, ntop=ntop, contrasts=NULL)
                })
                timings <- rbind(timings, c("connectivity", tt))
                
                db0 <- sub(".*/","",db)
                ngs$connectivity[[db0]] <- scores
                remove(scores)
            }
        }
        names(ngs$connectivity)        
    }

    ##------------------------------------------------------
    ## pretty collapse all timings
    ##------------------------------------------------------
    ##timings0 <- do.call(rbind, timings)
    timings <- as.matrix(timings)
    rownames(timings) <- timings[,1]
    timings0 <- apply(as.matrix(timings[,-1,drop=FALSE]),2,as.numeric)
    rownames(timings0) <- rownames(timings)
    timings0 <- apply( timings0, 2, function(x) tapply(x,rownames(timings0),sum))
    if(is.null(nrow(timings0))) {
        cn <- names(timings0)
        rn <- unique(rownames(timings))
        timings0 <- matrix(timings0, nrow=1)
        colnames(timings0) <- cn
        rownames(timings0) <- rn[1]
    }
    rownames(timings0) <- paste("[extra]",rownames(timings0))
    
    ngs$timings <- rbind(ngs$timings, timings0)
    
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


compute.drugActivityEnrichment <- function(ngs, lib.dir, combo=TRUE ) {
    ## -------------- drug enrichment
    ##source(file.path(RDIR,"pgx-drugs.R"))
    ##source(file.path(RDIR,"pgx-graph.R", local=TRUE)
    X <- readRDS(file=file.path(lib.dir,"l1000_es.rds"))
    x.drugs <- gsub("_.*$","",colnames(X))
    length(table(x.drugs))
    dim(X)

    res.mono = res.combo = NULL
    
    NPRUNE=-1
    NPRUNE=250
    res.mono <- pgx.computeDrugEnrichment(
        ngs, X, x.drugs, methods=c("GSEA","cor"),
        nprune=NPRUNE, contrast=NULL )

    if(is.null(res.mono)) {
        cat("[compute.drugActivityEnrichment] WARNING:: pgx.computeDrugEnrichment failed!\n")
        return(ngs)
    }

    ## attach annotation
    annot0 <- read.csv(file.path(lib.dir,"L1000_repurposing_drugs.txt"),
                       sep="\t", comment.char="#")
    annot0$drug <- annot0$pert_iname
    rownames(annot0) <- annot0$pert_iname
    head(annot0)
    
    annot1 <- NULL
    if(combo==TRUE) {

        res.combo <- pgx.computeComboEnrichment(
            ngs, X, x.drugs, res.mono=res.mono,
            contrasts = NULL,
            ntop=15, nsample=80, nprune=NPRUNE)
        names(res.combo)
                
        ## create combo annotation table from mono-drug
        combo <- rownames(res.combo$X)
        annot1 <- pgx.createComboDrugAnnot(combo, annot0)        
    }
    dim(res.mono[["GSEA"]]$X)

    ##ngs$drugs <- NULL
    ngs$drugs[["activity/L1000"]]  <- res.mono[["GSEA"]]
    ngs$drugs[["activity/L1000"]][["annot"]] <- annot0[,c("drug","moa","target")]
    
    if(!is.null(res.combo)) {
        ngs$drugs[["activity-combo/L1000"]] <- res.combo
        ngs$drugs[["activity-combo/L1000"]][["annot"]] <- annot1[,c("drug","moa","target")]
        names(ngs$drugs)
    }

    remove(X)
    remove(x.drugs)
    return(ngs)
}

##ref="CTRPv2";lib.dir="../lib";combo=FALSE
compute.drugSensitivityEnrichment <- function(ngs, lib.dir, combo=TRUE, ref.db=c("CTRPv2","GDSC") )
{
    ref <- ref.db[1]
    for(ref in ref.db) {
        ##X <- readRDS(file=file.path(lib.dir,"drugSX-GDSC-t25-g1000.rds"))
        ##X <- readRDS(file=file.path(lib.dir,"drugSX-CTRPv2-t25-g1000.rds"))
        X <- readRDS(file=file.path(lib.dir,paste0("drugSX-",ref,"-t25-g1000.rds")))
        x.drugs <- gsub("@.*$","",colnames(X))
        length(table(x.drugs))
        dim(X)
        
        res.mono = res.combo = NULL
        
        NPRUNE=-1
        NPRUNE=250
        res.mono <- pgx.computeDrugEnrichment(
            ngs, X, x.drugs, methods=c("GSEA","cor"),
            nprune=NPRUNE, contrast=NULL )
        
        if(is.null(res.mono)) {
            cat("[compute.drugActivityEnrichment] WARNING:: pgx.computeDrugEnrichment failed!\n")
            return(ngs)
        }
        
        ## attach annotation
        ##annot0 <- read.csv(file.path(lib.dir,"drugSX-GDSC-drugs.csv"))
        ##annot0 <- read.csv(file.path(lib.dir,"drugSX-CTRPv2-drugs.csv"))
        annot0 <- read.csv(file.path(lib.dir,paste0("drugSX-",ref,"-drugs.csv")))
        head(annot0)
        rownames(annot0) <- annot0$drug
        
        annot1 <- NULL
        if(combo==TRUE) {
            
            res.combo <- pgx.computeComboEnrichment(
                ngs, X, x.drugs, res.mono=res.mono,
                contrasts = NULL,
                ntop=15, nsample=80, nprune=NPRUNE)
            names(res.combo)
            
            ## create combo annotation table from mono-drug
            combo <- rownames(res.combo$X)
            annot1 <- pgx.createComboDrugAnnot(combo, annot0)
        }
        
        s1 <- paste0("sensitivity/",ref)
        s2 <- paste0("sensitivity-combo/",ref)
        ngs$drugs[[s1]] <- res.mono[["GSEA"]]
        ngs$drugs[[s1]][["annot"]] <- annot0[,c("moa","target")]
        
        if(!is.null(res.combo)) {
            ngs$drugs[[s2]] <- res.combo
            ngs$drugs[[s2]][["annot"]] <- annot1[,c("moa","target")]
        }
    } ## end of for rr
    
    names(ngs$drugs)
    
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

