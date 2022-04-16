##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## NOTE: needs global variables file


##-----------------------------------------------------------------------------
## GLOBAL variables
##-----------------------------------------------------------------------------

if(0) {
## Caching the init files
INIT.FILE <- file.path(OPG,"cache/global-init.rda") ## avoid rw permission
##unlink(INIT.FILE)
INIT.FILE

file.exists(INIT.FILE)

if(1 && file.exists(INIT.FILE)) {    

    message("[INIT] loading cached INIT file ",INIT.FILE)
    t0 <- Sys.time()
    load(INIT.FILE, verbose=1)
    message("Loading cache took: ", round(Sys.time() - t0), " seconds")

} else {
    
    message("[INIT] no INIT file! building INIT from scratch.")
    message("[INIT] INIT.FILE = ", INIT.FILE)    
    t0 <- Sys.time()

    oldvars <- ls()

    ## All gene families in Human UPPER CASE
    require(org.Hs.eg.db)
    GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
    GENE.SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = GENE.SYMBOL
    ##GSET.PREFIX.REGEX = paste(paste0("^",GSET.PREFIXES,"_"),collapse="|")
    GSET.PREFIX.REGEX="^BIOCARTA_|^C2_|^C3_|^C7_|^CHEA_|^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^KEA_|^KEGG_|^PID_|^REACTOME_|^ST_"
    GENE.SUMMARY = read.csv(file.path(FILES,"gene-summary.csv"),row.names=1)
    GENE.SUMMARY = array(GENE.SUMMARY[,1], dimnames=list(rownames(GENE.SUMMARY)))
    
    ## GENExGENE <- readRDS(file=file.path(FILES,"GENExGENE-cosSparseKNN500-XL.rds"))
    GSETxGENE <- readRDS(file.path(FILES,"gset-sparseG-XL.rds"))
    load(file.path(FILES,"gmt-all.rda"),verbose=1)
    GSETS = gmt.all;remove(gmt.all)

    message("[INIT] parsing gene families...")
    FAMILIES <- pgx.getGeneFamilies(GENE.SYMBOL, FILES=FILES, min.size=10, max.size=9999)
    fam.file <- file.path(FILES,"custom-families.gmt")
    if(file.exists(fam.file)) {
        custom.gmt = read.gmt(file.path(FILES,"custom-families.gmt"),add.source=TRUE)
        names(custom.gmt)
        FAMILIES= c(FAMILIES, custom.gmt)
    }
    FAMILIES[["<all>"]] <- GENE.SYMBOL
    f1 <- FAMILIES
    names(f1) <- paste0("FAMILY:",names(f1))
    names(f1) <- sub("FAMILY:<all>","<all>",names(f1))
    GSETS <- c(GSETS,f1)

    ## convert to integer list (more efficient)
    message("[INIT] converting GSETS to list of integers...")
    GSET.GENES <- sort(unique(unlist(GSETS)))  ## slow...
    iGSETS <- parallel::mclapply(GSETS, function(a) match(a,GSET.GENES))  ## slow...
    names(iGSETS) <- names(GSETS)
    getGSETS <- function(gs) {
        lapply(iGSETS[gs],function(i) GSET.GENES[i])
    }
        
    message("[INIT] parsing collections...")
    COLLECTIONS <- pgx.getGeneSetCollections(names(GSETS), min.size=10, max.size=99999)
    COLLECTIONS <- COLLECTIONS[order(names(COLLECTIONS))]

    remove(list=c("custom.gmt","f1","GSETS"))

    ##-----------------------------------------------------------------------------
    ## TISSUE/REFERENCE data sets
    ##-----------------------------------------------------------------------------
    load(file.path(FILES,"sig/rna_tissue.rda"))  ## TISSUE and TISSUE.grp

    ##-----------------------------------------------------------------------------
    ## Immune cell markers
    ##-----------------------------------------------------------------------------

    ## Really needed???
    IMMPROT <- read.csv(file.path(FILES,"sig/ImmProt-signature.csv"),row.names=1)
    IMMPROT_MARKERS <- rownames(read.csv(file.path(FILES,"sig/immprot-signature1000.csv"),row.names=1))
    DICE_MARKERS <- rownames(read.csv(file.path(FILES,"sig/DICE-signature1000.csv"),row.names=1))
    LM22 <- read.csv(file.path(FILES,"sig/LM22.txt"),sep="\t",row.names=1)
    LM22_MARKERS <- rownames(LM22)

    ##-----------------------------------------------------------------------------
    ## Meta MA-profiles (fold changes) of all experiments
    ##-----------------------------------------------------------------------------
    ##load( file.path(FILES,"allMA-pub.rda"), verbose=1)
    ##load(file.path(FILES,"allFoldChanges-pub-8k.rda"))
    ##PROFILES <- list(M=allM, A=allA, FC=allFC)
    ##allFC <- pgx.readDatasetProfiles(PGX.DIR, file="datasets-allFC.csv") 
    ##PROFILES <- list(FC=allFC)
    ##remove(allA)
    ##remove(allM)
    ##remove(allFC)

    ##-----------------------------------------------------------------------------
    ## Colors
    ##-----------------------------------------------------------------------------
        
    COLORS = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
    COLORS = rep(c(ggsci::pal_npg("nrc", alpha = 0.7)(10),
                   ggsci::pal_aaas("default", alpha = 0.7)(10),
                   ggsci::pal_d3("category10", alpha = 0.7)(10)),99)
##    BLUERED <- grDevices::colorRampPalette(
##        rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#EEEEEE",
##              "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
##    BLUERED <- function(n=64) suppressWarnings(gplots::colorpanel(n,low="royalblue3",mid="grey90",high="indianred3"))
    BLUERED <- colorRampPalette(c("royalblue3","grey90","indianred3"))
    PURPLEYELLOW <- colorRampPalette(c("purple","purple3","black","yellow3","yellow"))
    PURPLEYELLOW <- colorRampPalette(c("purple","purple4","black","yellow4","yellow"))

    newvars <- setdiff(ls(), oldvars)
    newvars

    message("Creating global init took: ", round(Sys.time() - t0), " seconds")
    message("[INIT] saving INIT file ", INIT.FILE)    
    save( list=newvars, file=INIT.FILE)
}    
}

pgx.initialize <- function(pgx) {

    ##---------------------------------------------------------------------
    ## This function must be called after creation of a PGX object
    ## and include some cleaning up and updating some internal
    ## structures to keep compatibility with new/old versions.
    ##---------------------------------------------------------------------
    message("[pgx.initialize] initializing pgx object")

    ##----------------- check object
    obj.needed <- c("genes", ## "deconv","collections", "families", "counts",
                    "GMT","gset.meta","gsetX","gx.meta","model.parameters",
                    "samples","tsne2d","X")
    all(obj.needed %in% names(pgx))
    if(!all(obj.needed %in% names(pgx))) {
        obj.missing <- setdiff(obj.needed, names(pgx))
        msg <- paste("invalid pgx object. missing parts in object: ",obj.missing)
        shiny::showNotification(msg,duration=NULL,type="error")
        ##stop(msg)
        return(NULL)
    }

    vars.needed <- c("group")
    if(FALSE && !all(vars.needed %in% colnames(pgx$samples))) {
        vars.missing <- setdiff(vars.needed, colnames(pgx$samples))
        msg <- paste("invalid pgx object. missing variables in object: ",vars.missing)
        shiny::showNotification(msg,duration=NULL,type="error")
        ##stop(msg)
        return(NULL)
    }
    
    ## for COMPATIBILITY: if no counts, estimate from X
    if(is.null(pgx$counts)) {
        cat("WARNING:: no counts table. estimating from X\n")
        ##pgx$counts <- (2**pgx$X-1) ##
        pgx$counts = pmax(2**pgx$X - 1,0)
        k = grep("lib.size|libsize",colnames(pgx$samples))[1]
        if(length(k)>0) {
            libsize = pgx$samples[colnames(pgx$counts),k]
            libsize
            pgx$counts = Matrix::t(Matrix::t(pgx$counts) * libsize)
        }
    }
    pgx$counts <- as.matrix(pgx$counts)
    if(!is.null(pgx$X)) pgx$X <- as.matrix(pgx$X)
    
    ##----------------------------------------------------------------
    ## model parameters
    ##----------------------------------------------------------------
    has.design <- !is.null(pgx$model.parameters$design)
    has.expmatrix <- !is.null(pgx$model.parameters$exp.matrix)
    if(!"group" %in% names(pgx$model.parameters) && has.design) {
        ii <- max.col(pgx$model.parameters$design)
        group <- colnames(pgx$model.parameters$design)[ii]
        names(group) <- rownames(pgx$model.parameters$design)
        pgx$model.parameters$group <- group
    }
    if(!"group" %in% names(pgx$model.parameters) && has.expmatrix) {
        group <- pgx.getConditions(pgx$model.parameters$exp.matrix,nmax=0)
        names(group) <- rownames(pgx$model.parameters$exp.matrix)
        pgx$model.parameters$group <- group
    }
    
    if(is.null(pgx$model.parameters$group)) {
        stop("[pgx.initialize] FATAL: group is null!!!")
    }

    ##----------------------------------------------------------------
    ## Convert to labeled contrast matrix (new style)
    ##----------------------------------------------------------------

    ## don't add if not exists for now...
    if(FALSE && !("contrasts" %in% names(pgx))) {
        design <- pgx$model.parameters$design
        expmat <- pgx$model.parameters$exp.matrix
        contr.mat <- pgx$model.parameters$contr.mat
        if(is.null(expmat)) {
            expmat <- design %*% contr.mat
        }
        pgx$contrasts <- contrastAsLabels(expmat)
    }

    is.numlev <- all(unique(pgx$contrasts) %in% c(NA,"",-1,0,1))
    is.samplewise <- all(rownames(pgx$contrasts) == rownames(pgx$samples))
    is.samplewise
    if("contrasts" %in% names(pgx) && (!is.samplewise || is.numlev)) {
        design <- pgx$model.parameters$design
        expmat <- pgx$model.parameters$exp.matrix
        contr.mat <- pgx$model.parameters$contr.matrix
        new.contr <- pgx$contrasts
        is.numlev <- all(unique(new.contr) %in% c(NA,"",-1,0,1))
        is.numlev <- is.numlev && (-1 %in% new.contr)  ## must have -1 !!
        if(is.numlev) {
            new.contr <- contrastAsLabels(new.contr)
        }
        is.groupwise <- all(rownames(new.contr) %in% pgx$samples$group)
        is.groupwise
        if(is.groupwise) {
            grp <- as.character(pgx$samples$group)
            new.contr <- new.contr[grp,,drop=FALSE]
            rownames(new.contr) <- rownames(pgx$samples)
        }
        pgx$contrasts <- new.contr
    }

    
    ##----------------------------------------------------------------
    ## Tidy up phenotype matrix (important!!!): get numbers/integers
    ## into numeric, categorical into factors....
    ##----------------------------------------------------------------
    ##pgx$samples <- tidy.dataframe(pgx$samples)  ## warning!! this converts all to CHR!!
    pgx$samples <- type.convert(pgx$samples)    ## autoconvert to datatypes
    pgx$samples <- pgx$samples[,which(colMeans(is.na(pgx$samples))<1),drop=FALSE]

    is.num  <- sapply(pgx$samples,class) %in% c('numeric','integer')
    numlev  <- apply(pgx$samples,2,function(x) length(unique(x[!is.na(x)])))
    is.numfac <- (is.num & numlev <= 3)
    is.numfac
    if(any(is.numfac)) {
        for(i in which(is.numfac)) pgx$samples[,i] <- as.character(pgx$samples[,i])
    }
    
    ## clean up: pgx$Y is a cleaned up pgx$samples
    kk = grep("batch|lib.size|norm.factor|repl|donor|clone|sample|barcode",
              colnames(pgx$samples),invert=TRUE,value=TRUE)
    kk = grep("lib.size|norm.factor|donor|clone|barcode",
              colnames(pgx$samples),invert=TRUE,value=TRUE)
    pgx$Y = pgx$samples[colnames(pgx$X),kk,drop=FALSE]
    pgx$Y <- type.convert(pgx$Y)   ## autoconvert to datatypes
    
    ## *****************************************************************
    ## ******************NEED RETHINK***********************************
    ## *****************************************************************
    ## ONLY categorical variables for the moment!!!
    ny1 <- nrow(pgx$Y)-1
    k1 = pgx.getCategoricalPhenotypes(pgx$Y, min.ncat=2, max.ncat=ny1)  ## exclude
    k2 = grep("OS.survival|cluster|condition|group",colnames(pgx$Y),value=TRUE) ## must include
    ##kk = sort(unique(c("group",k1,k2)))
    kk = sort(unique(c(k1,k2)))
    pgx$Y <- pgx$Y[,kk,drop=FALSE]
    colnames(pgx$Y)
    ## pgx$samples <- pgx$Y    ## REALLY? !!!!!!!!!!!!!!!!!!!!

    ##----------------------------------------------------------------
    ## Tidy up genes matrix
    ##----------------------------------------------------------------
    pgx$genes = pgx$genes[rownames(pgx$counts),,drop=FALSE]
    pgx$genes$gene_name = as.character(pgx$genes$gene_name)
    pgx$genes$gene_title = as.character(pgx$genes$gene_title)

    ## Add chromosome annotation if not
    if(!("chr" %in% names(pgx$genes))) {
        require(org.Hs.eg.db)
        symbol = sapply(as.list(org.Hs.egSYMBOL),"[",1)  ## some have multiple chroms..
        CHR = sapply(as.list(org.Hs.egCHR),"[",1)  ## some have multiple chroms..
        MAP <- sapply(as.list(org.Hs.egMAP),"[",1)  ## some have multiple chroms..
        names(CHR) = names(MAP) = symbol
        pgx$genes$chr <- CHR[pgx$genes$gene_name]
        pgx$genes$map <- MAP[pgx$genes$gene_name]
    }

    ##-----------------------------------------------------------------------------
    ## intersect and filter gene families (convert species to human gene sets)
    ##-----------------------------------------------------------------------------
    if("hgnc_symbol" %in% colnames(pgx$genes) ) {
        hgenes <- toupper(pgx$genes$hgnc_symbol)
        genes  <- pgx$genes$gene_name
        pgx$families <- lapply(FAMILIES, function(x) setdiff(genes[match(x,hgenes)],NA))
    } else {
        genes <- toupper(pgx$genes$gene_name)
        pgx$families <- lapply(FAMILIES, function(x) intersect(x,genes))
    }
    famsize <- sapply(pgx$families, length)
    pgx$families <- pgx$families[which(famsize>=10)]
    
    all.genes <- sort(rownames(pgx$genes))
    pgx$families[["<all>"]] <- all.genes
    ## rownames(pgx$GMT) <- toupper(rownames(pgx$GMT)) ## everything to human...
    
    ##-----------------------------------------------------------------------------
    ## Recompute geneset meta.fx as average fold-change of genes
    ##-----------------------------------------------------------------------------
    message("[pgx.initialize] Recomputing geneset fold-changes")
    nc <- length(pgx$gset.meta$meta)
    i=1
    for(i in 1:nc) {
        gs <- pgx$gset.meta$meta[[i]]
        fc <- pgx$gx.meta$meta[[i]]$meta.fx
        names(fc) <- rownames(pgx$gx.meta$meta[[i]])
        fc <- fc[which(toupper(names(fc)) %in% colnames(GSETxGENE))]
        ## G1 <- GSETxGENE[rownames(gs),toupper(names(fc))]
        G1 <- Matrix::t(pgx$GMT[names(fc),rownames(gs)])
        mx <- (G1 %*% fc)[,1]
        pgx$gset.meta$meta[[i]]$meta.fx <- mx
    }

    ##-----------------------------------------------------------------------------
    ## Recode survival
    ##-----------------------------------------------------------------------------
    pheno <- colnames(pgx$Y)
    ## DLBCL coding
    if(("OS.years" %in% pheno && "OS.status" %in% pheno)) {
        message("found OS survival data")
        event <- ( pgx$Y$OS.status %in% c("DECEASED","DEAD","1","yes","YES","dead"))
        pgx$Y$OS.survival <- ifelse(event, pgx$Y$OS.years, -pgx$Y$OS.years)            
    }

    ## cBioportal coding
    if(("OS_MONTHS" %in% pheno && "OS_STATUS" %in% pheno)) {
        message("[pgx.initialize] found OS survival data\n")
        event <- ( pgx$Y$OS_STATUS %in% c("DECEASED","DEAD","1","yes","YES","dead"))
        pgx$Y$OS.survival <- ifelse(event, pgx$Y$OS_MONTHS, -pgx$Y$OS_MONTHS)            
    }

    ##-----------------------------------------------------------------------------
    ## Check if clustering is done
    ##-----------------------------------------------------------------------------
    message("[pgx.initialize] Check if clustering is done...")
    if(!"cluster.genes" %in% names(pgx)) {
        message("[pgx.initialize] clustering genes...")
        pgx <- pgx.clusterGenes(pgx, methods='umap', dims=c(2), level='gene')
        pgx$cluster.genes$pos <- lapply( pgx$cluster.genes$pos, pos.compact )
    }
    if(!"cluster.gsets" %in% names(pgx)) {
        message("[pgx.initialize] clustering genesets...")
        pgx <- pgx.clusterGenes(pgx, methods='umap', dims=c(2), level='geneset')
        pgx$cluster.gsets$pos  <- lapply( pgx$cluster.gsets$pos, pos.compact )
    }

    ##-----------------------------------------------------------------------------
    ## Remove redundant???
    ##-----------------------------------------------------------------------------
    message("[pgx.initialize] Remove redundant phenotypes...")
    if(".gender" %in% colnames(pgx$Y) &&
        any(c("gender","sex") %in% tolower(colnames(pgx$Y)))) {
        pgx$Y$.gender <- NULL
    }
        
    ##-----------------------------------------------------------------------------
    ## Keep compatible with OLD formats
    ##-----------------------------------------------------------------------------
    message("[pgx.initialize] Keep compatible OLD formats...")
    if( any(c("mono","combo") %in% names(pgx$drugs)) ) {
        dd <- pgx$drugs[["mono"]]
        aa1 <- pgx$drugs[["annot"]]
        if(is.null(aa1)) {
            aa1 <- read.csv(file.path(FILES,"L1000_repurposing_drugs.txt"),
                            sep="\t", comment.char="#")
            aa1$drug <- aa1$pert_iname
            rownames(aa1) <- aa1$pert_iname
        }
        dd[["annot"]] <- aa1
        pgx$drugs[["activity/L1000"]] <- dd
        if("combo" %in% names(pgx$drugs)) {
            dd2 <- pgx$drugs[["combo"]]
            combo <- rownames(dd2$X)
            aa2 <- pgx.createComboDrugAnnot(combo, aa1)             
            dd2[["annot"]] <- aa2
            pgx$drugs[["activity-combo/L1000"]] <- dd2
        }
        pgx$drugs$mono  <- NULL
        pgx$drugs$annot <- NULL
        pgx$drugs$combo <- NULL        
    }
    
    ##-----------------------------------------------------------------------------
    ## remove large deprecated outputs from objects
    ##-----------------------------------------------------------------------------
    message("[pgx.initialize] Removing deprecated objects...")
    pgx$gx.meta$outputs <- NULL
    pgx$gset.meta$outputs <- NULL
    pgx$gmt.all <- NULL

    message("[pgx.initialize] done!")
    return(pgx)
}
