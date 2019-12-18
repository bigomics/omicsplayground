
library(survival)
library(knitr)
library(shiny)
library(shinyjs)
library(shinyBS)
library(shinyjqui)
library(rmarkdown)
library(shinycssloaders)
library(dragulaR)
library(shinyWidgets)

library(scatterD3)
library(fastcluster)
library(ComplexHeatmap)
library(plotly)
library(Matrix)
library(igraph)
library(DT)
library(ggplot2)
## library(Cairo)

##useShinyjs(rmd=TRUE)  
useShinyjs()
ht_global_opt(fast_hclust = TRUE)
options(shiny.maxRequestSize = 200*1024^2)  ## max 200Mb upload

source(file.path(RDIR,"gx-heatmap.r"))
source(file.path(RDIR,"gx-plot.r"))
source(file.path(RDIR,"gx-limma.r"))
source(file.path(RDIR,"gx-volcano.r"))
source(file.path(RDIR,"gx-combat.r"))
source(file.path(RDIR,"gx-util.r"))
       
source(file.path(RDIR,"gset-gsea.r"))
source(file.path(RDIR,"gset-fisher.r"))
source(file.path(RDIR,"gset-meta.r"))

source(file.path(RDIR,"ngs-cook.r"))
source(file.path(RDIR,"ngs-fit.r"))

source(file.path(RDIR,"pgx-functions.R"))
source(file.path(RDIR,"pgx-graph.R"))
source(file.path(RDIR,"pgx-deconv.R"))
source(file.path(RDIR,"pgx-cna.R"))
source(file.path(RDIR,"pgx-plotting.R"))
source(file.path(RDIR,"pgx-correct.R"))
source(file.path(RDIR,"pgx-predict.R"))
source(file.path(RDIR,"pgx-links.R"))
source(file.path(RDIR,"pgx-modules.R"))
source(file.path(RDIR,"pgx-upload.R"))
source(file.path(RDIR,"pgx-proteomics.R"))
source(file.path(RDIR,"pgx-drugs.R"))

source(file.path(RDIR,"xcr-graph.r"))
source(file.path(RDIR,"ui-code.R"))

source( file.path(RDIR,"compute2-genes.R"))
source( file.path(RDIR,"compute2-genesets.R"))
source( file.path(RDIR,"compute2-extra.R"))

## you need to override this!!!
PRO.VERSION=FALSE
DEV.VERSION=FALSE

## some custom code
##code.textInput

##-----------------------------------------------------------------------------
## Added GLOBAL info
##-----------------------------------------------------------------------------

## Caching the init files
INIT.FILE <- file.path(FILES,"global-init.rda")
INIT.FILE

if( file.exists(INIT.FILE)) {    
    
    load(INIT.FILE, verbose=1)

} else {

    oldvars <- ls()

    ## All gene families in Human UPPER CASE
    require(org.Hs.eg.db)
    GENE.TITLE  = unlist(as.list(org.Hs.egGENENAME))
    GENE.SYMBOL = unlist(as.list(org.Hs.egSYMBOL))
    names(GENE.TITLE) = GENE.SYMBOL
    ##GSET.PREFIX.REGEX = paste(paste0("^",GSET.PREFIXES,"_"),collapse="|")
    GSET.PREFIX.REGEX="^BIOCARTA_|^C2_|^C3_|^C7_|^CHEA_|^GOBP_|^GOCC_|^GOMF_|^HALLMARK_|^KEA_|^KEGG_|^PID_|^REACTOME_|^ST_"

    ##xGENExGENE <- readRDS(file=file.path(FILES,"GENExGENE-cosSparseKNN500-XL.rds"))
    GSETxGENE <- readRDS(file.path(FILES,"gset-sparseG-XL.rds"))
    load(file.path(FILES,"gmt-all.rda"),verbose=1)
    GSETS = gmt.all;remove(gmt.all)
    ##saveRDS(gmt.all, file.path(FILES,"gmt-all.rds"))
    ##GSETS <- readRDS(file.path(FILES,"gmt-all.rds"))

    cat("<init> parsing gene families...\n")
    FAMILIES <- pgx.getGeneFamilies(GENE.SYMBOL, FILES=FILES, min.size=10, max.size=9999)
    ##FAMILIES <- c(FAMILIES, list( "<LM22 markers>"=LM22_MARKERS,"<ImmProt markers>"=IMMPROT_MARKERS))
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

    cat("<init> parsing collections...\n")
    COLLECTIONS <- pgx.getGeneSetCollections(names(GSETS), min.size=10, max.size=99999)
    COLLECTIONS <- COLLECTIONS[order(names(COLLECTIONS))]

    remove(list=c("custom.gmt","f1"))

    ##-----------------------------------------------------------------------------
    ## TISSUE/REFERENCE data sets
    ##-----------------------------------------------------------------------------

    load(file.path(FILES,"rna_tissue.rda"))  ## TISSUE and TISSUE.grp
    IMMPROT <- read.csv(file.path(FILES,"ImmProt-signature.csv"),row.names=1)

    ##-----------------------------------------------------------------------------
    ## Immune cell markers
    ##-----------------------------------------------------------------------------

    IMMPROT_MARKERS <- rownames(read.csv(file.path(FILES,"immprot-signature1000.csv"),row.names=1))
    DICE_MARKERS <- rownames(read.csv(file.path(FILES,"DICE-signature1000.csv"),row.names=1))
    LM22 <- read.csv(file.path(FILES,"LM22.txt"),sep="\t",row.names=1)
    LM22_MARKERS <- rownames(LM22)

    ##-----------------------------------------------------------------------------
    ## Meta MA-profiles (fold changes) of all experiments
    ##-----------------------------------------------------------------------------
    ##load( file.path(FILES,"allMA-pub.rda"), verbose=1)
    load(file.path(FILES,"allFoldChanges-pub-8k.rda"))
    ##PROFILES <- list(M=allM, A=allA, FC=allFC)
    PROFILES <- list(FC=allFC)

    ##remove(allA)
    ##remove(allM)
    remove(allFC)

    ##-----------------------------------------------------------------------------
    ## Colors
    ##-----------------------------------------------------------------------------
    library(ggsci)
    library(RColorBrewer)
    COLORS = rep(brewer.pal(8,"Set2"),99)
    COLORS = rep(c(pal_npg("nrc", alpha = 0.7)(10),
                   pal_aaas("default", alpha = 0.7)(10),
                   pal_d3("category10", alpha = 0.7)(10)),99)
    BLUERED <- colorRampPalette(
        rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#EEEEEE",
              "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))
    PURPLEYELLOW <- colorRampPalette(c("purple","purple3","black","yellow3","yellow"))
    PURPLEYELLOW <- colorRampPalette(c("purple","purple4","black","yellow4","yellow"))

    newvars <- setdiff(ls(), oldvars)
    newvars
    save( list=newvars, file=file.path(FILES, INIT.FILE))
    
}

##================================================================================
##========================== FUNCTIONS ===========================================
##================================================================================

pgx.getFamilies <- function(ngs, nmin=10, extended=FALSE) {
    if(extended) {
        fam <- grep("^[<].*|^FAMILY|^TISSUE|^COMPARTMENT|^CELLTYPE|^GOCC|^DISEASE|^CUSTOM",
                    names(GSETS),value=TRUE)
        fam <- grep("^[<].*|^FAMILY|^COMPARTMENT|^CUSTOM",names(GSETS),value=TRUE)
    } else {
        fam <- grep("^[<].*|^FAMILY|^CUSTOM",names(GSETS),value=TRUE)
    }
    xgenes <- toupper(rownames(ngs$X))
    xgenes <- toupper(ngs$genes$gene_name)
    jj <- which(sapply(GSETS[fam],function(x) sum(x %in% xgenes)) >= nmin)
    sort(fam[jj])
}

pgx.initialize <- function(ngs) {

    cat("INFO <init:initialize> initializing ngs object for the Playground\n")

    ##----------------- check object
    obj.needed <- c("genes", ## "deconv","collections", "families", "counts",
                    "GMT","gset.meta","gsetX","gx.meta","model.parameters",
                    "samples","tsne2d","X")
    all(obj.needed %in% names(ngs))
    if(!all(obj.needed %in% names(ngs))) {
        obj.missing <- setdiff(obj.needed, names(ngs))
        msg <- paste("invalid ngs object. missing parts in object: ",obj.missing)
        showNotification(msg,duration=NULL,type="error")
        stop(msg)
        return(NULL)
    }

    ## for COMPATIBILITY: if no counts, estimate from X
    if(is.null(ngs$counts)) {
        cat("WARNING:: no counts table. estimating from X\n")
        ##ngs$counts <- (2**ngs$X-1) ##
        ngs$counts = pmax(2**ngs$X - 1,0)
        k = grep("lib.size|libsize",colnames(ngs$samples))[1]
        if(length(k)>0) {
            libsize = ngs$samples[colnames(ngs$counts),k]
            libsize
            ngs$counts = t(t(ngs$counts) * libsize)
        }
    }
    ngs$counts <- as.matrix(ngs$counts)
    if(!is.null(ngs$X)) ngs$X <- as.matrix(ngs$X)
    
    ##----------------------------------------------------------------
    ## Tidy up phenotype matrix (important!!!): get numbers/integers
    ## into numeric, categorical into factors....
    ##----------------------------------------------------------------
    ngs$samples <- tidy.dataframe(ngs$samples)  ## warning!! this converts all to CHR!!

    ## clean up: ngs$Y is a cleaned up ngs$samples
    ngs$samples$barcode <- NULL
    ngs$samples <- ngs$samples[,which(colMeans(is.na(ngs$samples))<1),drop=FALSE]
    kk = grep("group|batch|lib.size|norm.factor|repl|donor|clone|sample|barcode",
              colnames(ngs$samples),invert=TRUE)
    kk <- unique( c(grep("^group$",colnames(ngs$samples)),kk))
    ngs$Y = ngs$samples[colnames(ngs$X),kk,drop=FALSE]
    ngs$Y <- tidy.dataframe(ngs$Y) ## NEED CHECK!!!
    
    ##----------------------------------------------------------------
    ## Tidy up genes matrix
    ##----------------------------------------------------------------
    ngs$genes = ngs$genes[rownames(ngs$counts),,drop=FALSE]
    ngs$genes$gene_name = as.character(ngs$genes$gene_name)
    ngs$genes$gene_title = as.character(ngs$genes$gene_title)

    ## Add chromosome annotation if not
    if(!("chr" %in% names(ngs$genes))) {
        symbol = sapply(as.list(org.Hs.egSYMBOL),"[",1)  ## some have multiple chroms..
        CHR = sapply(as.list(org.Hs.egCHR),"[",1)  ## some have multiple chroms..
        MAP <- sapply(as.list(org.Hs.egMAP),"[",1)  ## some have multiple chroms..
        names(CHR) = names(MAP) = symbol
        ngs$genes$chr <- CHR[ngs$genes$gene_name]
        ngs$genes$map <- MAP[ngs$genes$gene_name]
    }

    ##-----------------------------------------------------------------------------
    ## intersect and filter gene families (convert species to human gene sets)
    ##-----------------------------------------------------------------------------
    if("hgnc_symbol" %in% colnames(ngs$genes) ) {
        hgenes <- toupper(ngs$genes$hgnc_symbol)
        genes  <- ngs$genes$gene_name
        ngs$families <- lapply(FAMILIES, function(x) setdiff(genes[match(x,hgenes)],NA))
    } else {
        genes <- toupper(ngs$genes$gene_name)
        ngs$families <- lapply(FAMILIES, function(x) intersect(x,genes))
    }
    famsize <- sapply(ngs$families, length)
    ngs$families <- ngs$families[which(famsize>=10)]
    
    all.genes <- sort(rownames(ngs$genes))
    ngs$families[["<all>"]] <- all.genes
    
    ##-----------------------------------------------------------------------------
    ## Recode survival
    ##-----------------------------------------------------------------------------
    pheno <- colnames(ngs$Y)
    ## DLBCL coding
    if(("OS.years" %in% pheno && "OS.status" %in% pheno)) {
        cat("found OS survival data\n")
        event <- ( ngs$Y$OS.status %in% c("DECEASED","DEAD","1","yes","YES","dead"))
        ngs$Y$OS.survival <- ifelse(event, ngs$Y$OS.years, -ngs$Y$OS.years)            
    }

    ## cBioportal coding
    if(("OS_MONTHS" %in% pheno && "OS_STATUS" %in% pheno)) {
        cat("found OS survival data\n")
        event <- ( ngs$Y$OS_STATUS %in% c("DECEASED","DEAD","1","yes","YES","dead"))
        ngs$Y$OS.survival <- ifelse(event, ngs$Y$OS_MONTHS, -ngs$Y$OS_MONTHS)            
    }

    ##-----------------------------------------------------------------------------
    ## Remove redundant???
    ##-----------------------------------------------------------------------------
    if(".gender" %in% colnames(ngs$Y) &&
        any(c("gender","sex") %in% tolower(colnames(ngs$Y)))) {
        ngs$Y$.gender <- NULL
    }
    
    ## *****************************************************************
    ## ******************NEED RETHINK***********************************
    ## *****************************************************************
    ## ONLY categorical variables for the moment!!!
    k1 = pgx.getCategoricalPhenotypes(ngs$Y, min.ncat=2, max.ncat=20)
    k2 = grep("OS.survival",colnames(ngs$Y),value=TRUE)
    kk = sort(unique(c("group",k1,k2)))
    ngs$Y <- ngs$Y[,kk]
    colnames(ngs$Y)
    ngs$samples <- ngs$Y    ## REALLY?
    
    ##-----------------------------------------------------------------------------
    ## remove large deprecated outputs from objects
    ##-----------------------------------------------------------------------------
    ngs$gx.meta$outputs <- NULL
    ngs$gset.meta$outputs <- NULL
    ngs$gmt.all <- NULL
    return(ngs)
}

dbg <- function(... ) {
    if(DEV.VERSION) {
        ##msg = paste0(ifelse(is.null(module),"",paste0("<",module,"> ")),msg)
        msg = sapply( list(...),paste,collapse=" ")
        cat(paste0("DBG ",sub("\n$","",paste(msg,collapse=" ")),"\n"))
    }
}
