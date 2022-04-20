##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

DECONV.METHODS = c("I-NNLS","CIBERSORT","DCQ","DeconRNAseq","EPIC","NNLM",
                   "cor","SingleR")

##low.th=0.01;add.unknown=1;collapse="no";min.prob=0.2
pgx.inferCellType <- function(counts, low.th=0.01, add.unknown=FALSE,
                              min.prob=0.2, min.count=3, scalex=FALSE,
                              normalize.mat = TRUE, method="NNLM",
                              collapse="max", markers=NULL)
{
    ## infer cell type from markers
    if(is.null(markers)) {
        ##M = read.csv(file.path(FILES,"sig/LM22.txt"),row.names=1,sep="\t",check.names=FALSE)
        ##colnames(M) <- gsub(" ","_",sub(" ",".",sub(" ","_",colnames(M))))
        ##colnames(M) <- sub("Macrophages_","Macrophages.",colnames(M))
        message("[pgx.inferCellType] using database: 'signature-immuneMeta.csv'")
        M = read.csv(file.path(FILES,"signature-immuneMeta.csv"),row.names=1,check.names=FALSE)
        M <- as.matrix(M)
    } else {
        message("[pgx.inferCellType] using provided markers list")
        marker.genes <- gsub("[-+]*$","",unlist(markers))
        marker.genes <- sort(unique(marker.genes))
        M <- matrix(0, nrow=length(marker.genes), ncol=length(markers))
        dimnames(M) <- list(marker.genes, names(markers))
        dim(M)
        k=7
        k=1
        for(k in 1:ncol(M)) {
            mm <- markers[[k]]
            ##m.pos <- sub("[+]$","",grep("[+]$",mm,value=TRUE))
            m.neg <- sub("[-]$","",grep("[-]$",mm,value=TRUE))
            m.pos <- setdiff(gsub("[+-]$","",mm),m.neg)
            if(length(m.pos)) M[match(m.pos,rownames(M)),k] <- +1
            if(length(m.neg)) M[match(m.neg,rownames(M)),k] <- -1
        }
        M <- pmax(M,0)  ## sorry.. no negative markers yet in this function
        M <- as.matrix(M)
    }

    ## Filter count matrix 
    X <- counts
    X <- X[(rowMeans(X >= min.count) > low.th),] ## OK???

    ## Match matrices
    rownames(X) <- toupper(rownames(X))
    rownames(M) <- toupper(rownames(M))
    gg <- intersect(rownames(M),rownames(X))
    X1 <- X[gg,]
    M1 <- M[gg,]
    M1 <- M1[,Matrix::colSums(M1!=0)>0,drop=FALSE]
    if(scalex) X1 <- X1 / (1 + rowMeans(X1))
    
    ## run deconvolution algorithm
    ##add.unknown=0;method="SingleR"
    out <- pgx.deconvolution(X1, ref=M1, methods=method,
                             add.unknown=add.unknown, normalize.mat=normalize.mat)
    names(out$results)
    P <- out$results[[method]]  ## choose specified method
    rownames(P) <- colnames(X1)
    P[is.na(P)] <- 0
    dim(P)
    
    ## Collapse to single cell.type
    if(collapse=="sum") {
        P <- tapply(1:ncol(P), colnames(P), function(i) rowSums(P[,i,drop=FALSE]))
    } else if(collapse=="mean") {
        P <- tapply(1:ncol(P), colnames(P), function(i) rowMeans(P[,i,drop=FALSE]))
    ## } else if(collapse=="max") {
    } else {
        P <- tapply(1:ncol(P), colnames(P), function(i) apply(P[,i,drop=FALSE],1,max))
    }
    P <- do.call(cbind,P)    
    dim(P)
    
    ## Get maximum probability cell.type
    P <- P / (1e-6 + rowSums(P,na.rm=TRUE))
    pmax <- apply(P,1,max,na.rm=TRUE)
    sel.dodgy <- (pmax < min.prob) ## not reliable
    celltype <- colnames(P)[max.col(P)]
    if(any(sel.dodgy)) celltype[which(sel.dodgy)] <- NA
    
    ## collapse small groups to 'other_cells'
    low.ct <- names(which(table(celltype) < low.th*length(celltype)))
    celltype[celltype %in% low.ct] <- "other_cells"
    names(celltype) <- colnames(counts)

    res <- list( celltype=celltype, probs=P)
    res
}

##low.th=0.01;add.unknown=0;celltype0=NULL;min.prob=0.2;min.count=3
pgx.inferCellTypeLM22 <- function(counts, low.th=0.01, add.unknown=FALSE,
                                  normalize.mat = TRUE, method="NNLM",
                                  min.count=3, celltype0=NULL, min.prob=0.2)
{
    ## Two-pass (2-level) cell type identification using LM22
    ##
    ##
    M = read.csv(file.path(FILES,"sig/LM22.txt"),row.names=1,sep="\t",check.names=FALSE)
    M <- as.matrix(M)
    colnames(M) <- gsub(" ","_",sub(" ",".",sub(" ","_",colnames(M))))
    colnames(M) <- sub("Macrophages_","Macrophages.",colnames(M))

    P0 <- NULL
    if(is.null(celltype0)) {
        res <- pgx.inferCellType(counts, low.th=0, add.unknown=add.unknown,
                                 normalize.mat = TRUE, method="NNLM",
                                 min.prob = min.prob)
        celltype0 <- res$celltype
        P0 <- res$probs
        table(celltype0)
    }

    ## Filter count matrix
    X <- counts
    X <- X[(rowMeans(X >= min.count) > low.th),] ## OK???
    dim(X)
    
    ## Match matrices
    rownames(X) <- toupper(rownames(X))
    rownames(M) <- toupper(rownames(M))
    gg <- intersect(rownames(X),rownames(M))
    X <- X[gg,]
    M <- M[gg,]
    dim(X)
    
    ## 2nd stage
    table(celltype0)
    celltype <- rep(NA, ncol(X))
    ct="Macrophages"
    ct="Monocytes"
    ct="T_cells"
    P1 <- matrix(NA,nrow=ncol(X),ncol=0)
    rownames(P1) <- colnames(X)
    for(ct in unique(celltype0)) {
        msel <- grep(ct,colnames(M))
        M1 <- M[,msel,drop=FALSE]
        jj <- which(celltype0==ct)
        jj <- 1:length(celltype0)
        ct3 <- celltype0[jj]
        dim(M1)
        if(ncol(M1)>0) {
            X1 <- X[,jj,drop=FALSE]
            X1 <- X1 / (1e-3 + rowMeans(X1)) ## center feature means??
            M1 <- M1 / (1e-3 + rowMeans(M1)) ## center feature means??
            res1 <- pgx.deconvolution(
                X1, ref=M1, methods="NNLM", normalize.mat = TRUE, 
                add.unknown=FALSE )
            C1 <- res1$results[["NNLM"]]
            C1 <- C1 / rowSums(C1)
            ct3 <- colnames(C1)[max.col(C1)]
            C0 <- P0[colnames(X1),ct]
            C1 <- C1 * C0
            C1 <- C1[match(rownames(P1),rownames(C1)),,drop=FALSE]
            C1[is.na(C1)] <- 0
            P1 <- cbind(P1, C1)
        } else {
            C1 <- P0[,ct]
            P1 <- cbind(P1, C1)
            colnames(P1)[ncol(P1)] <- ct
        }
        celltype[jj] <- ct3        
    }

    celltype <- colnames(P1)[max.col(P1)]
    table(celltype0)
    table(celltype)
    table(celltype0, celltype)
    dim(P1)
    
    ## collapse small groups to 'other_cells'
    low.ct <- names(which(table(celltype) < low.th*length(celltype)))
    celltype[celltype %in% low.ct] <- "other_cells"
    names(celltype) <- colnames(counts)

    res <- list(celltype=celltype, probs=P1)    
    res
}

pgx.checkCellTypeMarkers <- function(counts, min.count=3, markers=NULL)
{
    ## Checks if markers are expressed (positive markers) or
    ## not-expressed (negative markers).
    ##
    ##
    CANONICAL.MARKERS.SAVE = list(
        "B cells" = c("Ms4a1","CD79a","CD79b", "Fcmr","Ebr1"),
        "CD4 T cells" = c("Il17r","Ccr7","Cd3e","Cd3d","Cd3g"),
        "CD8 T cells" = c("Cd8a","Cd8b"),
        "NK cells" = c("Gnly","Nkg7","Gzma", "Klrb1c", "Klrk1", "Klra4"),
        "Dendritic cells" = c("Fcer1a","Cst3","Siglech","Fscn1","Ccl22"),
        "Macrophages" = c("C1qa","C1qb","C1qc","Lyz2"),
        "Monocytes" = c("Ly6c2","Lyz","Cd14","Fcgr3a","Ms4a7"),    
        "B16.melanoma" = c("Mlana", "Dct", "Tyrp1", "Mt1", "Mt2", "Pmel", "Pgk1")    
    )
    CANONICAL.MARKERS2 = list(
        ##"Bcells" = c("Ms4a1+","CD79a+","CD79b+","Fcmr+","Ebf1+"),
        "B_cells" = c("Ms4a1+","Cd79a+","Cd79b+"),
        "T_cells" = c("Cd3e+","Cd3d+","Cd3g+"),
        ##"TcellsCD4" = c("Il17r","Ccr7","Cd3e","Cd3d","Cd3g"),
        ##"TcellsCD8" = c("Cd8a","Cd8b1"),
        "NK_cells" = c("Gzma+", "Klrb1c+", "Klra4+"),
        ##"DendriticCells" = c("Fscn1","Ccl22","Cst3"),
        "Dendritic_cells" = c("Cst3+","H2-Ab1+","Adgre1-"),
        "Macrophages" = c("C1qa+","C1qb+","C1qc+"),
        "Monocytes" = c("Ly6c2+","Lyz2+","Ccr2+"),
        "Immune_cell" = c("Ptprc+")
        ## "B16.melanoma" = c("Mlana", "Dct", "Tyrp1", "Mt1", "Mt2", "Pmel", "Pgk1")    
    )
    if(is.null(markers)) {
        markers <- CANONICAL.MARKERS2
    }
    
    ##markers <- markers[order(names(markers))]
    marker.genes <- sort(unique(unlist(markers)))
    marker.genes <- gsub("[-+]$","",marker.genes)
    marker.genes
    
    M <- matrix(0, nrow=length(marker.genes), ncol=length(markers))
    dimnames(M) <- list(marker.genes,names(markers))
    dim(M)
    k=7
    k=1
    for(k in 1:ncol(M)) {
        mm <- markers[[k]]
        ##m.pos <- sub("[+]$","",grep("[+]$",mm,value=TRUE))
        m.neg <- sub("[-]$","",grep("[-]$",mm,value=TRUE))
        m.pos <- setdiff(gsub("[+-]$","",mm),m.neg)
        if(length(m.pos)) M[match(m.pos,rownames(M)),k] <- +1
        if(length(m.neg)) M[match(m.neg,rownames(M)),k] <- -1
    }
    ##M <- M[,Matrix::colSums(M!=0)>0]
    M
    X=counts
    ##X <- (counts / rowMeans(counts))**2
    rownames(M) <- toupper(rownames(M))
    rownames(X) <- toupper(rownames(X))
    gg <- intersect(rownames(M),rownames(X))
    counts0 <- counts[match(gg,toupper(rownames(counts))),,drop=FALSE]
    
    X1 <- 1*(X[gg,,drop=FALSE] >= min.count)
    M1 <- M[gg,,drop=FALSE]
    check.pos <-  t(pmax(M1,0)) %*% X1
    check.neg <-  t(pmin(M1,0)) %*% X1
    check.pos <-  t(M1 == +1) %*% X1 / Matrix::colSums(M1 == 1)    
    check.neg <-  t(M1 == -1) %*% X1 / Matrix::colSums(M1 == -1)
    table(check.pos)
    table(check.neg)
    check <- 1*t(check.pos & check.neg)

    list(check=check, marker.expr=counts0, markers=markers)
}


pgx.simplifyCellTypes <- function(ct, low.th=0.01)
{
    ## Simplifies cell types names from LM22, DICE, ImmProt and
    ## ImmunoStates to standardized classification names.
    ##
    
    ## LM22
    ct[grep("^B.cells",ct)] <- "B_cells"
    ct[grep("^T.cells",ct)] <- "T_cells"    
    ct[grep("^Macrophages",ct)] <- "Macrophages"
    ct[grep("^NK",ct)] <- "NK_cells"
    ct[grep("^Dendritic",ct)] <- "Dendritic_cells"
    ct[grep("^Mast",ct)] <- "Mast_cells"
    ct[grep("^Plasma",ct)] <- "Plasma_cells"
    ct[grep("Eosinophils|Neutrophils",ct)] <- "Granulocytes"

    ## DICE
    ct[grep("CD4|CD8|TH[12]|TFH|TREG|THSTAR",ct)] <- "T_cells"
    ct[grep("^M2$",ct)] <- "Macrophages"
    ct[grep("^B_CELL",ct)] <- "B_cells"
    ct[grep("^MONOCYTES$",ct)] <- "Monocytes"

    ## ImmProt
    ct[grep("^Bmem|^Bnav|^Bplasma",ct)] <- "B_cells"
    ct[grep("^T4|^Th17|^T8|Th[12]|Tregs",ct)] <- "T_cells"
    ct[grep("^NK",ct)] <- "NK_cells"
    ct[grep("^mDC|^pDC",ct)] <- "Dendritic_cells"
    ct[grep("^M2$",ct)] <- "Macrophages"
    ct[grep("^MOim|^MOnc|^MOcl",ct)] <- "Monocytes"
    ct[grep("^BSph|^ESph|^NTph",ct)] <- "Granulocytes"

    ## ImmunoStates
    ct[grep("B_cell",ct)] <- "B_cells"
    ct[grep("T_cell",ct)] <- "T_cells"    
    ct[grep("macrophage",ct)] <- "Macrophages"
    ct[grep("natural.killer",ct)] <- "NK_cells"
    ct[grep("dendritic",ct)] <- "Dendritic_cells"
    ct[grep("plasma",ct)] <- "Plasma_cells"
    ct[grep("MAST",ct)] <- "Mast_cells"
    ct[grep("monocyte",ct)] <- "Monocytes"
    ct[grep("eosinophil|neutrophil|basophil",ct)] <- "Granulocytes"
    ct[grep("hemato",ct)] <- "other_cells"

    ## otherCells (from EPIC)
    ct[grep("otherCells",ct)] <- "other_cells"
    
    ## collapse low frequency celltype to "other"
    low.ct <- names(which(table(ct) < low.th*length(ct)))
    ct[ct %in% low.ct] <- "other_cells"
    
    ##ct <- sub("cells","_cells",ct,ignore.case=TRUE)
    ## ct <- tolower(ct)
    ct
}

pgx.purify <- function( X, ref, k=3, method=2) {

    if(0) {
        X <- 2**ngs$X
        X <- Matrix::head(X[order(-apply(ngs$X,1,sd)),],1000)
        ref = "hepatocyte"
    }

    ##----------------------------------------------------------------------
    ## Using NNLM
    ##----------------------------------------------------------------------
    
    if(method==1) {
        ## contaminating profiles (e.g. stromal, normal cells)
        normalX <- cbind(const=1, X[,ref,drop=FALSE])  ## contamination (e.g. stromal, normal cells)
        colnames(normalX) <- paste0("N",1:ncol(normalX))

        ## compute proportion of tumour content using NNMF
        res.nmf <- NNLM::nnmf( X, k=k, init = list(W0 = normalX), check.k=FALSE)
        alpha <- with(res.nmf, Matrix::colSums(W[,1:k]%*%H[1:k,]) / Matrix::colSums(W %*% H))
        str(alpha)

        ## estimate "pure" matrix
        x.total <- res.nmf$W[,] %*% res.nmf$H[,]
        x.purified <- res.nmf$W[,1:k] %*% res.nmf$H[1:k,]
        Matrix::head(X)[,1:4]
        Matrix::head(x.purified)[,1:4]

        x.contaminant <- pmax(X - x.purified,0)

    } else if(method==2) {

        ##jj <- which(!(ngs$samples$group %in% ref)) ## samples of interest
        jj <- setdiff( 1:ncol(X), ref)
        tumorX <- cbind(const=1, X[,jj,drop=FALSE])
        colnames(tumorX) <- paste0("T",1:ncol(tumorX))

        ## compute proportion of contaminant content using NNMF
        res.nmf <- NNLM::nnmf( X, k=k, init = list(W0 = tumorX), check.k=FALSE)
        beta <- with(res.nmf, Matrix::colSums(W[,1:k]%*%H[1:k,]) / Matrix::colSums(W %*% H))
        str(beta)
        alpha <- (1 - beta)

        ## estimate "pure" matrix
        x.contaminant <- res.nmf$W[,1:k] %*% res.nmf$H[1:k,]
        x.purified <- pmax(X - x.contaminant,0)
    } else {
        cat("fatal error:: unknown method\n")
        stop()
    }

    res <- list( purified = x.purified,
                contaminant = x.contaminant,
                alpha = alpha )
    return(res)
}


##counts=ngs$counts
pgx.inferCellCyclePhase <- function(counts)
{
    

    ## List of cell cycle markers, from Tirosh et al, 2015
    ##
    ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
    cc.genes = strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",split=" ")[[1]]
    s_genes <- cc.genes[1:43]
    g2m_genes <- cc.genes[44:97]

    ## Create our Seurat object and complete the initalization steps
    rownames(counts) <- toupper(rownames(counts))  ## mouse...
    ##counts1 <- cbind(counts,counts,counts,counts,counts,counts)
    obj <- Seurat::CreateSeuratObject(counts)
    obj <- Seurat::NormalizeData(obj, verbose=0)
    suppressWarnings( obj <- Seurat::CellCycleScoring(obj, s.features=s_genes,
                                              g2m.features=g2m_genes, set.ident=TRUE))
    ## view cell cycle scores and phase assignments
    ##head(x = obj@meta.data)
    ##table(obj@meta.data$Phase)
    s.score <- obj@meta.data$S.Score
    g2m.score <- obj@meta.data$G2M.Score
    phase <- obj@meta.data$Phase
    if(is.null(phase) || length(phase)==0) return(NULL)
    return(phase)
}

##counts=ngs$counts
pgx.scoreCellCycle <- function(counts)
{
    

    ## List of cell cycle markers, from Tirosh et al, 2015
    ##
    ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
    cc.genes = strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",split=" ")[[1]]
    s_genes <- cc.genes[1:43]
    g2m_genes <- cc.genes[44:97]
    length(s_genes)
    length(g2m_genes)
    
    ## Create our Seurat object and complete the initalization steps
    rownames(counts) <- toupper(rownames(counts))  ## mouse...
    ##counts1 <- cbind(counts,counts,counts,counts,counts,counts)
    obj <- Seurat::CreateSeuratObject(counts)
    obj <- Seurat::NormalizeData(obj, verbose=0)
    suppressWarnings( obj <- Seurat::CellCycleScoring(obj, s_genes, g2m_genes, set.ident=TRUE))
    ## view cell cycle scores and phase assignments
    ##head(x = obj@meta.data)
    ##table(obj@meta.data$Phase)
    s.score <- obj@meta.data$S.Score
    g2m.score <- obj@meta.data$G2M.Score
    diff.score <- s.score - g2m.score
    phase <- obj@meta.data$Phase
    df <- data.frame(phase=phase, s_score=s.score, g2m_score=g2m.score,
                     diff_score=diff.score)
    if(nrow(df)==0) return(NULL)
    return(df)
}

pgx.inferGender <- function(X, gene_name=NULL) {
    ## List of cell cycle markers, from Tirosh et al, 2015
    ##
    ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
    if(is.null(gene_name)) gene_name <- toupper(sub(".*:","",rownames(X)))
    if(0) {
        X <- log2(1+ngs$counts)
        gene_name = rownames(X)
        x.genes <- rownames(ngs$genes)[grep("X",ngs$genes$chr)]
        y.genes <- rownames(ngs$genes)[grep("Y",ngs$genes$chr)]
        x.genes
        y.genes
    } else {
        y.genes = intersect(c("DDX3Y","RPS4Y1","USP9Y","KDM5D"),gene_name)
        y.genes
        x.genes = intersect(c("XIST"),gene_name)
        x.genes
    }
    if( length(y.genes)==0 && length(x.genes)==0 ) {
        cat("warning:: could not determine sex. missing some X/Y marker genes\n")
        sex <- rep(NA, ncol(X))
        return(sex)
    }
    sex <- rep(NA, ncol(X))
    if( length(y.genes)>0 && length(x.genes)>0 ) {
        x.expr <- colMeans(X[match(x.genes, gene_name),,drop=FALSE])
        y.expr <- colMeans(X[match(y.genes, gene_name),,drop=FALSE])
        x.expr
        y.expr
        mean.expr <- colMeans(X)
        ##plot(x.expr, y.expr)
        sex <- rep(NA, ncol(X))
        sex <- ifelse( x.expr > mean.expr & y.expr < mean.expr, "F", sex)
        sex <- ifelse( y.expr > mean.expr & x.expr < mean.expr, "M", sex)
        sex        
        return(sex)
    }
    return(sex)
}

##counts=ngs$counts
pgx.multipleDeconvolution <- function(counts, refmat, methods=DECONV.METHODS)
{
    methods
    timings <- c()
    results <- list()

    if(class(refmat)!="list") {
        refmat <- list("reference" = refmat)
    }

    refnames <- names(refmat)
    refnames
    i=1
    for(i in 1:length(refmat)) {
        message("[pgx.multipleDeconvolution] computing for ",refnames[i])
        ref <- refmat[[i]]
        dim(ref)
        res = pgx.deconvolution(counts, ref=ref, methods=methods)

        if(!is.null(res)) {
            names(res)
            names(res$results)
            m <- names(refmat)[i]
            results[[m]] <- res$results
            timings <- rbind(timings, res$timings)
        }
    }
    timings
    timings0 <- NULL
    if(!is.null(timings)) {
        timings0 <- apply(timings, 2, function(x) tapply(x, rownames(timings), sum))
    }
    res2 <- list(results=results, timings=timings0)
    return(res2)
}

##X=as.matrix(2**ngs$X);methods="NNLM"
pgx.deconvolution <- function(X, ref, methods=DECONV.METHODS,
                              add.unknown=FALSE, normalize.mat=TRUE)
{

    if(max(X)<50 || min(X)<0) {
        cat("WARNING:: pgx.deconvolution: is X really counts? (not logarithmic)\n")
    }

    dbg("[pgx.deconvolution] called")

    
    ## clean up matrix, remove duplicate names
    mat <- as.matrix(X)
    rownames(mat) <- gsub(".*:","",rownames(mat)) ## strip prefix
    rownames(mat) <- toupper(rownames(mat)) ## handle mouse??
    mat <- mat[order(-rowMeans(mat)),,drop=FALSE]
    mat <- as.matrix(mat[!duplicated(rownames(mat)),,drop=FALSE])

    ref <- as.matrix(ref)
    rownames(ref) <- toupper(rownames(ref))
    ref <- ref[order(-rowMeans(ref)),,drop=FALSE]
    ref <- as.matrix(ref[!duplicated(rownames(ref)),,drop=FALSE])

    dbg("[pgx.deconvolution] 1:")
    
    ## Add "unknown" class to reference matrix
    if(add.unknown) {
        dbg("[pgx.deconvolution] adding unknown group")        
        gg <- intersect(rownames(ref),rownames(mat))
        x1 <- log(1+ref[gg,,drop=FALSE])
        y1 <- log(1+rowMeans(mat[gg,,drop=FALSE]))
        x1 <- cbind(offset=1, x1)

        ## compute residual matrix by substracting all possible linear
        ## combinations of reference.
        
        x1 <- ref[gg,,drop=FALSE]
        y1 <- rowMeans(mat[gg,,drop=FALSE])
        cf <- NNLM::nnlm(x1,cbind(y1))$coefficients
        cf[is.na(cf)] <- 0
        resid <- pmax(y1 - x1 %*% cf,0) ## residual vector           
        resx <- rep(0,nrow(ref))
        names(resx) <- rownames(ref)
        resx[gg] <- resid
        ref <- cbind(ref, "other_cells"=resx)
    }
    
    ## normalize all matrices to CPM
    if(normalize.mat) {
        ref <- t(t(ref) / sqrt(1e-6 + Matrix::colSums(ref**2,na.rm=TRUE))) * 1e6
        mat <- t(t(mat) / sqrt(1e-6 + Matrix::colSums(mat**2,na.rm=TRUE))) * 1e6
    }
    
    ## add small noise, some methods need it...
    ref <- ref + 1e-2*matrix(rnorm(length(ref)),nrow(ref),ncol(ref))  
    mat <- mat + 1e-2*matrix(rnorm(length(mat)),nrow(mat),ncol(mat)) 
    ref <- pmax(ref,0)
    mat <- pmax(mat,0)

    dbg("[pgx.deconvolution] 2:")
    
    gg <- intersect(rownames(ref),rownames(mat))
    Matrix::head(gg)
    if(length(gg) < 10) {
        message("WARNING:: pgx.deconvolution: no enough marker genes")
        return(NULL)
    }

    ## conform??
    if(0) {
        
        length(gg)
        ##ref <- ref[gg,]
        ##mat <- mat[gg,]
        qx <- limma::normalizeQuantiles( cbind(ref[gg,],mat[gg,]))
        ref <- qx[,colnames(ref)]
        mat <- mat[,colnames(mat)]
    }
    
    timings <- list()
    results <- list()
    methods
    CIBERSORT.code <- "/opt/CIBERSORT/CIBERSORTmat.R"
    if("CIBERSORT" %in% methods && file.exists(CIBERSORT.code)) {
        ## CIBERSORT
        ##source(file.path(FILES,CIBERSORT.code))
        source(CIBERSORT.code)
        cat("starting deconvolution using CIBERSORT...\n")
        ciber.out <- NULL
        stime <- system.time(
            ##try( ciber.out <- CIBERSORT(ref, mat, perm=0, QN=TRUE) )
            try( ciber.out <- CIBERSORT(ref, mat, perm=0, QN=FALSE) )
        )
        if(!is.null(ciber.out)) {
            timings[["CIBERSORT"]] <- stime
            cat("deconvolution using CIBERSORT took",stime[3],"s\n")
            ciber.out <- ciber.out[,!(colnames(ciber.out) %in% c("P-value","Correlation","RMSE"))]
            results[["CIBERSORT"]] <- ciber.out
        } else {
            cat("WARNING:: CIBERSORT failed\n")
        }
    }

    if("EPIC" %in% methods) {
        ## EPIC
        ##devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
        dbg("[pgx.deconvolution] calculating EPIC...")
        
        out <- NULL
        gg = intersect(rownames(ref),rownames(mat))
        ref1 <- ref
        if("other_cells" %in% colnames(ref1)) {
            ref1 <- ref1[,setdiff(colnames(ref1),"other_cells")]
        }
        ref.list = list(refProfiles=ref1, sigGenes=gg)
        mat1 <- mat
        colnames(mat1) <- 1:ncol(mat)  ## EPIC doesnt like duplicated column names...
        stime <- system.time( try(
            out <- EPIC(bulk = mat1, reference = ref.list)
        ))
        remove(mat1)
        if(!is.null(out)) {
            message(paste0("deconvolution using EPIC took",stime[3],"s\n"))
            timings[["EPIC"]] <- stime
            out.mat <- out[["cellFractions"]]
            colnames(out.mat) <- sub("otherCells","other_cells",colnames(out.mat))
            rownames(out.mat) <- colnames(mat)
            results[["EPIC"]] <- out.mat
        } else {
            mssage("WARNING:: EPIC fail (no pun intended...)\n")
        }
    }

    if("DeconRNAseq" %in% methods) {
        ## DeconRNAseq
        ##BiocManager::install("DeconRNASeq", version = "3.8")
        ##BiocManager::install("DeconRNASeq")
        if("package:Seurat" %in% search()) detach("package:Seurat", unload=TRUE)
        dbg("[pgx.deconvolution] calculating DeconRNAseq...")
        
        ## uses psych::pca() from pcaMethods
        require(pcaMethods)  ## uses pcaMethods::prep
        drs <- NULL
        stime <- system.time(suppressMessages(suppressWarnings(
            drs <- try(DeconRNASeq::DeconRNASeq(data.frame(mat, check.names=FALSE),
                                   data.frame(ref, check.names=FALSE))$out.all)
        )))
        timings[["DeconRNAseq"]] <- stime
        class(drs)
        if(!is.null(drs) && class(drs)!="try-error") {
            cat("deconvolution using DeconRNAseq took",stime[3],"s\n")
            rownames(drs) <- colnames(mat)
            ##colnames(drs) <- colnames(ref)
            results[["DeconRNAseq"]] <- drs
        } else {
            cat("WARNING:: DeconRNAseq failed\n")
        }
    }

    ##--------- DCQ from ComICS ---------------------
    if("DCQ" %in% methods) {
        ## DCQ seems to work in logX, so we use log-transform
        ##install.packages("ComICS")
        dbg("[pgx.deconvolution] calculating DCQ...")
        
        res.dcq=NULL
        stime <- system.time(
            res.dcq <- try(
                ComICS::dcq(reference_data = log2(1+as.matrix(ref)),
                    mix_data = log2(1+as.matrix(mat)),  ## log data OK??
                    ##marker_set = matrix(rownames(ref),ncol=1),
                    marker_set = cbind(intersect(rownames(ref),rownames(mat))),
                    alpha_used=0.05, lambda_min=0.2, number_of_repeats=3,
                    precent_of_data=1.0)
            ))
        if(!is.null(res.dcq) && class(res.dcq)!="try-error") {
            timings[["DCQ"]] <- stime
            cat("deconvolution using DCQ took",stime[3],"s\n")
            results[["DCQ"]] <- res.dcq$average
        } else {
            cat("WARNING:: DCQ failed\n")
        }
    }

    if("I-NNLS" %in% methods) {
        ## !!!!!!!!!!!!!! WARNING:: needs more testing/validation !!!!!!!!!!!!!
        ##----- Constrained iterative non-negative least squares (Abbas et al. 2009) ----
        dbg("[pgx.deconvolution] calculating I-NNLS...")
        
        GetFractions.Abbas <- function(XX, y, w=NA){
            ## XX is immune expression data
            ## y is cancer expression data
            ss.remove=c()
            ss.names=colnames(XX)
            while(TRUE){
                if(length(ss.remove)==0) {
                    tmp.XX=XX
                } else{
                    if(is.null(ncol(tmp.XX))) return(rep(0, ncol(XX)))
                    tmp.XX=tmp.XX[,-ss.remove,drop=FALSE]
                }
                if(length(ss.remove)>0){
                    ss.names=ss.names[-ss.remove]
                    if(length(ss.names)==0) return(rep(0, ncol(XX)))
                }
                if(is.na(w[1])) {
                    tmp=lsfit(tmp.XX, y, intercept=FALSE)
                } else {
                    tmp=lsfit(tmp.XX, y, w, intercept=FALSE)
                }
                if(ncol(tmp.XX)==1) {
                    tmp.beta=tmp$coefficients[1]
                } else {
                    tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
                }
                if(min(tmp.beta>0)) break ## break if coefs are all positive
                ss.remove=which.min(tmp.beta)  ## removes most negative coeff
            }
            tmp.F=rep(0, ncol(XX))
            names(tmp.F)=colnames(XX)
            tmp.F[ss.names]=tmp.beta
            return(tmp.F)
        }
        gg <- intersect(rownames(ref),rownames(mat))
        XX = ref[gg,,drop=FALSE]
        YY = mat[gg,,drop=FALSE]
        res.abbas=NULL
        stime <- system.time(
            res.abbas <- try(apply(YY,2,function(y) GetFractions.Abbas(XX, y, w=NA)))
        )
        timings[["I-NNLS"]] <- stime
        cat("deconvolution using I-NNLS took",stime[3],"s\n")
        if(!is.null(res.abbas) && class(res.abbas)!="try-error") {
            rownames(res.abbas) <- colnames(ref)
            results[["I-NNLS"]] <- t(res.abbas)
        }
    }

    ## Own NNLM (non-negative linear modeling)...
    if("NNLM" %in% methods) {
        ## NNLM
        ##install.packages("NNLM")
        dbg("[pgx.deconvolution] calculating NNLM...")
        
        x1 <- log2(1+ref[gg,,drop=FALSE])
        x2 <- log2(1+mat[gg,,drop=FALSE])
        x1 <- cbind(offset=1, x1)
        stime <- system.time(
            cf <- NNLM::nnlm(x1,x2)$coefficients[-1,,drop=FALSE]
        )
        timings[["NNLM"]] <- stime
        cat("deconvolution using NNLM took",stime[3],"s\n")
        results[["NNLM"]] <- t(cf)

        ## very much the same as I-NNLS
        results[["NNLM.lin"]] <- t(NNLM::nnlm(ref[gg,,drop=FALSE], mat[gg,,drop=FALSE])$coefficients)

        r1 <- apply(ref[gg,,drop=FALSE],2,rank,na.last="keep")
        r2 <- apply(mat[gg,,drop=FALSE],2,rank,na.last="keep")
        r1 <- cbind(offset=1, r1)
        cf <- NNLM::nnlm(r1,r2)$coefficients[-1,,drop=FALSE]
        results[["NNLM.rnk"]] <- t(cf)

    }

    ## Simple (rank) correlation
    if("cor" %in% methods) {
        dbg("[pgx.deconvolution] calculating cor...")        
        ##results[["cor"]] <- stats::cor(log2(1+mat[gg,]), log2(1+ref[gg,]))
        r1 <- apply(mat[gg,,drop=FALSE],2,rank,na.last="keep")
        r2 <- apply(ref[gg,,drop=FALSE],2,rank,na.last="keep")
        stime <- system.time(
            cf <- stats::cor(r1,r2,use="pairwise")
        )
        timings[["cor"]] <- stime
        cat("deconvolution using COR took",stime[3],"s\n")
        results[["cor"]] <- cf
    }

    if("SingleR" %in% methods) {
        dbg("[pgx.deconvolution] calculating SingleR...")                
        stime <- system.time(
            sr1 <- SingleR(test=mat, ref=ref, labels=colnames(ref))
        )
        timings[["SingleR"]] <- stime
        cat("deconvolution using SingleR took",stime[3],"s\n")
        results[["SingleR"]] <- sr1$scores
    }
    ## clean up
    names(results)
    results <- results[which(!sapply(results,is.null))]
    results <- lapply(results, function(x) {x[is.na(x)]=0;x})

    dbg("[pgx.deconvolution] calculating meta values...")                
    
    ## meta
    if(length(results)>1) {
        jj <- colnames(ref)
        norm.results <- lapply( results, function(x) x[,jj,drop=FALSE] / (1e-8 + rowSums(x[,jj,drop=FALSE]))  )
        lognorm.results <- lapply( norm.results, function(x) log(0.001+pmax(x,0)))
        res.meta1 = Reduce( '+', norm.results) / length(norm.results)
        res.meta2 = exp(Reduce( '+', lognorm.results) / length(lognorm.results))
        results[["meta"]] <- res.meta1 / (1e-8 + rowSums(res.meta1,na.rm=TRUE))
        results[["meta.prod"]] <- res.meta2 / (1e-8+rowSums(res.meta2,na.rm=TRUE))
    }

    timings0 <- do.call(rbind, timings)
    res2 <- list(results=results, timings=timings0)

    dbg("[pgx.deconvolution] done!")
    
    return(res2)
}
