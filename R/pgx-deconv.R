DECONV.METHODS = c("I-NNLS","CIBERSORT","DCQ","DeconRNAseq","EPIC","NNLM","cor")

if(0) {
    ## WARNING: this snippet will overwrite your ngs object...
    ##
    ##

    load(file="../pgx/geiger2018-arginine-8k-LT.pgx",verbose=1)
    load(file="../pgx/sallusto2019-th1star-UC-12x-LT.pgx",verbose=1)
    imm.sig <- t(apply(ngs$count, 1, function(x) tapply(x, ngs$samples$group, mean)))
    imm.sig <- head(imm.sig[order(-apply(ngs$X,1,sd)),],1000)
    ##imm.sig <- head(imm.sig[order(-apply(imm.sig,1,sd)),],500)
    dim(imm.sig)
    colnames(imm.sig) <- sub("_S$",".resting",colnames(imm.sig))
    colnames(imm.sig) <- sub("_A$",".activated",colnames(imm.sig))
    write.csv(imm.sig, file="../files/immprot-signature1000.csv")
    remove(ngs)

    ##devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
    load(file="../pgx/schmiedel2018-DICE-mRNA-4k.pgx",verbose=1)
    dice.sig <- t(apply(ngs$count, 1, function(x) tapply(x, ngs$samples$group, mean)))
    dice.sig <- head(dice.sig[order(-apply(ngs$X,1,sd)),],1000)
    dice.sig <- head(dice.sig[order(-apply(dice.sig,1,sd)),],500)
    dim(dice.sig)
    write.csv(dice.sig, file="../files/DICE-signature1000.csv")
    remove(ngs)
}

pgx.testSignature <- function(ngs, sig) {



}

pgx.purify <- function( X, ref, k=3, method=2) {

    if(0) {
        X <- 2**ngs$X
        X <- head(X[order(-apply(ngs$X,1,sd)),],1000)
        ref = "hepatocyte"
    }

    ##----------------------------------------------------------------------
    ## Using NNLM
    ##----------------------------------------------------------------------
    require(NNLM)
    if(method==1) {
        ## contaminating profiles (e.g. stromal, normal cells)
        normalX <- cbind(const=1, X[,ref,drop=FALSE])  ## contamination (e.g. stromal, normal cells)
        colnames(normalX) <- paste0("N",1:ncol(normalX))

        ## compute proportion of tumour content using NNMF
        res.nmf <- nnmf( X, k=k, init = list(W0 = normalX), check.k=FALSE)
        alpha <- with(res.nmf, colSums(W[,1:k]%*%H[1:k,]) / colSums(W %*% H))
        str(alpha)

        ## estimate "pure" matrix
        x.total <- res.nmf$W[,] %*% res.nmf$H[,]
        x.purified <- res.nmf$W[,1:k] %*% res.nmf$H[1:k,]
        head(X)[,1:4]
        head(x.purified)[,1:4]

        x.contaminant <- pmax(X - x.purified,0)

    } else if(method==2) {

        ##jj <- which(!(ngs$samples$group %in% ref)) ## samples of interest
        jj <- setdiff( 1:ncol(X), ref)
        tumorX <- cbind(const=1, X[,jj,drop=FALSE])
        colnames(tumorX) <- paste0("T",1:ncol(tumorX))

        ## compute proportion of contaminant content using NNMF
        res.nmf <- nnmf( X, k=k, init = list(W0 = tumorX), check.k=FALSE)
        beta <- with(res.nmf, colSums(W[,1:k]%*%H[1:k,]) / colSums(W %*% H))
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
    require(Seurat)

    ## List of cell cycle markers, from Tirosh et al, 2015
    ##
    ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
    cc.genes = strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA",split=" ")[[1]]
    s_genes <- cc.genes[1:43]
    g2m_genes <- cc.genes[44:97]

    ## Create our Seurat object and complete the initalization steps
    rownames(counts) <- toupper(rownames(counts))  ## mouse...
    ##counts1 <- cbind(counts,counts,counts,counts,counts,counts)
    obj <- CreateSeuratObject(counts)
    obj <- NormalizeData(obj)
    obj <- CellCycleScoring(
        obj, s_genes, g2m_genes, set.ident = TRUE)

    ## view cell cycle scores and phase assignments
    ##head(x = obj@meta.data)
    ##table(obj@meta.data$Phase)
    phase <- obj@meta.data$Phase
    return(phase)
}

pgx.inferGender <- function(X, gene_name=NULL) {
    ## List of cell cycle markers, from Tirosh et al, 2015
    ##
    ##cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
    if(is.null(gene_name)) gene_name <- toupper(sub(".*:","",rownames(X)))
    y.genes = intersect(c("DDX3Y","RPS4Y1","USP9Y","KDM5D"),gene_name)
    y.genes
    x.genes = intersect(c("XIST"),gene_name)
    x.genes
    if( length(y.genes)==0 && length(x.genes)==0 ) {
        cat("warning:: could not determine sex. missing some X/Y marker genes\n")
        sex <- rep(NA, ncol(X))
        return(sex)
    }
    sex <- rep(NA, ncol(X))
    if( length(y.genes)>0 && length(x.genes)>0 ) {
        y.expr <- colMeans(X[match(y.genes, gene_name),,drop=FALSE])
        x.expr <- colMeans(X[match(x.genes, gene_name),,drop=FALSE])
        y.expr
        x.expr
        ##plot(x.expr, y.expr)
        sex <- rep(NA, ncol(X))
        sex <- ifelse( x.expr > y.expr, "F", "M")
        return(sex)
    }
    if( length(y.genes)>0 && length(x.genes)==0 ) {
        y.expr <- colMeans(X[match(y.genes, gene_name),,drop=FALSE])
        y.expr
        sex <- rep(NA, length(y.expr))
        sex <- ifelse(y.expr > log2(100), "M", sex)
        sex <- ifelse(y.expr < log2(10),  "F", sex)
        return(sex)
    }
    table(sex)
    return(sex)
}


##counts=ngs$counts
pgx.multiDeconvolution <- function(counts, refmat, methods=DECONV.METHODS)
{
    methods
    timings <- c()
    results <- list()
    i=1
    for(i in 1:length(refmat)) {
        ref <- refmat[[i]]
        dim(ref)
        res = pgx.deconvolution(counts[,], ref=ref, methods=methods)
        names(res)
        names(res$results)
        results[[i]] <- res$results
        names(results)[i] <- names(refmat)[i]
        timings <- rbind(timings, res$timings)
    }
    timings
    timings0 <- apply(timings, 2, function(x) tapply(x, rownames(timings), sum))
    res2 <- list(results=results, timings=timings0)
    return(res2)
}

##X=as.matrix(2**ngs$X);ref="LM22";methods=DECONV.METHODS
pgx.deconvolution <- function(X, ref, methods=DECONV.METHODS)
{

    if(max(X)<100 || min(X)<0) {
        cat("WARNING:: pgx.deconvolution: is X really counts? (not logarithmic)\n")
    }

    mat <- as.matrix(X)
    rownames(mat) <- toupper(gsub(".*:","",rownames(mat))) ## handle mouse??
    mat <- mat[order(-rowMeans(mat)),]
    mat <- as.matrix(mat[!duplicated(rownames(mat)),])
    head(mat)[,1:4]
    ref <- ref[order(-rowMeans(ref)),]
    ref <- as.matrix(ref[!duplicated(rownames(ref)),])
    head(ref)[,1:4]

    gg <- intersect(rownames(ref),rownames(mat))
    head(gg)
    if(length(gg)==0) {
        stop("ERROR:: pgx.deconvolution: no overlapping rownames")
    }

    ## conform??
    if(0) {
        require(limma)
        length(gg)
        ##ref <- ref[gg,]
        ##mat <- mat[gg,]
        qx <- normalizeQuantiles( cbind(ref[gg,],mat[gg,]))
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
        library(EPIC)
        out <- NULL
        gg = intersect(rownames(ref),rownames(mat))
        ref.list = list(refProfiles=ref, sigGenes=gg)
        mat1 <- mat
        colnames(mat1) <- 1:ncol(mat)  ## EPIC doesnt like duplicated column names...
        stime <- system.time( try(
            out <- EPIC(bulk = mat1[,], reference = ref.list)
        ))
        remove(mat1)
        if(!is.null(out)) {
            cat("deconvolution using EPIC took",stime[3],"s\n")
            timings[["EPIC"]] <- stime
            out.mat <- out[["cellFractions"]]
            rownames(out.mat) <- colnames(mat)
            results[["EPIC"]] <- out.mat
        } else {
            cat("WARNING:: EPIC fail (no pun intended...)\n")
        }
    }

    if("DeconRNAseq" %in% methods) {
        ## DeconRNAseq
        ##BiocManager::install("DeconRNASeq", version = "3.8")
        ##BiocManager::install("DeconRNASeq")
        if("package:Seurat" %in% search()) detach("package:Seurat", unload=TRUE)
        library(DeconRNASeq)
        ## uses pca() from pcaMethods
        drs <- NULL
        stime <- system.time(
            drs <- try(DeconRNASeq(data.frame(mat, check.names=FALSE),
                                   data.frame(ref, check.names=FALSE))$out.all)
        )
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
        require(ComICS)
        res.dcq=NULL
        stime <- system.time(
            res.dcq <- try(
                dcq(reference_data = log2(1+as.matrix(ref)),
                    mix_data = log2(1 + as.matrix(mat)),  ## log data OK??
                    ##marker_set = matrix(rownames(ref),ncol=1),
                    marker_set = cbind(rownames(ref)),
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
        ##----- Constrained iterative non-negative least squares (Abbas et al. 2009) ----
        GetFractions.Abbas <- function(XX, YY, w=NA){
            ## XX is immune expression data
            ## YY is cancer expression data
            ss.remove=c()
            ss.names=colnames(XX)
            while(T){
                if(length(ss.remove)==0) {
                    tmp.XX=XX
                } else{
                    if(is.null(ncol(tmp.XX))) return(rep(0, ncol(XX)))
                    tmp.XX=tmp.XX[, -ss.remove]
                }
                if(length(ss.remove)>0){
                    ss.names=ss.names[-ss.remove]
                    if(length(ss.names)==0) return(rep(0, ncol(XX)))
                }
                if(is.na(w[1])) {
                    tmp=lsfit(tmp.XX, YY, intercept=F)
                } else {
                    tmp=lsfit(tmp.XX, YY, w, intercept=F)
                }
                if(is.null(ncol(tmp.XX))) {
                    tmp.beta=tmp$coefficients[1]
                } else {
                    tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
                }
                if(min(tmp.beta>0)) break
                ss.remove=which.min(tmp.beta)
            }
            tmp.F=rep(0, ncol(XX))
            names(tmp.F)=colnames(XX)
            tmp.F[ss.names]=tmp.beta
            return(tmp.F)
        }
        gg <- intersect(rownames(ref),rownames(mat))
        XX=ref[gg,];YY=mat[gg,]
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
        require(NNLM)
        x1 <- log2(1+ref[gg,])
        x2 <- log2(1+mat[gg,])
        x1 <- cbind(offset=1, x1)
        stime <- system.time(
            cf <- nnlm(x1,x2)$coefficients[-1,]
        )
        timings[["NNLM"]] <- stime
        cat("deconvolution using NNLM took",stime[3],"s\n")
        results[["NNLM"]] <- t(cf)

        ## very much the same as I-NNLS
        results[["NNLM.lin"]] <- t(nnlm(ref[gg,], mat[gg,])$coefficients)

        r1 <- apply(ref[gg,],2,rank,na.last="keep")
        r2 <- apply(mat[gg,],2,rank,na.last="keep")
        r1 <- cbind(offset=1, r1)
        cf <- nnlm(r1,r2)$coefficients[-1,]
        results[["NNLM.rnk"]] <- t(cf)

    }

    ## Simple correlation
    if("cor" %in% methods) {
        ##results[["cor"]] <- cor(log2(1+mat[gg,]), log2(1+ref[gg,]))
        r1 <- apply(mat[gg,],2,rank,na.last="keep")
        r2 <- apply(ref[gg,],2,rank,na.last="keep")
        stime <- system.time(
            cf <- cor(r1,r2,use="pairwise")
        )
        timings[["cor"]] <- stime
        cat("deconvolution using COR took",stime[3],"s\n")
        results[["cor"]] <- cf
    }

    ## clean up
    names(results)
    results <- results[which(!sapply(results,is.null))]
    results <- lapply(results, function(x) {x[is.na(x)]=0;x})

    ## meta
    if(length(results)>1) {
        jj <- colnames(ref)
        norm.results <- lapply( results, function(x) x[,jj] / (1e-8 + rowSums(x[,jj]))  )
        lognorm.results <- lapply( norm.results, function(x) log(0.001+pmax(x,0)))
        res.meta1 = Reduce( '+', norm.results) / length(norm.results)
        res.meta2 = exp(Reduce( '+', lognorm.results) / length(lognorm.results))
        results[["meta"]] <- res.meta1 / (1e-8 + rowSums(res.meta1,na.rm=TRUE))
        results[["meta.prod"]] <- res.meta2 / (1e-8+rowSums(res.meta2,na.rm=TRUE))
    }

    timings0 <- do.call(rbind, timings)
    res2 <- list(results=results, timings=timings0)
    return(res2)
}
