if(0) {

    BiocManager::install("TCGAbiolinks")

    ## See: https://cran.r-project.org/web/packages/TCGAretriever
    BiocManager::install("TCGAretriever")
    library(TCGAretriever)
    ## Define a set of genes of interest
    q_genes <- c("TP53", "MDM2", "E2F1", "EZH2")
    q_cases <- "brca_tcga_pub_complete"
    rna_prf <- "brca_tcga_pub_mrna"
    mut_prf <- "brca_tcga_pub_mutations"
    brca_RNA <- TCGAretriever::get_profile_data(case_id = q_cases, gprofile_id = rna_prf, glist = q_genes)
    head(brca_RNA[, 1:5])

    BiocManager::install("curatedTCGAData")
    BiocManager::install("TCGAutils")
    library(curatedTCGAData)
    library(MultiAssayExperiment)
    library(TCGAutils)
    (mae <- curatedTCGAData("DLBC", c("RNASeq2GeneNorm"), FALSE))
    ##(mae <- curatedTCGAData("DLBC", c("RNASeqGene", "Mutation"), FALSE))
    counts <- assays(mae)[[1]]
    dim(counts)

}

pgx.testTCGAsurvival <- function(sig, matrix_file, lib.dir, ntop=100,
                                 sortby.p = FALSE, plot=TRUE, verbose=1 )
{                                          
    require(survival)
    require(rhdf5)

    ##matrix_file = file.path(ARCHS4.DIR,"tcga_matrix.h5")
    if(!file.exists(matrix_file)) {
        stop("cannot find TCGA H5 matrix file")
    }
    if(is.null(names(sig))) {
        stop("sig must have names")
    }

    ## get the top DE genes
    genes <- c(head(names(sig),ntop),tail(names(sig),ntop))
    head(genes)

    ## Read the H5 matrix file
    ## aa <- h5ls(matrix_file)[,1:2]
    ## aa
    ## ii <- which(aa[,1]=="/meta")[-1]
    ## aa.head <- lapply(ii,function(i) head(h5read(matrix_file, paste0("/meta/",aa[i,2]))))
    ## names(aa.head) <- aa[ii,2]
    ## aa.head

    if(verbose) cat("[pgx.TCGA.testSurvivalSignature] extracting expression from H5 matrix file\n")
    
    h5.samples = h5read(matrix_file, "/meta/gdc_cases.submitter_id")
    h5.genes = h5read(matrix_file, "/meta/genes")            
    h5.project = h5read(matrix_file, "/meta/gdc_cases.project.project_id")            
    
    sample_index <- 1:length(h5.samples)
    gene_index <- 1:length(h5.genes)            
    
    ##sample_index <- which(h5.samples %in% samples)
    gene_index <- which(h5.genes %in% genes)
    head(gene_index)
    
    expression = h5read(
        matrix_file, "data/expression",
        index = list(gene_index, sample_index)
    )
    colnames(expression) <- h5.samples[sample_index]
    rownames(expression) <- h5.genes[gene_index]
    dim(expression)

    ## Read the survival data
    if(verbose) cat("[pgx.TCGA.testSurvivalSignature] reading TCGA survival data...\n")
    
    surv.file <- file.path(lib.dir, "rtcga-survival.csv")
    surv <- read.csv(surv.file, row.names=1)
    head(surv)
    surv$months <- round(surv$times/365*12, 2)
    surv$status <- surv$patient.vital_status

    ## conform expression and surv matrices
    samples <- intersect(colnames(expression), rownames(surv))
    length(samples)    
    expression <- expression[,samples]
    surv <- surv[samples,]

    ## print KM survival plot for each study/cancertype
    all.studies <- sort(unique(surv$cancer_type))
    length(all.studies)
    study <- all.studies[1]

    if(verbose) cat("[pgx.TCGA.testSurvivalSignature] fitting survival probabilities...\n")
    
    surv.p <- rep(NA,length(all.studies))
    rho.list <- list()
    names(surv.p) <- all.studies
    for(study in head(all.studies,99)) {
        
        study

        ## calculate correlation with signature        
        sel <- which(surv$cancer_type == study)
        if(length(sel)<20) next()
        gg <- rownames(expression)
        rho <- cor( expression[,sel], sig[gg], use="pairwise")[,1]
        sel.data <- surv[sel,]
    
        ## fit survival curve on two groups
        poscor <- (rho > median(rho,na.rm=TRUE))
        table(poscor)

        sdf <- survdiff( Surv(months, status) ~ poscor, data = sel.data)
        p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

        surv.p[study] <- p.val
        rho.list[[study]] <- rho
    }

    if(plot) {
        if(verbose) cat("[pgx.TCGA.testSurvivalSignature] plotting KM curves...\n")    
        jj <- 1:length(rho.list)
        if(sortby.p) {
            ii <- order(surv.p)
            surv.p <- surv.p[ii]
            rho.list <- rho.list[ii]
        }
        surv.q <- p.adjust(surv.p)
        names(surv.q) <- names(surv.p)
        
        par(mfrow=c(5,7), mar=c(2,3,2,1))
        for(study in names(surv.p)) {
            
            study
            
            ## calculate correlation with signature        
            sel <- which(surv$cancer_type == study)
            if(length(sel)<20) next()
            rho <- rho.list[[study]]
            sel.data <- surv[sel,]
            
            ## fit survival curve on two groups
            poscor <- (rho > median(rho,na.rm=TRUE))
            table(poscor)
            library(survival)
            fit <- survfit( Surv(months, status) ~ poscor, data = sel.data )
            
            ##legend.labs <- paste(c("negative","positive"),"correlated")
            legend.labs <- paste(c("rho<0","rho>0"))
            if(1) {
                plot(fit, col=2:3, lwd=2, main=study, xlab="time   (days)", cex.main=1.1)
                legend("bottomleft", legend.labs, pch="__", lwd=2, col=2:3, cex=1.0)
                
                p.val <- round(surv.p[study], 3)
                q.val <- round(surv.q[study], 3)
                pq <- c(paste("p=",p.val), paste("q=",q.val))
                legend("bottomright", pq, bty='n', cex=1.0)
            } else {
                library(survminer)
                ggsurvplot(
                    fit, 
                    data = sel.data, 
                    size = 1,                 # change line size
                    ## palette = c("#E7B800", "#2E9FDF"),# custom color palettes
                    conf.int = TRUE,          # Add confidence interval
                    pval = TRUE,              # Add p-value
                    risk.table = TRUE,        # Add risk table
                    risk.table.col = "strata",# Risk table color by groups
                    ## legend.labs = c("Male", "Female"),    # Change legend labels
                    legend.labs = legend.labs,
                    risk.table.height = 0.20, # Useful to change when you have multiple groups
                    ggtheme = theme_bw()      # Change ggplot2 theme
                )
            } ## end of if
        } ## end of for
    }
    
    return(surv.p)
}        

cancertype="dlbc";variables="OS_"
cancertype="brca_tcga_pub"

pgx.selectTCGAstudies <- function(cancertype, variables)
{
    ## Scan the available TCGA studies for cancertype and clinical
    ## variables.
    ##
    ##
    ##
    library(cgdsr)
    mycgds <- CGDS("http://www.cbioportal.org/")
    all.studies <- sort(getCancerStudies(mycgds)[,1])
    studies <- grep(cancertype, all.studies, value=TRUE)    
    clin <- list()
    samples <- list()
    studies
    mystudy <-  studies[1]

    for(mystudy in studies) {

        mystudy
        myprofiles <- getGeneticProfiles(mycgds,mystudy)[,1]
        myprofiles

        ## mrna datatypes
        mrna.type <- "mrna"
        if(any(grepl("rna_seq_mrna$", myprofiles))) mrna.type <- "rna_seq_mrna"
        if(any(grepl("v2_mrna$", myprofiles))) mrna.type <- "rna_seq_v2_mrna"
        pr.mrna <- grep( paste0(mrna.type,"$"), myprofiles,value=TRUE)
        pr.mrna
        if(length(pr.mrna)==0) next()
        
        all.cases <- getCaseLists(mycgds,mystudy)[,1]
        all.cases
        ##if(!any(grepl("complete$",all.cases))) next        
        ##caselist <- grep("complete$",all.cases,value=TRUE)
        caselist <- grep(paste0(mrna.type,"$"),all.cases,value=TRUE)
        caselist
        clin0 <- getClinicalData(mycgds, caselist)
        head(clin0)[,1:4]
        rownames(clin0) <- gsub("[.]","-",rownames(clin0)) ## correct names...
        head(clin0)[,1:4]
        clin[[mystudy]] <- clin0
        samples[[mystudy]] <- rownames(clin0)
    }
    
    sel <- sapply(clin, function(v) any(grepl(variables,colnames(v))))
    sel
    sel.studies <- studies[sel]
    sel.clin    <- clin[sel]
    
    res <- list(
        studies = sel.studies,
        ##samples = sel.samples,
        clinicalData = sel.clin
    )
    return(res)
}

genes=NULL;study="brca_tcga_pub"
pgx.getTCGAdataset <- function(study, genes=NULL, matrix_file=NULL, from.h5=TRUE)
{
    ## For a specific TCGA study get the expression matrix and
    ## clinical data.
    ##

    ## check if H5 exists
    from.h5 <- (from.h5 && !is.null(matrix_file) && file.exists(matrix_file))
    
    ##BiocManager::install("cgdsr")    
    library(cgdsr)
    mycgds <- CGDS("http://www.cbioportal.org/")

    all.studies <- sort(getCancerStudies(mycgds)[,1])
    if(!all(study %in% all.studies)) {
        ss <- setdiff(study, all.studies)
        stop(ss,"is not in TCGA studies")
    }
    
    ## Gather data from all study
    X <- list()
    clin <- list()
    mystudy <-  study[1]
    for(mystudy in study) {

        cat("getting TCGA expression for",mystudy,"...\n")
        
        mystudy
        ##myprofiles = "ov_tcga_rna_seq_v2_mrna"        
        myprofiles <- getGeneticProfiles(mycgds,mystudy)[,1]
        myprofiles

        ## mrna datatypes
        mrna.type <- "_mrna"
        if(any(grepl("rna_seq_v2_mrna$", myprofiles))) mrna.type <- "rna_seq_v2_mrna"
        if(any(grepl("rna_seq_mrna$", myprofiles))) mrna.type <- "rna_seq_mrna"
        pr.mrna <- grep( paste0(mrna.type,"$"), myprofiles,value=TRUE)
        pr.mrna
        if(length(pr.mrna)==0) {
            cat("WARNING:: could not find mRNA for",mystudy,"...\n")
            next()
        }
        
        all.cases <- getCaseLists(mycgds,mystudy)[,1]
        all.cases
        ##if(!any(grepl("complete$",all.cases))) next        
        ##caselist <- grep("complete$",all.cases,value=TRUE)
        caselist <- grep(paste0(mrna.type,"$"),all.cases,value=TRUE)
        caselist
        samples <- NULL
        head(genes)
        if(!is.null(genes) && !from.h5) {
            cat("downloading...\n")
            ## If only a few genes, getProfileData is a faster way
            ##
            expression <- t(getProfileData(mycgds, genes, pr.mrna, caselist))
            samples <- gsub("[.]","-",colnames(expression))
            colnames(expression) <- samples            
            dim(expression)
        } else {
            cat("extracting from locally stored H5 matrix...\n")
            ## For all genes, getProfileData cannot do and we use
            ## locally stored H5 TCGA data file from Archs4.
            ##
            xx <- getProfileData(mycgds, "---", pr.mrna, caselist)
            samples <- gsub("[.]","-",colnames(xx))[3:ncol(xx)]
            head(samples)

            library("rhdf5")
            library("preprocessCore")
            h5closeAll()
            ##matrix_file = file.path(ARCHS4.DIR, "tcga_matrix.h5")
            has.h5 <- file.exists(matrix_file)
            has.h5
            
            if(!has.h5) {
                stop("FATAL: could not find tcga_matrix.h5 matrix. Please download from Archs4.")
            } else {
                
                ## Retrieve information from locally stored H5 compressed data            
                aa <- h5ls(matrix_file)[,1:2]
                aa
                ii <- which(aa[,1]=="/meta")[-1]
                lapply(ii,function(i) head(h5read(matrix_file, paste0("/meta/",aa[i,2]))))
                ##h5read(matrix_file, "/meta/gdc_cases.project.project_id")
                id1 = h5read(matrix_file, "/meta/gdc_cases.samples.portions.submitter_id")
                id2 = h5read(matrix_file, "/meta/gdc_cases.samples.submitter_id")
                id3 = h5read(matrix_file, "/meta/gdc_cases.submitter_id")    
                id2x <- substring(id2,1,15)
                
                h5.genes = h5read(matrix_file, "/meta/genes")            
                if(!is.null(genes)) h5.genes <- intersect(genes,h5.genes)
                samples = intersect(samples, id2x)
                sample_index <- which(id2x %in% samples)
                gene_index <- 1:length(h5.genes)            

                if(length(sample_index)==0 || length(gene_index)==0) {
                    return(list(X=NULL, clin=NULL))
                }
                
                expression = h5read(
                    matrix_file, "data/expression",
                    index = list(gene_index, sample_index)
                )
                H5close()
                dim(expression)
                colnames(expression) <- substring(id2[sample_index],1,15)
                rownames(expression) <- h5.genes
                expression <- expression[,order(-colSums(expression))]
                expression <- expression[,samples]
            }
            
        }
        dim(expression)
        this.clin <- getClinicalData(mycgds, caselist)
        rownames(this.clin) <- gsub("[.]","-",rownames(this.clin))
        this.clin <- this.clin[samples,,drop=FALSE]
        expression <- expression[,samples,drop=FALSE]
        X[[mystudy]] <- expression
        clin[[mystudy]] <- this.clin
    }

    res <- list(X=X, clin=clin)
    return(res)
}


pgx.getTCGA.multiomics.TOBEFINISHED <- function(studies, genes=NULL, batch.correct=TRUE,
                                                tcga.only=TRUE )
{
    ## Better use curatedTCGA bioconductor package!!!!
    ##
    
    ##BiocManager::install("cgdsr")    
    library(cgdsr)
    mycgds <- CGDS("http://www.cbioportal.org/")
    if(0) {
        all.studies <- sort(getCancerStudies(mycgds)[,1])
        tcga.studies <- grep("_tcga$",all.studies, value=TRUE)
        all.studies <- tcga.studies
        all.studies
        mystudy <- "ov_tcga"
        mystudy <- "brca_tcga"
        mystudy <- "thca_tcga"
    }
    ## all.profiles <- list()
    ## for(mystudy in all.studies) {    
    ##     myprofiles <- getGeneticProfiles(mycgds,mystudy)[,1]
    ##     all.profiles[[mystudy]] <- myprofiles
    ## }

    GENE = "CIITA"
    GENE = "NLRC5"

    ## Gather data from all cancers
    all.X <- list()
    mystudy <-  studies[1]
    for(mystudy in studies) {
        mystudy
        ##myprofiles = "ov_tcga_rna_seq_v2_mrna"        
        myprofiles <- getGeneticProfiles(mycgds,mystudy)[,1]
        myprofiles

        ## prioritize datatypes
        pr.mrna <- grep("rna_seq_v2_mrna$|rna_seq_mrna$",myprofiles,value=TRUE)[1]
        ## pr.prot <- paste0(mystudy,"_protein_quantification")
        pr.cna <- grep("_log2CNA$|_linear_CNA$",myprofiles,value=TRUE)[1]
        pr.gistic <- grep("_gistic$",myprofiles,value=TRUE)[1]
        pr.me  <- grep("_methylation_hm450|_methylation_hm27",myprofiles,value=TRUE)[1]
        pr.mut <- grep("_mutations",myprofiles,value=TRUE)[1]    

        all.cases <- getCaseLists(mycgds,mystudy)[,1]
        all.cases
        if(!any(grepl("complete$",all.cases))) next
        
        caselist <- grep("complete$",all.cases,value=TRUE) 
        cna=counts=cna.gistic=me=mut=gx=NULL    
        counts <- getProfileData(mycgds, genes, pr.mrna, caselist)
        cna    <- getProfileData(mycgds, GENE, pr.cna, caselist)
        ##prot <- getProfileData(mycgds, GENE, pr.prot, caselist)
        cna.gistic <- getProfileData(mycgds, GENE, pr.gistic, caselist)
        me   <- getProfileData(mycgds, GENE, pr.me, caselist)
        mut  <- getProfileData(mycgds, GENE, pr.mut, caselist)
        mut  <- 1*!is.na(mut)
        gx   <- log2(10 + as.matrix(counts))
        cna[is.na(cna)] <- NA
        if(grepl("linear",pr.cna)) cna <- log2(0.01 + 2 + cna)  ## assume diploid
    
        ##colnames(counts) <- paste0("GX:",colnames(counts))
        if(!is.null(cna)) colnames(cna) <- paste0("CN:",colnames(cna))
        if(!is.null(cna.gistic)) colnames(cna.gistic) <- paste0("CNA:",colnames(cna.gistic))
        if(!is.null(me))  colnames(me) <- paste0("ME:",colnames(me))
        if(!is.null(mut)) colnames(mut) <- paste0("MT:",colnames(mut))
        ##colnames(prot) <- paste0("PX:",colnames(prot))
        
        ##xx <- list(gx, cna, cna.gistic, me, mut)
        xx <- list(gx, cna.gistic, me, mut)
        xx <- xx[sapply(xx,nrow)>0]    
        X <- do.call(cbind, xx)
        dim(X)
        
        if(!is.null(X) && ncol(X)>=4 ) {
            X <- X[,colMeans(is.na(X)) < 0.5,drop=FALSE]
            X <- X[rowMeans(is.na(X)) < 0.5,,drop=FALSE]
            dim(X)
            all.X[[mystudy]] <- X
        }        
    }
}
