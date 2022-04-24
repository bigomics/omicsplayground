##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

SignatureBoard <- function(id, inputData, selected_gxmethods)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns ## NAMESPACE

    fullH = 800   ## full height of page
    tabH = '70vh'
    
infotext =
    "In the <strong>Signature Analysis module</strong>, users can test their gene signature by calculating an enrichment score. They can use a sample list provided on the platform or upload their own gene list. Instead of a short list, a profile can also be selected, which is a complete gene list resulted from one of the contrasts in the analysis.

<br><br>After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene.

<br><br>The <strong>Enrichment tab</strong> performs the enrichment analysis of the gene list against all contrasts by running the GSEA algorithm and plots enrichment outputs. The enrichment statistics can be found in the corresponding table

<br><br>Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, as well as the p and q values by the Fisher’s test for the overlap test.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=7' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
"

    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    IMMCHECK.GENES = "ADORA2A ARHGEF5 BTLA CD160 CD244 CD27 CD274 CD276 CD47 CD80 CEACAM1 CTLA4 GEM HAVCR2 ICOS IDO1 LAG3 PDCD1 TNFSF4 VISTA VTCN1 TIGIT PVR CD28 CD40 CD40LG ICOSLG TNFRSF9 TNFSF9 CD70 TNFRSF4 TNFRSF18 TNFSF18 SIRPA LGALS9 ARG1 CD86 IDO2 PDCD1LG2 KIR2DL3"
    APOPTOSIS.GENES = "BAD CRADD AGT FAS BCL2 PPIF S100A9 S100A8 BBC3 BCL2L11 FADD CTSH MLLT11 TRAF7 BCL2L1 HTRA2 BNIP3 BAK1 PMAIP1 LGALS9 BID"
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Signature Analysis Board</strong>"),
            shiny::HTML(infotext),
            easyClose = TRUE, size="l"))
    })

    ##------------------------ observe/reactive function  -----------------------------

    shiny::observeEvent(input$example1, { 
        shiny::updateTextAreaInput(session,"genelistUP", value=IMMCHECK.GENES)
    })
    shiny::observeEvent(input$example2, { 
        shiny::updateTextAreaInput(session,"genelistUP", value=APOPTOSIS.GENES)
    })
    shiny::observeEvent(input$example3, { 
        shiny::updateTextAreaInput(session,"genelistUP", value=CELLCYCLE.GENES)
    })

    shiny::observe({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        type="contrast"
        type <- input$type
        if(is.null(type)) type <- "<custom>"

        if(type=="contrast") {
            contr <- sort(names(ngs$gx.meta$meta))
            shiny::updateSelectInput(session, "feature", choices=contr, selected=contr[1])
        } else if(type=="hallmark") {
            ## collection
            gsets <- sort(grep("HALLMARK",names(iGSETS),value=TRUE))
            shiny::updateSelectInput(session, "feature", choices=gsets, selected=gsets[1])
        } else if(type=="KEGG") {
            ## collection
            gsets <- sort(grep("KEGG",names(iGSETS),value=TRUE))
            shiny::updateSelectInput(session, "feature", choices=gsets, selected=gsets[1])
        } else if(type=="geneset") {
            ## all genesets... this is a bit too much for selectInput (DO NOT USE!!)
            gsets <- sort(names(iGSETS))
            shiny::updateSelectizeInput(session, "feature", choices=gsets, selected=gsets[1], server=TRUE)
        } else {
            ## custom
            shiny::updateSelectInput(session, "feature", choices="<custom>", selected="<custom>")
        }
    })


    ##================================================================================
    ##======================= REACTIVE FUNCTIONS =====================================
    ##================================================================================

    input_genelistUP <- shiny::reactive({
        gg <- input$genelistUP
        if(is.null(gg)) return(NULL)
        gg <- strsplit(as.character(gg), split="[, \n\t]")[[1]]
        if(length(gg)==1 && gg[1]!="") gg <- c(gg,gg)  ## hack to allow single gene....
        return(gg)
    }) %>% shiny::debounce(1000)

    
    getCurrentMarkers <- shiny::reactive({
        ##
        ## Get current selection of markers/genes 
        ##
        ##
        
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        type="<custom>"
        type="contrast"
        type <- input$type
        ##if(is.null(type)) return(NULL)
        ##if(is.null(input$contrast)) return(NULL)
        ##if(is.null(input$feature)) return(NULL)
        shiny::req(input$type, input$feature)
        
        dbg("<signature:getCurrentMarkers> called\n")
        
        level = "gene"
        features = toupper(ngs$genes$gene_name)
        xfeatures = toupper(ngs$genes[rownames(ngs$X),"gene_name"])
        gset <- NULL
        if(input$feature=="<custom>") {
            gset <- input_genelistUP()
            if(is.null(gset) || length(gset)==0 || gset[1]=="") return(NULL)
            ##gset <- toupper(gset)        
            if(length(gset)==1) {
                gene <- sub("^[@#]","",gset[1])
                if(grepl("^@",gset[1]) && gene %in% xfeatures) {
                    ## most correlated with this genes
                    jj <- match(gene, xfeatures)  ## single gene
                    rho <- cor(t(ngs$X), ngs$X[jj,])[,1]
                    gset <- head(names(sort(abs(rho),decreasing=TRUE)),36)  ## how many?
                } else {
                    ## grep-like match
                    rx <- toupper(gset[1])
                    rx <- grep(rx, xfeatures, value=TRUE, ignore.case=TRUE)
                    gset <- rownames(ngs$X)[which(xfeatures %in% rx)]  ## all probes matching gene
                }
            }
        } else if(type=="contrast" &&
                  input$feature %in% names(ngs$gx.meta$meta) ) {
            contr=1
            contr <- input$feature
            fx <- ngs$gx.meta$meta[[contr]]$meta.fx
            probes <- rownames(ngs$gx.meta$meta[[contr]])
            genes <- toupper(ngs$genes[probes,"gene_name"])
            top.genes <- genes[order(-fx)]
            top.genes <- head(top.genes,100)
            top.genes0 <- paste(top.genes,collapse=" ")
            shiny::updateTextAreaInput(session,"genelistUP", value=top.genes0)
            gset <- top.genes
        } else if(input$feature %in% names(iGSETS)) {
            ##gset <- toupper(GSETS[[input$feature]])
            gset <- toupper(unlist(getGSETS(input$feature)))
            gset0 <- paste(gset, collapse=" ")
            shiny::updateTextAreaInput(session,"genelistUP", value=gset0)
        } else {
            return(NULL)
        }
        
        return(gset)
    })
    
    getSingleSampleEnrichment <- shiny::reactive({
        ##
        ## Calls calcSingleSampleValues() and calculates single-sample
        ## enrichment values for complete data matrix and reduced data by
        ## group (for currentmarkers)
        ##
        ##
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        ## select samples
        X = ngs$X
        sel = colnames(X)
        X <- X[,sel]
        
        ## get the signature
        gset <- strsplit(IMMCHECK.GENES,split=" ")[[1]]
        gset <- getCurrentMarkers()
        if(is.null(gset)) return(NULL)

        ##y = 1*(rownames(X) %in% gset)
        ##rownames(X)=toupper(rownames(X)); gset=toupper(gset)
        xgene <- ngs$genes[rownames(X),"gene_name"]
        y = 1*(toupper(xgene) %in% toupper(gset))
        names(y) <- rownames(X)
        table(y)
        
        ## expression by group
        ##grp = ngs$samples[colnames(X),"group"]
        grp <- ngs$model.parameters$group
        groups = unique(grp)
        gX <- sapply( groups, function(g) rowMeans(X[,which(grp==g),drop=FALSE]))
        colnames(gX) = groups
        dim(gX)
        dim(ngs$X)
        
        ## for large datasets pre-grouping is faster
        ss.bygroup  <- calcSingleSampleValues(gX, y, method=c("rho","gsva"))
        do.rho   = TRUE
        dbg("getSingleSampleEnrichment:: 2 : do.rho")
        ss1 <- calcSingleSampleValues(X[,], y, method=c("rho"))
        ##ss.bysample <- cbind(ss.bysample, rho=ss1)
        ss.bysample <- cbind(rho=ss1)        
        
        res <- list( by.sample=ss.bysample, by.group=ss.bygroup)
        return(res)
    })


    sigCalculateGSEA <- shiny::reactive({
        ## 
        ## Calculate fgsea for current marker selection and active
        ## datasets.
        ##
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        
        ## observe input list
        gset = head(rownames(ngs$X),100)
        gset <- getCurrentMarkers()
        if(is.null(gset)) return(NULL)
        ##if(is.null(input$enplotsdb)) return(NULL)
        
        ## get all logFC of this dataset
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        F <- meta$fc
        rownames(F) <- toupper(rownames(F))
                        
        ## cleanup matrix
        F = as.matrix(F)
        dim(F)
        F =  F[,which(!duplicated(colnames(F))),drop=FALSE]
        
        ## cleanup names and uppercase for mouse genes
        rownames(F) <- toupper(sub(".*:","",rownames(F)))
        gset <- toupper(sub(".*:","",gset))
        gset <- intersect(toupper(gset), rownames(F))
        length(gset)

        if(length(gset)==0) {
            cat("FATAL:: sigCalculateGSEA : gset empty!\n")
            return(NULL)
        }
        
        ## ------------ prioritize with quick correlation

        y = 1*(toupper(rownames(F)) %in% toupper(gset))
        ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)[,1]
        rho = cor(apply(F,2,ss.rank), y, use="pairwise")[,1]
        ##wt = c(mean(y==0),mean(y==1))[1+y]
        ##wt.rho = apply(F,2, function(x) weightedCorr((x), y, weights=wt, method="Pearson"))
        rho[is.na(rho)] <- 0
        names(rho) = colnames(F)
        
        ## ------- restrict to top 100 comparisons (fgsea is otherwise to
        ## ------- slow) but we do not expect many with so many
        ## ------- comparisons
        ntop = 100    
        jj <- head(order(-abs(rho)), ntop )
        F <- F[,jj,drop=FALSE]
        F <- F[!duplicated(rownames(F)),,drop=FALSE]
        F <- F + 1e-4*matrix(rnorm(length(F)),nrow(F),ncol(F))
        dim(F)
        
        ## ------------- do fast GSEA
        gmt = list("gset"=unique(gset))
        res <- NULL
        enrich_method="rcor"
        enrich_method="fgsea"
        ##enrich_method <- input$rankmethod
        
        if(enrich_method=="fgsea") {
            i=1
            dbg("sigCalculateGSEA:: starting fgsea...\n")            
            shiny::withProgress(message="computing GSEA ...", value=0.8, {
                res <- lapply(1:ncol(F), function(i) {
                    suppressWarnings( suppressMessages(
                        res <- fgsea::fgsea(gmt, stats=F[,i], nperm=1000)
                    ))
                    res <- as.data.frame(res[,c("pval","padj","ES","NES")])
                    rownames(res)[1] = colnames(F)[i]
                    return(res)
                })
            })
            dbg("sigCalculateGSEA:: fgsea done!\n")
            res1 <- data.frame(do.call(rbind, res))
            res1$ES <- NULL
        } else {
            i=1
            fx <- 1*(rownames(F) %in% gmt[[1]])
            rho <- cor(apply(F,2,rank,na.last="keep"), fx, use="pairwise")[,1]
            pv <- cor.pvalue(rho, nrow(F))
            qv <- p.adjust(pv,method="fdr")
            res1 <- data.frame(pval=pv, padj=qv, rho=rho, NES=NA)
            rownames(res1) <- names(pv)
        }

        ## columns are: NES, pval, fdr, contrast
        res1 <- as.matrix(res1)
        res1 <- res1[match(colnames(F),rownames(res1)),,drop=FALSE]
        
        if( nrow(res1) != ncol(F)) {
            cat("WARNING sigCalculateGSEA:: fgsea results are corrupted?\n")
            cat("WARNING sigCalculateGSEA:: got contrasts: ",res$contrast,"\n")
            cat("WARNING sigCalculateGSEA:: colnames.F= ",colnames(F),"\n")
        }
        
        ## make nice table
        ##nes   <- unlist(sapply(res, function(x) x$NES))
        ##pval  <- unlist(sapply(res, function(x) x$pval))
        nes <- res1[,"NES"]
        pval <- res1[,"pval"]
        qval <- p.adjust( pval, method="fdr")
        rho <- rho[colnames(F)]
        
        output <- as.matrix(cbind(NES=nes, p=pval, q=qval, rho=rho))
        rownames(output) <- colnames(F)    
        output <- output[order(-abs(output[,"NES"])),,drop=FALSE]
        F <- F[,rownames(output),drop=FALSE]    
        gsea <- list(F=as.matrix(F), gset=gset, output=output)
        dbg("sigCalculateGSEA:: done!\n")
        return(gsea)
    })

    ##X=gX;method=c("rho","gsva") 
    calcSingleSampleValues <- function(X, y, method=c("rho","gsva") ) {
        ##
        ## Calculates single-sample enrichment values for given matrix and
        ## binarized signature vector.
        ##
        ##
        ## very fast rank difference

        dbg("<signature:calcSingleSampleValues> called\n")
        
        if(is.null(names(y)) && length(y)!=nrow(X) ) {
            cat("<signature:calcSingleSampleValues> FATAL ERROR: y must be named if not matched\n")
            return(NULL)
        }
        
        if(!is.null(names(y)) && length(y)!=nrow(X) ) {
            y <- y[match(rownames(X),names(y))]
        }
        names(y) <- rownames(X)
        jj <- which(!is.na(y))
        X <- X[jj,]
        y <- y[jj]

        dbg("<signature:calcSingleSampleValues> 1\n")
        
        if(sum(y!=0)==0) {
            cat("<signature:calcSingleSampleValues> WARNING: y is all zero!\n")        
            matzero <- matrix(0, nrow=ncol(X), ncol=length(method))
            colnames(matzero) <- method
            rownames(matzero) <- colnames(X)
            return(matzero)
        }
        ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)[,1]
        
        S = list()
        if("rho" %in% method) {
            S[["rho"]] <- cor(apply(X, 2, ss.rank), y, use="pairwise")[,1]
            ##S$rho <- scale(S$rho)[,1]  ## should we scale??
        }

        dbg("<signature:calcSingleSampleValues> 2\n")
        
        ## calculate GSVA
        if("gsva" %in% method) {


            gset = names(y)[which(y!=0)]
            gmt <- list("gmt"=gset)
            res.gsva <- GSVA::gsva( X, gmt, method="gsva", parallel.sz=1) ## parallel=buggy
            res.colnames = colnames(res.gsva)
            fc = as.vector(res.gsva[1,])
            names(fc) = res.colnames
            S[["gsva"]] = fc[colnames(X)]
        }    
        s.names = names(S)
        if(length(S)>1) {
            S1 = do.call(cbind, S)
        } else {
            S1 <- S[[1]]
        }

        dbg("<signature:calcSingleSampleValues> done!\n")

        S1 = as.matrix(S1)
        rownames(S1) = colnames(X)
        colnames(S1) = s.names
        ##S1 <- S1[order(-S1[,"gsva"]),]
        return(S1)
    }


    ##================================================================================
    ## Enrichment {data-height=800}
    ##================================================================================
    
    enplots.RENDER <- shiny::reactive({
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        if(is.null(ngs)) return(NULL)
        

        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)
        

        ## filter with table selection/search
        ii  <- enrichmentContrastTable$rows_all()
        shiny::req(ii)
        ct <- rownames(gsea$output)[ii]
        F <- as.matrix(gsea$F[,ct,drop=FALSE])
        qv <- gsea$output[ct,"q"]
        gset <- gsea$gset
        

        cex.main=1.1
        nc=3
        par(mfrow=c(4,3), mar=c(0.3,3,3,0.5), mgp=c(1.9,0.7,0), oma=c(0,1,0,0) )
        if(ncol(F)>12) {
            par(mfrow=c(5,4), mar=c(0.2,2,3,0.6))
            cex.main=0.9
            nc=4
        }
        ## if(ncol(F)>24) par(mfrow=c(7,5), mar=c(1,2,2.5,0.6))
        for(i in 1:min(20,ncol(F))) {
            f <- colnames(F)[i]
            tt <- sub(".*\\]","",f)
            tt <- breakstring(substring(tt,1,50),28,force=TRUE)
            ylab <- ""
            if(i%%nc==1) ylab <- "rank metric"
            gsea.enplot(F[,i], gset, main=tt, cex.main=cex.main,
                        xlab="", ylab=ylab)
            qv1 <- paste("q=",round(qv[i],digits=3))
            legend("topright",qv1, cex=0.9, bty="n", adj=0)
            if(grepl("^\\[",f)) {
                db <- sub("\\].*","]",colnames(F)[i])
                legend("topleft",db, cex=0.9, bty="n", adj=0)
            }
        }
    })

    
    enplots_info = "<b>Enrichment plots.</b> Enrichment of the query signature in all constrasts. Positive enrichment means that this particular contrast shows similar expression changes as the query signature."

    enplots.opts = NULL
    shiny::callModule(
        plotModule,
        id = "enplots", 
        func = enplots.RENDER,
        func2 = enplots.RENDER,
        plotlib="base",
        info.text = enplots_info,
        options = enplots.opts,
        pdf.width=10, pdf.height=8,
        height = c(fullH-80,750),
        width = c('100%',1000),
        res=c(90,90),
        add.watermark = WATERMARK
    )

    ##================================================================================
    ## Volcano {data-height=800}
    ##================================================================================
    
    volcanoPlots.RENDER <- shiny::reactive({
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        if(is.null(ngs)) return(NULL)
        

        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)
        
        ## filter with table selection/search
        ii  <- enrichmentContrastTable$rows_all()
        shiny::req(ii)
        ct = colnames(ngs$model.parameters$contr.matrix)
        ct <- rownames(gsea$output)[ii]

        mm <- 'meta'
        mm <- selected_gxmethods()
        meta <- pgx.getMetaMatrix(ngs, methods=mm)
        F  <- meta$fc[,ct]
        qv <- meta$qv[,ct]
        score <- abs(F) * -log(qv)
        gset=head(rownames(F),100)
        gset <- intersect(gsea$gset,rownames(F))
        
        sel <- enrichmentGeneTable$rows_selected()
        sel.gene <- NULL
        if(length(sel)) {
            df <- getEnrichmentGeneTable()
            sel.gene <- df$gene[sel]
        }
                

        cex.main=1.2
        par(mfrow=c(2,2), mar=c(2,4,3,1), mgp=c(2.2,0.8,0) )
        if(ncol(F)>4) {        
            par(mfrow=c(3,3), mar=c(1,4,3,1), mgp=c(2.2,0.8,0) )
            cex.main=1            
        }
        if(ncol(F)>9) {
            par(mfrow=c(4,4), mar=c(0.2,2,3,0.6))
            cex.main=0.9
        }
        ## if(ncol(F)>24) par(mfrow=c(7,5), mar=c(1,2,2.5,0.6))
        i=1
        for(i in 1:min(16,length(ct))) {
            gset2 = head(gset[order(-score[gset,i])],30)
            cex2 = 0.8
            if(!is.null(sel.gene)) {
                gset2 <- sel.gene
                cex2 = 1.3
            }
            ##pgx.Volcano(ngs, ct[i], hilight=gset,
            ##            hilight2=gset2, cex=0.85,
            ##            cpal=c("grey80","grey80"), title='')
            xy <- cbind(fc=F[,i], z=-log10(qv[,i]))
            pgx.scatterPlotXY.BASE(
                xy, var=NULL, type="factor", title='',
                xlab = "differential expression (log2FC)",
                ylab = "significance (-log10q)",
                hilight = gset, hilight2 = gset2,
                cex = 0.9, cex.lab = cex2, cex.title = 1.0,
                legend = FALSE, col=c("grey80","grey80"),
                opacity = 1)
            title(ct[i], cex.main=cex.main, line=0.3)
        }

    })

    volcanoPlots_caption = "<b>Volcano plots.</b> Visualization of the query signature on the volcano plots of all constrasts. For positive enrichment, genes of the query signature would fall on the upper right of the volcano plot, for negative enrichment, on the upper left."
    volcanoPlots_info = volcanoPlots_caption

    volcanoPlots.opts = NULL
    shiny::callModule(
        plotModule,
        id = "volcanoPlots", 
        func = volcanoPlots.RENDER,
        func2 = volcanoPlots.RENDER,
        plotlib="base",
        info.text = volcanoPlots_info,
        options = volcanoPlots.opts,
        pdf.width=10, pdf.height=8,
        height = c(fullH-80,780),
        width = c('100%',1100),
        res = c(90,100),
        add.watermark = WATERMARK
    )
    
    ##================================================================================
    ## Overlap/similarity
    ##================================================================================

    getOverlapTable <- shiny::reactive({
        ##
        ##
        ##
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        markers <- head(rownames(ngs$X),100)
        markers <- getCurrentMarkers()
        if(is.null(markers)) return(NULL)
        
        ## fold change just for ranking of genes
        ##F <- sapply(ngs$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
        F <- sapply(ngs$gx.meta$meta, function(x) x$meta.fx)
        rownames(F) <- rownames(ngs$gx.meta$meta[[1]])
        fx <- rowMeans(F**2)
        
        ## fisher test
        ##ii <- setdiff(match(markers, colnames(GSETxGENE)),NA)
        ii <- setdiff(match(toupper(markers), colnames(GSETxGENE)),NA)
        N <- cbind(k1=Matrix::rowSums(GSETxGENE!=0), n1=ncol(GSETxGENE),
                   k2=Matrix::rowSums(GSETxGENE[,ii]!=0), n2=length(ii) )
        rownames(N) = rownames(GSETxGENE)
        ##N <- N[which(!(N[,1]==0 & N[,3]==0)), ]
        N <- N[which(N[,1]>0 | N[,3]>0), ]
        odds.ratio = ( N[,3]/ N[,4]) / ( N[,1]/ N[,2]) 
        dim(N)
        
        ## WOW THIS IS FAST!!!!!!!
        pv <- corpora::fisher.pval( N[,1], N[,2], N[,3], N[,4], log.p=FALSE)
        head(pv)
        names(pv) <- rownames(N)
        pv = pv[match(names(odds.ratio),names(pv))]
        ##qv = p.adjust(pv, method="fdr")
        qv = p.adjust(pv, method="bonferroni")
        A = data.frame( odds.ratio=odds.ratio, p.fisher=pv, q.fisher=qv)
        dim(A)

        ## limit the list??
        table( qv < 0.05)
        table( qv < 0.2)
        table( qv < 0.999)
        A <- A[which( A$q.fisher < 0.999),]
        ##A <- A[which( A$q.fisher < 0.05),]
        dim(A)

        ## get shared genes
        dbg("[getOverlapTable] determining shared genes...\n")
        aa = rownames(A)

        y <- 1*(colnames(GSETxGENE) %in% toupper(markers))
        names(y) <- colnames(GSETxGENE)
        ncommon <- Matrix::colSums(Matrix::t(GSETxGENE[aa,,drop=FALSE])*as.vector(y)!=0)
        ntotal  <- Matrix::rowSums(GSETxGENE[aa,,drop=FALSE]!=0)
        A$ratio <- ncommon / ntotal
        ratio.kk <- paste0(ncommon,"/",ntotal)    

        gg <- colnames(GSETxGENE)
        gset <- names(y)[which(y!=0)]
        G1 = GSETxGENE[aa,which(y!=0)]
        commongenes <- apply(G1, 1, function(x) colnames(G1)[which(x!=0)])
        ##commongenes <- lapply(commongenes, function(x) x[order(-fx[x])])
        ##commongenes <- parallel::mclapply(commongenes, function(x) x[order(-fx[x])])
        for(i in 1:length(commongenes)) {
            gg <- commongenes[[i]]
            gg <- gg[order(-abs(fx[gg]))]
            if(length(gg)>10) {
                others <- paste0("(+",length(gg)-10," others)")
                gg <- c(head(gg,10),others)
            }
            commongenes[[i]] <- paste(gg,collapse=",")        
        }
        ##commongenes <- sapply(commongenes,paste,collapse=",")        
        commongenes <- unlist(commongenes)
        
        ## construct results dataframe
        gset.names <- substring(rownames(A),1,72)    
        ##aa <- apply(A, 2, formatC, format="e", digits=3)
        A$ratio <- round(A$ratio, digits=3)
        A$log.OR <- round(log10(A$odds.ratio), digits=3)
        A$odds.ratio <- round(A$odds.ratio, digits=3)
        db = sub(":.*","",gset.names)
        score = (log10(A$odds.ratio) * -log10(A$q.fisher + 1e-40))**0.5
        score = round(score, digits=3)
        df <- cbind(db=db, geneset=gset.names, score=score, "k/K"=ratio.kk, A, common.genes=commongenes)
        
        if(DEV) {
            df <- df[,c("db","geneset","score","k/K","ratio","odds.ratio","log.OR","q.fisher","common.genes")]
        } else {
            df <- df[,c("db","geneset","score","k/K","odds.ratio","q.fisher","common.genes")]
        }
        
        ##df <- df[order(-df$odds.ratio),]
        df <- df[order(-df$score),]
        dbg("[getOverlapTable] done! \n")
        return(df)
    })

    overlapScorePlot.RENDER <- shiny::reactive({



        df <- getOverlapTable()
        sel <- 1:nrow(df)
        sel<- overlapTable$rows_all()
        shiny::req(df,sel)

        df1 <- df[sel,]
        df1$geneset = as.character(rownames(df1))
        df1$db = factor(df1$db)
        
        ntop = 1000
        ntop = as.integer(input$overlapScorePlot_ntop)
        df1 = df1[head(order(-df1$score),ntop),]
        jj = order(df1$db, -df1$score)
        df1 = df1[jj,]

        df1$idx = factor(1:nrow(df1), levels=1:nrow(df1))
        df1$idx <- as.integer(df1$idx)
        klr = rep(RColorBrewer::brewer.pal(8,"Set2"),10)[as.integer(df1$db)]
        
        plt <- plotly::plot_ly(
            df1, x = ~idx, y = ~score,
            type='bar',  ## orientation='v',
            ## text = ~geneset,
            hoverinfo = 'text',
            hovertemplate = paste0("%{text}<br>score: %{y}<extra>",df1$db,"</extra>"),
            ##hovertemplate = "%{y}",
            marker = list( color=klr ) ) %>%
            plotly::layout(
                showlegend = FALSE,
                dragmode= 'select',
                ##annotations = anntitle(colnames(rho)[i]),
                ##annotations = list(text="TITLE"),
                yaxis = list(##range = c(0,1),
                    titlefont = list(size=11),
                    tickfont = list(size=10),
                    showgrid = TRUE,
                    title = "overlap score" ),
                xaxis = list(
                    title = "",
                    showgrid = FALSE,
                    showline = FALSE,
                    showticklabels = FALSE,
                    showgrid = FALSE,
                    zeroline = FALSE)) 

        if( min(nrow(df1),ntop) < 100 && input$overlapScorePlot_shownames) {
            ## labeling the y-axis inside bars
            plt <- plt %>%
                plotly::add_annotations( yref='paper', xref = 'x',
                                x = ~idx, y=0.005, yanchor='bottom',
                                text = substring(df1$geneset,1,35),
                                textangle = -90,
                                font = list(size = 10),
                                showarrow = FALSE, align='right')
        }

        plt    
    })

    overlapTable.RENDER <- shiny::reactive({

        df <- getOverlapTable()
        shiny::req(df)    

        df$geneset <- wrapHyperLink(df$geneset, df$geneset)

        numeric.cols <- which(sapply(df, is.numeric))
        numeric.cols <- intersect(c("p.fisher","q.fisher"),colnames(df))
        
        DT::datatable(df, class='compact cell-border stripe',
                      rownames=FALSE, escape = c(-1,-2),
                      extensions = c('Scroller'),
                      selection='none',
                      fillContainer=TRUE,
                      options=list(
                          dom = 'frtip',
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, scrollY = tabH, scroller=TRUE ## deferRender=TRUE,
                      )  ## end of options.list 
                      ) %>%
            DT::formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("score",
                                background = color_from_middle( df$score, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    overlapScorePlot.opts = shiny::tagList(
        withTooltip(shiny::radioButtons(ns("overlapScorePlot_ntop"),
                            "Number of features",c(60,120,250),inline=TRUE),
               "Specify the number to top features to show.",
               placement="top", options = list(container = "body")),
        withTooltip(shiny::checkboxInput(ns("overlapScorePlot_shownames"),
                             "Show feature names",TRUE),
               "Select to show/hide the feature names in the plot.",
               placement="top", options = list(container = "body"))
    )

    shiny::callModule(
        plotModule,
        id = "overlapScorePlot", 
        func = overlapScorePlot.RENDER,
        plotlib = "plotly",
        title = "Signature overlap scores", label="a",
        info.text = "Top overlapping gene sets with selected signature. The vertical axis shows the overlap score of the gene set which combines the odds ratio and significance (q-value) of the Fisher's test.",
        options = overlapScorePlot.opts,
        pdf.width = 12, pdf.height = 6,
        height = 0.45*fullH, res=100,
        add.watermark = WATERMARK
    )
    
    overlapTable <- shiny::callModule(
        tableModule,
        id = "overlapTable",
        func = overlapTable.RENDER,
        title = "Overlap with other signatures", label="b",
        info.text = "Under the <strong>Overlap/similarity tab</strong>, users can find the similarity of their gene list with all the gene sets and pathways in the platform, including statistics such as the total number of genes in the gene set (K), the number of intersecting genes between the list and the gene set (k), the overlapping ratio of k/K, logarithm of the  odds ratio (log.OR), as well as the p and q values by the Fisher’s test for the overlap test.",
        height = 0.4*fullH
    )


    ##================================================================================
    ## Markers {data-height=800}
    ##================================================================================

    markers.RENDER <- shiny::reactive({
        ##if(!input$tsne.all) return(NULL)


        
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)

        dbg("<signature:markers.RENDER> called\n")        
        
        markers <- ngs$families[[2]]
        markers <- COLLECTIONS[[10]]
        markers <- getCurrentMarkers()
        if(is.null(markers)) return(NULL)
        
        level = "gene"
        ##markers <- intersect(markers,ngs$genes$gene_name)
        ##jj <- match(markers,ngs$genes$gene_name)
        xgene <- ngs$genes[rownames(ngs$X),]$gene_name
        jj <- match(toupper(markers), toupper(xgene))
        jj <- setdiff(jj,NA)
        gx <- ngs$X[jj,,drop=FALSE]
        
        if(nrow(gx)==0) {
            cat("WARNING:: Markers:: markers do not match!!\n")
            return(NULL)
        }
        
        ## get t-SNE positions of samples
        pos = ngs$tsne2d[colnames(gx),]
        gx = gx - min(gx,na.rm=TRUE) + 0.001 ## subtract background
        dim(gx)
        ##grp <- ngs$samples[colnames(gx),"group"]
        grp <- ngs$model.parameters$group
        zx <- t(apply(gx,1,function(x) tapply(x,as.character(grp),mean)))
        gx <- gx[order(-apply(zx,1,sd)),,drop=FALSE]
        rownames(gx) = sub(".*:","",rownames(gx))
        
        ## ---------------- get GSVA values
        res <- getSingleSampleEnrichment()
        if(is.null(res)) return(NULL)
        
        S <- res$by.sample
        if(NCOL(S)==1) {
            fc = S[,1]
        } else {
            fc = colMeans( t(S) / (1e-8+sqrt(colSums(S**2))) ) ## scaled mean
        }
        fc <- scale(fc)[,1]  ## scale??
        names(fc) = rownames(S)
        ##fc1 = sign(fc) * (fc/(1e-8+max(abs(fc))))**2
        fc1 = tanh(1.0*fc / (1e-4+sd(fc)))
        fc1 = fc1[rownames(pos)]
        
        cex1 = 1.2
        cex1 <- 0.7*c(1.6,1.2,0.8,0.5)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]    
        cex2 <- ifelse(level=="gene",1,0.8)
        klrpal = colorRampPalette(c("grey90", "grey60", "red3"))(16)

        nmax=NULL
        if(input$markers_layout=="6x6") {
            nmax = 35
            par(mfrow=c(6,6), mar=c(0,0.2,0.5,0.2), oma=c(2,1,2,1)*0.8 )
        }
        if(input$markers_layout=="4x4") {
            nmax = 15
            par(mfrow=c(4,4), mar=c(0,0.2,0.5,0.2), oma=c(2,1,2,1)*0.8 )
        }

        top.gx = head(gx,nmax)
        if(input$markers_sortby=="name") {
            top.gx = top.gx[order(rownames(top.gx)),,drop=FALSE]
        }
        if(input$markers_sortby=="probability") {
            top.gx = top.gx[order(-rowMeans(top.gx)),,drop=FALSE]
        }
        if(input$markers_sortby=="correlation") {
            rho <- cor(t(top.gx), fc1)[,1]
            top.gx = top.gx[order(-rho),,drop=FALSE]
        }
        
        i=1    
        for(i in 0:min(nmax,nrow(top.gx))) {
            jj <- 1:ncol(top.gx)
            if(i==0) {
                klr1 = BLUERED(16)[8 + round(7*fc1)]
                tt = "INPUT SIGNATURE"
                jj <- order(abs(fc1))
            } else {
                colvar = pmax(top.gx[i,],0) 
                colvar = 1+round(15*(colvar/(0.7*max(colvar)+0.3*max(top.gx))))
                klr1 = klrpal[colvar]
                gene <- substring(sub(".*:","",rownames(top.gx)[i]),1,80)
                tt <- breakstring(gene, n=20, force=TRUE)
                jj <- order(abs(top.gx[i,]))
            }
            klr1 = paste0(gplots::col2hex(klr1),"99")
            
            base::plot( pos[jj,], pch=19, cex=cex1, col=klr1[jj],
                 xlim=1.2*range(pos[,1]), ylim=1.2*range(pos[,2]),
                 fg = gray(ifelse(i==0,0.1,0.8)), bty = "o",
                 xaxt='n', yaxt='n', xlab="tSNE1", ylab="tSNE2")
            legend("topleft", tt, cex=cex2, col="grey30", text.font=ifelse(i==0,2,1),
                   inset=c(-0.1,-0.05), bty="n")

        }

        dbg("<signature:markers.RENDER> done!\n")        
        
    })


    markers_info = "After uploading a gene list, the <strong>Markers</strong> section produces a t-SNE plot of samples for each gene, where the samples are colored with respect to the upregulation (in red) or downregulation (in blue) of that particular gene."

    markers_caption = "<b>Markers t-SNE plot</b>. T-SNE plot for each gene, where the dot (corresponding to samples) are colored depending on the upregulation (in red) or downregulation (in blue) of that particular gene."
    
    markers.opts = shiny::tagList(
        withTooltip(shiny::radioButtons(ns("markers_sortby"),"Sort by:",
                            choices=c("correlation","probability","name"), inline=TRUE),
               "Sort by correlation, probability or name.", placement="top",
               options = list(container = "body")),
        withTooltip(shiny::radioButtons(ns("markers_layout"),"Layout:", choices=c("4x4","6x6"),
                            ## selected="6x6",
                            inline=TRUE),
               "Choose layout.", 
               placement="top", options = list(container = "body")),
    )

    shiny::callModule(
        plotModule,
        id = "markers",
        title = "Markers plot", 
        func = markers.RENDER,
        func2 = markers.RENDER,
        plotlib = "base",
        info.text = markers_info,
        options = markers.opts,
        pdf.width=8, pdf.height=8,
        height = c(fullH-100,750), res=c(100,95),
        add.watermark = WATERMARK
    )
    
    ##================================================================================
    ## Enrichment {data-height=800}
    ##================================================================================
    
    enrichmentContrastTable.RENDER <- shiny::reactive({
        
        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)

        dbg("enrichmentContrastTable.RENDER: reacted")
        
        output <- as.matrix(gsea$output)
        output <- round(output, digits=4)
        output <- data.frame( contrast=rownames(output), output)
        if(!DEV) {
            output$p <- NULL
            output$rho <- NULL
        }
        
        color_fx = as.numeric(output[,"NES"])
        color_fx[is.na(color_fx)] <- 0  ## yikes...        
        numeric.cols <- which(sapply(output, is.numeric))
        numeric.cols
        
        DT::datatable(output, class='compact cell-border stripe',
                      rownames=FALSE,
                      extensions = c('Scroller'),
                      ##selection='none',
                      selection = list(mode='single', target='row', selected=1),
                      ##selection = list(target='row', selected=1),
                      fillContainer=TRUE,                      
                      options = list(
                          dom = 'lrtip',
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, scrollY = tabH, scroller=TRUE,
                          deferRender=FALSE)
                      ) %>%  ## end of options.list 
            DT::formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("NES",
                                background = color_from_middle(color_fx,'lightblue','#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
        
    })

       
    getEnrichmentGeneTable <- shiny::reactive({
        
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        shiny::req(ngs)

        dbg("[getEnrichmentGeneTable] reacted!")
        
        gsea <- sigCalculateGSEA()
        if(is.null(gsea)) return(NULL)

        dbg("[getEnrichmentGeneTable] 1:")
        
        i=1
        i <- enrichmentContrastTable$rows_selected()
        if(is.null(i) || length(i)==0) return(NULL)

        dbg("[getEnrichmentGeneTable] 2:")
        
        meta <- pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        fc <- meta$fc
        qv <- meta$qv
        rownames(fc) <- toupper(rownames(fc))
        rownames(qv) <- toupper(rownames(qv))

        dbg("[getEnrichmentGeneTable] 3:")
        
        contr <- rownames(gsea$output)[i]
        fc <- fc[,contr,drop=FALSE]
        ##qv <- qv[,contr,drop=FALSE]

        dbg("[getEnrichmentGeneTable] 4:")
        
        gset <- getCurrentMarkers()
        if(is.null(gset)) return(NULL)

        gset <- setdiff(toupper(gset),c("",NA))
        genes <- intersect(gset,rownames(fc))
        dd1 <- setdiff(genes,rownames(fc))
        dd2 <- setdiff(genes,rownames(qv))        
        fc <- fc[genes,,drop=FALSE]
        qv <- qv[genes,,drop=FALSE]

        dbg("[getEnrichmentGeneTable] 5:")
        
        gene.tt <- substring(GENE.TITLE[toupper(rownames(fc))],1,40)
        names(gene.tt) <- rownames(fc)
        df <- data.frame(gene=rownames(fc), title=gene.tt, fc, check.names=FALSE)
        ##df <- df[order(-abs(df$FC)),]

        dbg("[getEnrichmentGeneTable] done!")
        
        df

    })

    enrichmentGeneTable.RENDER <- shiny::reactive({

        df <- getEnrichmentGeneTable()
        shiny::req(df)
        
        color_fx = as.numeric(df[,3:ncol(df)])
        color_fx[is.na(color_fx)] <- 0  ## yikes...

        numeric.cols <- colnames(df)[3:ncol(df)]
        numeric.cols
        
        DT::datatable(df, class='compact cell-border stripe',
                      rownames=FALSE,
                      extensions = c('Scroller'),
                      ## selection='none',
                      selection = list(mode='single', target='row', selected=NULL),
                      fillContainer=TRUE,
                      options=list(
                          dom = 'lrftip',
                          ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE, scrollY = tabH, scroller=TRUE,
                          deferRender=FALSE
                      )) %>%  ## end of options.list 
            DT::formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle(
                        numeric.cols,
                        background = color_from_middle(color_fx,'lightblue','#f5aeae'),
                        backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',
                        backgroundPosition = 'center')    
    })

    info.text1 = "<b>Enrichment by contrast.</b> Enrichment scores of query signature across all contrasts. The table summarizes the enrichment statistics of the gene list in all contrasts using the GSEA algorithm. The NES corresponds to the normalized enrichment score of the GSEA analysis.  "

    enrichmentContrastTable <- shiny::callModule(
        tableModule,
        id = "enrichmentContrastTable", 
        func = enrichmentContrastTable.RENDER,
        info.text = info.text1,
        caption2 = info.text1,                
        title = "Enrichment by contrasts", label="a",
        height = c(230,700)
    )

    info.text2 = "<b>Gene table.</b> Genes of the current signature corresponding to the selected contrast. Genes are sorted by decreasing (absolute) fold-change."
    enrichmentGeneTable <- shiny::callModule(
        tableModule,
        id = "enrichmentGeneTable", 
        func = enrichmentGeneTable.RENDER,
        info.text = info.text2,
        caption2 = info.text2,        
        title = "Genes in signature", label="b",
        height = c(360,700)
    )

  })
}  ## end-of-Board
