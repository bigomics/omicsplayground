##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

EnrichmentBoard <- function(id, inputData, selected_gxmethods)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns ## NAMESPACE
    
    fullH = 800
    rowH = 420  ## row height of panels
    imgH = 340  ## height of images
    tabV = "70vh"  ## height of tables
    tabH = 340  ## row height of panels
    tabH = "80vh"  ## height of tables
    
    gs_infotext = paste("Similar to the differential gene expression analysis, users can perform differential
        expression analysis on a geneset level in this page, which is also referred as gene set enrichment (GSE) analysis.
        The platform has more than 50.000 genesets (or pathways) in total that are divided into 30 geneset collections
        such as ",a_Hallmark,", ",a_MSigDB,", ",a_KEGG," and ",a_GO,". Users have to specify which comparison they want to
        visually analyze employing a certain geneset collection.<br><br>
        To ensure the statistical reliability, the platform performs Enrichment Analyses using multiple methods.
        The result from the statistical methods is displayed in <strong>Enrichment table</strong> panel.
        In the <strong>Top enriched</strong> panel, the top 10 differentially enriched geneses (pathways) are displayed. 
        In the <strong>Plots</strong> panel, a volcano plot of genes contained in the selected geneset and a barplot 
        of expressions per sample group are displayed. In the <strong>Compare</strong> panel, users can compare the
        differential expression status of that geneset for all other comparisons. Finally, volcano plots of genesets
        for all comparisons are displayed under the <strong>Volcano (all) </strong> tab. This allows users to have
        an overall picture across comparisons at the same time.<br><br>
        EXPERT MODE ONLY: To compare the different statistical methods, the <strong>Volcano (methods)</strong> 
        panel shows volcano plots of all methods. The <strong>FDR table</strong> panel reports the number of
        significant gene sets at different FDR thresholds for all contrasts.<br><br><br><br>
        <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=4' 
        frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' 
        allowfullscreen></iframe></center>"
    )
    
    GSET.DEFAULTMETHODS = c("gsva","camera","fgsea","fisher")
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$gs_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Enrichment Analysis Board</strong>"),
            shiny::HTML(gs_infotext),
            easyClose = TRUE, size="l" ))
    })

    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        meta <- ngs$gset.meta$meta
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons = sort(intersect(comparisons, names(meta)))
        shiny::updateSelectInput(session, "gs_contrast", choices=comparisons)

        ## get the computed geneset methods
        gset.methods = sort(colnames(meta[[1]]$fc))
        sel2 = c(intersect(GSET.DEFAULTMETHODS,gset.methods),gset.methods)
        sel2 = head(unique(sel2),3)

        shiny::updateCheckboxGroupInput(session, 'gs_statmethod',
                                 choices = sort(gset.methods),
                                 selected = sel2)
        
    })
    
    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        nn <- sapply(COLLECTIONS, function(k) sum(k %in% rownames(ngs$gsetX)))
        gsets.groups <- names(COLLECTIONS)[which(nn>=5)]
        gsets.groups <- c("<all>",sort(gsets.groups))
        sel = "<all>"
        hmark <- grep("^H$|hallmark|",gsets.groups,ignore.case=TRUE,value=TRUE)
        if(length(hmark)>0) sel <- hmark[1]
        shiny::updateSelectInput(session, "gs_features",choices=gsets.groups, selected=sel)
        
    })

    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================
    
    selected_gsetmethods <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        gset.methods0 = colnames(ngs$gset.meta$meta[[1]]$fc)
        ##test = head(intersect(GSET.DEFAULTMETHODS,gset.methods0),3) ## maximum three
        test = input$gs_statmethod
        test = intersect(test,gset.methods0) ## maximum three
        test
    })

    calcGsetMeta <- function(comparison, methods, ngs) {
        ##ngs <- inputData()
        mx = ngs$gset.meta$meta[[comparison]]
        if(is.null(mx)) return(NULL)
        mx.methods = colnames(unclass(mx$fc))
        mx.methods
        methods = intersect(methods, mx.methods)
        if(is.null(methods) || length(methods)==0) {
            cat("ERROR: calcGsetMeta:: no valid methods\n")
            return(NULL)
        }

        ## recalculate meta values
        pv = unclass(mx$p)[,methods,drop=FALSE]
        qv = unclass(mx$q)[,methods,drop=FALSE]
        fc = unclass(mx$fc)[,methods,drop=FALSE]

        ## !!!! Because the methods have all very difference "fold-change" !!!!
        ## estimators, we use the meta.fx (average of all genes in gset)
        fc = do.call(cbind,rep(list(mx$meta.fx),length(methods)))
        colnames(fc) <- methods
        
        pv[is.na(pv)] = 1
        qv[is.na(qv)] = 1
        fc[is.na(fc)] = 0
        score = fc * (-log10(qv))
        dim(pv)
        if(NCOL(pv)>1) {
            ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)
            ##fc = rowMeans(scale(fc,center=FALSE),na.rm=TRUE)  ## REALLY???
            fc = rowMeans(fc,na.rm=TRUE)  ## NEED RETHINK!!!
            ##pv = apply(pv,1,function(x) metap::allmetap(x,method="sumz")$p[[1]])
            ##pv = apply(pv,1,vec.combinePvalues,method="stouffer")
            ##qv = p.adjust(pv, method="fdr")
            pv = apply(pv,1,max,na.rm=TRUE)
            qv = apply(qv,1,max,na.rm=TRUE)
            ##score = rowMeans(scale(score,center=FALSE),na.rm=TRUE)
            score = rowMeans(apply(score, 2, ss.rank),na.rm=TRUE)
        }

        meta = cbind( score=score, fc=fc, pv=pv, qv=qv)
        rownames(meta) <- rownames(mx)
        colnames(meta) = c("score","fc","pv","qv")  ## need
        return(meta)
    }
    
    getFullGeneSetTable <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)
        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        if(!(comp %in% names(ngs$gset.meta$meta))) return(NULL)
        mx = ngs$gset.meta$meta[[comp]]
        dim(mx)
        
        outputs = NULL
        gsmethod = colnames(unclass(mx$fc))
        gsmethod <- input$gs_statmethod   
        if(is.null(gsmethod) || length(gsmethod)==0) return(NULL)

        lfc <- as.numeric(input$gs_lfc)
        fdr <- as.numeric(input$gs_fdr)
        
        ## filter gene sets for table
        gsfeatures="<all>"
        gsfeatures = input$gs_features
        if(is.null(input$gs_features)) return(NULL)
        if(1 && !(gsfeatures %in% c(NA,"","*","<all>"))  &&
           gsfeatures %in% names(COLLECTIONS)) {
            ##grp = paste(paste0("^",gsfeatures,":"),collapse="|")
            ##sel <- grep(grp,rownames(mx),ignore.case=TRUE)
            sel <- intersect(rownames(mx),COLLECTIONS[[gsfeatures]]) 
            mx = mx[sel,,drop=FALSE]
        }
        ##outputs = lapply(outputs, function(m) m[rownames(mx),])
        
        rpt = NULL
        ##length(gsmethod)==1 && any(grepl("gsea",gsmethod))
        if(is.null(outputs) || length(gsmethod)>1) {
            
            ## show meta-statistics table (multiple methods)        
            pv = unclass(mx$p)[,gsmethod,drop=FALSE]
            qv = unclass(mx$q)[,gsmethod,drop=FALSE]
            fx = unclass(mx$fc)[,gsmethod,drop=FALSE]

            ## !!!! Because the methods have all very difference "fold-change" !!!!
            ## estimators, we use the meta.fx (average of all genes in gset)
            fx = do.call(cbind,rep(list(mx$meta.fx),length(gsmethod)))
            colnames(fx) <- gsmethod

            pv[is.na(pv)] = 1
            qv[is.na(qv)] = 1
            fx[is.na(fx)] = 0
            
            is.sig <- (qv <= fdr & abs(fx) >= lfc)
            stars <- sapply(rowSums(is.sig,na.rm=TRUE), star.symbols, pch='\u2605')
            names(stars) <- rownames(mx)
            
            ##------------ calculate META parameters ----------------
            meta <- calcGsetMeta(comp, gsmethod, ngs=ngs)
            meta <- meta[rownames(mx),,drop=FALSE]        
            dim(meta)
            gset.size = Matrix::colSums(ngs$GMT[,rownames(mx),drop=FALSE]!=0)
            names(gset.size) <- rownames(mx)
            
            ## ---------- report *average* group expression FOLD CHANGE
            ## THIS SHOULD BETTER GO DIRECTLY WHEN CALCULATING GSET TESTS
            ##            
            s1 <- names(which(ngs$model.parameters$exp.matrix[,comp]>0))
            s0 <- names(which(ngs$model.parameters$exp.matrix[,comp]<0))
            jj <- colnames(ngs$GMT)
            jj <- rownames(mx)
            
            gsdiff.method <- "fc"  ## OLD default
            if(gsdiff.method=="gs") {
                AveExpr1 <- rowMeans(ngs$gsetX[jj,s1])
                AveExpr0 <- rowMeans(ngs$gsetX[jj,s0])
                meta.fc <- AveExpr1 - AveExpr0
            } else {
                ## WARNING!!! THIS STILL ASSUMES GENES AS rownames(ngs$X)
                ## and rownames(GMT)
                fc <- ngs$gx.meta$meta[[comp]]$meta.fx  ## stable
                names(fc) <- rownames(ngs$gx.meta$meta[[comp]])
                pp <- intersect(rownames(ngs$GMT),names(fc))
                rnaX <- ngs$X
                
                ## check if multi-omics
                is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]",names(fc)))
                is.multiomics
                if(is.multiomics) {
                    ii <- grep("\\[gx\\]|\\[mrna\\]",names(fc))
                    fc <- fc[ii]
                    rnaX <- ngs$X[names(fc),]
                    names(fc) <- sub(".*:|.*\\]","",names(fc))
                    rownames(rnaX) <- sub(".*:|.*\\]","",rownames(rnaX))
                    pp <- intersect(rownames(ngs$GMT),names(fc))
                    length(pp)
                }
                
                G <- Matrix::t(ngs$GMT[pp,jj] != 0)
                ngenes <- Matrix::rowSums(G)
                ## meta.fc <- as.vector(G %*% fc[pp] / ngenes)
                ## names(meta.fc) <- rownames(G)
                meta.fc <- ngs$gset.meta$meta[[comp]]$meta.fx
                names(meta.fc) <- rownames(ngs$gset.meta$meta[[comp]])
                
                AveExpr1 <- Matrix::rowMeans(G %*% rnaX[pp,s1]) / ngenes
                AveExpr0 <- Matrix::rowMeans(G %*% rnaX[pp,s0]) / ngenes
                remove(rnaX)
            }
            
            ## TWIDDLE means to reflect foldchange... 
            mean0 <- (AveExpr0 + AveExpr1)/2
            AveExpr1 <- mean0 + meta.fc/2
            AveExpr0 <- mean0 - meta.fc/2
            
            ##
            gs <- intersect(names(meta.fc),rownames(meta))
            length(gs)
            
            rpt = data.frame( size = gset.size[gs],
                             logFC = meta.fc[gs],
                             meta.q = meta[gs,"qv"],
                             stars  = stars[gs],
                             AveExpr0 = AveExpr0[gs],
                             AveExpr1 = AveExpr1[gs])
            
            ## add extra p/q value columns
            jj <- match(gs, rownames(mx))
            rpt <- cbind( rpt, q=qv[jj,])
            
            ##rownames(rpt) = gs
        }  else {
            ## show original table (single method)
            rpt = outputs[[gsmethod]]
        }

        ##rpt <- rpt[order(-abs(rpt$logFC)),]
        rpt <- rpt[order(-rpt$logFC),] ## positive        
        rpt = data.frame(rpt)
        
        return(rpt)
    })


    getFilteredGeneSetTable <- shiny::reactive({        
        
        if(is.null(input$gs_showall) || length(input$gs_showall)==0) return(NULL)
        if(is.null(input$gs_top10) || length(input$gs_top10)==0) return(NULL)
                
        res <- getFullGeneSetTable()
        
        ## just show significant genes
        if(!input$gs_showall && nrow(res)>0 ) {
            ##nmeth <- length(input$gs_statmethod)
            ##sel <- which(res$stars == star.symbols(nmeth))
            ##sel <- which(nchar(res$stars) == nmeth)
            lfc <- as.numeric(input$gs_lfc)
            fdr <- as.numeric(input$gs_fdr)
            dbg("[EnrichmentBoard::getFilteredGeneSetTable] lfc = ",lfc)
            dbg("[EnrichmentBoard::getFilteredGeneSetTable] fdr = ",fdr)            
            is.sig <- (abs(res$logFC) >= lfc & res$meta.q <= fdr)
            dbg("[EnrichmentBoard::getFilteredGeneSetTable] is.sig = ",table(is.sig))
            res <- res[is.sig,,drop=FALSE]
        }
        
        ## just show top 10        
        if(input$gs_top10 && nrow(res)>10 && length(input$gs_top10) ) {
            fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(res),value=TRUE)[1]
            dbg("[EnrichmentBoard::getFilteredGeneSetTable] fx.col = ",fx.col)
            fx  = as.numeric(res[,fx.col])
            names(fx) = rownames(res)
            pp <- unique(c(head(names(sort(-fx[which(fx>0)])),10),
                           head(names(sort(fx[which(fx<0)])),10)))
            res = res[pp,,drop=FALSE]
            fx  = as.numeric(res[,fx.col])
            res = res[order(-fx),,drop=FALSE]
        }

        ## limit to 1000 rows???
        ## rpt <- head(rpt, 1000)
        res <- data.frame(res)        
        
        if(nrow(res)==0) {
            shiny::validate(shiny::need(nrow(res) > 0, "warning. no genesets passed current filters."))
            return(NULL)
        }

        return(res)
    })


    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================
    
    plotTopEnriched <- function(ngs, rpt, comp, ntop, rowcol)
    {
        if(is.null(ngs)) return(NULL)
 
        gx.meta <- ngs$gx.meta$meta[[comp]]
        ##rnk0 <- gx.meta[,"fc"][,"trend.limma"]
        ##names(rnk0) = gx.meta[,"gene_name"]        
        rnk0 <- gx.meta$meta.fx
        names(rnk0) = ngs$genes[rownames(gx.meta),"gene_name"]
        rnk0 = rnk0 - mean(rnk0,na.rm=TRUE)  ## scaling/centering should be done in calculation...                
        fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(rpt))[1]
        qv.col = grep("meta.q|q$",colnames(rpt))[1]
        fx = rpt[,fx.col]
        qv = rpt[,qv.col]
        names(qv) <- names(fx) <- rownames(rpt)
        top <- rownames(rpt)
        top <- setdiff(top, c(NA,"NA"))
        if(is.null(top) || is.na(top[1])) return(NULL)
        
        par(mfrow=rowcol)
        if(ntop==1) {
            par(mar=c(1,6,2,6), mgp=c(1.6,0.6,0), oma=c(0.1,1,0,0.1))
        } else {
            par(mar=c(0.2,1.8,2.3,0.1), mgp=c(1.6,0.6,0), oma=c(0.1,1,0,0.1))
        }
        
        for(i in 1:ntop) {
            gs <- top[i]
            if(i > length(top) || is.na(gs) ) {
                frame()
            } else {
                genes = names(which(ngs$GMT[,gs]!=0))
                genes = toupper(genes)                
                names(rnk0) <- toupper(names(rnk0))
                ylab = ""
                ## if(i %in% c(1,6)) ylab = "Ranked list metric"
                if(i%%rowcol[2] == 1) ylab = "Rank metric"
                xlab=""
                gs1 = breakstring(gs,28,50,force=FALSE)                
                if(ntop==1) {
                    gs1 = breakstring(gs,100,200,force=FALSE)                
                    xlab = "Rank in ordered dataset"
                    ylab = "Rank metric"
                }
                gsea.enplot(rnk0, genes, names=NULL, ##main=gs,
                            main=gs1, xlab=xlab, ylab=ylab,
                            lab.line = c(0,1.8), cex.lab=0.75,
                            cex.main=0.78, len.main=200)
                qv1 = formatC(qv[gs],format="e", digits=2)
                legend("topright", paste("q=",qv1), bty="n",cex=0.85)
            }
        }
    }
    
    ## Top enriched    
    topEnriched.RENDER <- shiny::reactive({
        
        ngs <- inputData()
        shiny::req(ngs)       
        ##rpt <- getFullGeneSetTable()
        rpt <- getFilteredGeneSetTable()
        ##if(is.null(rpt)) return(NULL)
        shiny::req(rpt, input$gs_contrast)
        if(is.null(rpt)) return(NULL)
        
        comp=1
        comp = input$gs_contrast
        if(!(comp %in% names(ngs$gx.meta$meta))) return(NULL)

        ## selected
        sel = as.integer(gseatable$rows_selected())
        sel.gs <- NULL
        if(!is.null(sel) && length(sel)>0) sel.gs = rownames(rpt)[sel]
        
        ## filter on active rows (using search)
        ##ii  <- gseatable$rows_all()
        ##ii  <- gseatable$rows_current()        
        ##if(is.null(ii) || length(ii)==0) return(NULL)        
        ii  <- gseatable$rows_selected()
        jj  <- gseatable$rows_current()
        shiny::req(jj)

        dbg("[topEnriched.RENDER] dim.rpt", dim(rpt))
        if(nrow(rpt)==0) return(NULL)

        ## ENPLOT TYPE        
        if(length(ii)>0) {
            itop = ii[1]
        } else {
            itop <- head(jj,15)
        }
        if(length(itop)==1) {
            plotTopEnriched(ngs, rpt[itop,,drop=FALSE], comp=comp, ntop=1, rowcol=c(1,1))
        } else {
            plotTopEnriched(ngs, rpt[itop,,drop=FALSE], comp=comp, ntop=15, rowcol=c(3,5))
        }


    })
        
    topEnriched_text = "This plot shows the <strong>top enriched</strong> gene sets for the selected comparison in the <code>Contrast</code> settings. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group."

    shiny::callModule(
        plotModule,
        id = "topEnriched", label="a",
        func = topEnriched.RENDER,
        func2 = topEnriched.RENDER,
        info.text = topEnriched_text,
        height = c(imgH,720), width = c('auto',1500), 
        pdf.height = 6, pdf.width = 12, ## pdf.pointsize=16,
        res = c(90, 120),
        title = "Top enriched gene sets",
        ##caption = topEnriched_caption,
        add.watermark = WATERMARK
    )

            
    plotEnrichFreq <- function(ngs, rpt, ntop, ngenes, gset.weight, fcweight)
    {
        
        fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(rpt))[1]
        fx = rpt[,fx.col]
        names(fx) <- rownames(rpt)
        
        top <- rownames(rpt)        
        top <- head(top,ntop)
        if(!all(top %in% colnames(ngs$GMT))) return(NULL)
        
        F <- 1*(ngs$GMT[,top,drop=FALSE]>0)
        F <- as.matrix(F)
        wt = FALSE
        if(gset.weight) {
            F <- Matrix::t(Matrix::t(F)  / Matrix::colSums(F,na.rm=TRUE))
            wt = TRUE
        }
        F <- Matrix::t( Matrix::t(F) * sign(fx[top]))
        if(fcweight) {
            F <- Matrix::t( Matrix::t(F) * abs(fx[top]))
            wt = TRUE
        } 
        F <- head(F[order(-Matrix::rowSums(abs(F))),,drop=FALSE], ngenes)
        F <- F[order(-Matrix::rowSums(F)),,drop=FALSE]

        sel.zero <- which(Matrix::rowSums(abs(F)) < 1e-4)
        if(length(sel.zero)) rownames(F)[sel.zero] = ""
        
        par(mfrow=c(1,1), mar=c(6,4,2,0.5), mgp=c(2.2,0.8,0))
        col1 = grey.colors(ncol(F),start=0.15)
        ylab = ifelse(wt, "weighted frequency", "frequency")
        barplot(t(F), beside=FALSE, las=3, cex.names=0.90, col=col1,
                ylab=ylab)
    }

    topEnrichedFreq.RENDER <- shiny::reactive({

        ngs <- inputData()

        rpt <- getFilteredGeneSetTable()
        ##if(is.null(rpt)) return(NULL)
        shiny::req(ngs, rpt, input$gs_contrast)

        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        if(!(comp %in% names(ngs$gx.meta$meta))) return(NULL)
        
        ## filter on active rows (using search)
        ii <- gseatable$rows_current()
        rpt <- rpt[ii,,drop=FALSE]
        if(nrow(rpt)==0) return(NULL)
        ntop <- as.integer(input$gs_enrichfreq_ntop)
        gset.weight <- input$gs_enrichfreq_gsetweight
        fcweight <- input$gs_enrichfreq_fcweight

        plotEnrichFreq(ngs, rpt, ntop=ntop, ngenes=35, gset.weight, fcweight)

    })

    topEnrichedFreq.RENDER2 <- shiny::reactive({

        ngs <- inputData()

        rpt <- getFilteredGeneSetTable()
        ##if(is.null(rpt)) return(NULL)
        shiny::req(ngs, rpt, input$gs_contrast)

        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        if(!(comp %in% names(ngs$gx.meta$meta))) return(NULL)
        
        ## filter on active rows (using search)
        ##ii <- gseatable$rows_all()
        ii <- gseatable$rows_current()        
        rpt <- rpt[ii,,drop=FALSE]
        if(nrow(rpt)==0) return(NULL)
        ntop <- as.integer(input$gs_enrichfreq_ntop)
        gset.weight <- input$gs_enrichfreq_gsetweight
        fcweight <- input$gs_enrichfreq_fcweight

        plotEnrichFreq(ngs, rpt, ntop=ntop, ngenes=60, gset.weight, fcweight)

    })

    
    topEnrichedFreq_text = "<strong>Gene frequency.</strong> The plot shows the number of times a gene is present in the top-N genesets sorted by frequency. Genes that are frequently shared among the top enriched gene sets may suggest driver genes."
    
    topEnrichedFreq.opts = shiny::tagList(
        withTooltip( shiny::radioButtons(ns('gs_enrichfreq_ntop'),"Number of top sets",
                             c(5,10,15),inline=TRUE, selected=15),
               "Number of top genesets to consider for counting the gene frequency."),
        withTooltip( shiny::checkboxInput(ns('gs_enrichfreq_gsetweight'),
                              "Weight by geneset size", TRUE),
               "Weight by (inverse) gene set size."), 
        withTooltip( shiny::checkboxInput(ns('gs_enrichfreq_fcweight'),
                              "Weight by FC", TRUE),
               "Weight by fold-change of current contrast.")         
    )
    
    shiny::callModule(
        plotModule,
        id = "topEnrichedFreq", label="b",
        func = topEnrichedFreq.RENDER,
        func2 = topEnrichedFreq.RENDER2,
        options = topEnrichedFreq.opts,
        info.text = topEnrichedFreq_text,
        height = c(imgH,500), width = c('auto',1200),
        res = c(68,100),
        pdf.width = 10, pdf.height = 5, 
        title = "Frequency in top gene sets",
        ##caption = topEnrichedFreq_caption,
        add.watermark = WATERMARK
    )
    

    ##================================================================================
    ## Plots tab
    ##================================================================================

    
    ##comp0=colnames(ngs$model.parameters$contr.matrix)[1]
    getcolors <- function(ngs, comp0) {   
        ## get colors (what a mess...)
        contr.matrix <- ngs$model.parameters$contr.matrix

        exp.matrix <- ngs$model.parameters$exp.matrix        
        xgroup = as.character(ngs$Y$group)
        
        grp.name <- strsplit(comp0,split="[._ ]vs[._ ]")[[1]]
        grp.name <- c(grp.name, "other")
        xsign <- sign(exp.matrix[,comp0])
        xgroup = grp.name[1*(xsign>0) + 2*(xsign<0) + 1*(xsign==0)]
        table(xgroup)
        
        names(xgroup) = rownames(ngs$Y)
        table(xgroup)
        samples = names(which(exp.matrix[,comp0]!=0))
        
        xgroup1 <- xgroup[samples]
        table(xgroup1)
        ngrp <- length(unique(xgroup1))
        grp.klr = c("grey90",rep(RColorBrewer::brewer.pal(12,"Paired"),99)[1:ngrp])
        names(grp.klr) <- c("other",as.character(sort(unique(xgroup1))))
        grp.klr
        
        xgroup2 <- as.character(xgroup)
        xgroup2[which(!(xgroup %in% xgroup1))] <- "other"
        sample.klr = grp.klr[xgroup2]
        names(sample.klr) <- rownames(ngs$samples)
        table(sample.klr)
        list(samples=sample.klr, group=grp.klr)
    }    
    
    
    ##----------------------------------------------------------------------
    ## 0: Volcano plot in gene space
    ##----------------------------------------------------------------------
    subplot.MAR = c(3,3.5,1.5,0.5)
    subplot.MAR = c(2.8,4,4,0.8)
    
    subplot_volcano.RENDER <- shiny::reactive({
        
        par(mfrow=c(1,1), mgp=c(1.2,0.4,0), oma=c(0,0,0,0.4) )
        par(mar= subplot.MAR)
        ## par(mar= c(2.3,4,1.9,0))
        
        ngs <- inputData()    
        shiny::req(ngs)
        
        comp=1;gs=1
        comp = input$gs_contrast
        ngs <- inputData()
        shiny::req(ngs)
        
        gxmethods <- selected_gxmethods() ## from module-expression
        shiny::req(gxmethods)
        
        gx.meta <- ngs$gx.meta$meta[[comp]]
        ## limma1 = sapply(gx.meta[,c("fc","p","q")],function(x) x[,"trend.limma"])
        meta.q <- apply(gx.meta$q[,gxmethods,drop=FALSE],1,max)  ## max q-value
        limma1 = data.frame( meta.fx=gx.meta$meta.fx, meta.q=meta.q)
        gx.annot <- ngs$genes[rownames(gx.meta),c("gene_name","gene_title")]
        ##limma = cbind( gx.meta[,c("gene_name","gene_title")], limma1)
        limma = cbind(gx.annot, limma1)
    
        gs = gset_selected()
        ##if(is.null(gs)) return(NULL)
        if(is.null(gs) || length(gs)==0) {
            frame()
            text(0.5,0.5,"Please select a geneset",col="grey50")
            return()
        }
        gs <- gs[1]
        
        ##sel.genes = names(which(ngs$GMT[,gs]!=0))
        ##gset <- GSETS[[gs]]
        gset <- getGSETS(gs)[[1]]
        dbg("[subplot_volcano.RENDER] head.gset = ", head(gset,5) )        
        jj = match(toupper(gset), toupper(limma$gene_name))
        sel.genes <- setdiff(limma$gene_name[jj],c(NA,""," "))

        dbg("[subplot_volcano.RENDER] head.sel.genes = ", head(sel.genes,5) )

        fdr = 1
        fdr = as.numeric(input$gs_fdr)        
        fc.genes = as.character(limma[,grep("^gene$|gene_name",colnames(limma))])
        fx = limma[,grep("logFC|meta.fx|fc",colnames(limma))[1]]
        qval = limma[,grep("^q|adj.P.Val|meta.q|qval|padj",colnames(limma))[1]]
        ##sig.genes = fc.genes[which(qval <= fdr & abs(fx) > 0.1)]
        
        qval <- pmax(qval,1e-12)  ## prevent q=0
        qval[which(is.na(qval))] <- 1
        xlim = c(-1,1)*max(abs(fx),na.rm=TRUE)
        ylim = c(0,12)
        ylim = c(0, max(12, 1.1*max(-log10(qval),na.rm=TRUE)))
        ylim
        
        lfc=0.20
        lfc = as.numeric(input$gs_lfc)
        
        ##par(mar=c(4,3,3,1), mgp=c(2.0,0.8,0), oma=c(1,1.5,1,1.5) )
        gx.volcanoPlot.XY( x=fx, pv=qval, gene=fc.genes,
                          render="canvas", n=5000, nlab=10, 
                          xlim=xlim, ylim=ylim, ## hi.col="#222222",
                          use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                          cex=0.9, lab.cex=1.3,
                          cex.main=0.8, cex.axis=0.9,
                          xlab="fold change (log2)",
                          ylab="significance (log10q)",
                          highlight=sel.genes)
        gs = breakstring(gs,50)
        title(gs, cex.main=0.85)        
    })
    
    
    subplot_volcano.PLOTLY <- shiny::reactive({
        
        ##par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0.5,0.2)*2 )
        ##par(mar=subplot.MAR)
        
        ngs <- inputData()    
        shiny::req(ngs)
        
        comp=1;gs=1
        comp = input$gs_contrast
        ngs <- inputData()
        shiny::req(ngs)
        
        gxmethods <- selected_gxmethods() ## from module-expression
        shiny::req(gxmethods)

        gx.meta <- ngs$gx.meta$meta[[comp]]
        ##m1 <- intersect(c("trend.limma","notrend.limma","ttest"),colnames(gx.meta$p))[1]
        ##m1
        ##limma1 = sapply(gx.meta[,c("fc","p","q")],function(x) x[,m1])
        meta.q <- apply(gx.meta$q[,gxmethods,drop=FALSE],1,max,na.rm=TRUE)
        limma1 = data.frame( meta.fx=gx.meta$meta.fx, meta.q=meta.q)
        gx.annot <- ngs$genes[rownames(gx.meta),c("gene_name","gene_title")]
        ##limma = cbind( gx.meta[,c("gene_name","gene_title")], limma1)
        limma = cbind(gx.annot, limma1)
    
        gs = gset_selected()
        if(is.null(gs)) return(NULL)
        gs <- gs[1]        
        ##sel.genes = names(which(ngs$GMT[,gs]!=0))
        ##gset <- GSETS[[gs]]
        gset <- getGSETS(gs)[[1]]
        jj = match(toupper(gset), toupper(limma$gene_name))
        sel.genes <- setdiff(limma$gene_name[jj],c(NA,""," "))

        fdr = 1
        fdr = as.numeric(input$gs_fdr)
        
        fc.genes = as.character(limma[,grep("^gene$|gene_name",colnames(limma))])
        fx = limma[,grep("logFC|meta.fx|fc",colnames(limma))[1]]
        qval = limma[,grep("^q|adj.P.Val|meta.q|qval|padj",colnames(limma))[1]]
        ##sig.genes = fc.genes[which(qval <= fdr & abs(fx) > 0.1)]
        
        qval <- pmax(qval,1e-12)  ## prevent q=0
        qval[which(is.na(qval))] <- 1
        xlim = c(-1,1)*max(abs(fx),na.rm=TRUE)
        ylim = c(0,12)
        ylim = c(0, max(12, 1.1*max(-log10(qval),na.rm=TRUE)))
        ylim
        
        lfc=0.20
        lfc = as.numeric(input$gs_lfc)
        y <- -log10(qval+1e-20)
        scaled.fx <- scale(fx,center=FALSE)
        scaled.y <- scale(y,center=FALSE)

        impt <- function(g) {
            j = match(g, fc.genes)
            x1 = scaled.fx[j]
            y1 = scaled.y[j]
            x = sign(x1)*(x1**2 + 0.25*y1**2)
            names(x)=g
            x
        }
        lab.genes = c( head(sel.genes[order(impt(sel.genes))],10),
                      head(sel.genes[order(-impt(sel.genes))],10) )
        
        plotlyVolcano(
            x = fx, y = y, names=fc.genes,
            source = "plot1", marker.type = "scattergl",
            highlight = sel.genes, label = lab.genes,
            group.names = c("group1","group0"),
            ##xlim=xlim, ylim=ylim, ## hi.col="#222222",
            ##use.fdr=TRUE,
            psig = fdr, lfc = lfc,
            xlab = "effect size (log2FC)",
            ylab = "significance (-log10q)",
            marker.size = 4,
            displayModeBar = FALSE,
            showlegend = FALSE) %>%
            plotly::layout( margin = list(b=60) )        

    })

    ##----------------------------------------------------------------------
    ## 0: Single enrichment plot
    ##----------------------------------------------------------------------
    
    subplot_enplot.RENDER <- shiny::reactive({
    ##subplot_enplot.RENDER <- shiny::reactive({    

        dbg("[subplot_enplot.RENDER] reacted")        
        pgx <- inputData()
        shiny::req(pgx)
        
        comp=1;gs=100
        gs="H:HALLMARK_TNFA_SIGNALING_VIA_NFKB"        
        comp = input$gs_contrast

        ## !!!!!!!! SHOULD BE SELECTED GX METHODS ONLY???
        fc <- pgx$gx.meta$meta[[comp]]$meta.fx
        names(fc) <- rownames(pgx$gx.meta$meta[[comp]])

        gs = gset_selected()
        if(is.null(gs)) return(NULL)
        ##gs.genes = GSETS[[gs]]
        gs.genes = getGSETS(gs)[[1]]
        
        par(mfrow=c(1,1), mgp=c(1.95,0.8,0), oma=c(0,0,0.4,0.4) )
        par(mar =  c(2.8,4,3.8,0.8))        
        p1 <- NULL
        ##p1 <- ggenplot(fc, gs.genes, main=gs )
        gsea.enplot(fc, gs.genes, main=gs, cex.lab=0.9,
                    len.main=65, main.line=1.7)
        return(p1)
    })

    subplot_enplot.RENDER2 <- shiny::reactive({

        dbg("[subplot_enplot.RENDER] reacted")        
        pgx <- inputData()
        shiny::req(pgx)
        
        comp=1;gs=100
        gs="H:HALLMARK_TNFA_SIGNALING_VIA_NFKB"        
        comp = input$gs_contrast

        ## !!!!!!!! SHOULD BE SELECTED GX METHODS ONLY???
        fc <- pgx$gx.meta$meta[[comp]]$meta.fx
        names(fc) <- rownames(pgx$gx.meta$meta[[comp]])
        
        gs = gset_selected()
        if(is.null(gs)) return(NULL)
        ##gs.genes = GSETS[[gs]]
        gs.genes = getGSETS(gs)[[1]]

        p1 <- gsea.enplotly(fc, gs.genes, main=gs)
        return(p1)
    })

    
    ##----------------------------------------------------------------------
    ## 1: Gene set activation {data-width=200}
    ##----------------------------------------------------------------------
    subplot_barplot.RENDER <- shiny::reactive({
        
        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0,0.4) )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        shiny::req(ngs)
        

        gset = rownames(ngs$gsetX)[1]
        gset = gset_selected()
        if(is.null(gset) || length(gset)==0 ) return(NULL)
        gset <- gset[1]
        if(!gset %in% rownames(ngs$gsetX)) return(NULL)
        
        comp0 = colnames(ngs$model.parameters$contr.matrix)[1]
        comp0 = input$gs_contrast
        
        grouped=TRUE
        grouped=FALSE
        grouped <- !input$gs_ungroup1
        has.design <- !is.null(ngs$model.parameters$design)
        collapse.others <- ifelse(has.design, FALSE, TRUE)
        ##collapse.others=TRUE
        
        ngrp <- length(unique(ngs$samples$group))
        srt <- ifelse(!grouped || ngrp>4, 30, 0)
        if(!grouped && ncol(ngs$X) > 15) srt <- 60
        pgx.plotExpression(
            ngs, gset, comp=comp0, logscale=TRUE, level="geneset",
            collapse.others=collapse.others, grouped=grouped,
            cex=1.1, srt=srt, main="", ylab="enrichment (avg logFC)")
        title(breakstring(gset,42,80), cex.main=0.85)
        
    })

    ##----------------------------------------------------------------------
    ## 2: Gene expression {data-width=200}
    ##----------------------------------------------------------------------
    subplot_geneplot.RENDER <- shiny::reactive({

        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0,0.4) )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        shiny::req(ngs)

        comp0 = colnames(ngs$model.parameters$contr.matrix)[1]
        comp0 = input$gs_contrast

        has.design <- !is.null(ngs$model.parameters$design)
        collapse.others <- ifelse(has.design, FALSE, TRUE)
        ##collapse.others=TRUE

        sel  = gene_selected()
        if(is.null(sel) || is.na(sel) || length(sel)==0) {
            frame()
        } else {
            ##gene = sel$gene
            probe = sel$probe
            gene = sel$gene
            if(length(probe)>1) {
                probe <- grep("\\[gx|\\[mrna",probe,value=TRUE)
            }
            ngrp <- length(unique(ngs$samples$group))
            grouped=TRUE
            grouped <- !input$gs_ungroup2
            srt <- ifelse(!grouped || ngrp>4, 30, 0)
            if(!grouped && ncol(ngs$X) > 15) srt <- 60
            pgx.plotExpression(
                ngs, probe, comp=comp0, logscale=TRUE, level="gene",
                collapse.others=collapse.others, grouped=grouped,
                srt=srt, main="")
            title(gene, cex.main=0.9)
        }
    })


    ##----------------------------------------------------------------------
    ## 3: Gene - gene set correlation
    ##----------------------------------------------------------------------
    subplot_scatter.RENDER <- shiny::reactive({

        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0,0.4) )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        shiny::req(ngs)
        

        gene = rownames(ngs$X)[1]
        sel  = gene_selected()
        gset = gset_selected()
        if(is.null(sel)) return(NULL)
        if(is.null(gset)) return(NULL)

        if(is.null(sel) || length(sel)==0) {
            frame()
        } else {
            gene = sel$gene
            gset <- gset[1]    
            gx = ngs$X[sel$probe,]
            sx = ngs$gsetX[gset,]
            if( length(gx)==0 || length(sx)==0 ||
                length(gx)!=length(sx) ) {
                frame()
                return(NULL)
            }
            ## get colors
            comp0 = "Th17_mut_2h_VS_Th17_wt_2h_BLA"
            comp0 = input$gs_contrast
            klrs = getcolors(ngs, comp0)
            klr = klrs$samples[names(sx)]
            klr = paste0(gplots::col2hex(klr),"99")
            
            cex1 = c(1.4,0.8,0.3)[cut(length(gx),c(0,100,500,99999))]
            gset1 = breakstring(substring(gset,1,80),32)
                                        #tt = paste( breakstring(gset,40,80), "\nvs.", gene,"expression")
            tt = paste( breakstring(gset,40,80), " vs. ", gene)
            base::plot( gx, sx, col=klr, main=tt,
                 ylab = "gene set enrichment",
                 xlab = paste(gene,"expression"),
                 cex.lab=1, pch=19, cex=1.0*cex1, cex.main=0.85)
            abline( lm(sx ~ gx), lty=2, lwd=0.7, col="black" )
        }    
    })

    ##----------------------------------------------------------------------
    ## Calling Modules
    ##----------------------------------------------------------------------

    subplot_volcano_text = "</b>Volcano plot.<b> Volcano-plot showing significance versus fold-change on the y and x axes, respectively. Genes in the gene set that is selected from the enrichment analysis <b>Table I</b> are highlighted in blue."
    subplot_barplot_text = "An enrichment barplot per sample group for the gene set that is selected from the enrichment analysis Table <code>I</code>. Samples can be ungrouped in the barplot by selecting <code>ungroup samples</code> from the plot <i>Settings</i>."
    subplot_geneplot_text = "An expression barplot per sample group for the gene that is selected from the genes Table <code>II</code>. Samples can be ungrouped in the barplot by selecting <code>ungroup samples</code> from the plot <i>Settings</i>."
    subplot_scatter_text = "A scatter plot of enrichment scores versus expression values across the samples for the gene set selected from the enrichment analysis Table <code>I</code> and the gene selected from the genes Table <code>II</code>."
    

    shiny::callModule(
        plotModule,
        id = "subplot_volcano", 
        func = subplot_volcano.RENDER,
        plotlib = "base",
        func2 = subplot_volcano.PLOTLY,
        plotlib2 = "plotly",
        ##info.text = subplot_volcano_text,
        ##pdf.width=6, pdf.height=6,
        ##options = shiny::tagList(),
        ##res = c(72,100),
        height = c(imgH,750), width=c('auto',900),
        title="Volcano plot", label="b",
        add.watermark = WATERMARK
    )

    shiny::callModule(
        plotModule,
        id="subplot_barplot", 
        func = subplot_barplot.RENDER,
        func2 = subplot_barplot.RENDER,
        info.text = subplot_barplot_text,
        pdf.width=6, pdf.height=6,
        res = c(72,100),
        height = c(imgH,750),
        width = c('auto',900),        
        options = shiny::tagList(
            withTooltip( shiny::checkboxInput(
                ns('gs_ungroup1'),'ungroup samples',FALSE),
                "Ungroup samples in the plot", placement="top",
                options = list(container = "body"))),
        title="Enrichment barplot", label="c",
        add.watermark = WATERMARK
    )

    shiny::callModule(
        plotModule,
        id="subplot_geneplot", 
        func = subplot_geneplot.RENDER,
        func2 = subplot_geneplot.RENDER,
        info.text = subplot_geneplot_text,
        pdf.width=6, pdf.height=6,
        res = c(78,100),
        height = imgH,
        options = shiny::tagList(
            withTooltip( shiny::checkboxInput(ns('gs_ungroup2'),'ungroup samples',FALSE),
                   "Ungroup samples in the plot", placement="top",
                   options = list(container = "body"))),
        title = "Expression barplot", label="c",
        add.watermark = WATERMARK
    )

    shiny::callModule(
        plotModule,
        id="subplot_scatter",
        func = subplot_scatter.RENDER,
        func2 = subplot_scatter.RENDER,
        info.text = subplot_scatter_text,
        pdf.width=6, pdf.height=6, res=72,
        height = imgH,
        title = "Enrichment vs. expression", label="d",
        add.watermark = WATERMARK
    )

    shiny::callModule(
        plotModule,
        id = "subplot_enplot", 
        func = subplot_enplot.RENDER,
        func2 = subplot_enplot.RENDER2,         
        plotlib="base", plotlib2="plotly",
        info.text = '',
        pdf.width=6, pdf.height=4,
        res = c(68,110),
        height = imgH, 
        title="Enrichment plot", label="a",
        add.watermark = WATERMARK
    )
   

    ##================================================================================
    ## Compare
    ##================================================================================

    compare.RENDER <- shiny::reactive({
        
        ngs <- inputData()
        shiny::req(ngs,input$gs_contrast)
        
        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        
        gset = rownames(ngs$gsetX)[1]
        gset <- gset_selected()
        if(is.null(gset)) return(NULL)
        gset <- gset[1]

        score <- sapply(ngs$gset.meta$meta, function(x) x[gset,"meta.fx"])

        top.up   <- names(sort(score[which(score>0)],decreasing=TRUE))
        top.dn <- names(sort(score[which(score<0)]))
        genes    <- names(which(ngs$GMT[,gset]!=0))
        genes    <- toupper(sub(".*:","",genes))
        gx.meta  <- ngs$gx.meta$meta

        gsmethods <-  selected_gsetmethods()
        
        par(mfrow=c(2,5), mar=c(0.5,3.2,2.6,0.5), mgp=c(2,0.8,0))
        i=1
        for(i in 1:5) {
            if(i > length(top.up)) {
                frame()
            } else {
                cmp <- top.up[i]
                rnk0 <- gx.meta[[cmp]]$meta.fx
                names(rnk0) <- rownames(gx.meta[[1]])
                names(rnk0) <- toupper(sub(".*:","",names(rnk0)))

                ##qv0 <- ngs$gset.meta$meta[[cmp]][gset,"meta.q"]
                gs.meta <- ngs$gset.meta$meta[[cmp]]
                qv0 <- max(gs.meta[gset,"q"][,gsmethods],na.rm=TRUE)                
                
                gs1 = breakstring(gset,28,50,force=FALSE)
                cmp <- paste0(gset,"\n@",cmp)
                gsea.enplot(rnk0, genes, names=NULL, ##main=gs,
                            main=cmp, xlab="",
                            cex.main=0.80, len.main=72)
                qv1 = formatC(qv0,format="e", digits=2)
                legend("topright", paste("q=",qv1), bty="n",cex=0.85)
            }
        }
        for(i in 1:5) {
            if(i > length(top.dn)) {
                frame()
            } else {
                cmp <- top.dn[i]
                rnk0 <- gx.meta[[cmp]]$meta.fx
                names(rnk0) <- rownames(gx.meta[[1]])
                names(rnk0) <- toupper(sub(".*:","",names(rnk0)))

                ##qv0 <- ngs$gset.meta$meta[[cmp]][gset,"meta.q"]
                gs.meta <- ngs$gset.meta$meta[[cmp]]
                qv0 <- max(gs.meta[gset,"q"][,gsmethods],na.rm=TRUE)                

                gs1 = breakstring(gset,28,50,force=FALSE)
                cmp <- paste0(gset,"\n@",cmp)
                gsea.enplot(rnk0, genes, names=NULL, ##main=gs,
                            main=cmp, xlab="",
                            cex.main=0.80, len.main=72)
                qv1 = formatC(qv0,format="e", digits=2)
                legend("topright", paste("q=",qv1), bty="n",cex=0.85)
            }
        }
        
    })


    compare_text = "Under the <strong>Compare</strong> tab, enrichment profiles of the selected geneset in enrichment Table <code>I</code> can be visualised against all available contrasts."

    compare_module_opts = shiny::tagList()
    
    shiny::callModule(
        plotModule,
        id = "compare",
        func = compare.RENDER,
        func2 = compare.RENDER,
        options = compare_module_opts,
        height = c(imgH,450), width = c("auto",1500), res=95,
        pdf.width=14, pdf.height=4, 
        title = "Enrichment of geneset across multiple contrasts",
        info.text = compare_text,
        ##caption = compare_caption,
        add.watermark = WATERMARK
    )


    ##================================================================================
    ## Volcano (all)
    ##================================================================================

    volcanoAll.RENDER <- shiny::reactive({
        ##renderPlotly({

        ngs = inputData()
        shiny::req(ngs)
        if(is.null(input$gs_features)) return(NULL)
        
        meta = ngs$gset.meta$meta
        gsmethod = colnames(meta[[1]]$fc)
        gsmethod <- input$gs_statmethod
        if(is.null(gsmethod) || length(gsmethod)==0) return(NULL)

        fdr=1;lfc=1
        fdr = as.numeric(input$gs_fdr)    
        lfc = as.numeric(input$gs_lfc)
        sel.gsets <- NULL
        sel.gsets <- rownames(meta[[1]])
        sel.gsets = COLLECTIONS[[1]]
        sel.gsets = COLLECTIONS[[input$gs_features]]

        i=1
        mx.list <- list()
        for(i in 1:length(meta)) {
            mx.list[[i]] <- calcGsetMeta(i, gsmethod, ngs=ngs)
        }
        names(mx.list) <- names(meta)
        
        Q <- lapply(mx.list, function(mx) mx[,"qv"])
        names(Q) <- names(mx.list) 
        
        ## select maximum 24 comparisons (because of space...)
        q.score <- sapply(Q, function(q) mean(tail(sort(-log10(1e-99+q)),100)))
        sel <- head(names(sort(q.score, decreasing=TRUE)),20)
        ## sel <- sort(sel)
        Q <- Q[which(names(Q) %in% sel)]
        mx.list <- mx.list[names(Q)]
        ##ymax <- 1.2 * max(-log10(1e-99 + unlist(Q)), na.rm=TRUE)
        nlq <- -log10(1e-99 + unlist(Q))
        ymax <- max(3, 1.2 * quantile(nlq, probs=0.999, na.rm=TRUE)[1]) ## y-axis        
        
        ##------------- layout ----------------
        nplots <- length(mx.list)
        par(mfrow=c(1,1), mar=c(1,1,1,1)*0.2, mgp=c(2.6,1,0), oma=c(1,1,0,0)*2)
        if(nplots>24) {
            nc = max(ceiling(nplots/3),6)
            par(mfrow=c(3,nc))
        } else if(FALSE && nplots <= 3) {
            nc = 3
            par(mfrow=c(1,nc))
        } else {
            nc = max(ceiling(nplots/2),6)
            par(mfrow=c(2,nc))
        }        
        
        shiny::withProgress(message="computing volcano plots ...", value=0, {
            i=1
            for(i in 1:nplots) {

                mx <- mx.list[[i]]
                ## qv <- Q[[i]]
                is.sig <- (mx[,"qv"] <= fdr & abs(mx[,"fc"]) >= lfc)
                table(is.sig)
                sig.gs = rownames(mx)[which(is.sig)]
                if(!is.null(sel.gsets)) sig.gs <- intersect(sel.gsets, sig.gs)

                gx.volcanoPlot.XY(
                    x = mx[,"fc"], pv = mx[,"qv"],
                    use.fdr=TRUE, p.sig=fdr, lfc=lfc,                
                    gene = rownames(mx), 
                    xlab = "effect size (NES)", lab.cex=0, nlab=0,
                    render="canvas", n=1000, highlight=sig.gs,
                    cex=1, cex.axis=1.3, cex.main=1.4, axes=FALSE,
                    ylim=c(0,ymax), main="" )

                ## draw axis if first column or last row
                graphics::box(lwd=1, col="black", lty="solid")
                is.first = (i%%nc==1)
                last.row = ( (i-1)%/%nc == (nplots-1)%/%nc )
                if(is.first) axis(2, mgp=c(2,0.7,0), cex.axis=0.8)
                if(last.row) axis(1, mgp=c(2,0.7,0), cex.axis=0.8)
                legend("top", legend=names(mx.list)[i], box.lty=0,
                       x.intersp = 0.3, y.intersp = 0.5,
                       inset=c(0,0.01), cex=1.2, bg="white")
                shiny::incProgress( 1.0/nplots )
            }
            
        })
    })
    
    volcanoAll_text = "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for gene sets across all contrasts. This provides users an overview of the statistics across all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."

    shiny::callModule(
        plotModule,
        id = "volcanoAll",
        func = volcanoAll.RENDER,
        func2 = volcanoAll.RENDER,
        height = c(imgH,450), width = c("auto",1500), res=c(72,85),
        pdf.width=15, pdf.height=5, 
        title="Volcano plots for all contrasts",
        info.text = volcanoAll_text,
        ##caption = volcanoAll_caption,
        add.watermark = WATERMARK
    )
    

    ##================================================================================
    ## Volcano (methods)
    ##================================================================================

    volcanoMethods.RENDER <- shiny::reactive({
        ##renderPlotly({
        ngs <- inputData()    
        shiny::req(ngs, input$gs_features)
        ##if(is.null(input$gs_features)) return(NULL)
        
        cmp = 1
        cmp = input$gs_contrast
        mx = ngs$gset.meta$meta[[cmp]]
        fx = unclass(mx$fc)
        qv = unclass(mx$q)
        pv = unclass(mx$p)    

        fx[which(is.na(fx))] <- NA
        fx[which(is.infinite(fx))] <- NA
        qv[which(is.na(qv))] <- 1

        fdr=1;lfc=0
        fdr=0.05;lfc=1
        fdr = as.numeric(input$gs_fdr)    
        lfc = as.numeric(input$gs_lfc)
        sel.gsets <- rownames(mx)
        sel.gsets = COLLECTIONS[[1]]
        sel.gsets = COLLECTIONS[[input$gs_features]]
        
        nlq <- -log10(1e-99 + unlist(qv))
        ymax <- max(3, 1.2 * quantile(nlq, probs=0.999, na.rm=TRUE)[1]) ## y-axis        
        ##ymax <- 1.2 * max(-log10(1e-99 + qv), na.rm=TRUE)

        nplots = ncol(fx)
        nc = max(nplots/2,5)
        nn = c(2, nc)
        par(mfrow=nn, mar=c(1,1,1,1)*0.2, mgp=c(2.6,1,0), oma=c(1,1,0,0)*2)
        
        shiny::withProgress(message="computing volcano plots ...", value=0, {
            i=1
            for(i in 1:nplots) {
                
                is.sig <- ( qv[,i] <= fdr & abs(fx[,i]) >= lfc)
                sig.gs = rownames(mx)[which(is.sig)]
                sig.gs <- intersect(sel.gsets, sig.gs)
                
                method = colnames(fx)[i]
                gx.volcanoPlot.XY(
                    x = fx[,i], pv = qv[,i],
                    use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                    ##gene = substring(rownames(mx),1,35),
                    gene = rownames(mx),
                    xlab = "effect size (NES)", ylim=c(0,ymax), 
                    lab.cex=0, nlab=0, axes=FALSE, 
                    render="canvas", n=1000, highlight=sig.gs,
                    cex=1, cex.axis=1.3, main="")
                
                ##title(mt, line=-1.5, cex.main=1.4)
                graphics::box(lwd=1, col="black", lty="solid")                

                ##volcano_plot(limma, render="plotly", n=1000, cex=1, highlight=genes)
                ## draw axis if first column or last row
                is.first = (i%%nc==1)
                last.row = ( (i-1)%/%nc == (nplots-1)%/%nc )
                if(is.first) axis(2, mgp=c(2,0.7,0), cex.axis=0.8)
                if(last.row) axis(1, mgp=c(2,0.7,0), cex.axis=0.8)
                legend("top", legend=method, box.lty=0,
                       x.intersp = 0.3, y.intersp = 0.5,
                       inset=c(0,0.01), 
                       cex=1.2, bg="white")
                
                shiny::incProgress( 1/nplots )                
            }
        })


        
    })

    volcanoMethods_text = "The <strong>Volcano (methods)</strong> panel displays the volcano plots provided by different enrichment calculation methods. This provides users an quick overview of the sensitivity of the statistical methods at once. Methods showing better statistical significance will show volcano plots with 'higher' wings."

    shiny::callModule(
        plotModule,
        id="volcanoMethods",
        func = volcanoMethods.RENDER,
        func2 = volcanoMethods.RENDER,
        height = c(imgH,450), width = c('auto',1600), res=c(75,90),
        pdf.width=15, pdf.height=5, 
        title="Volcano plots for all methods",
        info.text = volcanoMethods_text,
        ##caption = volcanoMethods_caption,
        add.watermark = WATERMARK
    )
    

    ##================================================================================
    ## Enrichment table
    ##================================================================================
    
    gset_selected <- shiny::reactive({
        ##i = as.integer(input$gseatable_rows_selected)
        i = as.integer(gseatable$rows_selected())
        if(is.null(i) || length(i)==0) return(NULL)
        rpt = getFilteredGeneSetTable()
        gs = rownames(rpt)[i]
        return(gs)
    })

    geneDetails <- shiny::reactive({
        ## return details of the genes in the selected gene set
        ##
        
        ngs <- inputData()
        shiny::req(ngs,input$gs_contrast)
        gs=1;comp=1

        comp = input$gs_contrast
        gs = gset_selected()
        if(is.null(gs) || length(gs)==0) return(NULL)
        
        mx <- ngs$gx.meta$meta[[comp]]
        is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]",rownames(mx)))
        is.multiomics
        if(is.multiomics) {
            ii <- grep("\\[gx\\]|\\[mrna\\]",rownames(mx))
            mx <- mx[ii,]
            ##rownames(mx) <- sub(".*:|.*\\]","",rownames(mx))
        }            
        
        ##gxmethods <- c("trend.limma","ttest.welch")
        gxmethods <- selected_gxmethods() ## from module-expression
        shiny::req(gxmethods)
        ##limma1 = sapply(mx[,c("fc","p","q")], function(x) x[,"trend.limma"])
        ##limma1.fc <- rowMeans(mx$fc[,gxmethods,drop=FALSE],na.rm=TRUE)
        limma1.fc <- mx$meta.fx
        limma1.pq = sapply(mx[,c("p","q")], function(x) {
            apply(x[,gxmethods,drop=FALSE],1,max,na.rm=TRUE)
        })
        limma1 <- cbind( fc=limma1.fc, limma1.pq)
        ##limma  = cbind( ngs$gx.meta$meta[[comp]][,c("gene_name","gene_title")], limma1)
        rownames(limma1) <- rownames(mx)

        ## filter on significance?????
        if(FALSE && !input$gs_showall) {
            lfc=1;fdr=0.05
            lfc <- as.numeric(input$gs_lfc)
            fdr <- as.numeric(input$gs_fdr)
            is.sig <- abs(limma1[,"fc"]) >= lfc & limma1[,"q"] <= fdr
            table(is.sig)
            limma1 <- limma1[is.sig,,drop=FALSE]
        }
        
        ## in multi-mode we select *common* genes
        ns <- length(gs)
        gmt1 <- ngs$GMT[,gs,drop=FALSE]
        genes = rownames(gmt1)[which(Matrix::rowSums(gmt1!=0)==ns)]
        genes = intersect(genes, ngs$genes[rownames(limma1),"gene_name"])
        genes = setdiff(genes, c("",NA,"NA"," "))

        title <- rep(NA,length(genes))
        title = as.character(GENE.TITLE[genes])
        title[is.na(title)] <- " "
        
        rpt <- data.frame("gene_name"=genes, "gene_title"=as.character(title) )
        genes = rpt[,"gene_name"]
        genes1 <- ngs$genes[rownames(limma1),"gene_name"]
        limma1 = limma1[match(genes, genes1),,drop=FALSE ]  ## align limma1  
        ##avg.rho <- rowMeans(cor(t(ngs$X[rownames(limma1),,drop=FALSE]),
        ##                        t(ngs$gsetX[gs,,drop=FALSE])))
        ##rpt = cbind(rpt, limma1, gset.rho=avg.rho)
        rpt = cbind(rpt, limma1)
        rpt = rpt[which(!is.na(rpt$fc) & !is.na(rownames(rpt))),,drop=FALSE]
        ##rpt = data.frame(rpt, check.names=FALSE)

        if(nrow(rpt)>0) {
            rpt = rpt[order(-abs(rpt$fc)),,drop=FALSE]
        }
        return(rpt)
    })

    gene_selected <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)
        i = 1
        ##i = as.integer(input$genetable_rows_selected)
        i = as.integer(genetable$rows_selected())
        if(is.null(i) || is.na(i) || length(i)==0) i=1
        rpt <- geneDetails()
        if(is.null(rpt) || nrow(rpt)==0) {
            return(list(gene=NA, probe=NA))
        }
        sel.gene = rownames(rpt)[i]
        gene = as.character(rpt$gene_name[i])
        probe = rownames(ngs$genes)[match(gene, ngs$genes$gene_name)]
        return(list(gene=gene, probe=probe))
    })


    gseatable.RENDER <- shiny::reactive({

        ##rpt = getFullGeneSetTable()
        rpt = getFilteredGeneSetTable()
        if(is.null(rpt)) return(NULL)
        if(nrow(rpt)==0) return(NULL)
        ##cat("dim.rpt=",dim(rpt),"\n")
        
        if(!("GS" %in% colnames(rpt))) rpt = cbind(GS=rownames(rpt),rpt)
        if("GS" %in% colnames(rpt)) rpt$GS = shortstring(rpt$GS,72)
        if("size" %in% colnames(rpt)) rpt$size = as.integer(rpt$size)

        fx = NULL
        fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(rpt))[1]
        if(length(fx.col)>0) fx = rpt[,fx.col]
        
        jj = which(sapply(rpt,is.numeric))
        if(length(jj)>0) rpt[,jj] = round(rpt[,jj],digits=4)    
        jj = which( sapply(rpt,is.character) |  sapply(rpt,is.factor) )
        if(length(jj)>0) rpt[,jj] = apply(rpt[,jj,drop=FALSE],2,shortstring,100)

        if(!input$gs_showqvalues) {
            rpt <- rpt[,grep("^q[.]|^q$",colnames(rpt),invert=TRUE)]
        }
        
        ## wrap genesets names with known links.
        rpt$GS <- wrapHyperLink(rpt$GS, rownames(rpt))    
        selectmode = "single"
        selectmode

        is.numcol <- sapply(rpt, is.numeric)
        numcols <- which( is.numcol & !colnames(rpt) %in% c("size"))
        numcols <- colnames(rpt)[numcols]

        colnames(rpt) <- sub("GS","geneset",colnames(rpt))

        ##rpt = format(rpt, digits=4)
        DT::datatable(rpt,
                      class = 'compact cell-border stripe hover',
                      rownames=FALSE,
                      escape = c(-1,-5),
                      ##extensions = c('Buttons','Scroller'),
                      extensions = c('Scroller'),                  
                      fillContainer = TRUE,
                      selection = list(mode=selectmode, target='row', selected=NULL),
                      ## options=list(
                      ##     dom = 'lfrtip',
                      ##     ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                      ##     scrollX = TRUE,
                      ##     scrollY = tabH,
                      ##     scroller=TRUE,
                      ##     deferRender=TRUE
                      ## )
                      options=list(
                          dom = 'frtip',                          
                          paging = TRUE,
                          pageLength = 15, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = FALSE,                          
                          scroller = FALSE,
                          deferRender=TRUE,
                          search = list(
                              regex = TRUE,
                              caseInsensitive = TRUE
                            ##, search = 'HALLMARK'
                          )
                      )  ## end of options.list 
                      ) %>%
            DT::formatSignif(numcols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
            DT::formatStyle(fx.col, 
                            background = color_from_middle( fx, 'lightblue', '#f5aeae'))
    })

    genetable.RENDER <- shiny::reactive({

        rpt <- geneDetails()    
        ## if(is.null(rpt)) return(NULL)
        if(is.null(rpt) || nrow(rpt)==0) {
            shiny::validate(shiny::need(nrow(rpt) > 0, "warning. no genes."))
            return(NULL)
        }
        
        rpt$gene_title <- NULL    
        if(!is.null(rpt) && nrow(rpt)>0 ) {
            jj = which(sapply(rpt,is.numeric))
            rpt[,jj] = round(rpt[,jj],digits=4)
            jj = which( sapply(rpt,is.character) |  sapply(rpt,is.factor) )
            if(length(jj)>0) rpt[,jj] = apply(rpt[,jj,drop=FALSE],2,shortstring,60)
        } else {
            rpt <- data.frame('',0,0,0)[0,]
            colnames(rpt) <- c("gene_name","fc","p","q")
        }

        colnames(rpt) <- sub("^GS$","gene set",colnames(rpt))
        numeric.cols <- which(sapply(rpt, is.numeric))
        numeric.cols

        tbl <- DT::datatable(rpt,
                             class = 'compact cell-border stripe', rownames=FALSE,
                             extensions = c('Scroller'),
                             selection = list(mode="single", target='row', selected=1),
                             fillContainer = TRUE,
                             ## options=list(
                             ##     dom = 'lfrtip',
                             ##     ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                             ##     scrollX = TRUE,
                             ##     scrollY = tabH,
                             ##     scroller=TRUE,
                             ##     deferRender=TRUE
                             ## )
                             options=list(
                                 dom = 'frtip',                          
                                 paging = TRUE,
                                 pageLength = 15, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                                 scrollX = TRUE,
                                 scrollY = FALSE,                          
                                 scroller = FALSE,
                                 deferRender=TRUE,
                                 search = list(
                                     regex = TRUE,
                                     caseInsensitive = TRUE
                                     ##, search = 'HALLMARK'
                                 )
                      )  ## end of options.list 
                             ) %>%
            DT::formatSignif(numeric.cols,4)
        
        if( nrow(rpt)>0 && ("fc" %in% colnames(rpt)) ) {
            fx = rpt[,"fc"]
            tbl <- tbl %>%
                DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle("fc", background = color_from_middle( fx, 'lightblue', '#f5aeae'))
        }
        tbl
    })

    gseatable_text = paste("Similar to the differential gene expression analysis, users can perform differential expression analysis on a geneset level that is referred as gene set enrichment analysis. To ensure statistical reliability, the platform performs the gene set enrichment analysis using multiple methods, including",a_Spearman,", ",a_GSVA,", ",a_ssGSEA,", ",a_Fisher,", ",a_GSEA,", ",a_camera," and ",a_fry,".<br><br>The combined result from the methods is displayed in this table, where for each geneset the <code>meta.q</code> corresponds to the highest <code>q</code> value provided by the methods and the number of <code>stars</code> indicate how many methods identified the geneset as significant (<code>q < 0.05</code>). The table is interactive; users can sort it by <code>logFC</code>, <code>meta.q</code> and <code>starts</code>. Additionally, the list of genes in that geneset are displayed in the second table on the right. Users can filter top N = {10} differently enriched gene sets in the table by clicking the <code>top 10 gene sets</code> from the table <i>Settings</i>.")

    gseatable_opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns('gs_top10'),'top 10 gene sets',FALSE),
               "Display only top 10 differentially enirhced gene sets (positively and negatively) in the <b>enrihcment analysis</b> table.", placement="top", options = list(container = "body")),
        withTooltip(shiny::checkboxInput(ns('gs_showqvalues'),'show indivivual q-values',FALSE),
               "Show all q-values of each individual statistical method in the table.", 
               placement="top", options = list(container = "body"))    
    )

    gseatable <- shiny::callModule(
        tableModule, 
        id="gseatable",
        func = gseatable.RENDER,
        info.text = gseatable_text, label="I",
        options = gseatable_opts,
        title="Enrichment analysis",
        info.width="500px",
        height = c(285, 700)
    )

    genetable_text = "By clicking on a gene set in the table <code>I</code>, it is possible to see the gene list of that gene set in this table. By clicking on a gene in this table, users can check the expression status of the gene for the selected contrast in the <code>Expression</code> barplot and its correlation to the gene set in the <code>Gene to gene set correlation</code> scatter plot under the <code>Plots</code> section."
    
    genetable <- shiny::callModule(
        tableModule,
        id = "genetable",
        func=genetable.RENDER,
        info.text = genetable_text,
        title="Genes in gene set", label="II",
        height = c(285,700), width = c('auto',800)
    )
    

    ##================================================================================
    ## Enrichment (all)
    ##================================================================================

    fctable.RENDER <- shiny::reactive({
        
        ngs <- inputData()

        ## get all contrasts
        F <- sapply( ngs$gset.meta$meta, function(x) x[,"meta.fx"])
        colnames(F) <- gsub("_"," ",colnames(F))
        rownames(F) <- rownames(ngs$gset.meta$meta[[1]])
        fc.var <- round( rowMeans(F**2,na.rm=TRUE), digits=3)
        gs <- substring(rownames(F),1,60)
        F1 <- data.frame( geneset=gs, fc.var=fc.var, round(F,digits=3), check.names=FALSE)

        ## get current filtered geneset and extract names of gene sets
        ##rpt = getFullGeneSetTable()
        rpt = getFilteredGeneSetTable()
        F1 <- F1[intersect(rownames(rpt),rownames(F1)),,drop=FALSE]    
        F1$geneset <- wrapHyperLink(F1$geneset, rownames(F1))
        
        DT::datatable( F1, rownames=FALSE, escape=-1,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'frtip', 
                          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE,
                          deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle( "fc.var",
                                ##background = DT::styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(fc.var, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  %>%
                DT::formatStyle( colnames(F),
                                ##background = DT::styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(F[,], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    gx_fctable_text = "The <strong>Enrichment (all)</strong> panel reports the gene set enrichment for all contrasts in the selected dataset."
    
    gx_fctable_caption = "<b>Enrichment for all contrasts.</b> Table summarizing the enrichment for all gene sets across all contrasts. The column `fc.var` corresponds to the variance of the gene set across all contrasts."

    shiny::callModule(
        tableModule,
        id = "fctable",
        func = fctable.RENDER,
        title ="Gene set enrichment for all contrasts",
        info.text = gx_fctable_text,
        caption = gx_fctable_caption,
        height = c(295,750),
        width = c('100%',1600)
    )

    
    ##================================================================================
    ## FDR table
    ##================================================================================


    
    FDRtable.RENDER <- shiny::reactive({
        
        ngs <- inputData()    
        shiny::req(ngs, input$gs_statmethod)
        
        meta <- ngs$gset.meta
        test = GSET.DEFAULTMETHODS
        test <- input$gs_statmethod
        ##if(is.null(test)) return(NULL)

        if(length(test)==1) {
            sig.up   = meta$sig.counts[[test]][["up"]]
            sig.down = meta$sig.counts[[test]][["down"]]
            rownames(sig.up) = paste0(rownames(sig.up),"::",test[1])
            rownames(sig.down) = paste0(rownames(sig.down),"::",test[1])
        } else {
            sig.up = c()
            sig.down = c()
            for(i in 1:length(test)) {
                sig1 = meta$sig.counts[[test[i]]][["up"]]
                sig2 = meta$sig.counts[[test[i]]][["down"]]
                rownames(sig1) = paste0(rownames(sig1),"::",test[i])
                rownames(sig2) = paste0(rownames(sig2),"::",test[i])
                sig.up <- rbind(sig.up, sig1)
                sig.down <- rbind(sig.down, sig2)
            }
        }
        sig.up <- sig.up[order(rownames(sig.up)),,drop=FALSE]
        sig.down <- sig.down[order(rownames(sig.down)),,drop=FALSE]    
        pvals = sort( c(1e-16, 10**seq(-8,-2,2), 0.05, 0.1, 0.2, 0.5,1))
        kk = intersect(colnames(sig.up),pvals)
        sig.up = sig.up[,match(kk,colnames(sig.up)),drop=FALSE]
        sig.down = sig.down[,match(kk,colnames(sig.down)),drop=FALSE]
        
        colnames(sig.up)[1] = paste("UP   FDR = ",colnames(sig.up)[1])
        colnames(sig.down)[1] = paste("DOWN   FDR = ",colnames(sig.down)[1])
        colnames(sig.down) = paste0("  ",colnames(sig.down))
        sigcount = cbind( sig.down, sig.up[rownames(sig.down),] )
        dim(sigcount)    
        maxsig = 0.99 * max(sigcount,na.rm=TRUE)
        ##gs.up %>% kableExtra::kable("html") %>%
        ##    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
        ##                  font_size = 10)

        contr = sub("::.*","",rownames(sigcount))
        ##contr = rownames(sigcount)
        metd  = sub(".*::","",rownames(sigcount))
        D = data.frame( method=metd, contrast=contr, sigcount, check.names=FALSE)

        ##width  <- session$clientData$output_kegg_graph_width
        ##height <- session$clientData$output_kegg_graph_height    

        DT::datatable( D, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'frtip',
                          pageLength = 999,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE,
                          deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle(colnames(sig.up),
                                background = DT::styleColorBar(c(0,maxsig), '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  %>% 
                DT::formatStyle(colnames(sig.down),
                                background = DT::styleColorBar(c(0,maxsig), 'lightblue'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  
    })

    FDRtable_text = "The <strong>FDR table</strong> panel reports the number of significant gene sets at different FDR thresholds, for all contrasts and all methods. Using the table the user can determine which statistical methods perform better for a particular contrast."

    FDRtable_caption = "<b>FDR table.</b> Number of significant gene sets versus different FDR thresholds, for all contrasts and all methods. The blue color denote the number of downregulated genes, the red color for upregulated genes."
    
    shiny::callModule(
        tableModule,
        id = "FDRtable",
        func = FDRtable.RENDER,
        title = 'Number of significant gene sets',
        info.text = FDRtable_text,
        caption = FDRtable_caption,
        height = c(295,750),
        width = c('100%',1600)
    )

    ## reactive values to return to parent environment
    outx <- list(selected_gsetmethods=selected_gsetmethods)
    return(outx)

  }) ## end of moduleServer    
} ## end-of-Board
