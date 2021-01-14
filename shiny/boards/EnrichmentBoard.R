##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing EnrichmentBoard")

EnrichmentInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

EnrichmentUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        flex = c(1.4,1),
        height = 730,
        tabsetPanel(
            id = ns("tabs1"),
            tabPanel("Top enriched",uiOutput(ns("topEnriched_UI"))),
            tabPanel("Plots",uiOutput(ns("subplots_UI"))),
            tabPanel("Compare",uiOutput(ns("compare_UI"))),
            tabPanel("Volcano (all)",uiOutput(ns("volcanoAll_UI"))),
            tabPanel("Volcano (methods)",uiOutput(ns("volcanoMethods_UI"))),
            tabPanel("GeneMap",uiOutput(ns("genemap_UI")))
        ),
        tabsetPanel(
            id = ns("tabs2"),
            tabPanel("Table",uiOutput(ns("tables_UI"))),
            tabPanel("Foldchange (all)",uiOutput(ns("fctable_UI"))),
            tabPanel("FDR table",uiOutput(ns("FDRtable_UI")))                       
        )
    )
}


EnrichmentBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    ## reactive functions from shared environment
    inputData <- env[["load"]][["inputData"]]
    selected_gxmethods <- env[["expr"]][["selected_gxmethods"]]

    fullH = 730
    rowH = 345  ## row height of panels
    imgH = 265  ## height of images
    tabH = 160  ## height of tables
    tabH = "70vh" ## height of tables

    description = "<b>Geneset enrichment analysis.</b> Perform differential expression analysis on a geneset level, also called geneset enrichment analysis."
    output$description <- renderUI(HTML(description))
    
    gs_infotext = paste("Similar to the differential gene expression analysis, users can perform differential expression analysis on a geneset level in this page, which is also referred as gene set enrichment (GSE) analysis. The platform has more than 50.000 genesets (or pathways) in total that are divided into 30 geneset collections such as ",a_Hallmark,", ",a_MSigDB,", ",a_KEGG," and ",a_GO,". Users have to specify which comparison they want to visually analyze employing a certain geneset collection. 

<br><br>To ensure the statistical reliability, the platform performs Enrichment Analyses using multiple methods. The result from the statistical methods is displayed in <strong>Enrichment table</strong> panel. In the <strong>Top enriched</strong> panel, the top 10 differentially enriched geneses (pathways) are displayed. In the <strong>Plots</strong> panel, a volcano plot of genes contained in the selected geneset and a barplot of expressions per sample group are displayed. In the <strong>Compare</strong> panel, users can compare the differential expression status of that geneset for all other comparisons. Finally, volcano plots of genesets for all comparisons are displayed under the <strong>Volcano (all) </strong> tab. This allows users to have an overall picture across comparisons at the same time.

<br><br>EXPERT MODE ONLY: To compare the different statistical methods, the <strong>Volcano (methods)</strong> panel shows volcano plots of all methods. The <strong>FDR table</strong> panel reports the number of significant gene sets at different FDR thresholds for all contrasts.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=4' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
")
    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================
    
    FDR.VALUES2 = c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1)
    gs_testmethod_text = paste("Select a method for the statistical test. To ensure the statistical reliability, the platform performs the gene set enrichment analysis using multiple methods, including",a_Spearman,", ",a_GSVA,", ",a_ssGSEA,", ",a_Fisher,", ",a_GSEA,", ",a_camera," and ",a_fry,".")  ## does not work (tipify cannot handle HTML?)
    gs_testmethod_text = paste("Select a method or multiple methos for the statistical test.")
    
    GSET.DEFAULTMETHODS = c("fisher","gsva","camera","fgsea")
    GSET.DEFAULTMETHODS = c("fisher","gsva","camera")
    GSET.DEFAULTMETHODS = c("gsva","camera","fgsea","fisher")
    
    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("gs_info"), "Tutorial", icon = icon("youtube")),
                   "Show more information about this module."),
            hr(), br(),             
            tipify( selectInput(ns("gs_contrast"),"Contrast:", choices=NULL),
                   "Select a contrast of interest for the analysis.", placement="top"),
            tipify( selectInput(ns("gs_features"),"Gene set collection:", choices=NULL, multiple=FALSE),
                   "Choose a specific gene set collection for the analysis.", placement="top"),
            fillRow( flex=c(1,1),
                    tipify( selectInput(ns("gs_fdr"),"FDR", choices=FDR.VALUES2, selected=1),
                           "Set the false discovery rate (FDR) threshold.", placement="top"),
                    tipify( selectInput(ns("gs_lfc"),"logFC threshold", choices=c(0,0.2,0.5,1,2,5),
                                        selected=0),
                           "Set the logarithmic fold change (logFC) threshold.", placement="top")
                    ),
            br(),br(),br(),br(),
            tipify(actionLink(ns("gs_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            br(),br(),
            conditionalPanel(
                "input.gs_options % 2 == 1", ns=ns, 
                tagList(
                    tipify(checkboxGroupInput(ns('gs_method'),'Statistical methods:', choices=NULL),
                           gs_testmethod_text, placement="right", options = list(container="body"))
                )
            )
        )
        if(DEV) {
            uix <- tagList(
                hr(),h6("Developer options:"),
                tipify( radioButtons(ns('gs_lfcmethod'),'Score method: (dev)', choices=c("fc","gs")),
                       "Select reporting method for differential expression of gene sets. FC=average foldchange of genes in genest. GS= differential score derived from meta-score from single-sample geneset methods.", placement="right", options = list(container="body"))
            )
            ui <- c(ui, uix)
        }
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$gs_info, {
        showModal(modalDialog(
            title = HTML("<strong>Enrichment Analysis Board</strong>"),
            HTML(gs_infotext),
            easyClose = TRUE, size="l" ))
    })

    observe({
        ngs <- inputData()
        req(ngs)
        meta <- ngs$gset.meta$meta
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons = sort(intersect(comparisons, names(meta)))
        updateSelectInput(session, "gs_contrast", choices=comparisons)

        ## get the computed geneset methods
        gset.methods = sort(colnames(meta[[1]]$fc))
        sel2 = c(intersect(GSET.DEFAULTMETHODS,gset.methods),gset.methods)
        sel2 = head(unique(sel2),3)

        updateCheckboxGroupInput(session, 'gs_method',
                                 choices = sort(gset.methods),
                                 selected = sel2)
        
    })
    
    observe({
        ngs <- inputData()
        req(ngs)
        nn <- sapply(COLLECTIONS, function(k) sum(k %in% rownames(ngs$gsetX)))
        gsets.groups <- names(COLLECTIONS)[which(nn>=5)]
        gsets.groups <- c("<all>",sort(gsets.groups))
        sel = "<all>"
        hmark <- grep("^H$|hallmark|",gsets.groups,ignore.case=TRUE,value=TRUE)
        if(length(hmark)>0) sel <- hmark[1]
        updateSelectInput(session, "gs_features",choices=gsets.groups, selected=sel)
        
    })

    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================

    ## selected_gxmethods <- reactive({
    ##     sel <- SEL.GXMETHODS()
    ##     req(sel)
    ##     sel
    ## })
    
    selected_gsetmethods <- reactive({
        ngs <- inputData()
        req(ngs)
        gset.methods0 = colnames(ngs$gset.meta$meta[[1]]$fc)
        ##test = head(intersect(GSET.DEFAULTMETHODS,gset.methods0),3) ## maximum three
        test = input$gs_method
        test = intersect(test,gset.methods0) ## maximum three
        test
    })

    calculateMeta <- function(comparison, methods, ngs) {
        ##ngs <- inputData()
        mx = ngs$gset.meta$meta[[comparison]]
        if(is.null(mx)) return(NULL)
        mx.methods = colnames(unclass(mx$fc))
        mx.methods
        methods = intersect(methods, mx.methods)
        if(is.null(methods) || length(methods)==0) {
            cat("ERROR: calculateMeta:: no valid methods\n")
            return(NULL)
        }
        
        ## recalculate meta values
        pv = unclass(mx$p)[,methods,drop=FALSE]
        qv = unclass(mx$q)[,methods,drop=FALSE]
        fc = unclass(mx$fc)[,methods,drop=FALSE]
        pv[is.na(pv)] = 0.999
        qv[is.na(qv)] = 0.999
        fc[is.na(fc)] = 0
        score = fc * (-log10(qv))
        dim(pv)
        if(NCOL(pv)>1) {
            ss.rank <- function(x) scale(sign(x)*rank(abs(x)),center=FALSE)
            fc = rowMeans(scale(fc,center=FALSE),na.rm=TRUE)  ## REALLY???
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
    
    getGeneSetTable <- reactive({
        ngs <- inputData()
        req(ngs)
        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        if(!(comp %in% names(ngs$gset.meta$meta))) return(NULL)
        mx = ngs$gset.meta$meta[[comp]]
        dim(mx)
        
        outputs = NULL
        gsmethod = gsmethod0 = colnames(unclass(mx$fc))
        gsmethod <- input$gs_method   
        if(is.null(gsmethod) || length(gsmethod)==0) return(NULL)
        
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
            pv[is.na(pv)] = 0.999
            qv[is.na(qv)] = 0.999
            fx[is.na(fx)] = 0

            stars.symbols = sapply(1:20,function(i) paste(rep("\u2605",i),collapse=""))
            stars = c("",stars.symbols)[1+rowSums(qv < 0.05)]                
            names(stars) <- rownames(mx)
            
            ##------------ calculate META parameters ----------------
            meta <- calculateMeta(comp, gsmethod, ngs=ngs)
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
            if(DEV) gsdiff.method <- input$gs_lfcmethod
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
                
                G <- t(ngs$GMT[pp,jj] != 0)
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
                             stars = stars[gs],
                             AveExpr0=AveExpr0[gs],
                             AveExpr1=AveExpr1[gs])
            
            ## add extra p/q value columns
            jj <- match(gs, rownames(mx))
            rpt <- cbind( rpt, q=qv[jj,])
            
            ##rownames(rpt) = gs
        }  else {
            ## show original table (single method)
            rpt = outputs[[gsmethod]]
        }    
        rpt = data.frame(rpt)
        
        ## Get the meta q-value column, filter on q-value
        qv.col = grep("meta.q|fdr|adj.p.val|q.value|adjusted.p|padj|qv",
                      colnames(rpt),ignore.case=TRUE)[1]
        if(length(qv.col)==0) stop("could not parse q-value column")
        qv = rpt[,qv.col]
        fdr = 1
        fdr <- input$gs_fdr
        rpt = rpt[ which( qv <= as.numeric(fdr) ), ,drop=FALSE ]  ## important to cast numeric!

        ## filter on fold-change    
        fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(rpt))[1]
        fx  = as.numeric(rpt[,fx.col])
        names(fx) = rownames(rpt)
        lfc = as.numeric(input$gs_lfc)
        rpt = rpt[which(abs(fx) >= lfc),,drop=FALSE]
        fx  = rpt[,fx.col]
        rpt = rpt[order(-abs(fx)),]
        
        ## just show top 10
        if(length(input$gs_top10) && input$gs_top10==TRUE) {
            fx  = as.numeric(rpt[,fx.col])
            names(fx) = rownames(rpt)
            pp <- unique(c(head(names(sort(-fx[which(fx>0)])),10),
                           head(names(sort(fx[which(fx<0)])),10)))
            rpt = rpt[pp,,drop=FALSE]
            fx  = as.numeric(rpt[,fx.col])
            rpt = rpt[order(-fx),]
        }

        ## limit to 1000 rows???
        ## rpt <- head(rpt, 1000)
        rpt <- data.frame(rpt)
        
        return(rpt)
    })


    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================
    
    ## Top enriched    
    topEnriched.RENDER %<a-% reactive({

        ngs <- inputData()
        alertDataLoaded(session,ngs)
        rpt <- getGeneSetTable()
        ##if(is.null(rpt)) return(NULL)
        
        req(ngs, rpt, input$gs_contrast)

        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        if(!(comp %in% names(ngs$gx.meta$meta))) return(NULL)
        
        gx.meta <- ngs$gx.meta$meta[[comp]]
        ##rnk0 <- gx.meta[,"fc"][,"trend.limma"]
        ##names(rnk0) = gx.meta[,"gene_name"]        
        rnk0 <- gx.meta$meta.fx
        names(rnk0) = ngs$genes[rownames(gx.meta),"gene_name"]
        rnk0 = rnk0 - mean(rnk0,na.rm=TRUE)  ## scaling/centering should be done in calculation...
        
        ## filter on active rows (using search)
        ##ii <- input$gseatable_rows_all
        ii <- gseatable$rows_all()
        rpt <- rpt[ii,,drop=FALSE]
        if(nrow(rpt)==0) return(NULL)
        
        fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(rpt))[1]
        qv.col = grep("meta.q|q$",colnames(rpt))[1]
        fx = rpt[,fx.col]
        qv = rpt[,qv.col]
        names(qv) <- names(fx) <- rownames(rpt)

        ##top.up <- names(sort(fx[which(fx>0)],decreasing=TRUE))
        ##top.dn <- names(sort(fx[which(fx<0)]))
        top <- rownames(rpt)
        
        par(mfrow=c(2,5), mar=c(0.5,2.5,2.8,0), mgp=c(2,0.8,0))
        for(i in 1:10) {
            if(i > length(top)) {
                frame()
            } else {
                gs <- top[i]
                gs1 = breakstring(gs,28,50,force=FALSE)
                genes = toupper(names(which(ngs$GMT[,gs]!=0)))
                names(rnk0) <- toupper(names(rnk0))
                ylab = ""
                ## if(i %in% c(1,6)) ylab = "Ranked list metric"
                if(i%%5 == 1) ylab = "Ranked list metric"
                gsea.enplot(rnk0, genes, names=NULL, ##main=gs,
                            main=gs1, xlab="", ylab=ylab,
                            cex.main=0.78, len.main=80)
                qv1 = formatC(qv[gs],format="e", digits=2)
                legend("topright", paste("q=",qv1), bty="n",cex=0.85)
            }
        }
    })

    topEnriched_text = "The <strong>Top enriched</strong> section shows the enrichment plots for the top differentially (both positively and negatively) enriched gene sets for the selected comparison in the <code>Contrast</code> settings. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group."

    topEnriched_caption = "<b>Top enriched gene sets.</b> Enrichment plots of the top differentially enriched gene sets. Black vertical bars indicate the rank of genes in the gene set in the sorted list metric. The green curve corresponds to the 'running statistics' of the enrichment score."

    callModule(
        plotModule,
        id = "topEnriched", label="a",
        func = topEnriched.RENDER,
        func2 = topEnriched.RENDER,
        info.text = topEnriched_text,
        height = c(imgH,450), width = c('auto',1500), res=95,
        pdf.width = 14, pdf.height = 4, 
        title = "Top enriched gene sets"
        ##caption = topEnriched_caption
    )


    topEnrichedFreq.RENDER %<a-% reactive({

        ngs <- inputData()

        rpt <- getGeneSetTable()
        ##if(is.null(rpt)) return(NULL)
        req(ngs, rpt, input$gs_contrast)

        comp=1
        comp = input$gs_contrast
        if(is.null(comp)) return(NULL)
        if(!(comp %in% names(ngs$gx.meta$meta))) return(NULL)
        
        ## filter on active rows (using search)
        ##ii <- input$gseatable_rows_all
        ii <- gseatable$rows_all()
        rpt <- rpt[ii,,drop=FALSE]
        if(nrow(rpt)==0) return(NULL)
        
        fx.col = grep("score|fx|fc|sign|NES|logFC",colnames(rpt))[1]
        fx = rpt[,fx.col]
        names(fx) <- rownames(rpt)

        top <- rownames(rpt)        
        ntop <- as.integer(input$gs_enrichfreq_ntop)
        top <- head(top,ntop)
        if(!all(top %in% colnames(ngs$GMT))) return(NULL)
        
        F <- 1*(ngs$GMT[,top]>0)
        wt = FALSE
        if(input$gs_enrichfreq_gsetweight) {
            F <- t(t(F)  / colSums(F,na.rm=TRUE))
            wt = TRUE
        }
        F <- t(t(F) * sign(fx[top]))
        if(input$gs_enrichfreq_fcweight) {
            F <- t(t(F) * abs(fx[top]))
            wt = TRUE
        } 
        F <- head(F[order(-rowSums(abs(F))),,drop=FALSE], 32)
        F <- F[order(-rowSums(F)),,drop=FALSE]
        F <- as.matrix(F)
       
        par(mfrow=c(1,1), mar=c(6,4,2,0.5), mgp=c(2,0.8,0))
        col1 = grey.colors(ncol(F),start=0.15)
        ylab = ifelse(wt, "weighted frequency", "frequency")
        barplot(t(F), beside=FALSE, las=3, cex.names=0.80, col=col1,
                ylab=ylab)
        
    })

    topEnrichedFreq_text = "<strong>Gene frequency.</strong> The plot shows the number of times a gene is present in the top-N genesets sorted by frequency. Genes that are frequently shared among the top enriched gene sets may suggest driver genes."

    topEnrichedFreq_caption = "<strong>Gene frequency.</strong> The plot shows the number of times a gene is present in the top-N genesets sorted by frequency."
    
    topEnrichedFreq.opts = tagList(
        tipify( radioButtons(ns('gs_enrichfreq_ntop'),"Number of top sets",
                             c(5,10,25,100),inline=TRUE, selected=25),
               "Number of top genesets to consider for counting the gene frequency."),
        tipify( checkboxInput(ns('gs_enrichfreq_gsetweight'),
                              "Weight by geneset size", TRUE),
               "Weight by (inverse) gene set size."), 
        tipify( checkboxInput(ns('gs_enrichfreq_fcweight'),
                              "Weight by FC", TRUE),
               "Weight by fold-change of current contrast.")         
    )
    
    callModule(
        plotModule,
        id = "topEnrichedFreq", label="b",
        func = topEnrichedFreq.RENDER,
        func2 = topEnrichedFreq.RENDER,
        options = topEnrichedFreq.opts,
        info.text = topEnrichedFreq_text,
        height = c(imgH,450), width = c('auto',1000),
        res = c(72,100),
        pdf.width = 10, pdf.height = 5, 
        title = "Frequency in top gene sets"
        ##caption = topEnrichedFreq_caption
    )
    
    ## library(shinyjqui)
    topEnriched_captionALL <- paste(
        "<b>(a)</b>",topEnriched_caption,
        "<b>(b)</b>",topEnrichedFreq_caption)
    
    output$topEnriched_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1,NA),
            fillRow(
                flex = c(2.1,0.05,1),
                plotWidget(ns("topEnriched")),
                br(),
                plotWidget(ns("topEnrichedFreq"))
            ),
            div(HTML(topEnriched_captionALL), class="caption")
        )
    })

    ##================================================================================
    ## Plots
    ##================================================================================

    subplot1_text = "A volcano plot of genes contained in the gene set that is selected from the enrichment analysis Table <code>I</code>."
    subplot2_text = "An enrichment barplot per sample group for the gene set that is selected from the enrichment analysis Table <code>I</code>. Samples can be ungrouped in the barplot by selecting <code>ungroup samples</code> from the plot <i>Settings</i>."
    subplot3_text = "An expression barplot per sample group for the gene that is selected from the genes Table <code>II</code>. Samples can be ungrouped in the barplot by selecting <code>ungroup samples</code> from the plot <i>Settings</i>."
    subplot4_text = "A scatter plot of enrichment scores versus expression values across the samples for the gene set selected from the enrichment analysis Table <code>I</code> and the gene selected from the genes Table <code>II</code>."

    
    ##comp0=colnames(ngs$model.parameters$contr.matrix)[1]
    getcolors <- function(ngs, comp0) {   
        ## get colors (what a mess...)
        contr.matrix <- ngs$model.parameters$contr.matrix
        require(RColorBrewer)
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
        grp.klr = c("grey90",rep(brewer.pal(12,"Paired"),99)[1:ngrp])
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
    subplot.MAR = c(3.5,4,1.5,0.8)

    subplot1.RENDER %<a-% reactive({

        dbg("[subplot1.RENDER] reacted")
        
        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0.5,0.2)*2 )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        req(ngs)
        
        comp=1;gs=1
        comp = input$gs_contrast
        ngs <- inputData()
        req(ngs)

        dbg("[subplot1.RENDER] called")
        
        gx.meta <- ngs$gx.meta$meta[[comp]]
        ## limma1 = sapply(gx.meta[,c("fc","p","q")],function(x) x[,"trend.limma"])
        limma1 = data.frame( meta.fx=gx.meta$meta.fx, meta.q=gx.meta$meta.q)
        gx.annot <- ngs$genes[rownames(gx.meta),c("gene_name","gene_title")]
        ##limma = cbind( gx.meta[,c("gene_name","gene_title")], limma1)
        limma = cbind(gx.annot, limma1)
    
        gs = gset_selected()
        if(is.null(gs)) return(NULL)
        gs <- gs[1]

        dbg("[subplot1.RENDER] gs = ",gs)
        
        ##sel.genes = names(which(ngs$GMT[,gs]!=0))
        jj = match(toupper(GSETS[[gs]]), toupper(limma$gene_name))
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

        dbg("[subplot1.RENDER] drawing plot")
        
        ##par(mar=c(4,3,3,1), mgp=c(2.0,0.8,0), oma=c(1,1.5,1,1.5) )
        gx.volcanoPlot.XY( x=fx, pv=qval, gene=fc.genes,
                          render="canvas", n=5000, nlab=10, 
                          xlim=xlim, ylim=ylim, ## hi.col="#222222",
                          use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                          cex=0.9, lab.cex=1.3, cex.main=0.8,
                          xlab="fold change (log2)",
                          ylab="significance (log10q)",
                          highlight=sel.genes)
        gs = breakstring(gs,50)
        title(gs, cex.main=0.85)
        
    })

    subplot1.PLOTLY %<a-% reactive({
        
        ##par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0.5,0.2)*2 )
        ##par(mar=subplot.MAR)

        dbg("[subplot1.PLOTLY] reacted")
        
        ngs <- inputData()    
        req(ngs)
        
        comp=1;gs=1
        comp = input$gs_contrast
        ngs <- inputData()
        req(ngs)
        
        gx.meta <- ngs$gx.meta$meta[[comp]]
        ##m1 <- intersect(c("trend.limma","notrend.limma","ttest"),colnames(gx.meta$p))[1]
        ##m1
        ##limma1 = sapply(gx.meta[,c("fc","p","q")],function(x) x[,m1])
        limma1 = data.frame( meta.fx=gx.meta$meta.fx, meta.q=gx.meta$meta.q)
        gx.annot <- ngs$genes[rownames(gx.meta),c("gene_name","gene_title")]
        ##limma = cbind( gx.meta[,c("gene_name","gene_title")], limma1)
        limma = cbind(gx.annot, limma1)
    
        gs = gset_selected()
        if(is.null(gs)) return(NULL)
        gs <- gs[1]        
        ##sel.genes = names(which(ngs$GMT[,gs]!=0))
        jj = match(toupper(GSETS[[gs]]), toupper(limma$gene_name))
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
            source = "plot1",
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
            layout( margin = list(b=60) )        

    })
    
    ##----------------------------------------------------------------------
    ## 1: Gene set activation {data-width=200}
    ##----------------------------------------------------------------------
    subplot2.RENDER %<a-% reactive({
        
        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0.5,0.2)*2 )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        req(ngs)
        
        require(RColorBrewer)
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
        pgx.plotGeneExpression(
            ngs, gset, comp=comp0, logscale=TRUE, level="geneset",
            collapse.others=collapse.others, grouped=grouped,
            srt=srt, main="", ylab="gene set enrichment")
        title(breakstring(gset,42,80), cex.main=0.85)
        
    })

    ##----------------------------------------------------------------------
    ## 2: Gene expression {data-width=200}
    ##----------------------------------------------------------------------
    subplot3.RENDER %<a-% reactive({

        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0.5,0.2)*2 )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        req(ngs)

        comp0 = colnames(ngs$model.parameters$contr.matrix)[1]
        comp0 = input$gs_contrast

        has.design <- !is.null(ngs$model.parameters$design)
        collapse.others <- ifelse(has.design, FALSE, TRUE)
        ##collapse.others=TRUE

        require(RColorBrewer)
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
            pgx.plotGeneExpression(
                ngs, probe, comp=comp0, logscale=TRUE, level="gene",
                collapse.others=collapse.others, grouped=grouped,
                srt=srt, main="")
            title(gene, cex.main=0.9)
        }
    })


    ##----------------------------------------------------------------------
    ## 3: Gene - gene set correlation
    ##----------------------------------------------------------------------
    subplot4.RENDER %<a-% reactive({

        par(mfrow=c(1,1), mgp=c(1.8,0.8,0), oma=c(0,0,0.5,0.2)*2 )
        par(mar=subplot.MAR)
        
        ngs <- inputData()    
        req(ngs)
        
        require(RColorBrewer)
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
            klr = paste0(col2hex(klr),"99")
            
            cex1 = c(1.4,0.8,0.3)[cut(length(gx),c(0,100,500,99999))]
            gset1 = breakstring(substring(gset,1,80),32)
                                        #tt = paste( breakstring(gset,40,80), "\nvs.", gene,"expression")
            tt = paste( breakstring(gset,40,80), " vs. ", gene)
            plot( gx, sx, col=klr, main=tt,
                 ylab = "gene set enrichment",
                 xlab = paste(gene,"expression"),
                 cex.lab=1, pch=19, cex=1.0*cex1, cex.main=0.85)
            abline( lm(sx ~ gx), lty=2, lwd=0.7, col="black" )
        }    
    })


    callModule(
        plotModule,
        id = "subplot1", 
        func = subplot1.PLOTLY, plotlib="plotly",
        ##func = subplot1.RENDER,
        ##func2 = subplot1.RENDER,        
        info.text = subplot1_text,
        pdf.width=6, pdf.height=6,
        res = c(72,100),
        height = imgH, 
        title="Volcano plot", label="a"
    )

    callModule(
        plotModule,
        id="subplot2", 
        func = subplot2.RENDER,
        func2 = subplot2.RENDER,
        info.text = subplot2_text,
        pdf.width=6, pdf.height=6, res=72,
        height = imgH,
        options = tagList(
            tipify( checkboxInput(
                ns('gs_ungroup1'),'ungroup samples',FALSE),
                "Ungroup samples in the plot", placement="top",
                options = list(container = "body"))),
        title="Enrichment barplot", label="b"
    )

    callModule(
        plotModule,
        id="subplot3", 
        func = subplot3.RENDER,
        func2 = subplot3.RENDER,
        info.text = subplot3_text,
        pdf.width=6, pdf.height=6, res=72,
        height = imgH,
        options = tagList(
            tipify( checkboxInput(ns('gs_ungroup2'),'ungroup samples',FALSE),
                   "Ungroup samples in the plot", placement="top",
                   options = list(container = "body"))),
        title = "Expression barplot", label="c"
    )

    callModule(
        plotModule,
        id="subplot4",
        func = subplot4.RENDER,
        func2 = subplot4.RENDER,
        info.text = subplot4_text,
        pdf.width=6, pdf.height=6, res=72,
        height = imgH,
        title = "Enrichment vs. expression", label="d"
    )

    ## output <- attachModule(output, subplot1_module)
    ## output <- attachModule(output, subplot2_module)
    ## output <- attachModule(output, subplot3_module)
    ## output <- attachModule(output, subplot4_module)
    
    enrichplots_caption = "<b>Enrichment plots</b> associated with the gene set (selected in <b>Table I</b>) and gene (selected in <b>Table II</b>). <b>(a)</b> Volcano-plot showing significance versus fold-change on the y and x axes, respectively. Genes in the gene set are highlighted in blue. <b>(b)</b> Barplot of the gene set enrichment in the groups. <b>(c)</b> Barplot of the gene expression of the gene. <b>(d)</b> Scatter plot of the enrichment versus the expression of the selected geneset and gene, on the y and x axes, respectively."
   
    output$subplots_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1,0.05,NA),
            fillRow(
                id = ns("subplots"),
                height = imgH,
                flex=c(1,1,1,1), ##height = 370,
                plotWidget(ns("subplot1")),
                plotWidget(ns("subplot2")),
                plotWidget(ns("subplot3")),
                plotWidget(ns("subplot4"))
            ),
            br(),
            div(HTML(enrichplots_caption),class="caption")
        )
    })
    dragula(ns("subplots"))

    ##================================================================================
    ## Compare
    ##================================================================================

    compare.RENDER %<a-% reactive({
        
        ngs <- inputData()
        req(ngs,input$gs_contrast)
        
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
                qv0 <- ngs$gset.meta$meta[[cmp]][gset,"meta.q"]
                gs1 = breakstring(gset,28,50,force=FALSE)
                cmp <- paste0("@",cmp)
                if(i==1) cmp <- paste0(gset,"\n",cmp)
                gsea.enplot(rnk0, genes, names=NULL, ##main=gs,
                            main=cmp, xlab="",
                            cex.main=0.80, len.main=80)
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
                qv0 <- ngs$gset.meta$meta[[cmp]][gset,"meta.q"]
                gs1 = breakstring(gset,28,50,force=FALSE)
                cmp <- paste0("@",cmp)
                if(i==1) cmp <- paste0(gset,"\n",cmp)
                gsea.enplot(rnk0, genes, names=NULL, ##main=gs,
                            main=cmp, xlab="",
                            cex.main=0.80, len.main=80)
                qv1 = formatC(qv0,format="e", digits=2)
                legend("topright", paste("q=",qv1), bty="n",cex=0.85)
            }
        }
        
    })


    compare_text = "Under the <strong>Compare</strong> tab, enrichment profiles of the selected geneset in enrichment Table <code>I</code> can be visualised against all available contrasts."

    compare_caption = "<b>Enrichment across contrasts.</b> Enrichment plots for the selected gene set (in <b>Table I</b>) across multiple contrasts. The figure allows to quickly compare the enrichment of a certain gene set across all other comparisons."

    compare_module_opts = tagList()
    
    callModule(
        plotModule,
        id="compare",
        func = compare.RENDER,
        func2 = compare.RENDER,
        options = compare_module_opts,
        height = c(imgH,450), width = c("auto",1500), res=95,
        pdf.width=14, pdf.height=4, 
        title = "Enrichment of gene set across multiple contrasts",
        info.text = compare_text
        ##caption = compare_caption
    )
    ## output <- attachModule(output, compare_module)

    output$compare_UI <- renderUI({
        fillCol(
            height = rowH,
            flex=c(1,NA), ##height = 370,
            plotWidget( ns("compare")),
            div(HTML(compare_caption),class="caption")
        )
    })

    ##================================================================================
    ## Volcano (all)
    ##================================================================================

    volcanoAll.RENDER %<a-% reactive({
        ##renderPlotly({
        require(metap)
        ngs = inputData()
        req(ngs)
        if(is.null(input$gs_features)) return(NULL)
        
        meta = ngs$gset.meta$meta
        gsmethod0 = colnames(meta[[1]]$fc)
        gsmethod = intersect(gsmethod0, GSET.DEFAULTMETHODS)    
        gsmethod <- input$gs_method
        if(is.null(gsmethod) || length(gsmethod)==0) return(NULL)
        
        ng = length(meta)
        nn = c(2, max(ceiling(ng/2),5))
        ##if(ng>12) nn = c(3,8)
        par(mfrow=nn, mar=c(2,4,2.3,2)*0, mgp=c(2.6,1,0))
        n = ceiling(sqrt(ng))
        if(ng>24) {
            n = max(ceiling(ng/3),6)
            par(mfrow=c(3,n), mar=c(4,4,2,2)*0)
        } else if(FALSE && ng <= 3) {
            par(mfrow=c(1,3), mar=c(4,4,2,2)*0)
        } else {
            n = max(ceiling(ng/2),6)
            par(mfrow=c(2,n), mar=c(4,4,2,2)*0)
        }

        fdr = as.numeric(input$gs_fdr)    
        lfc = as.numeric(input$gs_lfc)
        sel.gsets = COLLECTIONS[[1]]
        sel.gsets = COLLECTIONS[[input$gs_features]]

        withProgress(message="computing volcano plots ...", value=0, {
            i=1
            for(i in 1:length(meta)) {
                mx <- calculateMeta(i, gsmethod, ngs=ngs)
                is.sig <- (mx[,"qv"] <= fdr & abs(mx[,"fc"]) >= lfc)
                sig.gs = rownames(mx)[which(is.sig)]
                sig.gs <- intersect(sel.gsets, sig.gs)
                gx.volcanoPlot.XY(
                    x = mx[,"fc"], pv = mx[,"qv"],
                    use.fdr=TRUE, p.sig=fdr, lfc=lfc,                
                    gene = substring(rownames(mx),1,35), 
                    xlab = "effect size (NES)", lab.cex=1.5, nlab=3,
                    render="canvas", n=1000, highlight=sig.gs,
                    cex=1, cex.axis=1.3, cex.main=1.4, axes=FALSE,
                    ylim=c(0,10), main="" )
                ##title(names(meta)[i],line=-1)
                legend("topright",names(meta)[i], cex=1.2, bg="white")

                ## draw axis if first column or last row
                ##n=nn[2]
                is.first = (i%%n==1)
                last.row = ( (i-1)%/%n == (length(meta)-1)%/%n )
                if(is.first) axis(2, tcl=0.5, mgp=c(-2,-1.5,0))
                if(last.row) axis(1, tcl=0.5, mgp=c(-2,-1.5,0))
                box()

                ##volcano_plot(limma, render="plotly", n=1000, cex=1, highlight=genes)
                incProgress( 1/length(meta) )
            }
        })
    })
    
    volcanoAll_text = "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for gene sets across all contrasts. This provides users an overview of the statistics across all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."


    volcanoAll_caption = "<b>Volcano plots for all contrasts.</b> Simultaneous visualisation of volcano plots of gene set enrichment across all contrasts. Volcano-plot are plotting enrichment score versus significance on the x and y axes, respectively. Experimental contrasts showing better statistical significance will show volcano plots with 'higher' wings."

    callModule(
        plotModule,
        id = "volcanoAll",
        func = volcanoAll.RENDER,
        func2 = volcanoAll.RENDER,
        height = c(imgH,450), width = c("auto",1500), res=c(72,85),
        pdf.width=15, pdf.height=5, 
        title="Volcano plots for all contrasts",
        info.text = volcanoAll_text
        ##caption = volcanoAll_caption        
    )
    ##output <- attachModule(output, volcanoAll_module)
    
    output$volcanoAll_UI <- renderUI({
        fillCol(
            height = rowH,
            flex=c(1,NA), ##height = 370,
            plotWidget(ns("volcanoAll")),
            div(HTML(volcanoAll_caption), class="caption")
        )
    })

    ##================================================================================
    ## Volcano (methods)
    ##================================================================================

    volcanoMethods.RENDER %<a-% reactive({
        ##renderPlotly({
        ngs <- inputData()    
        req(ngs, input$gs_features)
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
        
        ng = ncol(fx)
        nn = c(2, max(ng/2,5))
        par(mfrow=nn, mar=c(2,4,2.3,2)*0, mgp=c(2.6,1,0))

        i=1
        for(i in 1:ng) {
            
            is.sig <- ( qv[,i] <= fdr & abs(fx[,i]) >= lfc)
            sig.gs = rownames(mx)[which(is.sig)]
            sig.gs <- intersect(sel.gsets, sig.gs)
            
            method = colnames(fx)[i]
            gx.volcanoPlot.XY(
                x = fx[,i], pv = qv[,i],
                use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                ##gene = substring(rownames(mx),1,35),
                gene = rownames(mx),
                xlab = "effect size (NES)", ylim=c(0,10), 
                lab.cex=1.5, nlab=3, axes=FALSE, 
                render="canvas", n=1000, highlight=sig.gs,
                cex=1, cex.axis=1.3, main="")

            ##title(mt, line=-1.5, cex.main=1.4)
            legend("topright",method,bg="white", cex=1.2)
            
            ##volcano_plot(limma, render="plotly", n=1000, cex=1, highlight=genes)
            ## draw axis if first column or last row
            n=nn[2]
            is.first = (i%%n==1)
            last.row = ( (i-1)%/%n == (ng-1)%/%n )
            if(is.first) axis(2, tcl=0.5, mgp=c(-2,-1.5,0))
            if(last.row) axis(1, tcl=0.5, mgp=c(-2,-1.5,0))
            box()

        }

    })

    volcanoMethods_text = "The <strong>Volcano (methods)</strong> panel displays the volcano plots provided by different enrichment calculation methods. This provides users an quick overview of the sensitivity of the statistical methods at once. Methods showing better statistical significance will show volcano plots with 'higher' wings."

    volcanoMethods_caption = "<b>Volcano plots for all methods.</b> Simultaneous visualisation of volcano plots of gene sets for different enrichment methods. Methods showing better statistical significance will show volcano plots with 'higher' wings."

    callModule(
        plotModule,
        id="volcanoMethods",
        func = volcanoMethods.RENDER,
        func2 = volcanoMethods.RENDER,
        height = c(imgH,450), width = c("auto",1500), res=c(72,85),
        pdf.width=15, pdf.height=5, 
        title="Volcano plots for all methods",
        info.text = volcanoMethods_text
        ##caption = volcanoMethods_caption
    )
    ## output <- attachModule(output, volcanoMethods_module)

    output$volcanoMethods_UI <- renderUI({
        fillCol(
            height = rowH,
            flex=c(1,NA), ##height = 370,
            plotWidget(ns("volcanoMethods")),
            div(HTML(volcanoMethods_caption), class="caption")
        )
    })
    
    ##================================================================================
    ## GeneMap (dev)
    ##================================================================================

    genemap.RENDER %<a-% reactive({
        
        require(Matrix)

        ngs <- inputData()
        req(ngs, input$gs_contrast)
        
        ## -------------- get the gene-centric fold-changes (default LIMMA)
        comp1=1
        comp1 = input$gs_contrast
        gx.meta <- ngs$gx.meta$meta[[comp1]]
        ##limma1 = sapply(gx.meta[,c("fc","p","q")],function(x) x[,"trend.limma"])
        limma1 = data.frame( fc=gx.meta$meta.fx, meta.q=gx.meta$meta.q)        
        ##limma = cbind( ngs$gx.meta$meta[[comp1]][,c("gene_name","gene_title")], limma1)
        gx.annot <- ngs$genes[rownames(gx.meta),c("gene_name","gene_title")]
        ##limma = cbind( gx.meta[,c("gene_name","gene_title")], limma1)
        limma = cbind(gx.annot, limma1)
        
        gs = rownames(ngs$gsetX)[1]
        gs = gset_selected()
        if(is.null(gs)) return(NULL)
        gs <- gs[1]
        
        gsets = colnames(ngs$GMT)
        if(!(gs %in% gsets)) {
            cat("warning:: geneset",gs,"not in GSETS!!\n")
            return(NULL)
        }
        
        ## ------------------- compute closests neighbours (gene set)
        rpt = isolate( getGeneSetTable() )
        mx = ngs$gset.meta$meta[[comp1]]
        mx = mx[intersect(rownames(mx),rownames(rpt)),]  ## only those currently selected
        mx = mx[intersect(rownames(mx),gsets),]  ## only those that have GMT

        ##genes = rownames(ngs$X)
        ##genes = ngs$genes$gene_name
        kk <- intersect(rownames(ngs$gsetX), colnames(ngs$GMT))    
        G = ngs$GMT[,kk]    
        colnames(G) = kk
        gg = intersect(rownames(G),as.character(limma$gene_name))
        fc = abs(limma[match(gg,limma$gene_name),"fc"])    
        g1 = as.matrix(G[gg,rownames(mx)])*fc
        g2 = G[gg,gs]*fc
        suppressWarnings( rho1 <- cor(g1, g2)[,1] )
        
        ## ----------- create map with closest neigbours (what crappy coding...)
        gs.top = head(names(sort(-abs(rho1))),25)  ## how many gene sets
        fx = mx[gs.top,"meta.fx"]
        M = t(G[,gs.top]) * fx
        rownames(M) = gs.top
        M = M[,intersect(colnames(M),limma$gene_name)]
        fc = limma[match(colnames(M),limma$gene_name),"fc"]
        M = t(t(M) * abs(fc))
        jj = head(order(-Matrix::colSums(abs(M))), 80)  ## how many genes
        M = M[,jj]
        M = 1*(M!=0)
        fx = mx[rownames(M),"meta.fx"]
        fc = limma[match(colnames(M),limma$gene_name),"fc"]
        fc[is.na(fc)] = 0
        fx[is.na(fx)] = 0    
        M = t( t(M * rank(abs(fx)) ) * fc )
        M = sign(M) * abs(M)**0.33
        M = as.matrix(M)
        ##d3heatmap(M, scale="none", colors="Spectral", cexRow=0.7, cexCol=0.6 )
        hc <- gx.heatmap(M-1e-8, mar=c(6,38), cexRow=1.00, cexCol=0.95,
                         scale='none', keysize=0.3, key=FALSE)
        
    })

    genemap_text = "Co-activation heatmap of top N = {25} enriched gene sets and their common genes."

    genemap_caption = "<b>Co-activation heatmap.</b> Clustered heatmap of top most correlated gene sets and their shared genes. The top gene sets most correlated with the selected gene set (in Table I) are shown. Red corresponds to overexpression, blue to downregulation of the gene." 

    callModule(
        plotModule,
        id = "genemap",
        func = genemap.RENDER,
        func2 = genemap.RENDER,
        height = c(imgH,450), width = c("auto",1500), res=c(80,80),
        pdf.width=14, pdf.height=4, ## res=65,
        title = "Co-activation heatmap",
        info.text = genemap_text
        ##caption = genemap_caption
    )

    output$genemap_UI <- renderUI({
        fillCol(
            height = rowH,
            flex=c(1,NA), ##height = 370,
            plotWidget(ns("genemap")),
            div(HTML(genemap_caption), class="caption")
        )
    })

    ##================================================================================
    ## Enrichment table
    ##================================================================================
    
    gset_selected <- reactive({
        ##i = as.integer(input$gseatable_rows_selected)
        i = as.integer(gseatable$rows_selected())
        if(is.null(i) || length(i)==0) return(NULL)
        rpt = getGeneSetTable()
        gs = rownames(rpt)[i]
        return(gs)
    })

    geneDetails <- reactive({
        ## return details of the genes in the selected gene set
        ##
        
        ngs <- inputData()
        req(ngs,input$gs_contrast)
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
        req(gxmethods)
        ##limma1 = sapply(mx[,c("fc","p","q")], function(x) x[,"trend.limma"])
        ##limma1.fc <- rowMeans(mx$fc[,gxmethods,drop=FALSE],na.rm=TRUE)
        limma1.fc <- mx$meta.fx
        limma1.pq = sapply(mx[,c("p","q")], function(x) {
            apply(x[,gxmethods,drop=FALSE],1,max,na.rm=TRUE)
        })
        limma1 <- cbind( fc=limma1.fc, q=limma1.pq)
        ##limma  = cbind( ngs$gx.meta$meta[[comp]][,c("gene_name","gene_title")], limma1)
        rownames(limma1) <- rownames(mx)
        
        ## in multi-mode we select *common* genes
        ns <- length(gs)
        gmt1 <- ngs$GMT[,gs,drop=FALSE]
        genes = rownames(gmt1)[which(Matrix::rowSums(gmt1!=0)==ns)]
        genes = intersect(genes, ngs$genes[rownames(mx),"gene_name"])
        genes = setdiff(genes, c("",NA,"NA"," "))

        title = as.character(GENE.TITLE[genes])
        title[is.na(title)] <- " "
        
        rpt <- data.frame("gene_name"=genes, "gene_title"=as.character(title) )
        genes = rpt[,"gene_name"]
        genes1 <- ngs$genes[rownames(limma1),"gene_name"]
        limma1 = limma1[match(genes, genes1),,drop=FALSE ]    
        avg.rho <- rowMeans(cor(t(ngs$X[rownames(limma1),,drop=FALSE]),
                                t(ngs$gsetX[gs,,drop=FALSE])))
        
        rpt = cbind(rpt, limma1, gset.rho=avg.rho)
        rpt = rpt[which(!is.na(rpt$fc) & !is.na(rownames(rpt))),,drop=FALSE]
        ##rpt = data.frame(rpt, check.names=FALSE)

        if(nrow(rpt)>0) {
            rpt = rpt[order(-abs(rpt$fc)),,drop=FALSE]
        }
        return(rpt)
    })

    gene_selected <- reactive({

        ngs <- inputData()
        req(ngs)

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


    gseatable.RENDER <- reactive({

        rpt = getGeneSetTable()
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

        ## wrap with known links.
        rpt$GS <- wrapHyperLink(rpt$GS, rownames(rpt))    
        selectmode = "single"
        selectmode

        is.numcol <- sapply(rpt, is.numeric)
        numcols <- which( is.numcol & !colnames(rpt) %in% c("size"))
        numcols <- colnames(rpt)[numcols]

        ##rpt = format(rpt, digits=4)
        DT::datatable(rpt,
                      class = 'compact cell-border stripe hover',
                      rownames=FALSE, escape=-1,
                      ##extensions = c('Buttons','Scroller'),
                      extensions = c('Scroller'),                  
                      fillContainer = TRUE,
                      selection = list(mode=selectmode, target='row', selected=1),                      
                      options=list(
                          dom = 'lfrtip',
                          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      )) %>%
            formatSignif(numcols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle(fx.col, 
                                background = color_from_middle( fx, 'lightblue', '#f5aeae'))
    })

    genetable.RENDER <- reactive({

        rpt <- geneDetails()    
        if(is.null(rpt)) return(NULL)
        
        rpt$gene_title <- NULL    
        if(!is.null(rpt) && nrow(rpt)>0 ) {
            jj = which(sapply(rpt,is.numeric))
            rpt[,jj] = round(rpt[,jj],digits=4)
            jj = which( sapply(rpt,is.character) |  sapply(rpt,is.factor) )
            if(length(jj)>0) rpt[,jj] = apply(rpt[,jj,drop=FALSE],2,shortstring,60)
        } else {
            rpt <- data.frame(0,0,0,0,0)[0,]
            colnames(rpt) <- c("gene_name","fc","p","q","gset.rho")
        }
        ##rpt <- rpt[,c("gene_name","fc","p","q","gset.rho"),drop=FALSE]

        colnames(rpt) <- sub("^GS$","gene set",colnames(rpt))
        numeric.cols <- which(sapply(rpt, is.numeric))
        numeric.cols

        tbl <- DT::datatable(rpt,
                             class = 'compact cell-border stripe', rownames=FALSE,
                             extensions = c('Scroller'),
                             selection = list(mode="single", target='row', selected=1),
                             fillContainer = TRUE,
                             options=list(
                                 dom = 'lfrtip',
                                 ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                                 scrollX = TRUE, scrollY = tabH, scroller=TRUE, deferRender=TRUE
                             )) %>%
            formatSignif(numeric.cols,4)
        
        if( nrow(rpt)>0 && ("fc" %in% colnames(rpt)) ) {
            fx = rpt[,"fc"]
            tbl <- tbl %>%
                DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle("fc", background = color_from_middle( fx, 'lightblue', '#f5aeae'))
        }
        tbl
    })

    gseatable_text = paste("Similar to the differential gene expression analysis, users can perform differential expression analysis on a geneset level that is referred as gene set enrichment analysis. To ensure statistical reliability, the platform performs the gene set enrichment analysis using multiple methods, including",a_Spearman,", ",a_GSVA,", ",a_ssGSEA,", ",a_Fisher,", ",a_GSEA,", ",a_camera," and ",a_fry,".<br><br>The combined result from the methods is displayed in this table, where for each geneset the <code>meta.q</code> corresponds to the highest <code>q</code> value provided by the methods and the number of <code>stars</code> indicate how many methods identified the geneset as significant (<code>q < 0.05</code>). The table is interactive; users can sort it by <code>logFC</code>, <code>meta.q</code> and <code>starts</code>. Additionally, the list of genes in that geneset are displayed in the second table on the right. Users can filter top N = {10} differently enriched gene sets in the table by clicking the <code>top 10 gene sets</code> from the table <i>Settings</i>.")

    gseatable_opts = tagList(
        tipify( checkboxInput(ns('gs_top10'),'top 10 gene sets',FALSE),
               "Display only top 10 differentially enirhced gene sets (positively and negatively) in the <b>enrihcment analysis</b> table.", placement="top", options = list(container = "body")),
        tipify(checkboxInput(ns('gs_showqvalues'),'show q-values',FALSE),
               "Show q-values of all statistical methods in the table.", 
               placement="top", options = list(container = "body"))    
    )

    gseatable <- callModule(
        tableModule, 
        id="gseatable",
        func = gseatable.RENDER,
        info.text = gseatable_text, label="I",
        options = gseatable_opts,
        title="Enrichment analysis",
        info.width="500px",
        height = c(265, 700)
    )
    ## output <- attachModule(output, gseatable_module)

    genetable_text = "By clicking on a gene set in the table <code>I</code>, it is possible to see the gene list of that gene set in this table. By clicking on a gene in this table, users can check the expression status of the gene for the selected contrast in the <code>Expression</code> barplot and its correlation to the gene set in the <code>Gene to gene set correlation</code> scatter plot under the <code>Plots</code> section."
    
    genetable <- callModule(
        tableModule,
        id = "genetable",
        func=genetable.RENDER,
        info.text = genetable_text,
        title="Genes", label="II",
        height = c(265,700), width = c('100%',800)
    )
    ##output <- attachModule(output, genetable_module)
    
    tables_caption = "<b>Enrichment tables</b>. <b>(I)</b> Table summarizing the statistical results of the gene set enrichment analysis for selected contrast. The number of stars indicate how many methods identified the geneset significant. <b>(II)</b> Table showing the fold-change, statistics and correlation of the genes in the selected gene set."

    output$tables_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1, NA),
            fillRow(
                ## height = 200,
                flex = c(2,0.1,1), 
                tableWidget(ns("gseatable")),
                br(),
                tableWidget(ns("genetable"))        
            ),
            div(HTML(tables_caption),class="caption")
        )
    })

    ##================================================================================
    ## Enrichment (all)
    ##================================================================================

    fctable.RENDER <- reactive({
        
        ngs <- inputData()

        ## get all contrasts
        F <- sapply( ngs$gset.meta$meta, function(x) x[,"meta.fx"])
        colnames(F) <- gsub("_"," ",colnames(F))

        qv <- sapply( ngs$gset.meta$meta, function(x) x[,"meta.q"])
        rownames(qv) <- rownames(F) <- rownames(ngs$gset.meta$meta[[1]])
        fc.var <- round( rowMeans(F**2,na.rm=TRUE), digits=3)
        gs <- substring(rownames(F),1,60)
        F1 <- data.frame( geneset=gs, fc.var=fc.var, round(F,digits=3), check.names=FALSE)

        ## get current filtered geneset and extract names of gene sets
        rpt = getGeneSetTable()
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
                          scrollX = TRUE, scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle( "fc.var",
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(fc.var, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  %>%
                DT::formatStyle( colnames(F),
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(F[,], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    gx_fctable_text = "The <strong>Enrichment (all)</strong> panel reports the gene set enrichment for all contrasts in the selected dataset."
    
    gx_fctable_caption = "<b>Enrichment for all contrasts.</b> Table summarizing the enrichment for all gene sets across all contrasts. The column `fc.var` corresponds to the variance of the gene set across all contrasts."

    callModule(
        tableModule,
        id = "fctable",
        func = fctable.RENDER,
        title ="Gene set enrichment for all contrasts",
        info.text = gx_fctable_text,
        caption = gx_fctable_caption,
        height = c(280,700)
    )
    ##output <- attachModule(output, fctable_module)

    ## library(shinyjqui)
    output$fctable_UI <- renderUI({
        fillCol(
            height = rowH,
            plotWidget(ns("fctable"))
        )
    })

    
    ##================================================================================
    ## FDR table
    ##================================================================================

    require(kableExtra)
    
    FDRtable.RENDER <- reactive({
        
        ngs <- inputData()    
        req(ngs, input$gs_method)
        
        meta <- ngs$gset.meta
        test = GSET.DEFAULTMETHODS
        test <- input$gs_method
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
        ##gs.up %>% kable("html") %>%
        ##    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
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
                          scrollX = TRUE, scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle(colnames(sig.up),
                                background = styleColorBar(c(0,maxsig), '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  %>% 
                DT::formatStyle(colnames(sig.down),
                                background = styleColorBar(c(0,maxsig), 'lightblue'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  
    })

    FDRtable_text = "The <strong>FDR table</strong> panel reports the number of significant gene sets at different FDR thresholds, for all contrasts and all methods. Using the table the user can determine which statistical methods perform better for a particular contrast."

    FDRtable_caption = "<b>FDR table.</b> Number of significant gene sets versus different FDR thresholds, for all contrasts and all methods. The blue color denote the number of downregulated genes, the red color for upregulated genes."
    
    callModule(
        tableModule,
        id = "FDRtable",
        func = FDRtable.RENDER,
        title = 'Number of significant gene sets',
        info.text = FDRtable_text,
        caption = FDRtable_caption,
        height = c(280,700)
    )
    ##output <- attachModule(output, FDRtable_module)

    ## library(shinyjqui)
    output$FDRtable_UI <- renderUI({
        fillCol(
            height = rowH,
            tableWidget(ns("FDRtable"))
        )
    })


    ## reactive values to return to parent environment
    outx <- list(selected_gsetmethods=selected_gsetmethods)
    return(outx)

} ## end-of-Board


