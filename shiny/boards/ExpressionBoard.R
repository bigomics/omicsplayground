##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing ExpressionBoard")

ExpressionInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

ExpressionUI.test <- function(id) {
    ns <- NS(id)  ## namespace
    tabsetPanel(
        tabPanel("Table",uiOutput(ns("expr_tables_UI"))),
        tabPanel("Foldchange (all)",uiOutput(ns("expr_fctable_UI"))),
        tabPanel("FDR table",uiOutput(ns("expr_FDRtable_UI")))                       
    )
}

ExpressionUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        flex = c(1.4,1),
        height = 730,
        tabsetPanel(
            id = ns("tabs1"),
            tabPanel("Plot",uiOutput(ns("expr_plots_UI"))),
            tabPanel("Top genes",uiOutput(ns("expr_topgenesUI"))),
            tabPanel("Volcano (all)",uiOutput(ns("expr_volcanoAll_UI"))),
            tabPanel("Volcano (methods)",uiOutput(ns("expr_volcanoMethodsUI")))
        ),
        tabsetPanel(
            id = ns("tabs2"),
            tabPanel("Table",uiOutput(ns("expr_tables_UI"))),
            tabPanel("Foldchange (all)",uiOutput(ns("expr_fctable_UI"))),
            tabPanel("FDR table",uiOutput(ns("expr_FDRtable_UI")))                       
        )
    )
}

ExpressionBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    ## reactive functions from shared environment
    inputData <- env[["load"]][["inputData"]]

    fullH = 730
    rowH = 345  ## row height of panels
    imgH = 270  ## height of images
    tabH = 160  ## height of tables
    tabH = "70vh"  ## height of tables
    
    description = "<b>Differential Expression Analysis.</b> Compare expression between
two conditions. Determine which genes are significantly downregulated or overexpressed in one of the groups."
    output$description <- renderUI(HTML(description))

    gx_infotext ="The <strong>Differential Expression Analysis</strong> module compares expression between two conditions (i.e. tumor versus control), which is one of the fundamental analysis in the transcriptomics data analytics workflow. For each comparison of two conditions (also called \'contrast\'), the analysis identifies which genes are significantly downregulated or overexpressed in one of the groups.

<br><br>The <strong>Plots</strong> panel shows volcano and MA plots for the chosen contrast. It also shows the so-called \'signature\', i.e. the top downregulated and overexpressed genes, for that contrast. The <strong>Top genes</strong> panel shows the average expression plots across the samples for top differentially expressed genes within the selected comparison. A very useful feature of the platform is that it can display volcano plots for all comparisons simultaneously under the <strong>Volcano (all)</strong> panel. This provides users an overview of the statistics of all comparisons. The <strong>Table</strong> panel on the bottom shows the results of the statistical tests. The <strong>Foldchange (all)</strong> panel reports the gene fold changes for all contrasts.

<br><br>EXPERT MODE ONLY: To compare the different statistical methods, the <strong>Volcano (methods)</strong> panel shows volcano plots of all methods. The <strong>FDR table</strong> panel reports the number of significant genes at different FDR thresholds for all contrasts.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=3' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>"


    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    FDR.VALUES = c(1e-9,1e-6,1e-3,0.01,0.05,0.1,0.2,0.5,1)

    gx_testmethod_text = "Select a method for the statistical test. To increase the statistical reliability of the Omics Playground, we perform the DE analysis using commonly accepted methods in the literature, including t-test (standard, Welch), limma (no trend, trend, voom), edgeR (QLF, LRT), and DESeq2 (Wald, LRT), and merge the results."
    GX.DEFAULTTEST="trend.limma"
    GX.DEFAULTTEST=c("trend.limma","edger.qlf","deseq2.wald","edger.lrt")
    
    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("gx_info"), "Tutorial", icon = icon("youtube")),
                   "Show more information about this module."),
            hr(), br(),             
            tipify( selectInput(ns("gx_contrast"), "Contrast:", choices=NULL),
                   "Select a contrast of interest for the analysis.", placement="top"),
            tipify( selectInput(ns("gx_features"),"Gene family:", choices=NULL, multiple=FALSE),
                   "Choose a specific gene family for the analysis.", placement="top"),
            fillRow( flex=c(1,1),
                    tipify( selectInput(ns("gx_fdr"),"FDR", choices=FDR.VALUES, selected=1),
                           "Set the false discovery rate (FDR) threshold.", placement="top"),
                    tipify( selectInput(ns("gx_lfc"),"logFC threshold",
                                        choices=c(0,0.2,0.5,1,2,5), selected=0.5),
                           "Set the logarithmic fold change (logFC) threshold.", placement="top")
                    ),
            br(),br(),br(),br(),
            tipify( actionLink(ns("gx_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Toggle advanced options.", placement="top"),
            br(),br(),
            conditionalPanel(
                "input.gx_options % 2 == 1", ns=ns,
                tagList(
                    tipify( checkboxGroupInput(ns('gx_testmethod'),'Statistical methods:',
                                               choices=NULL, inline=TRUE),
                           gx_testmethod_text, placement="right", options = list(container = "body"))
                )
            )
        )
        return(ui)
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$gx_info, {
        showModal(modalDialog(
            title = HTML("<strong>Differential Expression Analysis Board</strong>"),
            HTML(gx_infotext),
            easyClose = TRUE, size="l"))
    })
    
    ## update choices upon change of data set 
    observe({
        ngs <- inputData()
        req(ngs)
        
        contr <- colnames(ngs$model.parameters$contr.matrix)
        updateSelectInput(session, "gx_contrast", choices=sort(contr))
        ##fam <- names(ngs$families)
        ##fam <- grep("^TISSUE|^COMPARTMENT|^CELLTYPE|^GOCC|^DISEASE|^CUSTOM",names(GSETS),value=TRUE)
        fam <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        updateSelectInput(session, "gx_features",choices=fam)

        ## available statistical methods
        gx.methods = colnames(ngs$gx.meta$meta[[1]]$fc) ## available
        sel1 = c(intersect(GX.DEFAULTTEST,gx.methods),gx.methods)
        sel1 = head(unique(sel1),3) ## maximum three!!

        updateCheckboxGroupInput(session, 'gx_testmethod',
                                 choices = sort(gx.methods),
                                 selected = sel1)

        updateCheckboxInput(session, "gx_ungroup", value= (ncol(ngs$X)<=8) )        
        
    })

    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================
    
    selected_gxmethods <- reactive({
        ngs <- inputData()
        gx.methods0 = colnames(ngs$gx.meta$meta[[1]]$fc)
        test = input$gx_testmethod
        test = intersect(test,gx.methods0) ## maximum three
        test
    })

    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================
    
    comparison=1;testmethods=c("trend.limma")
    getDEGtable <- function(ngs, testmethods, comparison, add.pq) {
        ##ngs = inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)
        
        if(is.null(testmethods)) return(NULL)
        if(is.null(comparison)) return(NULL)
        if(length(testmethods)==0 || testmethods=="") return(NULL)
        if(length(comparison)==0  || comparison=="") return(NULL)
        
        ## build meta table
        mx = ngs$gx.meta$meta[[comparison]]
        if(is.null(mx)) return(NULL)
        mm  = colnames(unclass(mx$p))
        testmethods = intersect(mm, testmethods)
        
        mx.p  = unclass(mx$p[,testmethods,drop=FALSE]) ## get rid of AsIs
        mx.q  = unclass(mx$q[,testmethods,drop=FALSE])
        mx.fc = unclass(mx$fc[,testmethods,drop=FALSE])
        ##mx$score = mx$fc * (-log10(1e-100+mx$q) )    

        mx.fc[ is.infinite(mx.fc) | is.nan(mx.fc) ] <- NA
        mx.p[ is.infinite(mx.p) | is.nan(mx.p) ] <- NA
        mx.q[ is.infinite(mx.q) | is.nan(mx.q) ] <- NA

        ## must recompute meta parameters (maxQ method)
        mx$meta.p = apply(mx.p,1,max,na.rm=TRUE)
        mx$meta.q = apply(mx.q,1,max,na.rm=TRUE)
        mx$meta.fx = rowMeans(mx.fc,na.rm=TRUE)  
        ##mx$meta.score = rowMeans(mx$score,na.rm=TRUE)  

        fdr = 0.05
        ##fdr = as.numeric(input$gx_fdr)
        star.symbols = sapply(1:20,function(i) paste(rep("\u2605",i),collapse=""))
        stars = c("",star.symbols)[ 1 + rowSums(mx.q <= fdr, na.rm=TRUE)]        
        
        ## recalculate group averages???
        y0 <- ngs$model.parameters$exp.matrix[,comparison]
        names(y0) <- rownames(ngs$model.parameters$exp.matrix)
        AveExpr1 <- rowMeans(ngs$X[rownames(mx),names(which(y0>0)),drop=FALSE])
        AveExpr0 <- rowMeans(ngs$X[rownames(mx),names(which(y0<0)),drop=FALSE])
        
        ##logFC <- unclass(ngs$gx.meta$meta[[comparison]][,"fc"])[,"trend.limma"]
        ##logFC <- ngs$gx.meta$meta[[comparison]][,"meta.fx"]
        logFC <- mx$meta.fx
        ##logFC <- (AveExpr1 - AveExpr0)  ## override ??? yes: see "contrast in R" Rose Maier 2015...
        ## [hack] adjust averages to match logFC...
        mean0 <- (AveExpr0+AveExpr1)/2
        AveExpr1 <- mean0 + logFC/2
        AveExpr0 <- mean0 - logFC/2

        ##gene.annot = mx[,grep("^gene|^chr",colnames(mx)),drop=FALSE]
        aa <- intersect(c("gene_name","gene_title","chr"), colnames(ngs$genes))
        gene.annot <- ngs$genes[rownames(mx),aa]
        metaq <- mx[,c("meta.q"),drop=FALSE]    
        res = data.frame( gene.annot, logFC=logFC, stars, metaq,
                         AveExpr0, AveExpr1, check.names=FALSE )
        
        if(add.pq) {
            ## add extra columns
            ##res <- cbind( res, q=mx$q, p=mx$p)
            colnames(mx.q) <- paste0("q.",colnames(mx.q))
            res <- cbind( res, mx.q)
        }
        rownames(res) = rownames(mx)
        return(res)
    }
    
    diffExprTable <- reactive({
        ## return the full DE table 
        ngs = inputData()
        if(is.null(ngs)) return(NULL)
        comp=1;test="trend.limma"
        comp = input$gx_contrast
        gx.methods = colnames(ngs$gx.meta$meta[[1]]$fc)
        gx.methods
        test = input$gx_testmethod
        if(is.null(comp)) return(NULL)
        if(is.null(test)) return(NULL)

        res = getDEGtable(ngs, testmethods=test, comparison=comp, add.pq=TRUE)
        return(res)
    })

    sigDiffExprTable <- reactive({
        ##
        ## DE table filtered by FDR and gene family
        ##
        ##
        dbg("[ExpressionBoard::sigDiffExprTable] reacted")

        ngs = inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs,input$gx_features,input$gx_fdr,input$gx_lfc)
        
        gx_features=1
        gx_features = input$gx_features
        ##if(is.null(gx_features)) return(NULL)
        ##res = getDEGtable(ngs, testmethods="trend.limma", comparison=1,add.pq=FALSE)
        res = diffExprTable()
        if(is.null(res)) return(NULL)

        fdr=1
        fdr = as.numeric(input$gx_fdr)    
        qv.col = grep("qval|adj.p|padj|fdr|meta.q",colnames(res),ignore.case=TRUE)[1]
        fx.col = grep("mean.diff|logfc|foldchange|meta.fx",colnames(res),ignore.case=TRUE)[1]        
        pval = res[,qv.col]

        if(is.null(ngs$families)) stop("FATAL:: no families in object")
        psel <- rownames(res)
        ##psel = filterProbes(ngs$genes, ngs$families[[gx_features]] )
        
        if(gx_features!="<all>") psel = filterProbes(ngs$genes, GSETS[[gx_features]] )
        res = res[which(pval <= fdr & rownames(res) %in% psel),,drop=FALSE]
        dim(res)
        
        ## filter on fold-change    
        fx.col = grep("fx|fc|sign|NES|logFC",colnames(res))[1]
        lfc = 1
        lfc = as.numeric(input$gx_lfc)
        fx  = as.numeric(res[,fx.col])
        names(fx) = rownames(res)
        res = res[which(abs(fx) >= lfc),,drop=FALSE]
        
        if(nrow(res)==0) {
            validate(need(nrow(res) > 0, "warning. no genes passed current filters."))
            return(NULL)
        }
        res = res[order(-abs(res[,fx.col])),]
        
        ## just show top 10
        if(length(input$gx_top10) && input$gx_top10==TRUE) {
            fx  = as.numeric(res[,fx.col])
            names(fx) = rownames(res)
            pp <- unique(c(head(names(sort(-fx[which(fx>0)])),10),
                           head(names(sort(fx[which(fx<0)])),10)))
            res = res[pp,,drop=FALSE]
            res = res[order(-res[,fx.col]),,drop=FALSE]
        }

        ## limit number of rows???
        ## res <- head(res, 1000)

        dbg("[ExpressionBoard::sigDiffExprTable] done!")
        
        return(res)
    })

    ##================================================================================
    ## Plots 
    ##================================================================================

    ## ------------------  Info messsages
    expr_plots_volcano_text = "A volcano plot of genes for the selected comparison under the <code>Contrast</code> settings plotting fold-change versus significance on the x and y axes, respectively."
    expr_plots_maplot_text = "An application of a Bland-Altman (MA) plot of genes for the selected comparison under the <code>Contrast</code> settings plotting mean intensity versus fold-change on the x and y axes, respectively."
    expr_plots_topgenesbarplot_text = "The top N = {12} differentially (both positively and negatively) expressed gene barplot for the selected comparison under the <code>Contrast</code> settings."
    expr_plots_topfoldchange_text = "The fold change summary barplot across all contrasts for a gene that is selected from the differential expression analysis table under the <code>Table</code> section."
    
    expr_plots_volcano.RENDER %<a-% reactive({
        
        comp1=1;fdr=0.10
        comp1 = input$gx_contrast
        if(length(comp1)==0) return(NULL)
        if(is.null(input$gx_features)) return(NULL)
        ngs <- inputData()
        req(ngs)
        
        fdr = 1
        fdr = as.numeric(input$gx_fdr)
        
        res = diffExprTable()
        if(is.null(res)) return(NULL)
        lfc=1
        lfc = as.numeric(input$gx_lfc)

        fam.genes = unique(unlist(GSETS[8]))
        fam.genes = res$gene_name
        ##fam.genes = unique(unlist(ngs$families[input$gx_features]))
        if(input$gx_features!="<all>") fam.genes = unique(unlist(GSETS[input$gx_features]))
        jj <- match(toupper(fam.genes),toupper(res$gene_name))
        sel.genes <- res$gene_name[setdiff(jj,NA)]
        
        fc.genes = as.character(res[,grep("^gene$|gene_name",colnames(res))])
        qval = res[,grep("adj.P.Val|meta.q|qval|padj",colnames(res))[1]]
        qval = pmax(qval, 1e-100)
        fx = res[,grep("logFC|meta.fx|fc",colnames(res))[1]]
        sig.genes = fc.genes[which(qval <= fdr & abs(fx) > lfc )]
        sel.genes = intersect(sig.genes, sel.genes)
        xlim = c(-1,1)*max(abs(fx),na.rm=TRUE)
        ylim = c(0, max(12, 1.1*max(-log10(qval),na.rm=TRUE)))

        par(mfrow=c(1,1), mar=c(4,3,2,1.5), mgp=c(2,0.8,0), oma=c(1,0,0.5,0))
        par(mfrow=c(1,1), mar=c(4,3,1,1.5), mgp=c(2,0.8,0), oma=c(0,0,0,0))
        gx.volcanoPlot.XY( x=fx, pv=qval, gene=fc.genes,
                          render="canvas", n=5000, nlab=10, 
                          xlim=xlim, ylim=ylim, ## hi.col="#222222",
                          use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                          cex=0.9, lab.cex=1.4, cex.main=1.0,
                          xlab="effect size (log2.FC)",
                          ylab="significance (-log10.q)",
                          ## main="Volcano plot",
                          highlight=sel.genes )


    })

    expr_plots_volcano.PLOTLY <- reactive({
        
        dbg("[expr_plots_volcano.PLOTLY] reacted")

        comp1=1;fdr=0.10
        comp1 = input$gx_contrast
        if(length(comp1)==0) return(NULL)
        if(is.null(input$gx_features)) return(NULL)
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        req(ngs)

        dbg("[expr_plots_volcano.PLOTLY] 1")
                
        fdr = 1
        fdr = as.numeric(input$gx_fdr)
        
        res = diffExprTable()
        if(is.null(res)) return(NULL)
        lfc=1
        lfc = as.numeric(input$gx_lfc)

        fam.genes = unique(unlist(GSETS[8]))
        fam.genes = res$gene_name
        ##fam.genes = unique(unlist(ngs$families[input$gx_features]))
        if(input$gx_features!="<all>") fam.genes = unique(unlist(GSETS[input$gx_features]))
        jj <- match(toupper(fam.genes),toupper(res$gene_name))
        sel.genes <- res$gene_name[setdiff(jj,NA)]
        
        fc.genes = as.character(res[,grep("^gene$|gene_name",colnames(res))])
        qval = res[,grep("adj.P.Val|meta.q|qval|padj",colnames(res))[1]]
        qval = pmax(qval, 1e-20)
        x = res[,grep("logFC|meta.fx|fc",colnames(res))[1]]
        y <- -log10(qval + 1e-12)
        
        sig.genes = fc.genes[which(qval <= fdr & abs(x) > lfc )]
        sel.genes = intersect(sig.genes, sel.genes)
        scaled.x <- scale(x,center=FALSE)
        scaled.y <- scale(y,center=FALSE)
        
        impt <- function(g) {
            j = match(g, fc.genes)
            x1 = scaled.x[j]
            y1 = scaled.y[j]
            x = sign(x1)*(x1**2 + 0.25*y1**2)
            names(x)=g
            x
        }
        lab.genes = c( head(sel.genes[order(impt(sel.genes))],10),
                      head(sel.genes[order(-impt(sel.genes))],10) )

        xlim = c(-1,1)*max(abs(x),na.rm=TRUE)
        ylim = c(0, max(12, 1.1*max(-log10(qval),na.rm=TRUE)))
        
        ## par(mfrow=c(1,1), mar=c(4,3,1,1.5), mgp=c(2,0.8,0), oma=c(0,0,0,0))
        plt <- plotlyVolcano(
            x=x, y=y, names=fc.genes,
            source = "plot1",
            highlight = sel.genes, label=lab.genes,
            group.names = c("group1","group0"),
            ##xlim=xlim, ylim=ylim, ## hi.col="#222222",
            ##use.fdr=TRUE,
            psig = fdr, lfc = lfc,
            xlab = "effect size (log2FC)",
            ylab = "significance (-log10q)",
            marker.size = 4,
            displayModeBar = FALSE,
            showlegend = FALSE) %>%
            layout( margin = list(b=65) )

        dbg("[expr_plots_volcano.PLOTLY] done!")
        return(plt)
    })
    
    ##expr_plots_volcano_module <- plotModule(
    callModule(
        plotModule,
        id = "expr_plots_volcano", ## ns=ns,
        ##func = expr_plots_volcano.RENDER,
        func = expr_plots_volcano.PLOTLY, plotlib = "plotly",       
        ##func2 = expr_plots_volcano.RENDER, 
        info.text = expr_plots_volcano_text, 
        title = "Volcano plot", label="a",
        height = imgH,  res=75,
        pdf.width=6, pdf.height=6
    )
    ## output <- attachModule(output, expr_plots_volcano_module)

    expr_plots_maplot.RENDER %<a-% reactive({
        comp1 = input$gx_contrast
        if(length(comp1)==0) return(NULL)

        ngs <- inputData()
        req(ngs)
        
        fdr=1;lfc=1
        fdr = as.numeric(input$gx_fdr)    
        lfc = as.numeric(input$gx_lfc)
        
        res = diffExprTable()
        if(is.null(res)) return(NULL)
        fc.genes = as.character(res[,grep("^gene$|gene_name",colnames(res))])
        ##pval = res$P.Value
        ##pval = res[,grep("P.Value|meta.p|pval|p.val",colnames(res))[1]]

        ## filter genes by gene family or gene set
        fam.genes = unique(unlist(ngs$families[10]))
        ##fam.genes = unique(unlist(ngs$families[input$gx_features]))
        fam.genes = res$gene_name
        if(input$gx_features!="<all>") fam.genes = unique(unlist(GSETS[input$gx_features]))
        jj <- match(toupper(fam.genes),toupper(res$gene_name))
        sel.genes <- res$gene_name[setdiff(jj,NA)]

        qval = res[,grep("adj.P.Val|meta.q|qval|padj",colnames(res))[1]]
        fx = res[,grep("logFC|meta.fx|fc",colnames(res))[1]]
        sig.genes = fc.genes[which(qval <= fdr & abs(fx) > lfc )]
        sel.genes = intersect(sig.genes, sel.genes)    
        xlim = c(-1,1)*max(abs(fx),na.rm=TRUE)
        ma = rowMeans(ngs$X[rownames(res),], na.rm=TRUE)

        par(mfrow=c(1,1), mar=c(4,3,2,1.5), mgp=c(2,0.8,0), oma=c(1,0,0.5,0))
        par(mfrow=c(1,1), mar=c(4,3,1,1.5), mgp=c(2,0.8,0), oma=c(0,0,0,0))
        gx.volcanoPlot.XY( x=fx, pv=qval, gene=fc.genes, lfc=lfc,
                          render="canvas", n=5000, nlab=10, 
                          xlim=xlim, ylim=c(0,15),
                          xlab="average expression (log2CPM)",
                          ylab="effect size (log2FC)",
                          ma_plot=TRUE, ma = ma, ## hi.col="#222222",
                          use.fdr=TRUE, p.sig=fdr, ##main=comp1,
                          ## main="MA plot",
                          cex=0.9, lab.cex=1.4, cex.main=1.0,
                          highlight=sel.genes)
    })

    expr_plots_maplot.PLOTLY %<a-% reactive({
        comp1 = input$gx_contrast
        if(length(comp1)==0) return(NULL)

        ngs <- inputData()
        req(ngs)

        dbg("[expr_plots_maplot.PLOTLY] reacted")
        
        fdr=1;lfc=1
        fdr = as.numeric(input$gx_fdr)    
        lfc = as.numeric(input$gx_lfc)
        
        res = diffExprTable()
        if(is.null(res)) return(NULL)
        fc.genes = as.character(res[,grep("^gene$|gene_name",colnames(res))])
        ##pval = res$P.Value
        ##pval = res[,grep("P.Value|meta.p|pval|p.val",colnames(res))[1]]

        ## filter genes by gene family or gene set
        fam.genes = unique(unlist(ngs$families[10]))
        ##fam.genes = unique(unlist(ngs$families[input$gx_features]))
        fam.genes = res$gene_name
        if(input$gx_features!="<all>") fam.genes = unique(unlist(GSETS[input$gx_features]))
        jj <- match(toupper(fam.genes),toupper(res$gene_name))
        sel.genes <- res$gene_name[setdiff(jj,NA)]

        qval = res[,grep("adj.P.Val|meta.q|qval|padj",colnames(res))[1]]
        y = res[,grep("logFC|meta.fx|fc",colnames(res))[1]]
        sig.genes = fc.genes[which(qval <= fdr & abs(y) > lfc )]
        sel.genes = intersect(sig.genes, sel.genes)    
        ylim = c(-1,1)*max(abs(y),na.rm=TRUE)
        x = rowMeans( ngs$X[rownames(res),], na.rm=TRUE)

        impt <- function(g) {
            j = match(g, fc.genes)
            x1 = scale(x,center=FALSE)[j]
            y1 = scale(y,center=FALSE)[j]
            x = sign(y1)*(1.0*x1**2 + 1.0*y1**2)
            names(x)=g
            x
        }
        lab.genes = c( head(sel.genes[order(impt(sel.genes))],10),
                      head(sel.genes[order(-impt(sel.genes))],10) )
        
        highlight=sel.genes;label=lab.genes;names=fc.genes
        plt <- plotlyMA(
            x=x, y=y, names=fc.genes,
            source = "plot1",
            highlight=sel.genes, label=lab.genes,
            group.names = c("group1","group0"),
            ##xlim=xlim, ylim=ylim, ## hi.col="#222222",
            ##use.fdr=TRUE,
            psig = fdr, lfc = lfc,
            xlab = "average expression (log2.CPM)",
            ylab = "effect size (log2.FC)",
            marker.size = 4,
            displayModeBar = FALSE,
            showlegend = FALSE) %>%
            layout( margin = list(b=65) )

        dbg("[expr_plots_maplot.PLOTLY] done!")
        
        return(plt)
    })
    
    callModule( plotModule,
        id="expr_plots_maplot", 
        ##func = expr_plots_maplot.RENDER,
        func = expr_plots_maplot.PLOTLY, plotlib="plotly",
        ##func2 = expr_plots_maplot.RENDER, 
        info.text = expr_plots_maplot_text, label="b",
        title = "MA plot",
        height = imgH,
        pdf.width=6, pdf.height=6, res=75
    )
    ## output <- attachModule(output, expr_plots_maplot_module)

    
    expr_plots_topgenesbarplot.RENDER %<a-% reactive({
        require(RColorBrewer)
        ngs = inputData()
        req(ngs)
        comp1 = input$gx_contrast

        dbg("expr_plots_topgenesbarplot.RENDER: reacted")
        
        if(length(comp1)==0) return(NULL)
        
        ## get table
        ##sel.row=1;pp=rownames(ngs$X)[1]
        ##sel.row = input$expr_genetable_rows_selected
        
        res = sigDiffExprTable()        
        if(is.null(res)) return(NULL)

        ##fc <- res$meta.fx
        fc <- res$logFC
        names(fc) <- rownames(res)
        top.up <- head(names(sort(fc[which(fc>0)],decreasing=TRUE)),12)
        top.dn <- head(names(sort(fc[which(fc<0)],decreasing=FALSE)),12)
        fc.top <- c(fc[top.up], fc[top.dn])
        klr.pal <- brewer.pal(4,"Paired")[2:1]
        klr <- c( rep(klr.pal[1],length(top.up)), rep(klr.pal[2],length(top.dn)) )
        names(fc.top) <- sub(".*:","",names(fc.top))
        
        ii <- order(fc.top)
        par(mfrow=c(1,1), mar=c(4,4,2,2)*1, mgp=c(2,0.8,0), oma=c(1,1,1,0.5)*0.2)
        par(mfrow=c(1,1), mar=c(5,3,1,1), mgp=c(2,0.8,0), oma=c(0,0,0,0))
        barplot(fc.top[ii], las=3, cex.names=0.75, ylab="fold change",
                col=klr[ii], ylim=c(-1.1,1.2)*max(abs(fc.top),na.rm=TRUE) )
        
        ## warning A_vs_B or B_vs_A not checked!!!
        groups <- strsplit(comp1,split="[._ ]vs[._ ]")[[1]]
        if(is.POSvsNEG(ngs)) groups <- rev(groups)
        tt <- c( paste("up in",groups[2]), paste("up in",groups[1]) )
        ##tt <- c( paste("up in",groups[1]), paste("down in",groups[1]) )
        legend("topleft", legend=tt, fill=klr.pal, cex=0.9, y.intersp=0.85, bty="n")
        ##title("top DE genes",cex.main=1)
        
        dbg("expr_plots_topgenesbarplot.RENDER: done\n")
        
    })

    callModule(
        plotModule,
        id="expr_plots_topgenesbarplot", ## ns=ns,
        func = expr_plots_topgenesbarplot.RENDER,
        func2 = expr_plots_topgenesbarplot.RENDER,        
        info.text = expr_plots_topgenesbarplot_text, label="c",
        title = "top DE genes",
        height = c(imgH,500), width=c('auto',800),
        pdf.width=6, pdf.height=6, res=75
    )
    ## output <- attachModule(output, expr_plots_topgenesbarplot_module)

    expr_plots_topfoldchange.RENDER %<a-% reactive({
        require(RColorBrewer)
        ngs = inputData()
        req(ngs)
        
        ## get table
        ##sel=1;pp=rownames(ngs$X)[1]
        ##sel = input$expr_genetable_rows_selected
        sel = expr_genetable$rows_selected()
        if(is.null(sel)) return(NULL)    

        res = sigDiffExprTable()
        if(is.null(res) || is.null(sel)) return(NULL)
        psel <- rownames(res)[sel]
        gene <- ngs$genes[psel,"gene_name"]
        
        ##fc <- res$meta.fx
        comp=1
        comp = input$gx_contrast
        if(is.null(comp) || length(comp)==0) return(NULL)
        fc <- sapply( ngs$gx.meta$meta, function(x) x[psel,"meta.fx"])

        top.up <- head(names(sort(fc[which(fc>0)],decreasing=TRUE)),10)
        top.dn <- head(names(sort(fc[which(fc<0)],decreasing=FALSE)),10)
        fc.top <- c(fc[top.up], fc[top.dn])
        fc.top <- fc.top[head(order(-abs(fc.top)),15)]
        fc.top <- sort(fc.top)
        fc.top <- head(c(fc.top, rep(NA,99)),15)
        klr.pal <- brewer.pal(4,"Paired")[2:1]
        ##klr.pal <- BLUERED(16)[c(3,14)]
        klr <- klr.pal[1 + 1*(sign(fc.top)<0)]
        
        par(mfrow=c(1,1), mar=c(4,4,2,2)*1, mgp=c(2,0.8,0), oma=c(1,1,1,0.5)*0.2)
        par(mfrow=c(1,1), mar=c(6,3,0,1), mgp=c(2,0.8,0), oma=c(1,0,0,0))
        nch <- max(nchar(names(fc.top)))
        m1 <- ifelse(nch > 12, 12, 8)
        m1 <- ifelse(nch > 30, 16, m1)

        ##par( mar=c(4,m1,2,0.5) )
        par( mar=c(3.2,m1-0.5,1,1) )
        cex1 <- 0.9
        nn <- sum(!is.na(fc.top))
        if(nn>15) cex1 <- 0.8
        barplot(fc.top, col=klr, horiz=TRUE, las=1,
                xlim=c(-1,1)*max(abs(fc.top),na.rm=TRUE),
                cex.names=cex1, xlab="fold change (log2)")
        title(gene, cex.main=1, line=-0.15)
        
    })

    callModule( plotModule,
        id = "expr_plots_topfoldchange", 
        func = expr_plots_topfoldchange.RENDER,
        func2 = expr_plots_topfoldchange.RENDER,
        info.text = expr_plots_topfoldchange_text,
        title = "Gene in contrasts", label = "d",
        height = c(imgH,500), width=c('auto',700),
        pdf.width=6, pdf.height=6, res=74
    )
    ##output <- attachModule(output, expr_plots_topfoldchange_module)
    
    expr_plots_boxplot.RENDER %<a-% reactive({
        require(RColorBrewer)
        ngs = inputData()
        req(ngs)
        
        ## get table
        ##sel=1
        sel = expr_genetable$rows_selected()
        if(is.null(sel)) return(NULL)    

        res = sigDiffExprTable()
        if(is.null(res) || is.null(sel)) return(NULL)

        psel <- rownames(res)[sel]
        ##gene = sub(".*:","",rownames(res)[sel])
        gene=ngs$genes[1,"gene_name"];comp=1;grouped=TRUE;logscale=TRUE;srt=45
        gene = ngs$genes[psel,"gene_name"]
        comp = input$gx_contrast
        grouped <- !input$gx_ungroup
        logscale <- input$gx_logscale
        srt <- ifelse(grouped, 0, 30)

        par(mfrow=c(1,1), mar=c(6,4,3,0)*1, mgp=c(2,0.8,0), oma=c(1,1,1,0.5)*0.2)
        par(mfrow=c(1,1), mar=c(4,3,1,1.5), mgp=c(2,0.8,0), oma=c(1,0,0,0))
        pgx.plotExpression(ngs, gene, comp=comp, grouped=grouped,
                           max.points=2000, names=TRUE,
                           logscale=logscale, srt=srt)    
    })

    callModule( plotModule,
        id = "expr_plots_boxplot", 
        func = expr_plots_boxplot.RENDER,
        func2 = expr_plots_boxplot.RENDER,
        info.text = "This is an example",
        height = imgH,
        pdf.width=6, pdf.height=6, res=75
    )

    expr_plots_caption = "<b>Expression plots</b> associated with the selected contrast. <b>(a)</b> Volcano-plot plotting fold-change versuson significance the x and y axes, respectively. <b>(b)</b> MA-plot plotting signal intensity versus fold-change on the x and y axes, respectively. <b>(c)</b> Sorted barplot of the top diffentially expressed genes with largest (absolute) fold-change for selected contrast. <b>(d)</b> Sorted barplot of the differential expression of the selected gene across all contrasts."

    ## library(shinyjqui)
    output$expr_plots_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1,0.25,NA),
            fillRow(
                id = "expr_plots",
                ##height = rowH,
                flex=c(1,1,1,1), ##height = 370,
                plotWidget(ns("expr_plots_volcano")),
                plotWidget(ns("expr_plots_maplot")),
                plotWidget(ns("expr_plots_topgenesbarplot")),
                plotWidget(ns("expr_plots_topfoldchange"))
            ),
            br(),
            div(HTML(expr_plots_caption), class="caption")
        )
    })
    ##dragula(ns("expr_plots"))
    
    ##================================================================================
    ## Top genes
    ##================================================================================

    expr_topgenes.RENDER %<a-% reactive({

        ngs <- inputData()
        req(ngs)
        
        ##res=getDEGtable(ngs, testmethods="trend.limma", comparison=6, add.pq=FALSE)
        res <- sigDiffExprTable()
        if(is.null(res) || nrow(res)==0) return(NULL)

        ## filter on active rows (using search)
        ##ii <- input$expr_genetable_rows_all
        ii  <- expr_genetable$rows_all()
        res <- res[ii,,drop=FALSE]
        if(nrow(res)==0) return(NULL)
        
        fx.col = grep("fc|fx|mean.diff|logfc|foldchange",tolower(colnames(res)))[1]
        fx.col
        fx = res[,fx.col]
        names(fx) <- rownames(res)
        top.up=top.down=rownames(ngs$X)
        top.up   <- names(sort(fx[fx>0],decreasing=TRUE))
        top.down <- names(sort(fx[fx<0]))
        head(top.up)
        head(top.down)
        
        y <- ngs$samples$group
        ngrp = length(unique(y))
        las = ifelse(ngrp>3, 3, 0)
        mar1 = ifelse(las==3, 4, 3)
        mar1 = ifelse(las==3, 3.5, 2.5)
        
        comp=1;grouped=0;logscale=1
        comp = input$gx_contrast
        grouped <- !input$gx_ungroup
        logscale <- input$gx_logscale
        showothers <- input$gx_showothers
        
        srt = 30
        mar1 = 3.5
        ylab = ifelse(logscale, "log2CPM", "CPM")
        show.names <- ifelse(!grouped & ngrp>25, FALSE, TRUE)
        ##nx = ifelse(grouped, ngrp, length(y))
        nx = ifelse(grouped, 3, length(y))
        nc = 4
        nc = 8
        if( nx <= 3) nc <- 10
        if( nx > 10) nc <- 5
        if( nx > 25) nc <- 4

        ##nc <- 10
        par(mfrow=c(2,nc), mar=c(mar1,3.5,1,1), mgp=c(2,0.8,0), oma=c(0.1,0.6,0,0.6) )
        i=1
        for(i in 1:nc) {

            if(i > length(top.up)) { frame() }
            gene = sub(".*:","",top.up[i])
            pgx.plotExpression(
                ngs, gene, comp=comp, grouped=grouped, max.points=1000,
                collapse.others=TRUE, showothers=showothers,
                ylab = ylab, xlab="",
                logscale=logscale, names=show.names, srt=srt, main="")
            title( gene, cex.main=1, line=-0.6)

        }

        for(i in 1:nc) {
            if(i > length(top.down)) { frame() }
            gene = sub(".*:","",top.down[i])
            pgx.plotExpression(
                ngs, gene, comp=comp, grouped=grouped,  max.points=1000,
                collapse.others=TRUE, showothers=showothers,
                ylab = ylab, xlab="",
                logscale=logscale, names=show.names, srt=srt, main="")
            title( gene, cex.main=1, line=-0.6)
            ##qv1 = formatC(qv[gs],format="e", digits=2)
            ##legend("topright", paste("q=",qv1), bty="n",cex=1)
        }
        
    })

    expr_topgenes_opts = tagList(
        tipify( checkboxInput(ns('gx_logscale'),'log scale',TRUE),
               "Logarithmic scale the counts (abundance levels).",
               placement="right", options = list(container = "body")),
        tipify( checkboxInput(ns('gx_ungroup'),'ungroup samples',FALSE),
               "Ungroup samples in the plot",
               placement="right", options = list(container = "body")),
        tipify( checkboxInput(ns('gx_showothers'),'show others',TRUE),
               "Show the 'others' class (if any)",
               placement="right", options = list(container = "body"))
        )

    expr_topgenes_text = "The <strong>Top genes</strong> section shows the average expression plots across the samples for the top differentially (both positively and negatively) expressed genes for the selected comparison from the <code>Contrast</code> settings. Under the plot <i>Settings</i>, users can scale the abundance levels (counts) or ungroup the samples in the plot from the <code>log scale</code> and <code>ungroup samples</code> settings, respectively."

    expr_topgenes_caption = "<b>Top differentially expressed genes.</b> Expression barplots of the top most differentially (both positively and negatively) expressed genes for the selected contrast."

    callModule( plotModule,
        id = "expr_topgenes", 
        func = expr_topgenes.RENDER,
        func2 = expr_topgenes.RENDER,
        options = expr_topgenes_opts,
        info.text = expr_topgenes_text,
        ##caption = expr_topgenes_caption,
        height = c(imgH,450), width = c('auto',1500), res=c(85,95),
        pdf.width=14, pdf.height=4, 
        title="Expression of top differentially expressed genes"
    )
    ##output <- attachModule(output, expr_topgenes_module)

    ## library(shinyjqui)
    output$expr_topgenesUI <- renderUI({
        fillCol(
            ## id = ns("expr_topgenes"),
            height = rowH,
            flex=c(1,NA,NA), ##height = 370,
            plotWidget(ns("expr_topgenes")),
            br(),
            div(HTML(expr_topgenes_caption),class="caption")
        )
    })
    outputOptions(output, "expr_topgenesUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ## Volcano (all contrasts)
    ##================================================================================

    expr_volcanoAll.RENDER %<a-% reactive({

        ngs = inputData()
        if( is.null(ngs)) return(NULL)
        
        comp = names(ngs$gx.meta$meta)
        if(length(comp)==0) return(NULL)
        if(is.null(input$gx_features)) return(NULL)
        
        fdr = 1
        fdr = as.numeric(input$gx_fdr)
        lfc = as.numeric(input$gx_lfc)
        
        sel.genes = GSETS[["<all>"]]
        sel.genes = rownames(ngs$X)
        if(input$gx_features!="<all>") sel.genes = unique(unlist(GSETS[input$gx_features]))

        ncomp <- length(comp)
        nr <- ceiling(sqrt(ncomp/3))
        NC <- 3*nr
        if(ncomp <= (nr-1)*NC) nr <- (nr-1)
        par(mfrow=c(nr,NC), mar=c(4,4,2,2)*0)

        ng = length(comp)
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

        
        ##comp <- head(comp,75)  ## maximum 75!!!
        i=1
        F <- list()
        Q <- list()
        for(i in 1:length(comp)) {
            
            test = GX.DEFAULTTEST
            test = input$gx_testmethod
            if(is.null(test)) return(NULL)
            res = getDEGtable(ngs, testmethods=test, comparison=i, add.pq=TRUE)
            
            fc.gene = res[,grep("^gene$|^gene_name$",colnames(res))]
            ##pv.col = grep("p.val|pval|meta.p",colnames(res),ignore.case=TRUE)[1]
            qv.col = grep("qval|adj.p|padj|fdr|meta.q",colnames(res),ignore.case=TRUE)[1]
            fx.col = grep("mean.diff|logfc|foldchange|meta.fx",colnames(res),ignore.case=TRUE)[1]
            qval = res[,qv.col]
            fx   = res[,fx.col]
            names(qval) <- names(fx) <- fc.gene
            F[[i]] <- fx
            Q[[i]] <- qval
        }

        ymax=15
        ymax <- 1.2 * max(-log10(1e-99 + unlist(Q)), na.rm=TRUE)        
        
        withProgress(message="computing volcano plots ...", value=0, {
            i=1
            for(i in 1:length(comp)) {
                qval <- Q[[i]]
                fx   <- F[[i]]
                fc.gene <- names(qval)
                sig.genes = fc.gene[which(qval <= fdr & abs(fx) >= lfc)]
                ##genes1 = intersect(sig.genes, sel.genes)
                genes1 = sig.genes[which(toupper(sig.genes) %in% toupper(sel.genes))]
                gx.volcanoPlot.XY( x=fx, pv=qval, gene=fc.gene,
                                  render="canvas", n=1000, nlab=5, 
                                  xlim=NULL, ylim=c(0,ymax), axes=FALSE, 
                                  use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                                  ##main=comp[i], 
                                  ## ma.plot=TRUE, use.rpkm=TRUE,
                                  cex=0.6, lab.cex=1.5, highlight=genes1)

                is.first = (i%%NC==1)
                last.row = ( (i-1)%/%NC == (length(comp)-1)%/%NC )
                if(is.first) axis(2, tcl=0.5, mgp=c(-2,-1.5,0))
                if(last.row) axis(1, tcl=0.5, mgp=c(-2,-1.5,0))
                box()
                legend("topright", legend=comp[i], cex=0.95,bg="white")

                incProgress( 1/length(comp) )
            }
        })  ## progress

    })

    expr_volcanoAll_text = "Under the <strong>Volcano (all)</strong> tab, the platform simultaneously displays multiple volcano plots for genes across all contrasts. This provides users an overview of the statistics for all comparisons. By comparing multiple volcano plots, the user can immediately see which comparison is statistically weak or strong."

    expr_volcanoAll_caption = "<b>Volcano plot for all contrasts.</b> Simultaneous visualisation of volcano plots of genes for all contrasts. Experimental contrasts with better statistical significance will show volcano plots with 'higher' wings."

    callModule( plotModule,
        id="expr_volcanoAll", 
        func = expr_volcanoAll.RENDER,
        func2 = expr_volcanoAll.RENDER,
        info.text = expr_volcanoAll_text,
        ##caption = expr_volcanoAll_caption,
        pdf.width=18, pdf.height=6,
        ##height = imgH, res=75,
        height = c(imgH,500), width = c('auto',1500), res=c(75,90),
        title="Volcano plots for all contrasts"
    )
    ##output <- attachModule(output, expr_volcanoAll_module)

    ## library(shinyjqui)
    output$expr_volcanoAll_UI <- renderUI({
        fillCol(
            ## id = ns("expr_topgenes"),
            height = rowH,
            flex=c(1,NA,NA), ##height = 370,
            plotWidget(ns("expr_volcanoAll")),
            br(),
            div(HTML(expr_volcanoAll_caption), class="caption")
        )
    })

    ##================================================================================
    ## Volcano (all methods)
    ##================================================================================

    expr_volcanoMethods.RENDER %<a-% reactive({
        comp = input$gx_contrast
        if(is.null(comp)) return(NULL)
        ngs = inputData()
        req(ngs)
        if(is.null(input$gx_features)) return(NULL)
        
        fdr = 1
        fdr = as.numeric(input$gx_fdr)
        lfc = as.numeric(input$gx_lfc)
        genes = NULL

        sel.genes = head(ngs$genes$gene_name,100)
        ##sel.genes = unique(unlist(ngs$families[input$gx_features]))
        sel.genes = unique(unlist(GSETS[input$gx_features]))
        
        ##methods = names(ngs$gx.meta$output)
        methods = colnames(ngs$gx.meta$meta[[1]]$fc)
        par(mfrow=c(2,6), mar=c(4,4,2,2)*0)
        if(length(methods)>12) {
            par(mfrow=c(3,8), mar=c(4,4,2,2)*0)
        }
        
        ## meta tables
        mx = ngs$gx.meta$meta[[comp]]
        fc = unclass(mx$fc)
        ##pv = unclass(mx$p)
        qv = unclass(mx$q)
        ## fc = cbind( meta=mx[,"meta.fx"], fc)
        ## qv = cbind( meta=mx[,"meta.q"], qv)

        ymax <- 1.2 * max(-log10(1e-99 + qv), na.rm=TRUE) ## y-axis       
        
        xlim = c(-1.1,1.1)*max(abs(fc))
        xlim = 1.3*c(-1,1) * quantile(abs(fc),probs=0.999)
        fc.genes = ngs$genes[rownames(mx),"gene_name"]
        nplots <- min(24,ncol(qv))
        
        withProgress(message="computing volcano plots ...", value=0, {
            i=1
            for(i in 1:nplots) {
                fx = fc[,i]
                ## pval = pv[,i]
                qval = qv[,i]
                sig.genes = fc.genes[which(qval <= fdr & abs(fx) >= lfc)]
                ##genes1 = intersect(sig.genes, sel.genes)
                genes1 = sig.genes[which(toupper(sig.genes) %in% toupper(sel.genes))]
                gx.volcanoPlot.XY(
                    x=fx, pv=qval, gene=fc.genes,
                    render="canvas", n=5000, nlab=5, 
                    xlim=xlim, ylim=c(0,ymax), axes=FALSE, 
                    use.fdr=TRUE, p.sig=fdr, lfc=lfc,
                    ##main=comp[i], 
                    ## ma.plot=TRUE, use.rpkm=TRUE,
                    cex=0.6, lab.cex=1.5, highlight=genes1)
                axis(2, tcl=0.5, mgp=c(-2,-1.5,0))
                axis(1, tcl=0.5, mgp=c(-2,-1.5,0))
                box()
                legend("topright", legend=colnames(fc)[i], cex=1.2, bg="white")

                incProgress( 1/length(nplots) )                
            }
        })

    })

    expr_volcanoMethods_text = "Under the <strong>Volcano (methods)</strong> tab, the platform displays the volcano plots provided by multiple differential expression calculation methods for the selected contrast. This provides users an overview of the statistics of all methods at the same time."

    expr_volcanoMethods_caption = "<b>Volcano plot for all statistical methods.</b> Simultaneous visualisation of volcano plots of genes by multiple differential expression methods for the selected contrast. Methods showing better statistical significance will show volcano plots with 'higher' wings."

    callModule( plotModule,
        id = "expr_volcanoMethods", 
        func = expr_volcanoMethods.RENDER,
        func2 = expr_volcanoMethods.RENDER,        
        title = "Volcano plots for all methods",
        info.text = expr_volcanoMethods_text, 
        ##caption = expr_volcanoMethods_caption,
        height = c(imgH,500), width = c('auto',1500), res=c(75,90),
        pdf.width=18, pdf.height=6
    )
    ##output <- attachModule(output, expr_volcanoMethods_module)

    ## library(shinyjqui)
    output$expr_volcanoMethodsUI <- renderUI({
        fillCol(
            height = rowH,
            flex=c(1,NA,NA), ##height = 370,
            plotWidget(ns("expr_volcanoMethods")),
            br(),
            div(HTML(expr_volcanoMethods_caption), class="caption")
        )
    })

    ##================================================================================
    ## Statistics Table
    ##================================================================================

    expr_genetable.RENDER <- reactive({
        
        dbg("[ExpressionBoard::expr_genetable.RENDER] reacted")
        
        res <- sigDiffExprTable()
        ##req(res)
        dbg("[ExpressionBoard::expr_genetable.RENDER] dim(res)=",dim(res),"\n")
        
        if(is.null(res) || nrow(res)==0) return(NULL)
        
        fx.col = grep("fc|fx|mean.diff|logfc|foldchange",tolower(colnames(res)))[1]
        fx.col
        fx = res[,fx.col]
        
        if("gene_title" %in% colnames(res)) res$gene_title = shortstring(res$gene_title,50)
        rownames(res) = sub(".*:","",rownames(res))

        if(!DEV) {
            kk = grep("meta.fx|meta.fc|meta.p",colnames(res),invert=TRUE)
            res <- res[,kk,drop=FALSE]
        }
        if(!input$gx_showqvalues) {
            kk = grep("^q[.]",colnames(res),invert=TRUE)
            res <- res[,kk,drop=FALSE]
        }

        numeric.cols <- which(sapply(res, is.numeric))
        numeric.cols <- colnames(res)[numeric.cols]
        dbg("[ExpressionBoard::expr_genetable.RENDER] numeric.cols=",numeric.cols)
        dbg("[ExpressionBoard::expr_genetable.RENDER] done!")
        
        DT::datatable( res, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip',
                          pageLength = 400,
                          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle(colnames(res)[fx.col],
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(fx, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
                                        #}, server=FALSE)
    })

    expr_genetable_text = "Table <strong>I</strong> shows the results of the statistical tests. To increase the statistical reliability of the Omics Playground, we perform the DE analysis using four commonly accepted methods in the literature, namely, <a href='https://en.wikipedia.org/wiki/Student%27s_t-test'>t-test</a> (standard, Welch), <a href='https://www.ncbi.nlm.nih.gov/pubmed/25605792'> limma</a> (no trend, trend, voom), <a href='https://www.ncbi.nlm.nih.gov/pubmed/19910308'> edgeR</a> (QLF, LRT), and <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4302049'> DESeq2</a> (Wald, LRT), and merge the results. 
<br><br>For a selected comparison under the <code>Contrast</code> setting, the results of the selected methods are combined and reported under the table, where <code>meta.q</code> for a gene represents the highest <code>q</code> value among the methods and the number of stars for a gene indicate how many methods identified significant <code>q</code> values (<code>q < 0.05</code>). The table is interactive (scrollable, clickable); users can sort genes by <code>logFC</code>, <code>meta.q</code>, or average expression in either conditions. Users can filter top N = {10} differently expressed genes in the table by clicking the <code>top 10 genes</code> from the table <i>Settings</i>."

    expr_genetable_opts = tagList(
        tipify(checkboxInput(ns("gx_top10"),"top 10 genes",FALSE),
               "Display only top 10 differentially (positively and negatively) expressed genes in the table.", 
               placement="top", options = list(container = "body")),
        tipify(checkboxInput(ns('gx_showqvalues'),'show q-values',FALSE),
               "Show each q-values of all statistical methods in the table.", 
               placement="top", options = list(container = "body"))    
    )

    expr_genetable <- callModule(
        tableModule,
        id = "expr_genetable",
        func = expr_genetable.RENDER, 
        info.text = expr_genetable_text,
        label="I", info.width="500px",
        options = expr_genetable_opts,
        title = "Differential expression analysis",
        height = c(265,700)
    )

    ##output$expr_genetable <- expr_genetable_module$render
    ##output <- attachModule(output, expr_genetable_module)

    gx_related_genesets <- reactive({

        dbg("[gx_related_genesets] reacted")
        
        ngs <- inputData()
        res <- sigDiffExprTable()
        if(is.null(res) || nrow(res)==0) return(NULL)
        contr <- input$gx_contrast
        if(is.null(contr)) return(NULL)

        dbg("[gx_related_genesets] 1")
        
        ## get table
        sel.row=1;
        ##sel.row = input$expr_genetable_rows_selected
        sel.row = expr_genetable$rows_selected()
        if(is.null(sel.row)) return(NULL)
        gene0 <- rownames(res)[sel.row]
        gene1 <- toupper(sub(".*:","",gene0))  ## always uppercase...

        j <- which(toupper(rownames(ngs$GMT))==gene1)
        gset <- names(which(ngs$GMT[j,] != 0))
        gset <- intersect(gset, rownames(ngs$gsetX))
        if(length(gset)==0) return(NULL)
        
        fx <- ngs$gset.meta$meta[[contr]]$meta.fx 
        names(fx) <- rownames(ngs$gset.meta$meta[[contr]])
        fx  <- round(fx[gset],digits=4)
        
        rho <- cor(t(ngs$gsetX[gset,]), ngs$X[gene0,])[,1]        
        rho <- round(rho, digits=3)
        gset1 <- substring(gset,1,60)
        
        df <- data.frame(geneset=gset1, rho=rho, fx=fx, check.names=FALSE)
        rownames(df) <- gset       
        df <- df[order(-abs(df$fx)),]

        dbg("[gx_related_genesets] done!")

        return(df)
    })

    expr_gsettable.RENDER <- reactive({

        df <- gx_related_genesets()        
        if(is.null(df)) return(NULL)

        df$geneset <- wrapHyperLink(df$geneset, rownames(df))
        
        DT::datatable(df,
                      class = 'compact cell-border stripe',
                      rownames=FALSE, escape = c(-1,-2),
                      extensions = c('Scroller'),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip', 
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      ),  ## end of options.list 
                      selection=list(mode='single', target='row', selected=1)) %>%
            ##formatSignif(1:ncol(df),4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle("fx", background = color_from_middle( df$fx, 'lightblue', '#f5aeae'))
                                        #}, server=FALSE)
    })

    expr_gsettable_text = "By clicking on a gene in the Table <code>I</code>, it is possible to see which genesets contain that gene in this table, and check the differential expression status in other comparisons from the <code>Gene in contrasts</code> plot under the <code>Plots</code> tab."

    callModule(
        tableModule,
        id = "expr_gsettable", 
        func = expr_gsettable.RENDER, 
        info.text = expr_gsettable_text, label="II",
        title="Gene sets",
        height = c(265,700), width = c('100%',800)        
    )

    expr_tablesUI_caption = "<b>Differential expression tables</b>. <b>(I)</b> Statistical results of the the differential expression analysis for selected contrast. The number of stars indicate how many statistical methods identified the gene significant. <b>(II)</b> Correlation and enrichment value of gene sets that contain the gene selected in Table I."
    
    output$expr_tables_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1,NA),
            fillRow(
                flex = c(2,0.1,1), 
                tableWidget(ns("expr_genetable")),
                br(),
                tableWidget(ns("expr_gsettable"))
            ),
            div(HTML(expr_tablesUI_caption),class="caption")
        )
    })

    ##================================================================================
    ## Foldchange (all)
    ##================================================================================

    expr_fctable.RENDER <- reactive({
        
        ngs <- inputData()
        res <- sigDiffExprTable()
        if(is.null(res) || nrow(res)==0) return(NULL)

        ##F <- sapply(ngs$gx.meta$meta, function(x) unclass(x$fc)[,"trend.limma"])
        F <- sapply(ngs$gx.meta$meta, function(x) x$meta.fx)
        rownames(F) <- rownames(ngs$gx.meta$meta[[1]])
        colnames(F) <- gsub("_"," ",colnames(F))
        dim(F)
        fc.var = F[,1]**2
        if(NCOL(F)>1) {
            fc.var <- round( apply(F,1,var), digits=4)
        }

        F1 <- data.frame( gene=rownames(F), fc.var=fc.var, round(F,digits=3), check.names=FALSE)
        F1 <- F1[order(-F1$fc.var),]
        F1 <- F1[intersect(rownames(F1),rownames(res)),]  ## take intersection of current comparison
        
        DT::datatable( F1, rownames=FALSE,
                      class = 'compact cell-border stripe hover',
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=c(1)),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip', 
                          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle( "fc.var",
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle( fc.var, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  %>%
                DT::formatStyle( colnames(F),
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle(F, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    })

    expr_fctable_text = "The <strong>Foldchange (all)</strong> tab reports the gene fold changes for all contrasts in the selected dataset."

    expr_fctable_caption = "<b>Differential expression (fold-change) across all contrasts.</b> The column `fc.var` corresponds to the variance of the fold-change across all contrasts."

    callModule(
        tableModule,
        id="expr_fctable",
        func = expr_fctable.RENDER, 
        title ="Gene fold changes for all contrasts",
        info.text = expr_fctable_text,
        caption = expr_fctable_caption,
        height = c(280,700)
    )

    ## library(shinyjqui)
    output$expr_fctable_UI <- renderUI({
        fillCol(
            height = rowH,
            tableWidget(ns("expr_fctable"))
        )
    })

    ##================================================================================
    ## FDR table
    ##================================================================================
    
    expr_FDRtable.RENDER <- reactive({

        dbg("expr_FDRtable.RENDER:: reacted")
        
        methods = GX.DEFAULTTEST
        methods = input$gx_testmethod
        ##methods = input$gx_testmethod
        if(is.null(methods)) return(NULL)

        dbg("expr_FDRtable.RENDER:: 1")
        
        ##comp <- input$gx_contrast
        ngs = inputData()

        kk = rownames(ngs$gx.meta$sig.counts[[1]][[1]])
        kk = intersect(methods, rownames(ngs$gx.meta$sig.counts[[1]][[1]]))
        counts.up <- ngs$gx.meta$sig.counts$up
        counts.down <- ngs$gx.meta$sig.counts$down
        counts.up <- lapply(counts.up, function(x) x[kk,,drop=FALSE])
        counts.down <- lapply(counts.down, function(x) x[kk,,drop=FALSE])
        for(i in 1:length(counts.up)) {
            rownames(counts.up[[i]]) = paste0(names(counts.up)[i],"::",rownames(counts.up[[i]]))
            rownames(counts.down[[i]]) = paste0(names(counts.down)[i],"::",rownames(counts.down[[i]]))
        }
        sig.up = do.call(rbind, counts.up)
        sig.down = do.call(rbind, counts.down)

        dbg("expr_FDRtable.RENDER:: 2")
        
        sig.up   <- sig.up[order(rownames(sig.up)),,drop=FALSE]
        sig.down <- sig.down[order(rownames(sig.down)),,drop=FALSE]    
        colnames(sig.up)[1] = paste("UP   FDR = ",colnames(sig.up)[1])
        colnames(sig.down)[1] = paste("DOWN   FDR = ",colnames(sig.down)[1])
        colnames(sig.down) = paste0("  ",colnames(sig.down))
        sigcount = cbind( sig.down, sig.up[rownames(sig.down),,drop=FALSE] )
        dim(sigcount)    
        maxsig = 0.99 * max(sigcount,na.rm=TRUE)

        dbg("expr_FDRtable.RENDER:: 3")
        
        contr = sub("::.*","",rownames(sigcount))
        ##contr = rownames(sigcount)
        metd  = sub(".*::","",rownames(sigcount))
        D = data.frame( method=metd, contrast=contr, sigcount, check.names=FALSE)
        
        DT::datatable( D, rownames=FALSE,
                      class = 'compact cell-border stripe hover',
                      fillContainer = TRUE,
                      extensions = c('Scroller'),                      
                      options=list(
                          dom = 'lfrtip',
                          pageLength = 999, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle(colnames(sig.up),
                                background = styleColorBar(c(0,maxsig), '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')  %>%
                DT::formatStyle(colnames(sig.down),
                                background = styleColorBar(c(0,maxsig), 'lightblue'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    })

    expr_FDRtable_text = "The <strong>FDR table</strong> tab reports the number of significant genes at different FDR thresholds for all contrasts within the dataset."

    expr_FDRtable_caption = "<b>Number of significant genes versus FDR.</b> This table reports the number of significant genes at different FDR thresholds for all contrasts and methods. This enables to quickly see which methods are more sensitive. The left part of the table (in blue) correspond to the number of significant down-regulated genes, the right part (in red) correspond to the number of significant overexpressed genes."
    
    callModule(
        tableModule,
        id="expr_FDRtable",
        func = expr_FDRtable.RENDER, 
        info.text = expr_FDRtable_text,
        title = 'Number of significant genes',
        caption = expr_FDRtable_caption,
        height = c(280, 700)
    )
    ##output <- attachModule(output, expr_FDRtable_module)

    ## library(shinyjqui)
    output$expr_FDRtable_UI <- renderUI({
        fillCol(
            height = rowH,
            tableWidget(ns("expr_FDRtable"))
        )
    })
    
    ## reactive values to return to parent environment
    outx <- list(selected_gxmethods=selected_gxmethods)
    return(outx)
    
} ## end-of-ExpressionBoard
