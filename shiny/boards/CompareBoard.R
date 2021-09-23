##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2021 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing CompareBoard")
##source("global.R")

CompareInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

CompareUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs1"),
            shiny::tabPanel("Compare", shiny::uiOutput(ns("compareScatter_UI"))),
            shiny::tabPanel("Foldchange", shiny::uiOutput(ns("FCcorrelation_UI"))),
            shiny::tabPanel("Gene Correlation", shiny::uiOutput(ns("GeneCorrelation_UI")))            
        )
    )
}

CompareBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    ## fullH = session$clientdata$output_foo_height
    fullH = 770       # row height of panel
    ## fullH = input$dimensions[2]
    tabH = '70vh'
    
    ## reactive functions from shared environment
    inputData <- env[["load"]][["inputData"]]
    ## selected_gxmethods <- env[["expr"]][["selected_gxmethods"]]
    ## selected_gsetmethods <- env[["enrich"]][["selected_gsetmethods"]]
    
    description =
    "<h3>Compare Datasets</h3> Compare expression and signatures between two datasets, from similar experiments or from different datatypes, e.g. transcriptomics and proteomics."
    
    output$description <- shiny::renderUI(shiny::HTML(description))

    infotext =
        "The <strong>Compare Datasets</strong> module enables users to compare their dataset to other datasets. This module allows side-by-side comparison of volcano, scatter or gene t-SNE plots. It provides pairwise correlation plots and/or enrichment plots with signatures from other data sets.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
"
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================
    
    output$inputsUI <- shiny::renderUI({
        ui <- shiny::tagList(
            shinyBS::tipify( shiny::actionLink(ns("info"), "Info", icon = shiny::icon("info-circle")),
                   "Show more information about this module"),
            shiny::hr(), shiny::br(),             
            shinyBS::tipify( shiny::selectInput(ns('contrast1'),'Dataset1:',
                                choices=NULL, multiple=TRUE),
                   "Select the contrast that you want to compare.",
                   placement="right", options = list(container = "body")
                   ),
            shiny::br(),            
            shinyBS::tipify( shiny::selectInput(ns('dataset2'),"Dataset2:", choices=NULL),
                   "Select second dataset to compare.",
                   placement="right", options = list(container = "body")),
            shinyBS::tipify( shiny::selectInput(ns('contrast2'),NULL, choices=NULL, multiple=TRUE),
                   "Select second contrast to compare.",
                   placement="right", options = list(container = "body")),
            shiny::br(),
            shinyBS::tipify( shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib="glyphicon")),
                   "Toggle advanced options.",
                   placement="right", options = list(container = "body")),
            shiny::br(),
            shiny::conditionalPanel(
                "input.options % 2 == 1", ns=ns,
                shiny::br(),
                shinyBS::tipify( shiny::radioButtons(ns('plottype'),"Plot type:",
                                     choices=c("volcano","MA","scatter","UMAP1","UMAP2","heatmap"),
                                     inline=TRUE),
                       "Select plot type.",
                       placement="right", options = list(container = "body")),
                shiny::br(),
                shinyBS::tipify( shiny::radioButtons(ns('hilighttype'),"Highlight genes:",
                                     choices=c("top scoring","custom"),
                                     inline=TRUE),
                       "Select highlight type.",
                       placement="right", options = list(container = "body")),
                shiny::conditionalPanel(
                    "input.hilighttype == 'custom'", ns=ns,
                    shinyBS::tipify( shiny::textAreaInput(ns("genelist"),NULL, value = NULL,
                                          height = "100px", width = "100%", 
                                          rows=5, placeholder="Paste your custom gene list"),
                           "Paste a custom list of genes to highlight.",
                           placement="right")
                ),
                shiny::br(),            
                shinyBS::tipify(
                    shiny::radioButtons( ns('ntop'),"ntop", choices=c(10,20,40,100),
                                 selected=20, inline=TRUE),
                    "number of top genes to show",
                    placement="right", options = list(container = "body"))
            )
        )
        ui
    })
    shiny::outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Compare Experiments</strong>"),
            shiny::HTML(infotext),
            easyClose = TRUE, size="l" ))
    })

    shiny::observe({
        ngs <- inputData()

        comparisons1 <- names(ngs$gx.meta$meta)
        sel1 <- comparisons1[1]
        shiny::updateSelectInput(session, "contrast1", choices=comparisons1, selected=sel1)        

        pgx.files <- sort(dir("../data",pattern="pgx$"))
        shiny::updateSelectInput(session, "dataset2", choices=c("<this>",pgx.files))        
       
        ##pheno <- colnames(ngs$samples)
        ##updateSelectInput(session, "colorby", choices=pheno, selected=pheno[1])                
    })

    shiny::observeEvent( input$contrast1, {
        ## all.comparisons <- names(ngs$gx.meta$meta)
        ct <- input$contrast1
        shiny::req(ct)
        shiny::updateSelectInput(session, "colorby", choices=ct, selected=ct[1])                
    })
    
    ## shiny::observeEvent({
    ##     input$contrast1
    ##     input$ntop
    ## }, {
    ##     ngs <- inputData()
    ##     shiny::req(ngs)
    ##     F <- pgx.getMetaMatrix(ngs)$fc
    ##     ct1 <- input$contrast1
    ##     shiny::req(ct1)
    ##     if(!all(ct1 %in% colnames(F))) return(NULL)
    ##     ntop <- as.integer(input$ntop)
    ##     shiny::req(ntop)
    ##     message("DBG [CompareBoard::observer] dimF = ",dim(F))        
    ##     message("DBG [CompareBoard::observer] ntop = ",ntop)
    ##     message("DBG [CompareBoard::observer] ct1 = ",ct1)
    ##     F <- F[,ct1,drop=FALSE]
    ##     message("DBG [CompareBoard::observer] dimF = ",dim(F))
    ##     higenes <- head(names(sort(-rowMeans(F**2))),ntop)
    ##     higenes <- paste(higenes, collapse=' ')
    ##     shiny::updateTextAreaInput(session, "genelist", value=higenes)
    ## })

    shiny::observe({
        df <- getOmicsScoreTable()
        if(is.null(df)) return(NULL)
        ## message("DBG*** [CompareBoard::observer] updating genelist...")        
        ntop <- as.integer(input$ntop)
        higenes <- rownames(df)[order(df$score,decreasing=TRUE)]        
        higenes <- head(higenes, ntop)
        higenes <- paste(higenes, collapse=' ')
        shiny::updateTextAreaInput(session, "genelist", value=higenes)
    })
    
    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================

    ##  example <- shiny::reactive({
    ##      123
    ##  })
    
    dataset2 <- shiny::reactive({
        shiny::req(input$dataset2)
        if(input$dataset2 == "<this>") {
            ngs <- inputData()
        } else {
            load(file.path("../data",input$dataset2))
        }
        comparisons2 <- names(ngs$gx.meta$meta)
        sel2 <- tail(head(comparisons2,2),1)
        shiny::updateSelectInput(session, "contrast2", choices=comparisons2, selected=sel2)
        ngs
    })

    getOmicsScoreTable <- shiny::reactive({
        ## scatter2.RENDER <- shiny::reactive({                    
        ngs1 <- inputData()
        ngs2 <- dataset2()
        shiny::req(ngs1)
        shiny::req(ngs2)
        
        ct1 <- head(names(ngs1$gx.meta$meta),2)
        ct2 <- head(names(ngs2$gx.meta$meta),2)        
        ct1 <- input.contrast1()
        ct2 <- input.contrast2()
        shiny::req(ct1)
        shiny::req(ct2)
        if(!all(ct1 %in% names(ngs1$gx.meta$meta))) return(NULL)
        if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)

        F1 <- pgx.getMetaMatrix(ngs1)$fc[,ct1,drop=FALSE]
        F2 <- pgx.getMetaMatrix(ngs2)$fc[,ct2,drop=FALSE]        

        ##gg <- igraph::union(rownames(F1),rownames(F2))
        gg <- intersect(rownames(ngs1$X),rownames(ngs2$X))
        F1 <- F1[match(gg,rownames(F1)),,drop=FALSE]
        F2 <- F2[match(gg,rownames(F2)),,drop=FALSE]
        rownames(F1) <- gg
        rownames(F2) <- gg        
        colnames(F1) <- paste0("1:",colnames(F1))
        colnames(F2) <- paste0("2:",colnames(F2))
        rho = 1
        
        kk <- intersect(colnames(ngs1$X),colnames(ngs2$X))
        if(length(kk) >= 10) {
            X1 <- scale(t(ngs1$X[gg,kk]))
            X2 <- scale(t(ngs2$X[gg,kk]))
            rho <- colSums(X1 * X2) / (nrow(X1)-1)
        }

        fc1 <- sqrt(rowMeans(F1**2))
        fc2 <- sqrt(rowMeans(F2**2))        
        score <- rho * fc1 * fc2

        title <- ngs1$genes[gg,"gene_title"]
        title <- substring(title,1,60)
            
        df <- data.frame( title, score, rho, F1, F2, check.names=FALSE )
        df <- df[order(-df$score),]        
        df
    })

    hilightgenes <- shiny::reactive({
        genes <- as.character(input$genelist)
        genes <- strsplit(genes, split='[\t, \n]')[[1]]
        gsub("[ ]","",genes)
    }) ## %>% shiny::debounce(1000)

    input.contrast1 <- shiny::reactive({
        input$contrast1
    }) %>% shiny::debounce(2500)

    input.contrast2 <- shiny::reactive({
        input$contrast2
    }) %>% shiny::debounce(2500)
    
    ##============================================================================
    ## ScatterPlot 1
    ##============================================================================

    createPlot <- function(ngs, ngs1, ngs2, ct, type, cex.lab, higenes, ntop) {
        p <- NULL
        if(type %in% c('UMAP1','UMAP2')) {
            if(type =='UMAP1') {
                pos = ngs1$cluster.genes$pos[['umap2d']]
            } else if(type=='UMAP2') {
                pos = ngs2$cluster.genes$pos[['umap2d']]
            }
            p <- pgx.plotGeneUMAP(
                ngs, contrast=ct, pos=pos,
                cex = 0.9, cex.lab = cex.lab,
                hilight = higenes, ntop=ntop,
                par.sq = TRUE, plotlib="base")
        } else if(type == 'heatmap') {
            gg <- intersect(higenes, rownames(ngs$X))
            if(length(gg)>1) {
                X1 <- ngs$X[gg,,drop=FALSE]
                Y1 <- ngs$samples
                gx.splitmap(X1, nmax=40, col.annot=Y1,
                            softmax=TRUE, show_legend=FALSE)
            }
        } else {
            p <- pgx.plotContrast(
                ngs, contrast=ct,
                hilight = higenes, ntop=ntop,
                cex.lab = cex.lab, ## dlim=0.06,
                par.sq = TRUE, 
                type=type, plotlib="base")
        }
        p
    }
    
    scatter1.RENDER %<a-% shiny::reactive({
        ## scatter1.RENDER <- shiny::reactive({                    
        ngs1 <- inputData()
        ngs2 <- dataset2()
        all.ct <- names(ngs1$gx.meta$meta)
        ct1 <- head(all.ct,2)
        ct1 <- input.contrast1()
        shiny::req(ct1)
        if(!all(ct1 %in% all.ct)) return(NULL)
        type <- input$plottype
        higenes <- hilightgenes()
        cex.lab = 1.3
        cex.lab = 1.0
        ntop = 9999
        ##ntop = as.integer(input$ntop)
        if(length(higenes) <= 3) cex.lab = 1.3
        createPlot(ngs1, ngs1, ngs2, ct1, type, cex.lab, higenes, ntop)
        
    })
    
    scatter1.opts <- shiny::tagList()
    
    scatter1_caption = "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."
    
    shiny::callModule(
        plotModule,
        "scatter1", label = "a",
        func = scatter1.RENDER,
        func2 = scatter1.RENDER,
        ##plotlib = 'ggplot',
        ## options = scatter1.opts,
        title = "DATASET 1",
        ##info.text = scatter1_info
        ##caption = scatter1_caption,
        pdf.height=6, pdf.width=6, 
        height = c(700,750), width=c("auto",900),
        res = c(90,110),
        add.watermark = WATERMARK
    )

    ##============================================================================
    ## ScatterPlot 2
    ##============================================================================
    
    scatter2.RENDER %<a-% shiny::reactive({
    ## scatter2.RENDER <- shiny::reactive({                    
        ngs1 <- inputData()
        ngs2 <- dataset2()
        ##alertDataLoaded(session,ngs)
        ct2 <- input.contrast2()
        shiny::req(ct2)
        if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)
        type <- input$plottype
        higenes <- hilightgenes()
        cex.lab = 1.3
        cex.lab = 1.0
        ntop = 9999
        ##ntop = as.integer(input$ntop)
        if(length(higenes) <= 3) cex.lab = 1.3
        p = NULL
        p = createPlot(ngs2, ngs1, ngs2, ct2, type, cex.lab, higenes, ntop)
        p
    })
    
    scatter2.opts <- shiny::tagList()
    
    scatter2_caption = "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."
    
    shiny::callModule(
        plotModule,
        "scatter2", label = "b",
        func = scatter2.RENDER,
        func2 = scatter2.RENDER,
        ##plotlib = 'ggplot',
        ## options = scatter2.opts,
        title = "DATASET 2",
        ##info.text = scatter2_info
        ##caption = scatter2_caption,
        pdf.height=6, pdf.width=6, 
        height = c(700,750), width=c("auto",900),
        res = c(90,110),
        add.watermark = WATERMARK
    )
    
    ##============================================================================
    ## FC correlation plot
    ##============================================================================
    
    fcfcplot.RENDER %<a-% shiny::reactive({
        ## scatter2.RENDER <- shiny::reactive({                    
        ngs1 <- inputData()
        ngs2 <- dataset2()
        ##alertDataLoaded(session,ngs)
        ct1 <- head(names(ngs1$gx.meta$meta),2)
        ct2 <- head(names(ngs2$gx.meta$meta),2)
        ct1 <- input.contrast1()
        ct2 <- input.contrast2()
        shiny::req(ct1)
        shiny::req(ct2)
        if(!all(ct1 %in% names(ngs1$gx.meta$meta))) return(NULL)
        if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)        

        F1 <- pgx.getMetaMatrix(ngs1)$fc[,ct1,drop=FALSE]
        F2 <- pgx.getMetaMatrix(ngs2)$fc[,ct2,drop=FALSE]        
        gg <- intersect(rownames(F1),rownames(F2))
        F1 <- F1[gg,,drop=FALSE]
        F2 <- F2[gg,,drop=FALSE]
        colnames(F1) <- paste0("1:",colnames(F1))
        colnames(F2) <- paste0("2:",colnames(F2))

        ##ff <- rowMeans(F1**2) * rowMeans(F2**2)
        ##higenes <- head(names(sort(ff, decreasing=TRUE)),12)
        higenes <- hilightgenes()
        
        p <- NULL
        ##p <- plot.ggsplom(F1, F2, title_cex=3, no.axes=FALSE)
        plot.SPLOM(F1, F2=F2, cex=0.3, cex.axis=0.95, hilight=higenes)            
        p
    })
    
    fcfcplot.opts <- shiny::tagList()
    
    fcfcplot_caption = "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."
    
    shiny::callModule(
        plotModule,
        "fcfcplot", label = "a",
        func = fcfcplot.RENDER,
        func2 = fcfcplot.RENDER,
        ##plotlib = 'ggplot',
        ## options = scatter2.opts,
        title = "FC CORRELATION",
        ##info.text = scatter2_info
        ##caption = scatter2_caption,
        pdf.height=6, pdf.width=6, 
        height = c(700,fullH), width=c("auto",900),
        res = c(85,100),
        add.watermark = WATERMARK
    )

    ##-----------------------------------------------------------------
    ## Cumulative foldchange
    ##-----------------------------------------------------------------

    cumfcplot.RENDER %<a-% shiny::reactive({
        ## scatter2.RENDER <- shiny::reactive({                    
        ngs1 <- inputData()
        ngs2 <- dataset2()

        ct1 <- head(names(ngs1$gx.meta$meta),2)
        ct2 <- head(names(ngs2$gx.meta$meta),2)        
        ct1 <- input.contrast1()
        ct2 <- input.contrast2()
        shiny::req(ct1)
        shiny::req(ct2)
        if(!all(ct1 %in% names(ngs1$gx.meta$meta))) return(NULL)
        if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)

        F1 <- pgx.getMetaMatrix(ngs1)$fc[,ct1,drop=FALSE]
        F2 <- pgx.getMetaMatrix(ngs2)$fc[,ct2,drop=FALSE]        

        gg <- igraph::union(rownames(F1),rownames(F2))
        gg <- intersect(rownames(F1),rownames(F2))
        F1 <- F1[match(gg,rownames(F1)),,drop=FALSE]
        F2 <- F2[match(gg,rownames(F2)),,drop=FALSE]
        ##colnames(F1) <- paste0("1:",colnames(F1))
        ##colnames(F2) <- paste0("2:",colnames(F2))
        rownames(F1) <- rownames(F2) <- gg
        
        F <- cbind(F1, F2)
        F[is.na(F)] <- 0
        
        if(1) {
            sel <- head(order(-rowMeans(F**2)),50)
            sel <- rownames(F)[sel]
        } else {
            ## get top genes from score table
            df <- getOmicsScoreTable()
            if(is.null(df)) return(NULL)
            sel <- score_table$rows_all()      ## from module  
            shiny::req(sel)
            sel <- head(rownames(df)[sel],50)        
            if(length(sel)==0) return(NULL)
        }

        sel <- sel[order(rowMeans(F[sel,]))]
        F  <- F[sel,,drop=FALSE]
        F1 <- F1[sel,,drop=FALSE]
        F2 <- F2[sel,,drop=FALSE]

        if(0) {
            par(mfrow=c(1,1), mar=c(4.5,10,0,2))
            barplot( t(F), beside=FALSE, las=1, horiz=TRUE,
                    cex.names = 0.9,
                    xlab = "cumulative foldchange", ylab = "" )
            legend("bottomright", colnames(F), fill=grey.colors(ncol(F)),
                   cex=0.85, y.intersp=0.9, inset=c(0.01,0.02) )
        }
        
        par(mfrow=c(1,1), mar=c(4.5,0,1,2), mgp=c(2.2,0.8,0))
        plotly::layout(matrix(c(1, 2, 3), nrow=1, byrow=T),widths=c(0.5,1,1))

        frame()
        mtext( rownames(F), cex=0.80, side=2, at=(1:nrow(F)-0.5)/nrow(F),
              las=1, line=-12)        
        ##barplot( t(F1), beside=FALSE, las=1, horiz=TRUE, cex.names = 0.01,
        ##        xlab = "cumulative foldchange", ylab = "" )
        col1 <- grey.colors(ncol(F1))
        if(ncol(F1)==1) col1 <- "grey50"
        pgx.stackedBarplot( F1, hz=TRUE, las=1, col=col1,
                           cex.names = 0.01, cex.lab=1.4, space=0.25,
                           xlab = "cumulative foldchange", ylab = "" )
        legend("bottomright", colnames(F1), fill=grey.colors(ncol(F1)),
               cex=0.9, y.intersp=0.9, inset=c(-0.03,0.02), xpd=TRUE )
        title("DATASET1", line=-0.35, cex.main=1.2)
        
        ##barplot( t(F2), beside=FALSE, las=1, horiz=TRUE, cex.names = 0.01,
        ##        xlab = "cumulative foldchange", ylab = "" )
        col2 <- grey.colors(ncol(F2))
        if(ncol(F2)==1) col2 <- "grey50"
        pgx.stackedBarplot( F2, hz=TRUE, las=1, col=col2,
                           cex.names = 0.01, cex.lab=1.4, space=0.25,
                           xlab = "cumulative foldchange", ylab = "" )
        legend("bottomright", colnames(F2), fill=grey.colors(ncol(F2)),
               cex=0.9, y.intersp=0.9, inset=c(-0.03,0.02), xpd=TRUE )
        title("DATASET2", line=-0.35, cex.main=1.2)
        
        
    })
    
    cumfcplot.opts <- shiny::tagList()    
    cumfcplot_caption = "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."
    
    shiny::callModule(
        plotModule,
        "cumfcplot", label = "b",
        func = cumfcplot.RENDER,
        func2 = cumfcplot.RENDER,
        ## plotlib = 'ggplot',
        ## options = scatter2.opts,
        title = "CUMULATIVE FOLDCHANGE",
        ## info.text = scatter2_info
        ## caption = scatter2_caption,
        pdf.height=8, pdf.width=8, 
        height = c(700,750), width=c("auto",900),
        res = c(80,98),
        add.watermark = WATERMARK
    )

    ##============================================================================
    ## Gene correlation plot
    ##============================================================================
    
    genecorr.RENDER %<a-% shiny::reactive({
        ## scatter2.RENDER <- shiny::reactive({                    
        ngs1 <- inputData()
        ngs2 <- dataset2()
        ##alertDataLoaded(session,ngs)
        ct1 <- head(names(ngs1$gx.meta$meta),2)
        ct2 <- head(names(ngs2$gx.meta$meta),2)
        ct1 <- input.contrast1()
        ct2 <- input.contrast2()
        shiny::req(ct1)
        shiny::req(ct2)
        if(!all(ct1 %in% names(ngs1$gx.meta$meta))) return(NULL)
        if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)        

        dbg("[genecorr.RENDER] 1:")
        
        gg <- intersect(rownames(ngs1$X),rownames(ngs2$X))
        kk <- intersect(colnames(ngs1$X),colnames(ngs2$X))

        dbg("[genecorr.RENDER] 2:")
        
        if(length(kk)==0) {
            par(mfrow=c(1,1))
            frame()
            text(0.5,0.6, "Warning: no common samples", col='black')
            text(0.5,0.5, "To compute gene correlation both datasets\nneed to have common samples",
                 col='black')
            return()
        }
        if(length(kk) < 10) {
            par(mfrow=c(1,1))
            frame()
            text(0.5,0.6, "Error: too few samples", col='red3')
            text(0.5,0.5, "For gene correlation we need at least 10 common samples", col='red3')
            return()
        }

        dbg("[genecorr.RENDER] 3:")
        
        ## conform matrices
        X1 <- ngs1$X[gg,kk]
        X2 <- ngs2$X[gg,kk]        
        Y1 <- ngs1$samples[kk,]
        Y2 <- ngs2$samples[kk,]    
        
        dbg("[genecorr.RENDER] 4:")

        dset1 <- paste0("[dataset1]  expression")
        dset2 <- paste0("[dataset2]  expression")
        dset1 <- paste0("1: expression")
        dset2 <- paste0("2: expression")
        ##dset2 <- paste0("[",input$dataset2,"] expression")

        if(0) {
            F <- pgx.getMetaMatrix(ngs1)$fc
            higenes <- names(head(sort(-rowMeans(F**2)),16))
            higenes <- hilightgenes()
            higenes <- intersect(higenes,rownames(X1))
            higenes <- head(higenes, 16)
        }

        dbg("[genecorr.RENDER] 5:")
        
        df <- getOmicsScoreTable()
        if(is.null(df)) return(NULL)
        
        ##sel <- input$score_table_rows_all
        sel <- score_table$rows_all()      ## from module  
        shiny::req(sel)
        if(is.null(sel)) return(NULL)        

        dbg("[genecorr.RENDER] 6: dim.df = ",dim(df))
        
        higenes <- head(rownames(df)[sel],16)        
        if(length(higenes)==0) return(NULL)

        dbg("[genecorr.RENDER] 7: higenes = ",higenes)
        
        ## Set color for points
        ##C1 <- ngs1$model.parameters$exp.matrix
        ##klr1 <- 2+as.integer(C1[,ct1[1]])
        klrpal <- rep(1:7,99)
        klrpal <- rep(RColorBrewer::brewer.pal(12,"Paired"),99)
        
        colorby="ER_STATUS"
        colorby = ct1[1]
        ## colorby <- input$colorby
        if(0) {
            grp <- factor(Y1[,colorby])
            klr1 <- klrpal[as.integer(grp)]
        } else {
            grp <- pgx.getContrastGroups(ngs1, colorby, as.factor=TRUE)
            ## grp <- factor(grp)
            grp <- grp[colnames(X1)]
            klr1 <- klrpal[as.integer(grp)]            
        }
        
        nc <- ceiling(sqrt(length(higenes)))
        nr <- (length(higenes)-1) %/% nc
        par(mfrow=c(nc,nc), mar=c(2.6,2.3,1.5,0.5)*1, oma=c(3,3,0,0), mgp=c(2.4,0.7,0) )
        i=1
        for(i in 1:length(higenes)) {
            j <- match(higenes[i],rownames(X1))
            KEGGgraph::plot( X1[j,], X2[j,],
                 ##xlab=dset1, ylab=dset2,
                 xlab="", ylab="",
                 pch=20, col=klr1, cex=1.2)
            title(higenes[i], line=0.4, cex.main=1.1)
            if( i%%nc == 1) {
                mtext(dset2, 2, line=2, cex=0.8)
            }
            if((i-1)%/%nc==nr) {
                mtext(dset1, 1, line=2, cex=0.8)
            }
            ##if(i==1) {
            if(i%%nc == 1) {
                tt <- c("   ",levels(grp))
                legend("topleft", legend=tt,
                       fill = c(NA,klrpal),
                       border = c(NA,"black","black"), bty='n',
                       cex=0.92, box.lwd=0, pt.lwd=0,
                       x.intersp=0.5, y.intersp=0.8)
                legend("topleft", colorby, x.intersp=-0.2,
                       cex=0.92, y.intersp=0.45, bty='n')
            }
        }

    })
    
    genecorr.opts <- shiny::tagList(
        shinyBS::tipify( shiny::selectInput(ns('colorby'),'Color by:',choices=NULL, multiple=FALSE),
               "Color samples by phenotype.",
               placement="right", options = list(container = "body")
               )
    )
    
    genecorr_caption = "<b>FC scatter plots.</b> Scatter plots of gene expression scatter values between two contrasts. Scatters that are similar show high correlation, i.e. are close to the diagonal."
    
    shiny::callModule(
        plotModule,
        "genecorr", label = "c",
        func = genecorr.RENDER,
        func2 = genecorr.RENDER,
        ##plotlib = 'ggplot',
        options = genecorr.opts,
        title = "GENE CORRELATION",
        ##info.text = scatter2_info
        ##caption = scatter2_caption,
        pdf.height=6, pdf.width=6, 
        height = c(740,750), width=c('auto',900),
        res = c(80,90),
        add.watermark = WATERMARK
    )
    
    ##============================================================================
    ## barplots
    ##============================================================================

    multibarplot.RENDER %<a-% shiny::reactive({

        ngs1 <- inputData()
        ngs2 <- dataset2()

        ct1 <- head(names(ngs1$gx.meta$meta),3)
        ct2 <- head(names(ngs2$gx.meta$meta),3)
        ct1 <- input.contrast1()
        ct2 <- input.contrast2()
        shiny::req(ct1)
        shiny::req(ct2)
        if(!all(ct1 %in% names(ngs1$gx.meta$meta))) return(NULL)
        if(!all(ct2 %in% names(ngs2$gx.meta$meta))) return(NULL)

        gene="ESR1"
        gene="ERBB2"
        genes = c("ERBB2","ESR1","FOXA1","PGR","ERBB2","ESR1","FOXA1","PGR")
        genes = rownames(ngs1$X)
        genes <- hilightgenes()

        df <- getOmicsScoreTable()
        if(is.null(df)) return(NULL)        
        sel <- score_table$rows_all()      ## from module  
        shiny::req(sel)
        genes <- rownames(df)[sel]
                
        xgenes <- intersect(rownames(ngs1$X), rownames(ngs2$X))
        genes <- head(intersect(genes, xgenes),8)
        
        par(mfrow=c(2,4), mar=c(6,4,1,0), mgp=c(2.0,0.7,0), oma=c(0,0,0,0))
        for(gene in genes) {
            x1 <- ngs1$X[gene,]
            x2 <- ngs2$X[gene,]
            e1 <- contrastAsLabels(ngs1$model.parameters$exp.matrix[,ct1,drop=FALSE])
            e2 <- contrastAsLabels(ngs2$model.parameters$exp.matrix[,ct2,drop=FALSE])
            m1 <- lapply(e1, function(y) tapply(x1, y, mean))
            m2 <- lapply(e2, function(y) tapply(x2, y, mean))

            ##gx.barplot(m1, srt=45, main=gene);gx.barplot(m2, srt=45, main=gene)
            grp1 <- paste0("1:",sub(":.*","",names(m1)))
            grp2 <- paste0("2:",sub(":.*","",names(m2)))
            grp.names <- c(grp1, grp2)

            b1 <- as.vector(sapply(m1,names))
            b2 <- as.vector(sapply(m2,names))
            bar.names <- c(b1,b2)
            bar.names <- toupper(substring(bar.names,1,1))
            
            srt=10
            srt=0
            srt <- ifelse(max(nchar(grp.names)) <= 5, 0, 25)
            
            mm <- cbind(do.call(cbind, m1), do.call(cbind, m2))
            mm.group <- c(rep(1,length(m1)), rep(2,length(m2)) )
            ##gx.barplot(mm, legend=FALSE, srt=15)
            gx.barplot(mm,
                       srt=srt, main=gene, cex.main=1.0,
                       group=mm.group, cex.names=0.85,
                       group.names = grp.names,
                       bar.names = bar.names, voff=3.5,
                       legend=FALSE, cex.legend=0.9,                       
                       ylab="expression")
            
        }

    })
    
    multibarplot.opts <- shiny::tagList()    
    multibarplot_info = "<b>Multi barplots.</b> Barplots of expression values for multiple comparisons in the two datasets (blue and green). "
    
    shiny::callModule(
        plotModule,
        "multibarplot", label = "a",
        func = multibarplot.RENDER,
        func2 = multibarplot.RENDER,
        ##plotlib = 'ggplot',
        ## options = scatter2.opts,
        title = "EXPRESSION",
        info.text = multibarplot_info,
        ##caption = scatter2_caption,
        pdf.height=6, pdf.width=8, 
        height = c(440,700),        
        width=c("auto",1280),
        res = c(95,130),
        add.watermark = WATERMARK
    )

    
    ##============================================================================
    ## Score table
    ##============================================================================

    score_table.RENDER <- shiny::reactive({

        df <- getOmicsScoreTable()
        if(is.null(df)) return(NULL)
        shiny::req(df)
        numeric.cols <- 2:ncol(df)
        
        DT::datatable(
                df, rownames=TRUE, ## escape = c(-1,-2),
                extensions = c('Buttons','Scroller'),
                selection=list(mode='single', target='row', selected=NULL),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', ##buttons = c('copy','csv','pdf'),
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller = TRUE,
                    deferRender = TRUE
                )  ## end of options.list 
            ) %>%
            DT::formatSignif(numeric.cols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 

    })

    score_table_info = "In this table, users can check mean expression values of features across the conditions for the selected genes."
    
    score_table <- shiny::callModule(
        tableModule, id = "score_table",
        func = score_table.RENDER, ## ns=ns,
        info.text = score_table_info,
        title = "CORRELATION SCORE",
        label = "b",
        height = c(235,750),
        width = c("auto",1600)
        ## caption = parcoord_caption
    )
    
    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================

    ## ------------------------------------------------------
    ## --------------------- tab1 ---------------------------
    ## ------------------------------------------------------
    
    compareScatter_title = "<h4>Compare Expression</h4>"
    
    output$compareScatter_UI <- shiny::renderUI({
        ##req(input$dimensions)
        shiny::fillCol(
            height = fullH,
            ## height = input$dimensions[2], ## dynamics with JS  
            flex = c(NA,0.015,1),
            shiny::div(shiny::HTML(compareScatter_title), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.1,1),
                plotWidget(ns("scatter1")),
                shiny::br(),
                ##plotWidget(ns("connectivityMap"))
                plotWidget(ns("scatter2"))
            )
        )
    })
    shiny::outputOptions(output, "compareScatter_UI", suspendWhenHidden=FALSE) ## important!!!

    ## ------------------------------------------------------
    ## --------------------- tab2 ---------------------------
    ## ------------------------------------------------------
    
    FCcorrelation_title = "<h4>Compare Foldchange</h4>"
    
    output$FCcorrelation_UI <- shiny::renderUI({
        ##req(input$dimensions)
        shiny::fillCol(
            height = fullH,
            ## height = input$dimensions[2], ## dynamics with JS  
            flex = c(NA,0.015,1),
            shiny::div(shiny::HTML(FCcorrelation_title), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.05,0.9),
                plotWidget(ns("fcfcplot")),
                shiny::br(),
                plotWidget(ns("cumfcplot"))
            )
        )
    })
    shiny::outputOptions(output, "FCcorrelation_UI", suspendWhenHidden=FALSE) ## important!!!

    ## ------------------------------------------------------
    ## --------------------- tab3 ---------------------------
    ## ------------------------------------------------------
    
    GeneCorrelation_title = "<h4>Gene Correlation</h4>"
    
    output$GeneCorrelation_UI <- shiny::renderUI({
        ##req(input$dimensions)
        shiny::fillCol(
            height = fullH,
            ## height = input$dimensions[2], ## dynamics with JS  
            flex = c(NA,0.015,1),
            shiny::div(shiny::HTML(GeneCorrelation_title), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1.0,0.05,1),
                shiny::fillCol(
                    flex = c(2.2,0.01,1),
                    plotWidget(ns("multibarplot")),                    
                    shiny::br(),
                    tableWidget(ns("score_table"))
                ),
                shiny::br(),
                plotWidget(ns("genecorr"))
            )
        )
    })
    shiny::outputOptions(output, "GeneCorrelation_UI", suspendWhenHidden=FALSE) ## important!!!         
    
} ## end-of-Board 


if(0) {
    source(file.path(RDIR,"pgx-include.R"))    ## lots of libraries and source()
    
    load("../data/tcga-brca2-gx.pgx")
    ngs1 = pgx.initialize(ngs)
    load("../data/tcga-brca2-px.pgx")
    ngs2 = pgx.initialize(ngs)
    
    gg <- intersect(rownames(ngs1$X),rownames(ngs2$X))
    kk <- intersect(colnames(ngs1$X),colnames(ngs2$X))
    X1 <- ngs1$X[gg,kk]
    X2 <- ngs2$X[gg,kk]

    Y <- ngs1$samples[kk,]
    i <- grep("ERBB2",rownames(X1))
    klr1 <- 1+as.integer(Y$HER2_STATUS)
    KEGGgraph::plot(X1[i,], X2[i,], pch=20, col=klr1, cex=1)

    pgx.plotExpression(ngs1, "ERBB2", comp="HER2_STATUS:1_vs_0", grouped=TRUE)
    pgx.plotExpression(ngs1, "ESR1", comp="ER_STATUS:1_vs_0", grouped=TRUE)
    pgx.plotExpression(ngs, "ESR1", comp="ER_STATUS:1_vs_0", grouped=TRUE)    

    all.ct = names(ngs$gx.meta$meta)
    pgx.Volcano(ngs, contrast=all.ct[2], plotlib='ggplot')
    
    gene="ESR1"
    gene="ERBB2"
    genes = c("ERBB2","ESR1","FOXA1","PGR","CA12","FABP7")
    par(mfrow=c(2,4), mar=c(10,4,4,2))
    for(gene in genes) {
        x1 <- ngs1$X[gene,]
        x2 <- ngs2$X[gene,]
        e1 <- ngs1$model.parameters$exp.matrix[,c(2,5,7,8)]
        e2 <- ngs2$model.parameters$exp.matrix[,c(2,5,7,8)]        
        m1 <- apply(e1, 2, function(y) tapply(x1, sign(y), mean)[c("-1","1")])
        m2 <- apply(e2, 2, function(y) tapply(x2, sign(y), mean)[c("-1","1")])    
    
        ##gx.barplot(m1, srt=45, main=gene);gx.barplot(m2, srt=45, main=gene)
        colnames(m1) <- paste0("[1] ",colnames(m1))
        colnames(m2) <- paste0("[2] ",colnames(m2))
        mm <- cbind(m1, m2)
        mm.group <- c( rep(1,ncol(m1)), rep(2,ncol(m2)) )
        gx.barplot(mm, srt=45, main=gene, group=mm.group,
                   cex.names=0.8, voffset=0.06, offset=NULL)
    }
    

    plot.ggbarplot(m1[,], srt=45, beside=TRUE, main=gene,
              legend.pos = c(0.016,1), 
              bar_width=0.6, legend.cex=1.5)

    plot.ggbarplot(m2[,], srt=45, beside=TRUE, main=gene,
              bar_width=0.6, legend.cex=1.5)

    legend.cex=2
    p + ggplot2::theme(
            legend.title = ggplot2::element_blank(), 
            ##legend.justification = legend.pos,
            legend.text = ggplot2::element_text(size=9*legend.cex),
            ##legend.position = legend.pos,
            legend.key.size = ggplot2::unit(9*legend.cex, "pt"),
            legend.key.height = ggplot2::unit(7*legend.cex, "pt"))

}
