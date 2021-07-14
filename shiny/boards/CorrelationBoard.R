##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing CorrelationBoard")

CorrelationInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

CorrelationUI <- function(id) {
    ns <- NS(id)  ## namespace
    ui <- fillCol(
        height = 750,
        tabsetPanel(
            id = ns("tabs"),
            tabPanel("Correlation",uiOutput(ns("corAnalysis_UI"))),
            ## tabPanel("Functional",uiOutput(ns("corFunctional_UI"))),
            tabPanel("Graph",uiOutput(ns("corGraph_UI"))),            
            tabPanel("Differential",uiOutput(ns("corDiff_UI")))
        )
    )
    ui
}

CorrelationBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    inputData <- env[["load"]][["inputData"]]

    fullH = 800  ## full height of page
    rowH  = 340  ## full height of page
    
    description = "<b>Correlation Analysis.</b> Compute the correlation between genes and find coregulated modules."
    output$description <- renderUI(HTML(description))


    cor_infotext ="The <strong>Correlation Analysis Board</strong> provides statistical correlation analysis on gene level with visualisations. During the visual analysis, users can filter out some samples or collapse the samples by predetermined groups. The dark shaded area in the barplot estimates the partial correlation."
    
    require(RColorBrewer)
    COL <- brewer.pal(12,"Paired")[seq(1,12,2)]
    COL <- brewer.pal(9,"Set1")[c(2,1,3:9)]
    COL2 <- rev(grey.colors(2))
    COL2 <- brewer.pal(2,"Paired")[1:2]
    COL2 <- COL[1:2]
    
    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================

    corAnalysis_caption = "<h3>Gene Correlation Analysis</h3><b>(a)</b> <b>Top-ranked correlation.</b> Top correlated features with respect to selected gene. <b>(b)</b> <b>Correlation table</b> of correlation and partial correlation with respect to selected gene. <b>(c)</b> <b>Scatter plots</b> of gene expression of top correlated genes."
    
    output$corAnalysis_UI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,
            div(HTML(corAnalysis_caption), class="caption"),
            br(),
            fillRow(
                flex = c(1,0.01,1.2),
                fillCol(
                    flex = c(1,1),
                    height = fullH-80,
                    plotWidget(ns('cor_barplot')),
                    ##plotWidget(ns('cor_graph'))
                    tableWidget(ns('cor_table'))
                ),	
                br(), ## spacer
                fillCol(
                    flex = c(NA,1),
                    height = fullH-80,
                    plotWidget(ns('cor_scatter'))
                )
            )
        )
    })
    outputOptions(output, "corAnalysis_UI", suspendWhenHidden=FALSE) ## important!!!

    
    corfunctional_caption ="<b>(a)</b> <b>Correlation GSEA.</b> Top enriched gene sets using the correlation as rank metric. The black bars denote the genes in the gene set and their position in the sorted rank metric. <b>(b)</b> <b>Enrichment table.</b> Statistical results from GSEA analysis. <b>(c)</b> <b>Gene frequency.</b> Frequency of leading edge genes in top correlated genesets. <b>(d)</b> <b>Leading edge table.</b> Leading edge genes and rank statistics (rho) of the selected geneset."

    output$corFunctional_UI <- renderUI({
        fillCol(
            flex = c(NA,0.025,1),
            height = fullH,
            div(HTML(corfunctional_caption), class="caption"),
            br(),
            fillRow(
                flex = c(1.8,0.11,1),
                height = fullH - 60,
                fillCol(
                    flex = c(1,0.07,0.6),
                    plotWidget(ns("corGSEA_plots")),
                    br(),
                    tableWidget(ns("corGSEA_table"))
                ),
                br(),
                fillCol(
                    flex = c(1,0.04,1.3),
                    ##flex = c(1),
                    plotWidget(ns("corGSEA_cumFC")),
                    br(),
                    tableWidget(ns("corGSEA_LeadingEdgeTable"))
                    ##div(HTML("caption"), class="caption")
                )
            )
        )
    })
    ##outputOptions(output, "hm_annotateUI", suspendWhenHidden=FALSE) ## important!!!

    corDiff_caption = "<h3>Differential Gene Correlation Analysis (DGCA)</h3>Compute and analyze differential correlations between gene pairs across multiple conditions."
    
    output$corDiff_UI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,
            div(HTML(corDiff_caption), class="caption"),
            br(),
            fillRow(
                flex = c(1,0.01,1.2),
                fillCol(
                    flex = c(1.3,1),
                    height = fullH-80,
                    plotWidget(ns('dgca_barplot')),
                    tableWidget(ns('dgca_table'))
                ),	
                br(), ## spacer
                fillCol(
                    flex = c(NA,1),
                    height = fullH-80,
                    plotWidget(ns('dgca_scatter'))
                )
            )
        )
    })
    outputOptions(output, "corDiff_UI", suspendWhenHidden=FALSE) ## important!!!    

    corGraph_caption = "<h3>Gene Correlation Network</h3>Visualization of gene correlation as network or UMAP. <b>(a)</b> <b>Partial correlation network</b> around the selected gene. <b>(b)</b> <b>Correlation UMAP</b>."
    
    output$corGraph_UI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,
            div(HTML(corGraph_caption), class="caption"),
            br(),
            fillRow(
                flex = c(1,0.05,1),
                plotWidget(ns('cor_graph')),
                br(),
                plotWidget(ns('cor_umap'))
            )
        )
    })
    outputOptions(output, "corGraph_UI", suspendWhenHidden=FALSE) ## important!!!    
    
    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================
    require(htmltools)
    
    output$inputsUI <- renderUI({
        ui <- tagList(
            actionLink(ns("cor_info"), "Info", icon=icon("info-circle")),
            hr(), br(),             

            ## data set parameters
            tipify( selectInput(ns("cor_gene"),"Gene:", choices=NULL),
                   "Choose a gene for the correlation analysis.", placement="top"),
            br(),
            tipify( selectInput(ns("cor_group"),"Color by:", choices=NULL, multiple=FALSE),
                   "Variable to split and color by groups.", placement="top"),
            br(),
            actionLink(ns("cor_options"), "Options", icon=icon("cog", lib = "glyphicon")),
            br(),br(),
            conditionalPanel(
                "input.cor_options % 2 == 1", ns=ns,
                tagList(
                    tipify( selectInput(ns("cor_features"),"Filter genes:", choices=NULL, multiple=FALSE),
                           "Filter gene features.", placement="top"),
                    conditionalPanel(
                        "input.cor_features == '<custom>'", ns=ns,
                        tipify( textAreaInput(ns("cor_customfeatures"),
                                              NULL, value = NULL,
                                              height = "100px", width = "100%", 
                                              rows=5, placeholder="Paste your custom gene list"),
                               "Paste a custom list of genes to be used as features.",
                               placement="top")
                    ),                    
                    ##tipify( selectInput(ns("cor_samplefilter"),"Filter samples",
                    ##                    choices=NULL, multiple=TRUE),
                    ##       "Filter (include) samples for the analysis", placement="top"),
                    tipify( checkboxInput(ns("dgca.allpairs"),"All pairs for DGCA"),
                           "Compute all DGCA pairs", placement="top")
                )
            )
        )

    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================

    observeEvent( input$cor_info, {
        showModal(modalDialog(
            title = HTML("<strong>Correlation Analysis Board</strong>"),
            HTML(cor_infotext),
            easyClose = TRUE ))
    })
    
    ## update filter choices upon change of data set 
    observe({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        
        ## levels for sample filter
        ##levels = getLevels(ngs$Y)
        ##updateSelectInput(session, "cor_samplefilter", choices=levels)

        genes <- sort(ngs$genes[rownames(ngs$X),]$gene_name)
        sel = genes[1]  ## most var gene
        sel = names(head(sort(-rowMeans(pgx.getMetaMatrix(ngs)$fc**2)),1))
        updateSelectizeInput(session,'cor_gene', choices=genes, selected=sel, server=TRUE)

        fam <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        fam <- sort(c("<custom>",fam))
        updateSelectInput(session, "cor_features", choices=fam)
        
        px <- colnames(ngs$Y)
        updateSelectInput(session, "cor_group", choices=px)
        
    })


    getFilteredExpression <- reactive({

        ngs <- inputData()
        req(ngs,input$cor_gene)
        X <- ngs$X    
        gene <- rownames(X)[1]
        gene <- input$cor_gene
        
        ## filter genes
        ft <- input$cor_features
        if(ft=="<custom>" && input$cor_customfeatures!="") {
            genes <- toupper(ngs$genes$gene_name)
            gg1 = strsplit(input$cor_customfeatures,split="[, ;\n\t]")[[1]]
            if(length(gg1)==1) gg1 <- paste0(gg1,"*")
            gg1 = gsub("[ \n\t]","",gg1)
            starred = grep("[*]",gg1)
            if(length(starred)>0) {
                gg2 = lapply(gg1[starred], function(a)
                    genes[grep(paste0("^",sub("[*]","",a)),genes,ignore.case=TRUE)])
                gg1 = unique(c(gg1,unlist(gg2)))
            }
            gg1 <- gg1[which(toupper(gg1) %in% toupper(ngs$genes$gene_name))]
            psel <- filterProbes(ngs$genes, c(gg1,gene))
            psel <- intersect(psel,rownames(X))
            X = X[psel,,drop=FALSE]
            
        } else if(ft!="<all>" && ft %in% names(GSETS)) {
            ft <- input$cor_features
            psel = filterProbes(ngs$genes, c(gene,GSETS[[ft]]) )
            ##psel = unique(c(gene, psel))
            psel <- intersect(psel,rownames(X))
            X = X[psel,,drop=FALSE]
        }

        if(0) {
            TISSUE            
        }       
        
        X <- X + 1e-3*matrix(rnorm(length(X)),nrow(X),ncol(X))
        X
    })
    
    getPartialCorrelationMatrix <- reactive({
        ngs <- inputData()
        req(ngs,input$cor_gene)

        gene = rownames(ngs$X)[1]
        gene <- input$cor_gene

        ## filter gene expression matrix
        X <- getFilteredExpression()        
        
        showNotification(paste("computing correlation...\n"))        
        NTOP = 50
        NTOP = as.integer(input$pcor_ntop)
        ##res <- pgx.computePartialCorrelationAroundGene(
        ##    X, gene, method=methods, nmax=NTOP, fast=FALSE)
        res <- pgx.computeGlassoAroundGene(X, gene, nmax=NTOP)    
        res$meta.pcor <- res$pcor
        
        j <- which(rownames(res$pcor)==gene)
        P <- res$pcor
        diag(P) <- 0
        rho1 <- min(head(sort(P,decreasing=TRUE),200))
        max1 <- round(max(P),digits=3)
        updateSliderInput(session, "cor_graph_threshold", value=rho1, max=max1)
        updateSliderInput(session, "dcga_graph_threshold", value=rho1, max=max1)        

        res
    })
    
    getPartialCorrelation <- reactive({
        res <- getPartialCorrelationMatrix()
        gene <- rownames(res$cor)[1]
        gene <- input$cor_gene
        rho <- res$cor[gene,]
        prho <- res$pcor[gene,]
        df <- data.frame(cor=rho, pcor=prho)
        df
    })

    getGeneCorr <- reactive({
        require(RColorBrewer)
        ngs <- inputData()
        req(ngs)	
        
        dbg("[getGeneCorr] reacted!")

        samples=colnames(ngs$X);gene="CD4"
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        gene <- input$cor_gene
        dbg("[getGeneCorr] 1a: gene = ",gene)        
        if(is.null(gene)) return(NULL)    

        dbg("[getGeneCorr] 1b:")
        
        ## corr always in log.scale and restricted to selected samples subset
        zx <- ngs$X
        zx <- getFilteredExpression()        
        dim(zx)
        zx.genes0 <- rownames(zx)
        ##rownames(zx) <- toupper(sub(".*:","",rownames(zx)))  ## NEED RETHINK!
        zx.genes <- as.character(ngs$genes[rownames(zx),]$gene_name)
        rownames(zx) <- toupper(zx.genes)        
        xref <- list("cor" = 2**zx,
                     "cor.HPA" = as.matrix(TISSUE),
                     "cor.ImmProt" = as.matrix(IMMPROT))
        gene0 <- toupper(gene)  ## uppercase mouse

        dbg("[getGeneCorr] 2: gene0 = ", gene0)
        
        R <- pgx.getGeneCorrelation(gene0, xref=xref)    
        if(is.null(R)) return(NULL)        
        R <- R[rownames(zx),]
        
        zx <- zx - rowMeans(zx, na.rm=TRUE)
        sdx <- sqrt(rowMeans(zx**2))
        R <- cbind(R, cov=R[,"cor"] * sdx * sdx[gene])
        
        dbg("[getGeneCorr] 3: dim.R = ",dim(R))
        
        rho.genes = rownames(zx)
        if("hgnc_symbol" %in% colnames(ngs$genes)) {
            rho.genes = as.character(ngs$genes[zx.genes0,]$hgnc_symbol)
        }
        R <- R[match(rho.genes,rownames(R)),,drop=FALSE]
        rownames(R) <- zx.genes0

        dbg("[getGeneCorr] 4: dim.R = ",dim(R))

        R <- R[order(R[,"cor"],decreasing=TRUE),]
        
        R
    })
    

    ##================================================================================
    ##======================= PLOTTING FUNCTIONS =====================================
    ##================================================================================


    ##-----------------------------------------------------------
    ## Correlation barplot
    ##-----------------------------------------------------------
    
    cor_barplot.PLOTFUN %<a-% reactive({

        df <- getPartialCorrelation()
        
        ##rho <- df$cor
        ##names(rho) <- rownames(df)
        ##rho <- rho[order(-rho)]

        R <- getGeneCorr()            
        dbg("[cor_barplot.PLOTFUN] 1: dim.R = ",dim(R))
        sel <- cor_table$rows_all()
        req(sel)
        NTOP=50
        sel <- head(sel,NTOP)
        rho <- R[sel,"cor"]
        if(length(sel)==1) names(rho) <- rownames(R)[sel]
        
        prho <- df$pcor
        names(prho) <- rownames(df)
        prho <- prho[match(names(rho),names(prho))]
        names(prho) <- names(rho)
        ##prho[is.na(prho)] <- 0            

        ylim0 <- c(-1,1)*max(abs(rho))*1.05
            
        par(mfrow=c(1,1), mar=c(10,4,1,0.5))
        barplot(rho, beside=FALSE, las=3,
                ylim = ylim0,
                ylab = "correlation",
                cex.names=0.85 )        
        barplot(prho, beside=FALSE, add=TRUE,                
                col="grey40", names.arg="")
        legend("topright", cex=0.85, y.intersp=0.85,
               inset = c(0.035,0),
               c("correlation","partial correlation"),
               fill=c("grey70","grey40"))   
        
    })

    cor_barplot.opts <- tagList(
        ##radioButtons(ns('cor_partialpc'),'partial correlation',
        ##             c("none","fast","all methods"), selected="none",
        ##             inline=TRUE)
        radioButtons(ns('pcor_ntop'),'nr of top genes to compute partial correlation.',
                     c(50,100,250), selected=100, inline=TRUE)
    )

    cor_barplot.info = "<b>Top correlated genes.</b> Highest correlated genes in respect to the selected gene. The height of the bars correspond to the Pearson correlation value. The dark grey bars correspond to the 'partial correlation' which essentially corrects the correlation value for indirect effects and tries to estimate the amount of direct interaction."
   
    callModule(
        plotModule,
        id = "cor_barplot", 
        func = cor_barplot.PLOTFUN,
        func2 = cor_barplot.PLOTFUN,        
        csvFunc = getPartialCorrelation,
        info.text = cor_barplot.info,
        options = cor_barplot.opts,
        title = "Top correlated genes", label = "a",
        caption2 = cor_barplot.info,
        pdf.width = 10, pdf.height = 5, 
        height = c(0.45*fullH,700),
        width = c('auto',1200),
        res = c(63,100),
        add.watermark = WATERMARK        
    )

    ##-----------------------------------------------------------------------------
    ## Correlation scatter plots
    ##-----------------------------------------------------------------------------
    
    cor_scatter.PLOTFUN %<a-% reactive({

        ngs <- inputData()
        req(input$cor_gene)
        X <- getFilteredExpression()

        this.gene=rownames(ngs$X)[1]
        this.gene <- input$cor_gene

        dbg("[cor_scatter.PLOTFUN] 0:")
        
        NTOP = 25
        if(0) {
            res <- getPartialCorrelationMatrix()
            j <- which(rownames(res$cor) == this.gene)
            rho <- res$cor[j,-j]
            rho <- head(rho[order(-abs(rho))],NTOP)
            rho <- rho[order(-rho)]
        } else {
            R <- getGeneCorr()            
            dbg("[cor_scatter.PLOTFUN] 1: dim.R = ",dim(R))
            sel <- cor_table$rows_all()
            req(sel)
            rho <- head(R[sel,"cor"],NTOP)
            if(length(sel)==1) names(rho) <- rownames(R)[sel]
        }

        dbg("[cor_scatter.PLOTFUN] 2: len.rho = ",length(rho))
        if(length(rho)==0) return(NULL)

        dbg("[cor_scatter.PLOTFUN] 2: head.names.rho = ",head(names(rho)))
        
        colorby=1
        colorby <- input$cor_group
        req(colorby)
        
        ph <- factor(ngs$samples[,colorby])
        klrpal <- rep(COL,99)
        klr <- klrpal[as.integer(ph)]

        dbg("[cor_scatter.PLOTFUN] 2: levels.ph = ",paste(levels(ph),collapse=' '))
        
        ndim <- ncol(ngs$X)
        cex = 1.2
        cex = ifelse( ndim > 40, 0.8, 1.2)
        cex = ifelse( ndim > 100, 0.5, cex)
        cex = ifelse( ndim > 200, 0.2, cex)                

        par(mfrow=c(5,5), mar=c(4,3.5,0.3,0),
            mgp=c(1.9,0.7,0), oma=c(0,0,0.5,0.5))
        par(mfrow=c(5,5), mar=c(3,3.5,0.5,0),
            mgp=c(1.7,0.6,0), oma=c(0,0,0.5,0.5))

        dbg("[cor_scatter.PLOTFUN] 3:")
        
        swapaxis = FALSE
        swapaxis = TRUE
        swapaxis <- input$corscatter.swapaxis
        ylab = gene2
        xlab = this.gene

        dbg("[cor_scatter.PLOTFUN] 4: this.gene = ", this.gene)
        
        i=1
        for(i in 1:min(25,length(rho))) {
            gene2 <- names(rho)[i]                       
            ##dbg("[cor_scatter.PLOTFUN] 4: gene2 = ", gene2)
            ##dbg("[cor_scatter.PLOTFUN] 4: gene2.in.X = ", gene2 %in% rownames(ngs$X))
            if(swapaxis) {
                x <- ngs$X[gene2,]
                y <- ngs$X[this.gene,]
                xlab = gene2
                ylab = this.gene
            } else {
                y <- ngs$X[gene2,]
                x <- ngs$X[this.gene,]
                ylab = gene2
                xlab = this.gene
            }
            plot(x, y, pch=19, cex=cex, col = klr,
                 ylab=ylab, xlab=xlab)
            
            y <- y + 1e-3*rnorm(length(y))
            x <- x + 1e-3*rnorm(length(x))            
            abline(lm(y ~ x), col="black", lty=2)

            if(i%%5==1) {
                ##legend("topleft", legend=levels(ph),
                ##       fill=klrpal, cex=0.85,
                ##       x.intersp=0.8, y.intersp=0.8)
                tt <- c("   ",levels(ph))
                legend("topleft", legend=tt,
                       fill = c(NA,klrpal), inset=c(0.02,0.02),
                       border = c(NA,rep("black",99)),
                       cex=0.9, box.lwd=0, pt.lwd=0,
                       x.intersp=0.5, y.intersp=0.8)
                legend("topleft", colorby, x.intersp=-0.2,
                       cex=0.9, y.intersp=0.45, bty='n')
            }

            
        }

    })

    cor_scatter.opts <- tagList(
        checkboxInput(ns("corscatter.swapaxis"),"swap axes")
        ## selectInput(ns("corscatter.colorby"),"color by:", choices=NULL),        
    )

    cor_scatter.info = "<b>Correlation scatter plots.</b> Pairwise scatter plots for the co-expression of correlated gene pairs across the samples. The straight line correspond to the (linear) regression fit."
    
    callModule(
        plotModule,
        id = "cor_scatter", label = "c",
        func = cor_scatter.PLOTFUN,
        func2 = cor_scatter.PLOTFUN,        
        info.text = cor_scatter.info,
        caption2 = cor_scatter.info,        
        options = cor_scatter.opts,
        title = "Correlation scatter plots",
        pdf.width = 12, pdf.height = 12, 
        height = c(fullH-80,760),
        width = c('auto',900),
        res = c(80,95),
        add.watermark = WATERMARK        
    )


    ##----------------------------------------------------------------------
    ##  Cumulative correlation plot (stacked)
    ##----------------------------------------------------------------------

    cum_corplot_data <- reactive({
        
        ngs <- inputData()
        req(ngs)	
        R <- getGeneCorr()
        
        ## get top correlated genes
        ##jj = head(order(rowSums(R),decreasing=FALSE),35)
        rsum <- rowSums(R,na.rm=TRUE)
        jj = head(order(abs(rsum),decreasing=TRUE),35)
        jj = head(order(-abs(rsum)),30)
        jj <- c( head(order(rsum),15), head(order(-rsum),15))
        jj <- jj[order(-rsum[jj])]
        head(rsum[jj])
        Rtop = R[jj,,drop=FALSE]
        rownames(Rtop) = sub(".*:","",rownames(Rtop))
        offset = min(Rtop, na.rm=TRUE)*0.95
        offset=0
        klr <- grey.colors(ncol(Rtop),start=0.3,end=0.7)
        
        ## --- color test -----##
        klr <- grey.colors(ncol(Rtop),start=0.3,end=0.7)
        klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.2)), alpha = TRUE)(ncol(Rtop))
        
        dbg("[cum_corplot_data()] done!")

        res = list(Rtop=Rtop, offset=offset, klr=klr)
        return(res)
    })
    
    cum_corplot.RENDER %<a-% reactive({

        dbg("[cum_corplot.RENDER] reacted")
        
        res <- cum_corplot_data()
        if(is.null(res)) return(NULL)
        
        par(mar=c(6,4,2,1), mgp=c(2.2,0.8,0))        
        mar=MARGINS1
        par(mar=mar, mgp=c(1.5,0.5,0))
        
        barplot( t(res$Rtop) - res$offset, col=res$klr, border=NA, ##horiz=TRUE, 
                las=3, cex.names=0.73, ##names.arg=rep(NA,nrow(R)),
                offset = res$offset, ylab="cumulative correlation (r)"            
                ##cex.main=1.2, main="cumulative correlation\nwith other data sets"
                )
        if(!is.null(colnames(res$Rtop))) {
            legend("topright", legend=rev(colnames(res$Rtop)), fill=rev(res$klr),
                   cex=0.8, y.intersp=0.8)
        }
        dbg("[cum_corplot.RENDER] done!")
        
    })

    cum_corplot_text = paste0('Top cumulative positively and negatively correlated genes with the selected gene in the current dataset as well as in public datasets such as ',a_ImmProt,' and ',a_HPA,'. The correlations of genes are colored by dataset.')
    
    callModule(
        plotModule, "cum_corplot",
        func = cum_corplot.RENDER,
        func2 = cum_corplot.RENDER,
        info.text = cum_corplot_text,
        height = 400,
        pdf.width=8, pdf.height=6,
        label="f",
        title="Cumulative correlation",
        add.watermark = WATERMARK
    )
    
    ##--------------------------------------------------------------------------------
    ## Correlation table
    ##--------------------------------------------------------------------------------

    cor_table.RENDER <- reactive({

        ngs <- inputData()
        req(ngs)

        dbg("[cor_table.RENDER] reacted!")

        R <- getGeneCorr()
        dbg("[cor_table.RENDER] nrow.R = ", nrow(R))        
        if(is.null(R)) return(NULL)
        
        P <- getPartialCorrelation()        
        pcor <- P[match(rownames(R),rownames(P)),"pcor"]
        
        title <- ngs$genes[rownames(R),"gene_title"]
        title <- substring(title, 1, 80)
        df <- data.frame(gene=rownames(R), title=title, cor=R[,"cor"], pcor=pcor)
        
        numeric.cols = colnames(df)[3:ncol(df)]
        ##selectmode <- ifelse(input$corGSEAtable_multiselect,'multiple','single')

        dbg("[cor_table.RENDER] dim.df = ",dim(df))
        dbg("[cor_table.RENDER] colnames.df = ",colnames(df))
        dbg("[cor_table.RENDER] rendering...")
        
        DT::datatable(
                df, rownames=FALSE, ## escape = c(-1),
                extensions = c('Buttons','Scroller'),
                ##selection=list(mode='multiple', target='row', selected=c(1)),
                selection=list(mode='single', target='row', selected=c(1)),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', 
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE,
                    deferRender=TRUE
                )  ## end of options.list 
            )  %>%
            formatSignif(numeric.cols,4)  %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
            DT::formatStyle("cor",
                            background = color_from_middle(
                                df[,"cor"],'lightblue', '#f5aeae'),
                            backgroundSize = '98% 88%',
                            backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center')
        
    })

    cor_table.info = "<b>DGCA table.</b> Statistical results from the DGCA computation for differentially correlated gene pairs."

    cor_table.opts <- tagList(
        ##checkboxInput(ns("cor_table.full"),"full table")
    )
    
    cor_table <- callModule(
        tableModule, 
        id = "cor_table", 
        func = cor_table.RENDER,
        options = cor_table.opts,
        info.text = cor_table.info,
        caption2 = cor_table.info,
        title = "Correlation table", label="b",
        height = c(360,700), width=c('auto',1400)
        ##caption = dgca_caption
    )    


    ##-----------------------------------------------------------
    ## Correlation network
    ##-----------------------------------------------------------
    
    getCorGraph <- reactive({

        req(input$cor_gene)

        dbg("[getCorGraph] reacted!")
        
        res <- getPartialCorrelationMatrix()
        gene="XIST";rho.min=0.3;layout="kk"
        gene <- input$cor_gene        
        dbg("[getCorGraph] gene = ",gene)

        rho.min <- input$cor_graph_threshold
        dbg("[getCorGraph] rho.min = ",rho.min)        
        
        layout <- input$cor_graph_layout
        ##fixed <- input$cor_graph_fixed
        numnodes <- nrow(res$cor)
        vsize = ifelse(numnodes > 50, 10, 12)
        vsize = ifelse(numnodes > 100, 8, vsize)
        
        radius = as.integer(input$cor_graph_radius)

        dbg("[getCorGraph] compiling...")        
        gr <- pgx.plotPartialCorrelationGraph(
            res, gene, ##what="graph", ## degree=deg,
            plot = FALSE,
            rho.min = rho.min, nsize = -1,
            layout = layout, radius = radius,
            vsize = vsize, edge.width = 10)
        gr
    })

    cor_graph.PLOTFUN %<a-% reactive({
        gr <- getCorGraph()
        par(mar=c(1,1,1,1)*0)
        plot(gr,
             ##layout = L1,
             ##vertex.label.cex = label.cex,
             ##vertex.size = 12,
             vertex.color = 'lightskyblue1',
             vertex.frame.color = 'skyblue'
             )
    })

    cor_graph.VISNETWORK <- reactive({
        require(visNetwork)
        gr <- getCorGraph()
        if(is.null(gr)) return(NULL)                
        
        visdata <- toVisNetworkData(gr, idToLabel=FALSE)
        visdata$edges$width <- 2*visdata$edges$width 
        
        graph <- visNetwork::visNetwork(
                                 nodes = visdata$nodes,
                                 edges = visdata$edges) %>%
            visNodes(font = list(size = 12))  %>%
            visEdges(hidden=FALSE, width=4, color=list(opacity=0.9))  %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            ## visPhysics(enabled=TRUE)  %>%
            ## visInteraction(hideEdgesOnDrag = TRUE) %>%
            ## visIgraphLayout(layout="layout.norm", layoutMatrix=pos)
            visIgraphLayout(layout="layout_nicely")
        graph

    })
    
    GRAPH.LAYOUTS = c("Fruchterman-Reingold"="fr", "Kamada-Kawai"="kk",
                      "graphopt"="graphopt","tree layout"="tree")

    cor_graph.opts <- tagList(
        sliderInput(ns('cor_graph_radius'),'radius:', 1, 8, 3, 1),
        sliderInput(ns('cor_graph_threshold'),'pcor threshold:', 0, 1, 0.90),        
        selectInput(ns('cor_graph_layout'),'layout:', choices=GRAPH.LAYOUTS)
    )

    cor_graph_info <- "<b>Gene correlation network.</b> Correlation graph centered on selected gene with top most correlated features. Red edges correspond to negative (marginal) correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation of the gene pair."
    cor_graph_info <- "<b>Partial correlation network.</b> Partial correlation graph centered on selected gene with top most correlated features. Red edges correspond to negative correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation value of the gene pair."
    
    callModule(
        plotModule,
        id = "cor_graph", label = "a",
        ##func = cor_graph.PLOTFUN,
        func = cor_graph.VISNETWORK,        
        ##func2 = cor_graph.PLOTFUN,        
        plotlib = 'visnetwork',
        info.text = cor_graph_info,
        options = cor_graph.opts,
        title = "PARTIAL CORRELATION NETWORK",
        ## caption = topEnriched_caption
        caption2 = cor_graph_info,
        ##pdf.width = 14, pdf.height = 4, 
        height = c(700,720),
        width = c('auto',1000),
        res = c(72,80),
        add.watermark = WATERMARK
    )

    ##-----------------------------------------------------------
    ## Correlation UMAP
    ##-----------------------------------------------------------
    
    ##cor_umap.PLOTFUN %<a-% reactive({
    cor_umap.PLOTFUN <- reactive({    

        ngs <- inputData()
        req(ngs)
        req(input$cor_gene)
        
        dbg("[cor_umap.PLOTFUN] reacted!")
        if(!"cluster.genes" %in% names(ngs)) {
            par(mfrow=c(1,1))
            frame()
            text(0.5,0.6, "Error: gene cluster position in PGX object", col='red3')
            return(NULL)
        }
        
        R <- getGeneCorr()
        dbg("[cor_table.RENDER] nrow.R = ", nrow(R))        
        if(is.null(R)) return(NULL)
        gene='ESR1'
        gene <- input$cor_gene
        pos <- ngs$cluster.genes$pos[['umap2d']]
        if(input$umap_param=='cov') {
            rho <- R[,'cov']
        } else {
            rho <- R[,"cor"]
        }
        rho <- rho[match(rownames(pos),names(rho))]
        names(rho) <- rownames(pos)
        ##rho[is.na(rho)] <- 0
        
        higenes <- c(gene)
        higenes <- names(tail(sort(rho**2),20))
        higenes <- unique(names(c(head(sort(rho),10),tail(sort(rho),10))))
        cexlab <- ifelse(length(higenes)==1, 2.5, 1.6)

        if(0) {

            rho = R[,1]
            rho = c("neg","yes")[1+1*(sign(rho)==1)]
            names(rho) <- rownames(R)

            rho <- rho[match(rownames(pos),names(rho))]
            names(rho) <- rownames(pos)
            rho[is.na(rho)] <- 0
            
            pgx.scatterPlotXY(
                pos, var=rho, ##type="factor", col=c("blue","red"),
                hilight = higenes, hilight.cex=1.2,
                zsym = TRUE, softmax=1,
                ## hilight2 = hilight2,
                cex = 1, cex.lab = 2.0,
                legend = TRUE, plotlib="plotly")
            
        }
       
        dbg("[cor_umap.PLOTFUN] dim.pos =",dim(pos))
        dbg("[cor_umap.PLOTFUN] len.rho =",length(rho))        
        
        p <- pgx.plotGeneUMAP(
            ngs, pos = pos, ##contrast=ct,
            value = rho, title='',
            cex = 0.9, cex.lab = cexlab, 
            hilight = higenes, ntop = 20,
            plotlib = "plotly")

        ##p
        if(!is.null(p)) return(p)
    })

    cor_umap.opts <- tagList(
        radioButtons(ns('umap_param'),'parameter:', choices=c("cor","cov"), inline=TRUE)
    )

    cor_umap_info <- "<b>Gene UMAP.</b> Partial correlation graph centered on selected gene with top most correlated features. Red edges correspond to negative correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation value of the gene pair."
    
    callModule(
        plotModule,
        id = "cor_umap", label = "b",
        func = cor_umap.PLOTFUN,
        ##func2 = cor_umap.PLOTFUN,
        plotlib = 'plotly',
        ##plotlib = 'ggplot',
        info.text = cor_umap_info,
        options = cor_umap.opts,
        title = "CORRELATION UMAP",
        ## caption = topEnriched_caption
        caption2 = cor_umap_info,
        ##pdf.width = 14, pdf.height = 4, 
        height = c(700,750),
        width = c('auto',900),
        res = c(72,80),
        add.watermark = WATERMARK
    )
    
    ##--------------------------------------------------------------------------------
    ## Correlation GSEA
    ##--------------------------------------------------------------------------------

    getCorrelationGSEA <- reactive({

        ngs <- inputData()
        alertDataLoaded(session,ngs)
        req(ngs)

        pgx.showSmallModal("Calculating GSEA...<br>please wait")
        
        gene = "CD4"
        gene = rownames(ngs$X)[1]
        gene <- input$cor_gene

        ## single gene correlation as rank metric
        gx <- ngs$X[gene,]
        rho <- cor(t(ngs$X), gx, use="pairwise")[,1]         
        names(rho) <- toupper(names(rho))
        
        gmt <- GSETS[colnames(ngs$GMT)]
        ## gmt <- GSETS  ## all???
        gsea <- fgsea(gmt, rho, minSize=15, maxSize=1000)
        gsea <- gsea[order(-gsea$NES),]
        head(gsea)

        removeModal()
        
        res <- list(gsea=gsea, rho=rho)
        return(res)
    })

    
    corGSEA_plots.RENDER %<a-% reactive({
        require(RColorBrewer)
        res = getCorrelationGSEA()
        ##if(is.null(rho)) return(NULL)

        ii <- corGSEA_table$rows_all()
        req(ii)
        gsea <- res$gsea[ii,,drop=FALSE]
        
        ## ENPLOT TYPE
        NTOP = 20
        par(oma=c(0,2,0,0))
        par(mfrow=c(4,5), mar=c(1,1.5,2.3,1))
        i=1
        for(i in 1:min(NTOP,nrow(gsea))) {
            gs <- gsea$pathway[i]
            gs
            gmt <- GSETS[[gs]]
            length(gmt)
            ##if(length(gmtdx) < 3) { frame(); next }
            ylab = ""
            if( i%%5 == 1) ylab = "correlation (rho)"            
            gsea.enplot( res$rho, gmt, xlab="", ylab=ylab,
                        main=substring(gs,1,58),
                        cex.main=0.92, len.main=30 )
            nes <- round(gsea$NES[i],2)
            qv  <- round(gsea$padj[i],3)
            tt <- c( paste("NES=",nes), paste("q=",qv) )
            legend("topright", legend=tt, cex=0.9)
        }
        
    })

    corGSEA_plots_opts = tagList(
        tipify( selectInput( ns("xann.refset"), "Reference set:", choices="", width='80%'),
               "Specify a reference set to be used in the annotation.",
               placement="left",options = list(container = "body"))
    )

    corGSEA_plots_info = "<b>Correlation GSEA.</b> Functional GSEA enrichment of correlated genes. Black vertical bars indicate the rank of genes in the gene set in the sorted correlation metric. The green curve corresponds to the 'running statistics' of the enrichment score (ES). The more the green ES curve is shifted to the upper left of the graph, the more the gene set is enriched in the first group. Conversely, a shift of the ES curve to the lower right, corresponds to more enrichment in the second group."

    
    ##corGSEA_plots_module <- plotModule(
    callModule(
        plotModule, 
        id = "corGSEA_plots", ##ns=ns,
        func = corGSEA_plots.RENDER,
        func2 = corGSEA_plots.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = corGSEA_plots_opts,
        info.text = corGSEA_plots_info,        
        title="Correlation GSEA", label="a",
        height = c(0.57*fullH,700), width = c('auto',1400),
        pdf.width=12, pdf.height=7, res=c(72,100),
        add.watermark = WATERMARK
    )
    ## output <- attachModule(output, corGSEA_plots_module)
    
    corGSEA_table.RENDER <- reactive({
        
        res = getCorrelationGSEA()
        
        ##rho = data.frame(cbind( name=rho.name, rho))
        gs <- res$gsea$pathway
        ##link <- wrapHyperLink(rep("link",length(gs)), gs)
        link <- wrapHyperLink(rep("&#x1F517;",length(gs)), gs)        
        df = data.frame( pathway=gs, link=link,
                        res$gsea[,c("NES","pval","padj","size")] )
        df$pathway <- shortstring(df$pathway,80)
        numeric.cols = c("pval","padj","NES")
        ##selectmode <- ifelse(input$corGSEAtable_multiselect,'multiple','single')
        
        DT::datatable(
                df, rownames=FALSE, escape = c(-1,-2),
                extensions = c('Buttons','Scroller'),
                ##selection=list(mode='multiple', target='row', selected=c(1)),
                selection=list(mode='single', target='row', selected=c(1)),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', 
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE,
                    deferRender=TRUE
                )  ## end of options.list 
            )  %>%
            formatSignif(numeric.cols,4)  %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("NES",
                                background = color_from_middle(df[,"NES"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    corGSEA_table_info = "<b>Enrichment table.</b> Statistical results from the GSEA computation for functional enrichment of correlated genes. The column 'pval' and 'padj' correspond to the p-value and (multiple testing) adjusted p-value of the GSEA test, respectively. The 'NES' column reports the normalized enrichment score."

    corGSEA_table_opts <- tagList(
        checkboxInput(ns("corGSEAtable_multiselect"),"enable multi-select")
    )
    
    corGSEA_table <- callModule(
        tableModule, 
        id = "corGSEA_table", 
        func = corGSEA_table.RENDER,
        ## options = corGSEA_table_opts,
        info.text = corGSEA_table_info,
        title = "Correlation GSEA table", label="b",
        height = c(220,700), width=c('auto',1000)
        ##caption = corGSEA_caption
    )    

    corGSEA_cumFC.RENDER %<a-% reactive({

        require(RColorBrewer)
        res = getCorrelationGSEA()
        ##if(is.null(rho)) return(NULL)

        ii <- 1:nrow(res$gsea)
        ii <- corGSEA_table$rows_all()
        req(ii)
        ii <- head(ii,20)
        le.genes <- res$gsea[ii,]$leadingEdge
        names(le.genes) <- res$gsea[ii,]$pathway
        all.le <- unique(unlist(le.genes))
        
        F <- 1*sapply(le.genes, function(g) all.le %in% g)
        rownames(F) <- all.le
        colnames(F) <- names(le.genes)
        if(0) {
            if(input$gs_enrichfreq_gsetweight) {
                F <- t(t(F)  / colSums(F,na.rm=TRUE))
            }
            F <- t(t(F) * sign(fx[top]))
            if(input$gs_enrichfreq_fcweight) {
                F <- t(t(F) * abs(fx[top]))
            }
        }

        F <- head(F[order(-rowSums(F**2)),,drop=FALSE],30)
        F <- F[order(-rowSums(F)),,drop=FALSE]
        F <- as.matrix(F)
       
        par(mfrow=c(1,1), mar=c(6,4,2,0.5), mgp=c(2,0.8,0))
        col1 = grey.colors(ncol(F),start=0.15)
        barplot(t(F), beside=FALSE, las=3, cex.names=0.80, col=col1,
                ylab="frequency")
        
    })

    corGSEA_cumFC.opts = tagList()
    corGSEA_cumFC.info = "<b>Leading-edge gene frequency.</b> Number of genesets that include the gene in their leading-edge. Genes with higher frequency are shared more often among the top genesets and may indicate a driver gene."    

    ##corGSEA_cumFC_module <- plotModule(
    callModule(
        plotModule, 
        id = "corGSEA_cumFC", ##ns=ns,
        func = corGSEA_cumFC.RENDER,
        func2 = corGSEA_cumFC.RENDER, 
        download.fmt = c("png","pdf"),
        ## options = corGSEA_cumFC.opts,
        info.text = corGSEA_cumFC.info,        
        caption2 = corGSEA_cumFC.info,        
        title = "Leading-edge gene frequency", label="c",
        height = c(280,650), width = c('auto',1000),
        pdf.width=10, pdf.height=6, res=c(72,90),
        add.watermark = WATERMARK
    )

    
    corGSEA_LeadingEdgeTable.RENDER <- reactive({
        
        res = getCorrelationGSEA()
        sel=1
        sel <- corGSEA_table$rows_selected()
        req(sel)
        le.genes <- res$gsea[sel[1],]$leadingEdge[[1]]
        if(length(sel)>1) {
            for(i in 2:length(sel)) {
                le.genes2 <- res$gsea[sel[i],]$leadingEdge[[1]]
                le.genes <- intersect(le.genes, le.genes2)
            }
        }
        
        ##rho = data.frame(cbind( name=rho.name, rho))
        rho1 <- res$rho[le.genes]
        title <- shortstring(GENE.TITLE[le.genes],50)
        df = data.frame( gene = le.genes, rho = rho1, title=title)
        rownames(df) = le.genes
        numeric.cols = c("rho")
        
        DT::datatable(
                df, rownames=FALSE, ## escape = c(-1,-2),
                extensions = c('Buttons','Scroller'),
                selection=list(mode='single', target='row', selected=c(1)),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', 
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE,
                    deferRender=TRUE
                )  ## end of options.list 
            )  %>%
            formatSignif(numeric.cols,4)  %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("rho",
                                background = color_from_middle(df[,"rho"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    corGSEA_LeadingEdgeTable.info = "<b>Leading-edge table</b> Leading edge genes as reported by GSEA corresponding to the selected geneset. The 'rho' column reports the correlation with respect to the query gene."
    
    corGSEA_LeadingEdgeTable <- callModule(
        tableModule, 
        id = "corGSEA_LeadingEdgeTable", 
        func = corGSEA_LeadingEdgeTable.RENDER,
        info.text = corGSEA_LeadingEdgeTable.info,
        title = "Leading edge genes", label="d",
        height = c(378,700), width=c('auto',1000)
        ##height = c(665,700), width=c('auto',1000)
        ##caption = corGSEA_caption
    )
   
    ##-----------------------------------------------------------
    ## Function for Differential Gene Correlation Analysis
    ##-----------------------------------------------------------

    dgca.output <- reactive({
        
        ngs <- inputData()
        gene=NULL
        gene="IRF4";ph='activated'
        gene="KCNN4";ph='activated'
        gene="ERBB2";ph='HER2_STATUS'
        gene="ESR1";ph='ER_STATUS'
        gene=NULL;ph='ER_STATUS'
        gene=NULL;ph='state'        
        
        gene <- input$cor_gene
        ph   <- input$cor_group
        req(gene)
        req(ph)        
        if(input$dgca.allpairs) gene <- NULL
        
        message("[dgca.output] gene= ",gene)
        message("[dgca.output] ph= ",ph)
        
        grp = factor(ngs$samples[,ph])
        ii <- which(!is.na(grp))
        grp <- grp[ii]
        D = model.matrix( ~ 0 + grp)
        colnames(D) <- sub("^grp","",colnames(D))
        rownames(D) <- rownames(ngs$samples)[ii]
        dim(D)

        X <- ngs$X
        X <- getFilteredExpression()        
        X <- X[,ii,drop=FALSE]
        sdx1 <- apply(X[,D[,1]==1],1,sd)
        sdx2 <- apply(X[,D[,2]==1],1,sd)        
        jj <- order(-sdx1*sdx2)
        if(!is.null(gene)) jj <- unique(c(match(gene,rownames(X)),jj))

        message("[dgca.output] length(jj)= ",length(jj))
        X <- X[jj,,drop=FALSE]        
        if(!is.null(gene))  {
            ii <- head(unique(c(gene,rownames(X))),1000)
            X <- X[ii,]
        } else {
            X <- head(X, 100)
        }
        message("[dgca.output] dimX= ",paste(dim(X),collapse='x'))
        message("[dgca.output] gene.inX= ",gene %in% rownames(X))
        message("[dgca.output] dimD= ",paste(dim(D),collapse='x'))
        
        require(DGCA)
        res = DGCA::ddcorAll(inputMat = X,
                             splitSet = gene,
                             design = D,
                             compare = levels(grp),
                             nPairs = "all",
                             adjust = "fdr",
                             nPerms = 0
                             )
        dim(res)
        message("[dgca.output] 1 : dim.res= ",paste(dim(res),collapse='x'))
        
        ## remove non-valid
        res <- res[rowSums(is.na(res))==0,]
        
        message("[dgca.output] 2 : dim.res= ",paste(dim(res),collapse='x'))
        dim(res)
        
        ## filter significant only??
        if(0) {
            sel <- (res$pValDiff_adj < 0.05)
            res <- res[sel,,drop=FALSE]
            message("[dgca.output] 3 : dim.res= ",paste(dim(res),collapse='x'))
        }
        
        dim(res)
        res$avg_cor <- rowMeans(res[,c(3,5)])
        
        head(res,20)
        message("[dgca.output] done!")
        res

    })

    
    dgca_barplot.PLOTFUN %<a-% reactive({

        res <- dgca.output()        

        ii  <- dgca_table$rows_all()
        if(is.null(ii) || length(ii)==0) return(NULL)
        res <- res[ii,,drop=FALSE]
        ##res <- res[rowSums(is.na(res))==0,]
        
        rho <- res[,c(3,5)]
        pv  <- res[,c(4,6)]        
        rho[is.na(rho)] <- 0
        pv[is.na(pv)] <- 1
        rho <- as.matrix(rho * exp(-pv/0.05))
        rownames(rho) <- res$Gene1
        if(input$dgca.allpairs) {
            rownames(rho) <- paste0(res$Gene1," x ",res$Gene2)
        }

        zerorho <- matrix(0,nrow=20,ncol=ncol(rho))
        rownames(zerorho) <- NULL
        rho <- rbind(rho, zerorho)
        
        NTOP=40
        rho <- head(rho,NTOP)
        ##rho <- head(rho[order(-rowMeans(rho**2)),],NTOP)
        if(input$dgca_barplot.meansort) {
            rho <- rho[order(-rowMeans(rho)),]
        }

        ##klr.pal = rev(grey.colors(2)),
        klr.pal <- brewer.pal(4,"Paired")[1:2]
        klr.pal <- COL2
        
        par(mfrow=c(1,1), mar=c(12,4,1,1))
        barplot(t(rho),
                col = klr.pal, ylim=c(-1,1)*1.1,
                beside=TRUE, las=3,
                xlab = "",
                ylab = "p-weighted correlation",
                cex.names = 0.95 )

        grp=""
        grp <- input$cor_group
        ## title(paste("by",grp), cex=1)

        vv <- sub("","",colnames(rho))
        legend("topright", legend=vv, fill=klr.pal,
               cex=0.9, y.intersp=0.9, inset=c(0.03,0) )
        
    })

    dgca_barplot.opts <- tagList(
        checkboxInput(ns("dgca_barplot.meansort"),"sort on mean correlation")
    )

    dgca_barplot.info = "<b>Differentially correlated gene pairs.</b> The height of the bars correspond to the Pearson correlation value of a gene pair in each group. Differentially correlated gene pairs show different correlation value in each group.."
   
    callModule(
        plotModule,
        id = "dgca_barplot", 
        func = dgca_barplot.PLOTFUN,
        func2 = dgca_barplot.PLOTFUN,        
        info.text = dgca_barplot.info,
        options = dgca_barplot.opts,
        title = "Top correlated genes", label = "a",
        caption2 = dgca_barplot.info,
        pdf.width = 10, pdf.height = 5, 
        height = c(0.45*fullH,700),
        width = c('auto',1200),
        res = c(63,100),
        add.watermark = WATERMARK
    )

    ##-----------------------------------------------------------
    ## DGCA Table
    ##-----------------------------------------------------------

    dgca_table.RENDER <- reactive({
        
        res <- dgca.output()        
        
        df <- res
        ## Gene1,Gene2,Classes,0_cor,1_cor,avg_cor,zScoreDiff
        ## sel <- c("Gene1","Gene2","activated_cor","activated_pVal","resting_cor","resting_pVal","zScoreDiff","pValDiff","pValDiff_adj","Classes")
        if(input$dgca_table.full) {
            df <- res[,c(1,2,10,3,5,11,7,4,6,8,9)]
        } else {
            df <- res[,c(1,2,10,3,5,11,7)]
        }

        numeric.cols = colnames(df)[4:ncol(df)]
        ##selectmode <- ifelse(input$corGSEAtable_multiselect,'multiple','single')
        
        DT::datatable(
                df, rownames=FALSE, escape = c(-1,-2),
                extensions = c('Buttons','Scroller'),
                ##selection=list(mode='multiple', target='row', selected=c(1)),
                selection=list(mode='single', target='row', selected=c(1)),
                class = 'compact cell-border stripe hover',
                fillContainer = TRUE,
                options=list(
                    dom = 'lfrtip', 
                    ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                    scrollX = TRUE, ##scrollY = TRUE,
                    ##scrollY = 170,
                    scrollY = '70vh',
                    scroller=TRUE,
                    deferRender=TRUE
                )  ## end of options.list 
            )  %>%
            formatSignif(numeric.cols,4)  %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle("zScoreDiff",
                                background = color_from_middle(df[,"zScoreDiff"],
                                                               'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    dgca_table.info = "<b>DGCA table.</b> Statistical results from the DGCA computation for differentially correlated gene pairs."

    dgca_table.opts <- tagList(
        checkboxInput(ns("dgca_table.full"),"full table")
    )
    
    dgca_table <- callModule(
        tableModule, 
        id = "dgca_table", 
        func = dgca_table.RENDER,
        options = dgca_table.opts,
        info.text = dgca_table.info,
        caption2 = dgca_table.info,
        title = "DGCA table", label="b",
        height = c(305,700), width=c('auto',1400)
        ##caption = dgca_caption
    )    

    ##-----------------------------------------------------------------------------
    ## DGCA Correlation scatter plots
    ##-----------------------------------------------------------------------------
            
    gene1="MKI67"
    gene2="NCAPD2"
    
    dgca.scatterplot <- function(X, gene1, gene2, grp, cex=1, key=TRUE,
                          rho=TRUE, col=c("grey40","red2"))
    {

        grp.levels = sort(unique(setdiff(as.character(grp),NA)))
        x1 <- X[gene1,which(grp==grp.levels[1])]
        y1 <- X[gene2,which(grp==grp.levels[1])]
        x2 <- X[gene1,which(grp==grp.levels[2])]
        y2 <- X[gene2,which(grp==grp.levels[2])]

        col1 <- paste0(col2hex(col),"AA") ## add opacity
        rgb2col <- function(cc) rgb(cc[1],cc[2],cc[3],maxColorValue=255)
        col2 <- apply( 0.66*col2rgb(col),2,rgb2col)
        
        ##par(mfrow=c(2,2), mar=c(4,4,2,2), oma=c(0,0,0,0))
        par(mgp=c(1.6,0.6,0))
        plot(x1, y1, col=col1[1], pch=19, cex=cex,
             xlim=range(c(x1,x2)), ylim=range(c(y1,y2)),
             ##main = 'differential correlation',
             xlab=gene1, ylab=gene2)

        y1 <- y1 + 1e-3*rnorm(length(y1))
        x1 <- x1 + 1e-3*rnorm(length(x1))            
        y2 <- y2 + 1e-3*rnorm(length(y2))
        x2 <- x2 + 1e-3*rnorm(length(x2))            

        abline(lm(y1 ~ x1), col=col2[1], lty=1, lwd=1.5)
        points(x2, y2, col=col1[2], pch=19, cex=cex)
        abline(lm(y2 ~ x2), col=col2[2], lty=1, lwd=1.5)

        if(key) {
            legend("topleft", legend=grp.levels,
                   fill=col, bty='n',
                   cex=0.9, y.intersp=0.8)
        }
        if(rho) {
            r1 <- round( cor(x1,y1), 3)
            r2 <- round( cor(x2,y2), 3)
            rr <- paste0("R_",grp.levels," = ", c(r1,r2))
            legend("bottomright",
                   legend = rr,
                   ## fill=c('red','blue'),
                   bty='n',
                   cex=0.9, y.intersp=0.85)
        }
    }
    
    dgca_scatter.PLOTFUN %<a-% reactive({
        
        message("[dgca_scatter.PLOTFUN] reacted!")

        req(input$cor_gene)
        res <- dgca.output()

        ii  <- dgca_table$rows_all()
        if(is.null(ii) || length(ii)==0) return(NULL)
        res1 <- res[ii,,drop=FALSE]
        ##res1 <- res1[rowSums(is.na(res1))==0,]
        
        ph='ER_STATUS';
        ngs <- inputData()
        ph  <- input$cor_group
        req(ph)
        
        grp  <- factor(ngs$samples[,ph])
        ndim <- nrow(ngs$samples)
        sel1 <- which(as.integer(grp)==1)
        sel2 <- which(as.integer(grp)==2)
        NTOP = 25        
        
        if(0) {
            sdx1 <- apply(ngs$X[res1$Gene1,sel1],1,sd) * apply(ngs$X[res1$Gene2,sel1],1,sd)
            sdx2 <- apply(ngs$X[res1$Gene1,sel2],1,sd) * apply(ngs$X[res1$Gene2,sel2],1,sd)        
            score <- res1$zScoreDiff * (sdx1 * sdx2)**0.1
            res1 <- res1[order(-abs(score)),]
        }
        
        message("[dgca_scatter.PLOTFUN] start render")
        nplots <- min(NTOP,nrow(res1))
        nr = ceiling(sqrt(nplots))
        par(mfrow=c(nr,nr), mar=c(3,3.5,0.5,0),
            mgp=c(1.9,0.7,0), oma=c(0,0,0.5,0.5))

        cex = 1.2
        cex = ifelse( ndim > 40, 0.8, cex)
        cex = ifelse( ndim > 100, 0.5, cex)
        cex = ifelse( ndim > 200, 0.2, cex)                
        klrpal = c("grey30","red2")
        klrpal = COL2
        klrpal = COL  ## more colors        
        if(nr <= 2) cex = cex*2
        
        i=1
        for(i in 1:nplots) {
            gene1 <- res1$Gene2[i]                       
            gene2 <- res1$Gene1[i]
            X1 <- ngs$X[c(gene1,gene2),]
            dgca.scatterplot(X1, gene1, gene2, grp=grp, cex=cex,
                             key=0, rho=1, col=klrpal) 
            if(i==1) {
                tt <- c("   ",levels(grp))
                legend("topleft", legend=tt,
                       fill = c(NA,klrpal), inset=c(0.01,0.01),
                       border = c(NA,"black","black"),
                       cex=0.9, box.lwd=0, pt.lwd=0,
                       x.intersp=0.5, y.intersp=0.8)
                legend("topleft", ph, x.intersp=-0.2, inset=c(0.01,0.01),
                       cex=0.9, y.intersp=0.45, bty='n')
            }
        }

        message("[dgca_scatter.PLOTFUN] done!")
        
    })

    dgca_scatter.opts <- tagList(
        ##checkboxInput(ns("dgcascatter.swapaxis"),"swap axes"),
        ##selectInput(ns("dgcascatter.colorby"),"color by:", choices=NULL),        
    )

    dgca_scatter.info = "<b>DGCA scatter plots.</b> Pairwise scatter plots for the co-expression of the gene pairs in two different conditions. Differentially correlated gene pairs will show different correlation values measured by the difference in their z-score ('zScoreDiff'). The straight lines correspond to the linear regression fits. "
    
    callModule(
        plotModule,
        id = "dgca_scatter", label = "c",
        func = dgca_scatter.PLOTFUN,
        func2 = dgca_scatter.PLOTFUN,        
        info.text = dgca_scatter.info,
        options = dgca_scatter.opts,
        title = "DGCA correlation scatter plots",
        ##pdf.width = 14, pdf.height = 4, 
        height = c(fullH-80,760),
        width = c('auto',900),
        res = c(80,95),
        add.watermark = WATERMARK
    )
    
} ## end of Board




if(0) {
    source(file.path(RDIR,"pgx-include.R"))    ## lots of libraries and source()

    load("../data/tcga-brca2-gx.pgx")
    load("../data/geiger2016-arginine.pgx")
    colnames(ngs$samples)
    grp = ngs$samples[,"ER_STATUS"]
    table(grp)
    
}
