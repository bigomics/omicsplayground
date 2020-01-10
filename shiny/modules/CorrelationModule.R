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
            id = ns("tabs1"),
            tabPanel("Correlation",uiOutput(ns("corAnalysis_UI")))
        )
    )    
    if(DEV.VERSION) {
        ui <- fillCol(
            height = 750,
            tabsetPanel(
                id = ns("tabs1"),
                tabPanel("Correlation",uiOutput(ns("corAnalysis_UI"))),
                tabPanel("Functional",uiOutput(ns("corFunctional_UI")))
            )
        )
    }
    ui
}

CorrelationModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    inputData <- env[["load"]][["inputData"]]
    usermode  <- env[["load"]][["usermode"]]

    fullH = 750  ## full height of page
    rowH  = 340  ## full height of page
    
    description = "<b>Correlation Analysis.</b> Statistical correlation analysis on a
gene or gene set level with visualisations. Compute the correlation
between genes and find coregulated modules."
    output$description <- renderUI(HTML(description))


    cor_infotext ="The <strong>Correlation Analysis Module</strong> provides statistical correlation analysis on gene level or gene set level with visualisations. During the visual analysis, users can filter out some samples or collapse the samples by predetermined groups. It also correlates the gene to the expressions of other genes across datasets such as ImmProt and HPA, and plots the cumulative correlation. Furthermore, it displays the tissue expression for a selected gene using the genotype-tissue expression (GTEx) dataset."
    

    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================

    corAnalysis_caption = "<b>Correlation analysis.</b> <b>(a)</b> Top correlated features with selected gene. <b>(b)</b> Correlation network around the selected gene. <b>(c)</b> Scatter plots of top correlated expression features."
    
    output$corAnalysis_UI <- renderUI({
        fillCol(
            flex = c(1,NA),
            height = fullH,
            fillRow(
                flex = c(1,0.1,1.2),
                fillCol(
                    flex = c(1,1),
                    height = fullH-80,
                    plotWidget(ns('cor_barplot')),
                    plotWidget(ns('cor_graph'))
                ),
                br(), ## spacer
                fillCol(
                    flex = c(NA,1),
                    height = fullH-80,
                    plotWidget(ns('cor_scatter'))
                )
            ),
            div(HTML(corAnalysis_caption), class="caption")
        )
    })


    corfunctional_caption = "<b>Correlation GSEA.</b> Functional annotation of the correlated genes as defined by Pearson correlation. <b>(a)</b> Top enriched gene sets using the correlation as rank metric. The black bars denote the genes in the gene set and their position in the sorted rank metric. <b>(b)</b> GSEA statistics table. <b>(c)</b> Leading edge table of the genes in the selected gene set."

    output$corFunctional_UI <- renderUI({

        fillCol(
            flex = c(1,NA),
            height = fullH,
            fillRow(
                flex = c(2,0.1,1),
                height = fullH - 60,
                fillCol(
                    flex = c(1,0.05,0.6),
                    plotWidget(ns("corGSEA_plots")),
                    br(),
                    tableWidget(ns("corGSEA_table"))
                ),
                br(),
                fillCol(
                    ##flex = c(1,0.05,0.65,NA),
                    flex = c(1),
                    ##plotWidget(ns("corGSEA_plots")),
                    ##br(),
                    ##br(),
                    tableWidget(ns("corGSEA_LeadingEdgeTable"))
                    ##div(HTML("caption"), class="caption")
                )
            ),
            div(HTML(corfunctional_caption), class="caption")
        )
    })
    ##outputOptions(output, "hm_annotateUI", suspendWhenHidden=FALSE) ## important!!!
    
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
            actionLink(ns("cor_options"), "Options", icon=icon("cog", lib = "glyphicon")),
            br(),br(),
            conditionalPanel(
                "input.cor_options % 2 == 1", ns=ns,
                tagList(
                    tipify( selectInput(ns("cor_features"),"Gene family:", choices=NULL, multiple=FALSE),
                           "Filter for a specific gene family.", placement="top"),
                    tipify( selectInput(ns("cor_samplefilter"),"Filter samples",
                                        choices=NULL, multiple=TRUE),
                           "Filter (include) samples for the analysis")                    
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
            title = HTML("<strong>Correlation Analysis Module</strong>"),
            HTML(cor_infotext),
            easyClose = TRUE ))
    })
    
    ## update filter choices upon change of data set 
    observe({
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        
        ## levels for sample filter
        levels = getLevels(ngs$Y)
        updateSelectInput(session, "cor_samplefilter", choices=levels)

        genes <- sort(ngs$genes[rownames(ngs$X),]$gene_name)
        sel = genes[1]  ## most var gene
        updateSelectInput(session,'cor_gene', choices=genes, selected=sel)

        fam <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        updateSelectInput(session, "cor_features",choices=fam)

    })

    ##================================================================================
    ##======================= PLOTTING FUNCTIONS =====================================
    ##================================================================================


    ##-----------------------------------------------------------
    ## Correlation barplot
    ##-----------------------------------------------------------

    getGenePartialCorrelation <- reactive({
        ngs <- inputData()
        req(ngs,input$cor_gene)
        X <- ngs$X    
        gene <- rownames(X)[1]
        gene
        gene <- input$cor_gene
        methods = c("cor")
        if(input$cor_partialpc == "fast") {
            methods = PCOR.FAST
        }
        if(input$cor_partialpc == "all methods") {
            methods = PCOR.METHODS
        }

        if(input$cor_features!="<all>") {
            ft <- input$cor_features
            psel = filterProbes(ngs$genes, GSETS[[ft]] )
            psel = unique(c(gene, psel))
            psel <- intersect(psel,rownames(X))
            X = X[psel,,drop=FALSE]
        }
        
        showNotification(paste("computing correlation...\n"))        
        res <- pgx.computePartialCorrelationAroundGene(
            X, gene, method=methods, nmax=100, fast=FALSE)    

        j <- which(rownames(res$cor)==gene)
        rho1 <- min(abs(head(res$cor[j,-j],20)))
        updateSliderInput(session, "cor_graph_threshold", value=rho1)

        res
    })
    
    cor_barplot.PLOTFUN %<a-% reactive({

        req(input$cor_gene)
        res <- getGenePartialCorrelation()

        gene <- rownames(res$cor)[1]
        gene <- input$cor_gene
        gene
        NTOP = 40
        prho <- head( res$meta.pcor[gene,], NTOP)
        rho  <- head( res$cor[gene,], NTOP)
        
        par(mfrow=c(1,1), mar=c(10,4,1,1))
        barplot(rho, beside=FALSE, las=3,
                ylab = "correlation" )
                
        if(!is.null(prho)) {
       
            prho[prho>0] <- pmin(prho[prho>0],rho[prho>0])
            prho[prho<0] <- pmax(prho[prho<0],rho[prho<0])
            barplot(prho, beside=FALSE, add=TRUE,
                    col="grey40", names.arg="")
            legend("topright",
                   c("correlation","partial correlation"),
                   fill=c("grey70","grey40"))
        }

    })

    cor_barplot.opts <- tagList(
        ##checkboxInput(ns('cor_partialpc'),'partial correlation', value=FALSE)
        radioButtons(ns('cor_partialpc'),'partial correlation',
                     c("none","fast","all methods"), selected="fast",
                     inline=TRUE )
    )

    cor_barplot.info = "<b>Top correlated genes.</b> Highest correlated genes in respect to the selected gene. The height of the bars correspond to the Pearson correlation value. The dark grey bars correspond to the 'partial correlation' which essentially corrects the correlation value for indirect effects and tries to estimate the amount of direct interaction."
    
    callModule(
        plotModule,
        id = "cor_barplot", 
        func = cor_barplot.PLOTFUN,
        func2 = cor_barplot.PLOTFUN,        
        info.text = cor_barplot.info,
        options = cor_barplot.opts,
        title = "Top correlated genes", label = "a",
        caption2 = "Top correlated genes.",
        ##pdf.width = 14, pdf.height = 4, 
        height = c(0.45*fullH,500),
        width = c('auto',800),
        res = c(65,80)
    )

    ##-----------------------------------------------------------
    ## Correlation network
    ##-----------------------------------------------------------
    
    cor_graph.PLOTFUN %<a-% reactive({

        req(input$cor_gene)
        res <- getGenePartialCorrelation()
        gene="XIST";rho.min=0.3;layout="kk"
        gene <- input$cor_gene
        rho.min <- input$cor_graph_threshold
        layout <- input$cor_graph_layout
        ##fixed <- input$cor_graph_fixed

        par(mfrow=c(1,1))
        gr <- pgx.plotPartialCorrelationAroundGene(
            res, gene, what="graph", ## degree=deg,
            rho.min = rho.min, nsize = 20, layout = layout, 
            edge.width = 3, main= "correlation graph")
    })

    GRAPH.LAYOUTS = c("Fruchterman-Reingold"="fr", "Kamada-Kawai"="kk",
                      "graphopt"="graphopt","tree layout"="tree")

    cor_graph.opts <- tagList(
        sliderInput(ns('cor_graph_threshold'),'rho:', 0, 1, 0.90),
        selectInput(ns('cor_graph_layout'),'layout:', choices=GRAPH.LAYOUTS)
        )

    cor_graph_info <- "</b>Correlation network.</b> Correlation graph centered on selected gene with top most correlated features. Red edges correspond to negative (marginal) correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation of the gene pair."
    
    callModule(
        plotModule,
        id = "cor_graph", label = "b",
        func = cor_graph.PLOTFUN,
        func2 = cor_graph.PLOTFUN,        
        info.text = cor_graph_info,
        options = cor_graph.opts,
        title = "Correlation network",
        ## caption = topEnriched_caption
        caption2 = cor_graph_info,
        ##pdf.width = 14, pdf.height = 4, 
        height = c(0.45*fullH,720),
        width = c('auto',1000),
        res = c(72,80)
    )

    ##--------------------------------------------------------------------------------
    ## Correlation scatter plots
    ##--------------------------------------------------------------------------------
    
    cor_scatter.PLOTFUN %<a-% reactive({

        req(input$cor_gene)
        res <- getGenePartialCorrelation()
        this.gene <- input$cor_gene

        NTOP = 25
        j <- which(rownames(res$cor) == this.gene)
        rho  <- head(res$cor[j,-j],NTOP)
        
        ngs <- inputData()
        par(mfrow=c(5,5), mar=c(4,3.5,0.5,0.2),
            mgp=c(2.1,0.8,0), oma=c(0,0,1,0))
        for(i in 1:min(25,length(rho))) {
            gene2 <- names(rho)[i]
            x <- ngs$X[gene2,]
            y <- ngs$X[this.gene,]
            plot(x, y, pch=20, cex=0.95, xlab=gene2, ylab=this.gene)
            abline(lm(y ~ x), col="red")
            
        }
    })

    cor_scatter.opts <- tagList()

    cor_scatter.info = "<b>Correlation scatter plots.</b> The plots show the pairwise scatter plots for the expression values of the gene pairs across the samples. The red line correspond to the (linear) regression fit."
    
    callModule(
        plotModule,
        id = "cor_scatter", label = "c",
        func = cor_scatter.PLOTFUN,
        func2 = cor_scatter.PLOTFUN,        
        info.text = cor_scatter.info,
        options = cor_scatter.opts,
        title = "Correlation scatter plots",
        ##caption = topEnriched_caption
        ##pdf.width = 14, pdf.height = 4, 
        height = c(fullH-80,720),
        width = c('auto',900),
        res = c(80,85)
    )


    ##--------------------------------------------------------------------------------
    ## Correlation GSEA
    ##--------------------------------------------------------------------------------
    getCorrelationGSEA <- reactive({
        ngs <- inputData()
        gene = "CD4"
        gene <- input$cor_gene
        gx <- ngs$X[gene,]
        rho <- cor( t(ngs$X), gx, use="pairwise")[,1]         

        gmt <- GSETS[colnames(ngs$GMT)]
        ## gmt <- GSETS  ## all???
        gsea <- fgsea(gmt, rho, nperm=1000, minSize=15, maxSize=1000)
        gsea <- gsea[order(gsea$pval),]

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
        NTOP = 16
        par(oma=c(0,1,0,0))
        par(mfrow=c(4,4), mar=c(2,1.5,4,1))
        par(mfrow=c(4,4), mar=c(0.5,1.5,2.8,1))
        i=1
        for(i in 1:min(NTOP,nrow(gsea))) {
            gs <- gsea$pathway[i]
            gs
            gmt <- GSETS[[gs]]
            length(gmt)
            ##if(length(gmtdx) < 3) { frame(); next }
            gsea.enplot( res$rho, gmt, main=gs, cex.main=0.9,
                        xlab="" )
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

    ##corGSEA_plots_module <- plotModule(
    callModule(
        plotModule, 
        id = "corGSEA_plots", ##ns=ns,
        func = corGSEA_plots.RENDER,
        func2 = corGSEA_plots.RENDER, 
        download.fmt = c("png","pdf"),
        options = corGSEA_plots_opts,
        ##info.text = corGSEA_plots_text,        
        title="Correlation GSEA", label="a",
        height = c(0.5*fullH,650), width = c('auto',1200),
        pdf.width=8, pdf.height=5, res=c(72,85)
    )
    ## output <- attachModule(output, corGSEA_plots_module)
    
    corGSEA_table.RENDER <- reactive({
        
        res = getCorrelationGSEA()
        
        ##rho = data.frame(cbind( name=rho.name, rho))
        gs <- res$gsea$pathway
        link <- wrapHyperLink(rep("link",length(gs)), gs)
        df = data.frame( pathway=res$gsea$pathway, link=link,
                        res$gsea[,c("pval","padj","NES","size")] )
        rownames(df) = gs
        numeric.cols = c("pval","padj","NES")
        
        DT::datatable(
                df, rownames=FALSE, escape = c(-1,-2),
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
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    corGSEA_table_info = "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings."

    corGSEA_table <- callModule(
        tableModule, 
        id = "corGSEA_table", 
        func = corGSEA_table.RENDER,
        info.text = corGSEA_table_info,
        title = "Correlation GSEA table", label="b",
        height = c(220,700), width=c('auto',1000)
        ##caption = corGSEA_caption
    )
    
    corGSEA_LeadingEdgeTable.RENDER <- reactive({
        
        res = getCorrelationGSEA()

        sel=1
        sel <- corGSEA_table$rows_selected()
        req(sel)
        
        ##rho = data.frame(cbind( name=rho.name, rho))
        le.genes <- res$gsea[sel,]$leadingEdge[[1]]
        rho1 <- res$rho[le.genes]
        title <- substring(GENE.TITLE[le.genes],1,40)
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
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') 
    })

    corGSEA_LeadingEdgeTable_info = "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings."

    corGSEA_LeadingEdgeTable <- callModule(
        tableModule, 
        id = "corGSEA_LeadingEdgeTable", 
        func = corGSEA_LeadingEdgeTable.RENDER,
        info.text = corGSEA_LeadingEdgeTable_info,
        title = "Leading edge genes", label="c",
        ##height = c(230,700), width=c('auto',1000)
        height = c(655,700), width=c('auto',1000)
        ##caption = corGSEA_caption
    )

    ##--------------------------------------------------------------------------------
    ## WGCNA
    ##--------------------------------------------------------------------------------

    getCalculateWGCNA <- reactive({
        ngs <- inputData()
        X <- ngs$X
        pheno <- ngs$samples
        res <- calculateWGCNA(X, pheno)
        return(res)
    })
    
    calculateWGCNA <- function(X, pheno) {
        require(WGCNA)
        
        ## Re-cluster samples
        X1 <- head(X[order(-apply(X,1,sd,na.rm=TRUE)),],1000)
        sampleTree2 = hclust(dist(X1), method = "average")
        ## Convert traits to a color representation: white means low, red means high, grey means missing entry

        
        traitColors = numbers2colors(pheno, signed = FALSE)
        ## Plot the sample dendrogram and the colors underneath.
        plotDendroAndColors(sampleTree2, traitColors,
                            groupLabels = names(datTraits),
                            main = "Sample dendrogram and trait heatmap")

    }

    


} ## end of module
