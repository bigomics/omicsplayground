CorrelationInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

CorrelationUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        height = 750,
        tabsetPanel(
            id = ns("tabs1"),
            tabPanel("Correlation",uiOutput(ns("corAnalysis_UI")))
        )
    )
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

        res <- pgx.computePartialCorrelationAroundGene(
            X, gene, method=methods, nmax=100, fast=FALSE)    

        rho1 <- min(abs(head(res$cor[gene,-1],18)))
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

    callModule(
        plotModule,
        id = "cor_barplot", 
        func = cor_barplot.PLOTFUN,
        func2 = cor_barplot.PLOTFUN,        
        ##info.text = topEnriched_text,
        options = cor_barplot.opts,
        title = "Top correlation", label = "a",
        caption2 = "Top correlated features.",
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

    callModule(
        plotModule,
        id = "cor_graph", label = "b",
        func = cor_graph.PLOTFUN,
        func2 = cor_graph.PLOTFUN,        
        ##info.text = topEnriched_text,
        options = cor_graph.opts,
        title = "Correlation network",
        ##caption = topEnriched_caption
        caption2 = "</b>Correlation network.</b> Correlation graph centered on selected gene with top most correlated features. Red edges correspond to negative (marginal) correlation, grey edges to positive correlation. Width of the edges is proportional to the absolute partial correlation of the gene pair.",
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
        gene <- input$cor_gene

        NTOP = 25
        rho  <- head(res$cor[gene,-1],NTOP)
        
        ngs <- inputData()
        par(mfrow=c(5,5), mar=c(4,3.5,0.5,0.2),
            mgp=c(2.1,0.8,0), oma=c(0,0,1,0))
        for(i in 1:min(25,length(rho))) {
            gene2 <- names(rho)[i]
            x <- ngs$X[gene,]
            y <- ngs$X[gene2,]
            plot(x, y, pch=20, cex=0.95, xlab=gene, ylab=gene2)
            abline(lm(y ~ x), col="red")
            
        }
    })

    cor_scatter.opts <- tagList()

    callModule(
        plotModule,
        id = "cor_scatter", label = "c",
        func = cor_scatter.PLOTFUN,
        func2 = cor_scatter.PLOTFUN,        
        ##info.text = topEnriched_text,
        options = cor_scatter.opts,
        title = "Correlation scatter plots",
        ##caption = topEnriched_caption
        ##pdf.width = 14, pdf.height = 4, 
        height = c(fullH-80,720),
        width = c('auto',900),
        res = c(80,85)
    )

    ##--------------------------------------------------------------------------------
    ## WGCNA
    ##--------------------------------------------------------------------------------

    require(WGCNA)
    HTML("To be implemented...")



} ## end of module
