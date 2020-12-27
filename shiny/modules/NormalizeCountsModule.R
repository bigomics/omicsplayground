##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


##=====================================================================================
##===================== NORMALIZE GADGET UI ===========================================
##=====================================================================================

if(0) {
    load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")    

    NormalizeCountsGadget(X=ngs$X, pheno=ngs$samples)

    out <- gadgetize2(
        NormalizeCountsUI, NormalizeCountsServer,
        title = "UploadGadget", height=640, size="l", 
        X = ngs$X, pheno=ngs$samples )
    names(out)
    
}

NormalizeCountsGadget <- function(X, pheno, height=720) {
    gadgetize(NormalizeCountsUI, NormalizeCountsServer,
              title="NormalizeCounts",
              X=X, pheno=pheno, height=height)
}

NormalizeCountsUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("UI"))
}

NormalizeCountsServerRT <- function(id, counts, height=720) {
    ##
    ## counts:  reactive input object
    ##
    moduleServer(
        id,
        function(input, output, session) {
            
            all.methods <- c("none","scale","quantile",
                             "CPM","TMM","RLE","upperquartile")    
            
            output$UI <- renderUI({
                ##ns <- NS(id)  ## namespace
                ##ns <- function(id) id
                ns <- session$ns
                sidebarLayout(
                    sidebarPanel(
                        ##helpText(bc_info), 
                        ##radioButtons(ns("selectmethod"),"Select normalization:",
                        ##choices=all.methods, selected="CPM"),
                        ##br(),
                        tipify2(
                            checkboxInput(ns("postqn"),"Post quantile normalization"),
                            "Apply additional quantile normalization after scaling method."
                        ),
                        tipify2(
                            checkboxInput(ns("addnoise"),"Simulate unnormalized"),
                            "Simulated unnormalized data by adding random scaling to raw data"
                        ),
                        width = 2
                    ),
                    mainPanel(
                        plotOutput(ns("canvas"), width="100%", height=height) %>% withSpinner(),
                        width = 10
                    )
                )
            })

            pgx <- reactive({
                if(is.null(input$addnoise)) return(NULL)
                pgx <- list(counts=counts())                
                ## Simulate raw signals
                if(input$addnoise) {
                    mult <- 1e3 * runif(ncol(pgx$counts))
                    pgx$counts <- t( t(pgx$counts) * mult )
                }
                pgx
            })
            
            output$canvas <- renderPlot({
                req(counts())
                ## Show all methods
                postqn <- input$postqn                
                viz.NormalizeCounts(
                    pgx(),
                    methods = all.methods,
                    post.qn = postqn,
                    title=NULL, subtitle=NULL, caption=NULL)
            })
            
            normalized_counts <- reactive({
                ##req(input$selectmethod)
                ##method <- input$selectmethod
                method = "CPM"
                norm.counts <- NULL
                pgx <- pgx()
                if(method %in% all.methods) {
                    norm.counts <- pgx.countNormalization(pgx$counts, method)
                } else {
                    norm.counts <- pgx$counts  ## no-normalization
                }
                norm.counts
            })

            return(normalized_counts)  ## pointing to reactive
        } ## end-of-server
    )
}

