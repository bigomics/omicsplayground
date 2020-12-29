##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


##=============================================================================
##==================== BATCHCORRECT GADGET UI =================================
##=============================================================================

if(0) {
    load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")    

    SuperBatchCorrectGadget(X=ngs$X, pheno=ngs$samples)

    out <- gadgetize2(
        SuperBatchCorrectUI, SuperBatchCorrectServer,
        title = "UploadGadget", height=640, size="l", 
        X = ngs$X, pheno=ngs$samples )
    names(out)
    
}

SuperBatchCorrectGadget <- function(X, pheno, height=720) {
    gadgetize(SuperBatchCorrectUI, SuperBatchCorrectServer,
              title="SuperBatchCorrect",
              X=reactive(X), pheno=reactive(pheno), height=height)
}

SuperBatchCorrectUI <- function(id, height=720) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            uiOutput(ns("inputsUI")),
            width = 2
        ),
        mainPanel(
            plotOutput(ns("canvas"), width="100%", height=height) %>% withSpinner(),
            width = 10
        )
    )
}

SuperBatchCorrectInputsUI <- function(id) {
    ns <- NS(id)
    uiOutput(ns("inputsUI"))
}

SuperBatchCorrectCanvas <- function(id, height=720) {
    ns <- NS(id)
    plotOutput(ns("canvas"), width="100%", height=height) %>% withSpinner()
}

SuperBatchCorrectServer <- function(id, X, pheno, is.count=FALSE, height=720) {
    moduleServer(
        id,
        function(input, output, session) {

            observeEvent( input$bc_strength, {
                if(input$bc_strength=="low") {
                    updateSelectInput(session,"bc_batchpar",selected="*")
                    updateCheckboxGroupInput(session, "bc_methods", selected="")
                }
                if(input$bc_strength=="medium") {
                    sel <- c("*","<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
                    updateSelectInput(session,"bc_batchpar",selected=sel)
                    updateCheckboxGroupInput(session,"bc_methods", selected=c("PCA","HC"))
                }
                if(input$bc_strength=="strong") {
                    sel <- c("*","<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
                    updateSelectInput(session,"bc_batchpar",selected=sel)
                    updateCheckboxGroupInput(session, "bc_methods", selected="SVA")
                }

            })
            
            outobj <- eventReactive(input$bc_compute, {
                
                req(input$bc_modelpar)
                ## input$bc_compute
                
                mp="";bp="Chemotherapy"
                mp="dlbcl.type";bp="*"
                mp <- isolate(input$bc_modelpar)
                bp <- isolate(input$bc_batchpar)
                bc <- isolate(input$bc_methods)
                ##bv <- NULL
                
                if(is.null(mp) && is.null(bp)) return(NULL)
                cat("outobj: 1\n")

                lib.correct <- FALSE
                bio.correct <- c()
                if("<libsize>" %in% bp)  lib.correct <- TRUE
                if("<mito/ribo>" %in% bp)  bio.correct <- c(bio.correct, c("mito","ribo"))
                if("<cell_cycle>" %in% bp) bio.correct <- c(bio.correct, c("cc.score"))
                if("<gender>" %in% bp)     bio.correct <- c(bio.correct, c("gender"))

                cat("outobj: 2\n")                
                hc.correct  <- "HC"  %in% bc
                pca.correct <- "PCA" %in% bc
                sva.correct <- "SVA" %in% bc
                ##mnn.correct <- "MNN" %in% bc
                mnn.correct  <- NULL
                nnm.correct <- "NNM" %in% bc
                if(mp[1]=="") mp <- NULL
                
                X0 <- X()
                cat("[SuperBatchCorrectServer] is.count=",is.count,"\n")
                if(is.count) {
                    X0 <- log2(1 + X0)  ## X0: normalized counts (e.g. CPM)
                }
                
                cat("Performing SuperBatchCorrection\n")
                out <- pgx.superBatchCorrect(
                    X = X0,
                    pheno = pheno(),
                    model.par = mp,
                    batch.par = bp,
                    ## batch.cov = bv,
                    lib.correct = lib.correct,                    
                    bio.correct = bio.correct,
                    hc.correct = hc.correct,
                    pca.correct = pca.correct,
                    sva.correct = sva.correct,
                    mnn.correct = mnn.correct,
                    nnm.correct = nnm.correct
                )

                out$params <- list(
                    model.par = mp,
                    batch.par = bp,
                    bio.correct = bio.correct,
                    hc.correct = hc.correct,
                    pca.correct = pca.correct,
                    sva.correct = sva.correct,
                    mnn.correct = mnn.correct,
                    nnm.correct = nnm.correct
                )
                
                out
            }, ignoreInit=FALSE, ignoreNULL=FALSE )
            
            output$canvas <- renderPlot({
                out <- outobj()
                if(is.null(out$X) || length(out$X)==0) return(NULL)

                is.valid <- all(dim(out$X)==dim(X()) &&
                                nrow(out$Y)==ncol(X()))
                if(!is.valid) {
                    cat("[SuperBatchCorrectServer] dim(out$X)=",dim(out$X),"\n")
                    cat("[SuperBatchCorrectServer] dim(out$Y)=",dim(out$Y),"\n")
                    cat("[SuperBatchCorrectServer] dim(X())=",dim(X()),"\n")                    
                    cat("[SuperBatchCorrectServer] WARNING: matrices not valid!!\n")
                    return(NULL)
                }
                
                mp <- isolate(input$bc_modelpar)
                p1 <- NULL
                if(!is.null(mp)) p1 <- mp[1]

                do.pca <- (input$bc_maptype == "PCA")                    
                ##viz.BatchCorrection( ngs, out$X, cX2=NULL, phenotype=p1,
                ##title=NULL, subtitle=NULL, caption=NULL)

                X0 <- X()
                if(is.count) {
                    X0 <- log2(1 + X0)  ## X0: normalized counts (e.g. CPM)
                }
                nmax <- as.integer(input$bc_nmax)
                
                viz.BatchCorrectionMatrix(
                    X0=X0, pheno=out$Y, cX=out$X, phenotype=p1,
                    pca.heatmap=do.pca, nmax=nmax, cex=1.8, npca=3, 
                    title=NULL, subtitle=NULL, caption=NULL)
            })
            
            output$inputsUI <- renderUI({

                ns <- session$ns
                bc.options = c("PCA","HC","SVA","NNM")
                bc.selected = ""
                ##bc.selected = c("mito/ribo","cell.cycle","gender")
                
                bc_info <- NULL
                bc_info <- "Batch correction can clean your data from 'unwanted variables'. Please specify your parameters of interest.\n"
                
                ##pheno.par <- colnames(ngs$samples)
                pheno.par <- sort(colnames(pheno()))
                sel.par   <- c(grep("^[.<]",pheno.par,invert=TRUE,value=TRUE),pheno.par)[1]
                batch.par <- c("*",pheno.par,"<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
                
                tagList(
                    ##helpText(bc_info),
                    p(bc_info),
                    ## shinyWidgets::prettySwitch(ns("on_off"),"on/off", inputId="s1"),
                    tipify2(
                        selectInput(ns("bc_modelpar"), "Model parameters:", pheno.par,
                                    selected=sel.par, multiple=TRUE),
                        "Please specify your model parameters. These are the parameters of interest that will determine your groupings.", placement="left"),
                    tipify(                    
                        radioButtons(ns("bc_strength"), NULL,
                                     c("low","medium","strong"), inline=TRUE),
                        "Choose the strength of batch correction: <b>low</b> corrects only for explicit batch parameters, <b>medium</b> corrects for additional unwanted biological effects (inferred from your data), <b>strong</b> applies SVA or NNM (nearest neighbour matching). You can tune the selection under the advanced options.", 
                        placement="left", options=list(container="body")),                    
                    actionButton(ns("bc_compute"),"Batch correct",
                                 ## icon=icon("exclamation-triangle"),
                                 class="run-button"),
                    br(), br(),
                    tipify( actionLink(ns("bc_options"), "Advanced",
                                       icon=icon("cog", lib = "glyphicon")),
                           "Toggle advanced options.", 
                           placement="left", options=list(container="body")),
                    br(),
                    conditionalPanel(
                        "input.bc_options%2 == 1", ns=ns,
                        tagList(
                            tipify(
                                selectInput(ns("bc_batchpar"), "Batch parameters:", batch.par,
                                            selected="*", multiple=TRUE),
                                "Specifiy the batch parameters that you want to correct for. The <b>star</b> stands for all remaining (unwanted) parameters not specified as parameter of interest. Bracketed parameters are technical/biological summaries computed from your data.",
                                placement="left", options=list(container="body")
                            ),
                            tipify(
                                checkboxGroupInput(
                                    ns('bc_methods'),'Correction methods:',
                                    choices=bc.options, bc.selected, inline=FALSE),
                                "Advanced batch correction methods. <b>PCA</b> corrects PC components not correlated to any model parameters; <b>HC</b> iteratively corrects hierarchical clustering; <b>SVA</b> applies surrogate variable analysis (Leek et al.); <b>NNM</b> applies nearest neighbour matching, a quasi-pairing approach for incomplete matched data.",
                                placement="left", options=list(container="body")
                            ),
                            tipify( radioButtons(ns("bc_nmax"), "Nmax:",c(40,200,1000),
                                                 inline=TRUE),
                                   "Maximum number of genes in heatmap",
                                   placement="left", options=list(container = "body")),
                            tipify( radioButtons(ns("bc_maptype"), "Heatmap type:",
                                                 c("topSD","PCA"), inline=TRUE),
                                   "Type of heatmap: top SD or PCA.",
                                   placement="left", options=list(container="body"))
                        )
                    )
                )                
            })
                

            return(outobj)  ## pointing to reactive
        } ## end-of-server
    )
}

