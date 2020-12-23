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
                
                bio.correct <- c()
                if("mito/ribo" %in% bc)  bio.correct <- c(bio.correct, c("mito","ribo"))
                if("cell.cycle" %in% bc) bio.correct <- c(bio.correct, c("cc.phase","cc.score"))
                if("gender" %in% bc)     bio.correct <- c(bio.correct, c("gender"))

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
                
                cat("[SuperBatchCorrectServer] libsizes=",colSums(2**X0),"\n")
                cat("[SuperBatchCorrectServer] dim(pheno)=",dim(pheno()),"\n")
                cat("[SuperBatchCorrectServer] dim(X0)=",dim(X0),"\n")                

                cat("Performing SuperBatchCorrection\n")
                out <- pgx.superBatchCorrect(
                    X = X0,
                    pheno = pheno(),
                    model.par = mp,
                    batch.par = bp,
                    ## batch.cov = bv,
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
                cat("[SuperBatchCorrectServer] max(X) = ",max(X()),"\n")
                cat("[SuperBatchCorrectServer] max(cX) = ",max(out$X),"\n")
                cat("[SuperBatchCorrectServer] dim(X) = ",dim(X()),"\n")
                cat("[SuperBatchCorrectServer] dim(cX) = ",dim(out$X),"\n")

                X0 <- X()
                if(is.count) {
                    X0 <- log2(1 + X0)  ## X0: normalized counts (e.g. CPM)
                }

                viz.BatchCorrectionMatrix(
                    X0=X0, pheno=out$Y, cX=out$X, phenotype=p1,
                    pca.heatmap=do.pca, nmax=40, cex=1.8, npca=3, 
                    title=NULL, subtitle=NULL, caption=NULL)
            })
            
            output$inputsUI <- renderUI({

                ns <- session$ns
                bc.options = c("mito/ribo","cell.cycle","gender","PCA","HC","SVA","NNM")
                bc.selected = ""
                ##bc.selected = c("mito/ribo","cell.cycle","gender")
                
                bc_info <- NULL
                bc_info <- "Batch correction. Specify your parameters of interest ('model parameters') and batch parameters you want to correct for.\n\n"
                bc_info <- "Batch correction can clean your data from 'unwanted variables'. Please specify your parameters of interest.\n"
                
                ##pheno.par <- colnames(ngs$samples)
                pheno.par <- colnames(pheno())
                sel.par   <- c(grep("^[.]",pheno.par,invert=TRUE,value=TRUE),pheno.par)[1]

                tagList(
                    ##helpText(bc_info),
                    p(bc_info),
                    ## shinyWidgets::prettySwitch(ns("on_off"),"on/off", inputId="s1"),
                    selectInput(ns("bc_modelpar"), "Model parameters:", pheno.par,
                                selected=sel.par, multiple=TRUE),
                    actionButton(ns("bc_compute"),"Apply correction",
                                 icon=icon("exclamation-triangle"), class="run-button"),
                    br(), br(),
                    tipify( actionLink(ns("bc_options"), "Advanced",
                                       icon=icon("cog", lib = "glyphicon")),
                           "Toggle advanced options.", placement="top"),
                    br(),
                    conditionalPanel(
                        "input.bc_options%2 == 1", ns=ns,
                        tagList(
                            selectInput(ns("bc_batchpar"), "Batch parameters:", c("*",pheno.par),
                                        selected="*", multiple=TRUE),                       
                            tipify( checkboxGroupInput(
                                ns('bc_methods'),'Correction methods:',
                                choices=bc.options, bc.selected, inline=FALSE),
                                "Advanced options", placement="right",
                                options = list(container = "body")),
                            radioButtons(ns("bc_maptype"), "Heatmap type:",
                                         c("topSD","PCA"), inline=TRUE)
                        )
                    )
                )

            })
                

            return(outobj)  ## pointing to reactive
        } ## end-of-server
    )
}

