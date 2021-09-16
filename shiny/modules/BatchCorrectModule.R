##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##


##=============================================================================
##==================== BATCHCORRECT GADGET UI =================================
##=============================================================================

if(0) {
    load("~/Playground/omicsplayground/data/GSE10846-dlbcl-nc.pgx")    

    BatchCorrectGadget(X=ngs$X, pheno=ngs$samples)

    out <- gadgetize2(
        BatchCorrectUI, BatchCorrectServer,
        title = "UploadGadget", height=640, size="l", 
        X = ngs$X, pheno=ngs$samples )
    names(out)
    
}

BatchCorrectGadget <- function(X, pheno, height=720) {
    gadgetize(BatchCorrectUI, BatchCorrectServer,
              title="BatchCorrect",
              X=reactive(X), pheno=reactive(pheno), height=height)
}

BatchCorrectUI <- function(id, height=720) {
    ns <- shiny::NS(id)
    shiny::sidebarLayout(
        shiny::sidebarPanel(
            shiny::uiOutput(ns("inputsUI")),
            width = 2
        ),
        shiny::mainPanel(
            shiny::plotOutput(ns("canvas"), width="100%", height=height) %>% shinycssloaders::withSpinner(),
            width = 10
        )
    )
}

BatchCorrectInputsUI <- function(id) {
    ns <- shiny::NS(id)
    shiny::uiOutput(ns("inputsUI"))
}

BatchCorrectCanvas <- function(id, height=720) {
    ns <- shiny::NS(id)
    shiny::plotOutput(ns("canvas"), width="100%", height=height) %>% shinycssloaders::withSpinner()
}

BatchCorrectServer <- function(id, X, pheno, is.count=FALSE, height=720) {
    shiny::moduleServer(
        id,
        function(input, output, session) {

            shiny::observeEvent( input$bc_strength, {
                if(input$bc_strength=="low") {
                    shiny::updateSelectInput(session,"bc_batchpar",selected="*")
                    shiny::updateCheckboxGroupInput(session, "bc_methods", selected="")
                }
                if(input$bc_strength=="medium") {
                    sel <- c("*","<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
                    shiny::updateSelectInput(session,"bc_batchpar",selected=sel)
                    shiny::updateCheckboxGroupInput(session,"bc_methods", selected=c("PCA","HC"))
                }
                if(input$bc_strength=="strong") {
                    sel <- c("*","<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
                    shiny::updateSelectInput(session,"bc_batchpar",selected=sel)
                    shiny::updateCheckboxGroupInput(session, "bc_methods", selected="SVA")
                }

            })
            
            outobj <- shiny::eventReactive(input$bc_compute, {
                
                message("[event:bc_compute] reacted...")
                shiny::req(X())
                ## input$bc_compute
                message("[event:bc_compute] shiny::req(X) OK..")
                
                mp="";bp="Chemotherapy"
                mp="dlbcl.type";bp="*"
                mp <- input$bc_modelpar
                bp <- input$bc_batchpar
                bc <- input$bc_methods
                ##bv <- NULL
                message("[event:bc_compute] 0")
                message("[event:bc_compute] 0 : mp = ", mp)
                message("[event:bc_compute] 0 : bp = ", bp)
                message("[event:bc_compute] 0 : bc = ", bc)
                
                if(is.null(mp) && is.null(bp)) return(NULL)
                message("[event:bc_compute] 1")

                lib.correct <- FALSE
                bio.correct <- c()
                if("<libsize>" %in% bp)  lib.correct <- TRUE
                if("<mito/ribo>" %in% bp)  bio.correct <- c(bio.correct, c("mito","ribo"))
                if("<cell_cycle>" %in% bp) bio.correct <- c(bio.correct, c("cell_cycle"))
                if("<gender>" %in% bp)     bio.correct <- c(bio.correct, c("gender"))

                message("[event:bc_compute] 2 : bc = ", bc)
                
                hc.correct  <- "HC"  %in% bc
                pca.correct <- "PCA" %in% bc
                sva.correct <- "SVA" %in% bc
                ##mnn.correct <- "MNN" %in% bc
                mnn.correct  <- NULL
                nnm.correct <- "NNM" %in% bc
                
                message("[event:bc_compute] 3 : mp = ",mp)
                
                ## if(length(mp) == 0) mp <- NULL
                ## if(mp=="") mp <- NULL
                
                message("[event:bc_compute] 4" )
                
                X0 <- X()
                message("[event:bc_compute] 4a : dim(X0) =",paste(dim(X0),collapse="x") )
                if(is.count) {
                    X0 <- log2(1 + X0)  ## X0: normalized counts (e.g. CPM)
                }
                
                message("[event:bc_compute] 4b : dim(X0) =",paste(dim(X0),collapse="x") )
                
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

                message("[event:bc_compute] done!")
                
                out
            }, ignoreInit=FALSE, ignoreNULL=FALSE )
            
            output$canvas <- shiny::renderPlot({

                out <- outobj()
                if(is.null(out$X) || length(out$X)==0) return(NULL)

                ## is.valid <- all(dim(out$X)==dim(X()) &&
                ##                 nrow(out$Y)==ncol(X()))
                ## if(!is.valid) {
                ##     message("[BatchCorrectServer] dim(out$X)=",dim(out$X))
                ##     message("[BatchCorrectServer] dim(out$Y)=",dim(out$Y))
                ##     message("[BatchCorrectServer] dim(X())=",dim(X()))                    
                ##     message("[BatchCorrectServer] WARNING: matrices not valid!!")
                ##     return(NULL)
                ## }
                
                mp <- shiny::isolate(input$bc_modelpar)
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

                cX <- out$X
                rownames(X0) <- sub("[;|,].*","",rownames(X0))
                rownames(cX) <- sub("[;|,].*","",rownames(cX))

                message("[BatchCorrect] canvas.renderPlot : dim(X0) =",paste(dim(X0),collapse="x") )
                message("[BatchCorrect] canvas.renderPlot : dim(cX) =",paste(dim(cX),collapse="x") )
                
                viz.BatchCorrectionMatrix(
                    X0=X0, pheno=out$Y, cX=cX, phenotype=p1,
                    pca.heatmap=do.pca, nmax=nmax, cex=1.8, npca=3, 
                    title=NULL, subtitle=NULL, caption=NULL)

            })
            
            output$inputsUI <- shiny::renderUI({
                
                ns <- session$ns
                bc.options = c("PCA","HC","SVA","NNM")
                bc.selected = ""
                ##bc.selected = c("mito/ribo","cell.cycle","gender")
                
                bc_info <- NULL
                bc_info <- "Batch correction can clean your data from 'unwanted variables'. Please specify your parameters of interest.\n"
                
                ##pheno.par <- colnames(ngs$samples)
                px <- pheno()
                if(is.null(px)) return(NULL)
                sel <- apply(px,2,function(x) length(unique(x[!is.na(x)])))
                pheno.par <- sort(colnames(px)[sel>1])
                sel.par   <- c(grep("^[.<]",pheno.par,invert=TRUE,value=TRUE),pheno.par)[1]
                batch.par <- c("*",pheno.par,"<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
                
                shiny::tagList(
                    ##helpText(bc_info),
                    shiny::p(bc_info),
                    ## shinyWidgets::prettySwitch(ns("on_off"),"on/off", inputId="s1"),
                    tipify2(
                        shiny::selectInput(ns("bc_modelpar"), "Model parameters:", pheno.par,
                                    selected=sel.par, multiple=TRUE),
                        "Please specify your model parameters. These are the parameters of interest that will determine your groupings."),
                    shinyBS::tipify(                    
                        shiny::radioButtons(ns("bc_strength"), NULL,
                                     c("low","medium","strong"), inline=TRUE),
                        "Choose the strength of batch correction: <b>low</b> corrects only for explicit batch parameters, <b>medium</b> corrects for additional unwanted biological effects (inferred from your data), <b>strong</b> applies SVA or NNM (nearest neighbour matching). You can tune the selection under the advanced options.", 
                        placement="left", options=list(container="body")),                    
                    shiny::actionButton(ns("bc_compute"),"Batch correct",
                                 ## icon=icon("exclamation-triangle"),
                                 class="run-button"),
                    shiny::br(), shiny::br(),
                    shinyBS::tipify( shiny::actionLink(ns("bc_options"), "Advanced",
                                       icon=icon("cog", lib = "glyphicon")),
                           "Toggle advanced options.", 
                           placement="left", options=list(container="body")),
                    shiny::br(),
                    shiny::conditionalPanel(
                        "input.bc_options%2 == 1", ns=ns,
                        shiny::tagList(
                            shinyBS::tipify(
                                shiny::selectInput(ns("bc_batchpar"), "Batch parameters:", batch.par,
                                            selected="*", multiple=TRUE),
                                "Specifiy the batch parameters that you want to correct for. The <b>star</b> stands for all remaining (unwanted) parameters not specified as parameter of interest. Bracketed parameters are technical/biological summaries computed from your data.",
                                placement="left", options=list(container="body")
                            ),
                            shinyBS::tipify(
                                shiny::checkboxGroupInput(
                                    ns('bc_methods'),'Correction methods:',
                                    choices=bc.options, bc.selected, inline=FALSE),
                                "Advanced batch correction methods. <b>PCA</b> corrects PC components not correlated to any model parameters; <b>HC</b> iteratively corrects hierarchical clustering; <b>SVA</b> applies surrogate variable analysis (Leek et al.); <b>NNM</b> applies nearest neighbour matching, a quasi-pairing approach for incomplete matched data.",
                                placement="left", options=list(container="body")
                            ),
                            shinyBS::tipify( shiny::radioButtons(ns("bc_nmax"), "Nmax:",c(40,200,1000),
                                                 inline=TRUE),
                                   "Maximum number of genes in heatmap",
                                   placement="left", options=list(container = "body")),
                            shinyBS::tipify( shiny::radioButtons(ns("bc_maptype"), "Heatmap type:",
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

