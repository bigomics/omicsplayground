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
      

      if(0) {
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
      }


      max.rho=0.5
      getNotCorrelatedBatchPars <- function(pheno, model.par, max.rho=0.5) {
        
        getModelMatrix <- function(v) {
          y <- as.character(pheno[,v])
          y[is.na(y)] <- "NA"  ## or impute???
          m1 <- model.matrix( ~ y)[,-1,drop=FALSE]
          colnames(m1) <- sub("^y",paste0(v,"="),colnames(m1))
          m1
        }

        model.par <- intersect(model.par, colnames(pheno))
        mod1 <- do.call(cbind,lapply(model.par, getModelMatrix))
        rownames(mod1) <- rownames(pheno)

        batch.prm <- setdiff(colnames(pheno),model.par)
        mod0 <- do.call(cbind,lapply(batch.prm, getModelMatrix))        
        rho <- stats::cor(mod0,mod1)
        rho
        rho[is.na(rho)] <- 0
        max(abs(rho),na.rm=TRUE)
        if(max(abs(rho),na.rm=TRUE) > max.rho) {
          idx <- which(abs(rho) > max.rho, arr.ind=TRUE)
          idx
          for(i in 1:nrow(idx)) {
            v0 <- colnames(mod0)[idx[i,1]]
            v1 <- colnames(mod1)[idx[i,2]]
            cat(paste0("WARNING:: '",v0,"' is confounded with '",v1,"' ",
                       ": rho= ",round(rho[idx[i,1],idx[i,2]],3),"\n"))
          }
          confounding.pars <- colnames(mod0)[idx[,1]]
          confounding.pars <- unique(gsub("=.*","",confounding.pars))
          cat("WARNING:: removing confounding batch factors:",confounding.pars,"\n")
          batch.prm <- setdiff(batch.prm, confounding.pars)
        }
        batch.prm        
      }
             
      ## reset batch parameters choices if modelparamters change.
      observeEvent( input$bc_modelpar, {

        ##pheno.par <- colnames(ngs$samples)
        px <- pheno()
        if(is.null(px)) return(NULL)
        sel <- apply(px,2,function(x) length(unique(x[!is.na(x)])))
        pheno.par <- sort(colnames(px)[sel>1])
        sel.par   <- c(grep("^[.<]|batch",pheno.par,invert=TRUE,value=TRUE),pheno.par)[1]
        ##batch.par <- c("*",pheno.par,"<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
        batch.par <- c(pheno.par,"<cell_cycle>","<gender>","<libsize>","<mito/ribo>")        

        bc <- setdiff(batch.par, input$bc_modelpar)
        shiny::updateSelectInput(session,"bc_batchpar",choices=bc, selected=bc)        

      })
      
      uiOK <- reactive({
        px <- pheno()
        mp <- input$bc_modelpar
        ok <- all(mp %in% colnames(px)) && length(mp)>0
        ok
      })
      
      output$canvas <- shiny::renderPlot({

        dbg("[BatchCorrect::canvas] reacted!")
        
        out <- outobj()
        
        dbg("[BatchCorrect::canvas] 0: dim.outX = ", dim(out$X))
        dbg("[BatchCorrect::canvas] 0: names.X = ", names(out))        

        if(is.null(out$X) || length(out$X)==0) {
          dbg("[BatchCorrect::canvas] WARNING :: X=NULL")
          return(NULL)
        }

        dbg("[BatchCorrect::canvas] 0: dim.X = ", dim(X()))
        
        shiny::req(X())

        dbg("[BatchCorrect::canvas] 1: ")
        
        
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

        if(is.null(mp) || length(mp)==0) {
          dbg("[BatchCorrect::canvas] WARNING :: mp = ",mp)
        }        

        p1 <- NULL
        if(!is.null(mp)) p1 <- mp[1]

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

        if(ncol(X0) != ncol(cX)) {
            message("[BatchCorrect] canvas.renderPlot : ***ERROR*** dimension mismatch!")
            ## return(NULL)
        }
        req(ncol(X0) == ncol(cX))

        do.pca <- (input$bc_maptype == "PCA")                    
        show_row <- (nmax < 50)
        
        viz.BatchCorrectionMatrix(
          X0=X0, pheno=out$Y, cX=cX, phenotype=p1,
          pca.heatmap=do.pca, npca=3,
          nmax=nmax, cex=1.8, show_rownames=show_row,
          title=NULL, subtitle=NULL, caption=NULL)

      })
      
      output$inputsUI <- shiny::renderUI({
        
        ns <- session$ns
        bc.options = c("none","PCA","HC","SVA","NNM")
        bc.options = c("none","PCA","SVA","NNM")        
        bc.selected = "none"
        ##bc.selected = c("mito/ribo","cell.cycle","gender")
        
        bc_info <- NULL
        bc_info <- "Batch correction can clean your data from 'unwanted variables'. Please specify your parameters of interest.\n"
        
        ##pheno.par <- colnames(ngs$samples)
        px <- pheno()
        if(is.null(px)) return(NULL)
        sel <- apply(px,2,function(x) length(unique(x[!is.na(x)])))
        pheno.par <- sort(colnames(px)[sel>1])
        sel.par   <- c(grep("^[.<]|batch",pheno.par,invert=TRUE,value=TRUE),pheno.par)[1]
        ##batch.par <- c("*",pheno.par,"<cell_cycle>","<gender>","<libsize>","<mito/ribo>")
        batch.par <- c(pheno.par,"<cell_cycle>","<gender>","<libsize>","<mito/ribo>")        
        
        shiny::tagList(
          ##helpText(bc_info),
          shiny::p(bc_info),
          ## shinyWidgets::prettySwitch(ns("on_off"),"on/off", inputId="s1"),

          shiny::br(),          
          shiny::actionButton(ns("bc_compute_button"),"Batch correct",
                              ## icon=icon("exclamation-triangle"),
                              class="run-button"),
          shiny::br(),
          shiny::br(),

          tipify2(
            shiny::selectInput(ns("bc_modelpar"), "Model parameters:", pheno.par,
                               selected=sel.par, multiple=TRUE),
            "Please specify <b>all</b> your model parameters. These are the parameters of interest that will determine your groupings."),
          
          ## withTooltip(                    
          ##   shiny::radioButtons(ns("bc_strength"), NULL,
          ##                       c("low","medium","strong"), inline=TRUE),
          ##   "Choose the strength of batch correction: <b>low</b> corrects only for explicit batch parameters, <b>medium</b> corrects for additional unwanted biological effects (inferred from your data), <b>strong</b> applies SVA or NNM (nearest neighbour matching). You can tune the selection under the advanced options.", 
          ##   placement="left", options=list(container="body")),                    

##          withTooltip( shiny::actionLink(ns("bc_options"), "Advanced",
##                                             icon=icon("cog", lib = "glyphicon")),
##                          "Toggle advanced options.", 
##                          placement="left", options=list(container="body")),
##          shiny::br(),

##          shiny::conditionalPanel(
##            "input.bc_options%2 == 1", ns=ns,
##            shiny::tagList(

##          shiny::conditionalPanel(
##            "input.bc_methods == 'limma'", ns=ns,
            withTooltip(
              shiny::selectInput(ns("bc_batchpar"), "Batch parameters:", batch.par,
                                 selected=batch.par, multiple=TRUE),
              "Specifiy the batch/unwanted parameters that you want to correct for. Bracketed parameters are technical/biological summaries computed from your data. These factors will be subtracted from your data using linear modelling (limma).",
              placement="left", options=list(container="body")
##            )
              ),

          
          withTooltip(
            ##shiny::checkboxGroupInput(
            shiny::radioButtons(                
              ns('bc_methods'),'Unsupervised correction:',
              choices=bc.options, selected=bc.selected, inline=FALSE),
            "Unsupervised correction methods. Correction will be performed additional to the (supervised) corrections above. <b>PCA</b> iteratively corrects low rank PC components not correlated to any model parameters; <b>SVA</b> applies surrogate variable analysis (Leek et al.); <b>NNM</b> applies nearest neighbour matching, a quasi-pairing approach for incomplete matched data (unpublished).",
            placement="left", options=list(container="body")
          ),
                    
          withTooltip( shiny::radioButtons(ns("bc_nmax"), "Nmax:",c(40,200,1000),
                                               sel=200, inline=TRUE),
                          "Maximum number of genes in heatmap",
                          placement="left", options=list(container = "body")),
          
          withTooltip( shiny::radioButtons(ns("bc_maptype"), "Heatmap type:",
                                               c("topSD","PCA"), inline=TRUE),
                              "Type of heatmap: top SD or PCA.",
                          placement="left", options=list(container="body"))
          
        )
##      )
##        )                
      })                

      ## ------------------------------------------------------------
      ## Reactive return object
      ## ------------------------------------------------------------
      outobj <- shiny::eventReactive( {
        list(input$bc_compute_button, X(), pheno())
      } , {
        
        dbg("[BatchCorrectionModule::outobj] reacted...")
        dbg("[BatchCorrectionModule::outobj] is.null(X) = ",is.null(X()))
        
        shiny::req(X())

        dbg("[BatchCorrectionModule::outobj] 0 : uiOK = ",uiOK())
        req(uiOK())
        
        mp="";bp="Chemotherapy"
        mp="dlbcl.type";bp="*"
        mp <- input$bc_modelpar
        bp <- input$bc_batchpar
        bc <- input$bc_methods

        dbg("[BatchCorrectionModule::outobj] 1 : mp = ",mp)
        dbg("[BatchCorrectionModule::outobj] 1 : bp = ",bp)
        dbg("[BatchCorrectionModule::outobj] 1 : bc = ",bc)
        dbg("[BatchCorrectionModule::outobj] 1 : colnames.pheno = ",colnames(pheno()))

        if(0) {
          mp <- intersect(mp, colnames(pheno()))  ## check        
          ##shiny::req(mp)
          validate( need(length(mp)>0, "need valid model parameters"))
          if(is.null(mp) || length(mp)==0) {
            dbg("[event:outobj] ***WARNING*** no model parameter mp = ",mp)
            return(NULL)
          }
        }
        
        lib.correct <- FALSE
        bio.correct <- c()
        if("<libsize>" %in% bp)    lib.correct <- TRUE
        if("<mito/ribo>" %in% bp)  bio.correct <- c(bio.correct, c("mito","ribo"))
        if("<cell_cycle>" %in% bp) bio.correct <- c(bio.correct, c("cell_cycle"))
        if("<gender>" %in% bp)     bio.correct <- c(bio.correct, c("gender"))
        
        hc.correct  <- "HC"  %in% bc
        pca.correct <- "PCA" %in% bc
        sva.correct <- "SVA" %in% bc
        ##mnn.correct <- "MNN" %in% bc
        mnn.correct  <- NULL
        nnm.correct <- "NNM" %in% bc

        if(0) {
            ## disable parameter correction if doing fancy stuff??
            if(pca.correct || sva.correct || hc.correct || nnm.correct ) {
                bio.correct <- c()
                bp <- c()
            }
        }
        
        X0 <- X()
        if(is.count) {
          X0 <- log2(1 + X0)  ## X0: normalized counts (e.g. CPM)
        }

        dbg("[event:outobj] dim(X0) = ",dim(X0))
        dbg("[event:outobj] dim(pheno) = ",dim(pheno()))
        dbg("[event:outobj] mp = ",mp)
        dbg("[event:outobj] bp = ",bp)              
        
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

        message("[event:outobj] done!")
        
        out
      }, ignoreInit=FALSE, ignoreNULL=FALSE )

      return(outobj)  ## pointing to reactive
    } ## end-of-server
  )
}

