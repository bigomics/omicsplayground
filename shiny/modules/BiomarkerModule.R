##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

dbg(">>> sourcing BiomarkerModule \n")

BiomarkerInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

BiomarkerUI <- function(id) {
    ns <- NS(id)  ## namespace
    ui <- fillCol(
        flex = c(1),
        height = 760,
        tabsetPanel(
            id=ns("tabs"),
            tabPanel("Importance", uiOutput(ns("pdx_biomarker_UI")))
        )
    )
    ui
}

BiomarkerModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    inputData <- env[["load"]][["inputData"]]
    fullH = 750  ## full height of panel
    rowH  = 320  ## row height of panel
    imgH  = 260
        
    description = "<b>Biomarker Module.</b> Select biomarkers that can be used for
classification or prediction purposes. The phenotype of interest can
be multiple categories (classes) or patient survival data."
    output$description <- renderUI(HTML(description))

    pdx_infotext =
        "The <strong>Biomarker Module</strong> performs the biomarker selection that can be used for classification or prediction purposes.

<br><br>To better understand which genes, mutations, or gene sets influence the final phenotype the most, Playground calculates a variable importance score for each feature using state-of-the-art machine learning algorithms, including LASSO, elastic nets, random forests, and extreme gradient boosting, and provides the top 50 features according to cumulative ranking by the algorithms. By combining several methods, the platform aims to select the best possible biomarkers.

<br><br>The phenotype of interest can be multi-categorical classes or patient survival data. Instead of choosing a phenotype, users can also specify a particular contrast from the analysis and perform biomarker selection. The platform also provides a heatmap of samples based on identified top features. 

<br><br>In addition, it generates a classification tree using top features and provides expression boxplots by phenotype classes for features present in the tree. The platform can also provide a survival tree analysis using top features and provides expression boxplots by phenotype classes for features present in the tree."


    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("pdx_info"), "Info", icon = icon("info-circle")),
                   "Show more information about this module."),
            hr(), br(),             
            tipify(selectInput(ns("pdx_predicted"),"Predicted target:", choices=NULL),
                   "Select the target variable for biomarker selection.",placement="top"),
            ##tipify( selectInput(ns("pdx_level"),"Feature level:", choices=c("gene","geneset")),
            ##       "Select feature level: gene or geneset", placement="top"),
            tipify( selectInput(ns("pdx_filter"),"Feature filter:", choices=NULL),
                   "Select a filter for the features.", placement="top"),
            conditionalPanel(
                "input.pdx_filter == '<custom>'", ns=ns,
                tipify(
                    ##selectizeInput("pdx_select","Custom features:", choices=NULL, multiple=TRUE),
                    div(class='gene-list',
                        textAreaInput(ns("pdx_select"), "Custom features:", value = NULL,
                                      rows=8, placeholder="Paste your gene list")),
                    "Paste a custom gene list to be used as features.", placement="top")
            ),
            tipify(actionButton(ns("pdx_runbutton"), label="Compute", class="run-button"),
                   "Click to start biomarker computation.", placement="right")
        )
        if(DEV.VERSION) {
            uix <- tagList(
                hr(),br(),
                h6("Developer options:"),
                checkboxInput(ns('pdx_multiomics'),'multi-omics weighting',FALSE)
            )
            ui <- c(ui, uix)
        }
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ##======================= REACTIVE/OBSERVE FUNCTIONS =============================
    ##================================================================================
    
    observeEvent( input$pdx_info, {
        dbg("<module-biomarker::observe pdxinfo> reacted")
        showModal(modalDialog(
            title = HTML("<strong>Biomarker Module</strong>"),
            HTML(pdx_infotext),
            easyClose = TRUE, size="l"))
    })

    ##input_pdx_select <- reactive({
    ##    input$pdx_select
    ##}) %>% debounce(3000)

    input_pdx_select <- reactive({
        dbg("[BiomarkerModule:<input_pdx_select>]  reacted")
        gg <- input$pdx_select
        if(is.null(gg)) return(NULL)
        ##req(gg)
        gg <- strsplit(as.character(gg), split="[, \n\t]")[[1]]
        if(length(gg)==0) return(NULL)
        if(length(gg)==1 && gg[1]!="") gg <- c(gg,gg)  ## hack to allow single gene....
        return(gg)
    }) %>% debounce(1000)

    observe({
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        req(ngs)
        dbg("[BiomarkerModule::observe1] reacted")
        ct <- colnames(ngs$Y)
        ##ct <- grep("group|sample|patient|donor",ct,value=TRUE,invert=TRUE)
        ct <- grep("sample|patient|donor",ct,value=TRUE,invert=TRUE)
        updateSelectInput(session, "pdx_predicted", choices=ct )
    })

    observe({
        ngs <- inputData()
        req(ngs)
        ## input$pdx_runbutton
        dbg("[BiomarkerModule::observe2] reacted")

        if(FALSE && isolate(input$pdx_level=="geneset")) {
            ft <- names(COLLECTIONS)
            nn <- sapply(COLLECTIONS, function(x) sum(x %in% rownames(ngs$gsetX)))
            ft <- ft[nn >= 10]
        } else {
            ## gene level
            ft <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
        }
        ft <- sort(ft)
        ##if(input$pdx_level == "gene") ft = sort(c("<custom>",ft))
        ft = sort(c("<custom>",ft))
        updateSelectInput(session, "pdx_filter", choices=ft, selected="<all>")    
    })

    calcVariableImportance <- eventReactive( input$pdx_runbutton, {
        ## 
        ## This code also features a progress indicator.
        ##
        
        ## input$pdx_runbutton
        dbg("[BiomarkerModule::calcVariableImportance] reacted on runbutton")
        
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        req(ngs, input$pdx_predicted)
        
        ct=2
        ct=12
        colnames(ngs$Y)    
        isolate(ct <- input$pdx_predicted)
        
        dbg("  calcVariableImportance","predicting = ",ct,"\n")
        do.survival <- grepl("survival",ct,ignore.case=TRUE)
        dbg("  calcVariableImportance","do.survival = ",do.survival,"\n")
        
        if(is.null(ct)) return(NULL)
        
        dbg("  calcVariableImportance","2")

        NFEATURES=50
        NFEATURES=60

        ## Create a Progress object
        progress <- shiny::Progress$new()
        ## Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Variable importance", value = 0)
        
        if(!(ct %in% colnames(ngs$Y))) return(NULL)
        y0 <- ngs$Y[,ct]
        names(y0) <- rownames(ngs$Y)
        y <- y0[!is.na(y0)]
        
        table(y)
        if(length(y)<40) y <- head(rep(y,10),100)  ## augment to 100 samples
        table(y)
        
        ##-------------------------------------------
        ## select features
        ##-------------------------------------------
        ## group prediction
        if(FALSE && isolate(input$pdx_level)=="geneset") {
            X <- ngs$gsetX[,names(y)]
        } else {
            X <- ngs$X[,names(y)]
        }
        dim(X)
        X0 <- X
        length(y)
        
        ## ----------- filter with selected features
        progress$inc(1/10, detail = "Filtering features")
        
        ft="<all>"
        ft <- isolate(input$pdx_filter)
        if(is.null(ft)) return(NULL)
        
        pp <- rownames(X)
        if(FALSE && isolate(input$pdx_level=="geneset")) {
            if(!(ft %in% names(COLLECTIONS))) return(NULL)
            pp <- COLLECTIONS[[ft]]
        } else {
            xfam <- c(names(ngs$families),names(GSETS))
            if(!ft %in% xfam) return(NULL)
            gg <- GSETS[[8]]
            gg <- NULL
            if(ft %in% names(ngs$families)) {
                gg <- ngs$families[[ft]]
            } else if(ft %in% names(GSETS)) {
                gg <- GSETS[[ft]]
            }
            pp = filterProbes(ngs$genes, gg)
        }
        pp <- intersect(pp,rownames(X))
        X <- X[pp,,drop=FALSE]
        dim(X)
        
        ## ------------- filter with user selection
        ##sel <- input$pdx_select
        isolate(sel <- input_pdx_select())
        if(ft=='<custom>' && !is.null(sel) && length(sel)>0) {
            if(sel[1]!="") {
                ##pp <- intersect(rownames(X),sel)
                pp <- rownames(X)[(toupper(rownames(X)) %in% toupper(sel))]
                X <- X[pp,,drop=FALSE]
            }
        }
        
        ## ----------- restrict to top 100
        dim(X)
        X <- head(X[order(-apply(X,1,sd)),,drop=FALSE], 10*NFEATURES)  ## top 100
        sdx <- mean(apply(X,1,sd))
        X <- X + 0.25*sdx*matrix(rnorm(length(X)),nrow(X),ncol(X))  ## add some noise
        dim(X)
        
        progress$inc(4/10, detail = "computing scores")
        
        ##-------------------------------------------
        ## compute importance values
        ##-------------------------------------------
        if(do.survival) {
            time <- abs(y)
            status <- (y > 0) ## dead is positive time
            methods=c("glmnet","randomforest","boruta","xgboost","pls")
            methods=c("glmnet","randomforest","xgboost","pls")
            P <- pgx.survivalVariableImportance(
                X, time=time, status=status, methods=methods)    
        } else {
            methods=c("glmnet","randomforest","boruta","xgboost","pls")
            methods=c("glmnet","randomforest","xgboost","pls")
            X1=X;y1=y
            names(y1) = colnames(X1) = paste0("x",1:ncol(X))
            P <- pgx.multiclassVariableImportance(X1, y1, methods=methods)    
            ##P <- pgx.variableImportance(X1, y1, methods=methods)
        }
        P <- abs(P)
        head(P)
        
        P[is.na(P)] <- 0
        P[is.nan(P)] <- 0
        P <- t( t(P) / (1e-3+apply(P,2,max,na.rm=TRUE)))
        ##P <- pmax(P,0.1)
        P <- P[order(-rowSums(P,na.rm=TRUE)),,drop=FALSE]
        head(P)
        
        R <- P
        if(nrow(R)>1) {
            R <- (apply(P,2,rank)/nrow(P))**4
            R <- R[order(-rowSums(R)),,drop=FALSE]
            head(R)
        }
        
        if(DEV.VERSION) {
            is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]",rownames(R)))
            is.multiomics    
            do.multiomics <- (is.multiomics && isolate(input$pdx_multiomics))
            if(do.multiomics) {
                dbg("calcVariableImportance:: 5: EXPERIMENTAL: multi-omics weighting")
                ## EXPERIMENTAL: multi-omics weighting 
                rr <- rowMeans(R)
                dtype <- ngs$genes[rownames(R),"data_type"]
                gene  <- ngs$genes[rownames(R),"gene_name"]            
                gg <- sort(unique(gene))
                dtypes <- unique(dtype)
                dt <- "gx"
                rho <- c()
                for(dt in dtypes) {
                    r1 <- rr[which(dtype == dt)]
                    g1 <- gene[which(dtype == dt)]
                    r1 <- r1[match(gg,g1)]
                    rho <- cbind(rho, r1)
                }
                rownames(rho) <- gg
                colnames(rho) <- dtypes
                rho[is.na(rho)] <- 0
                rho <- rho[order(-rowSums(rho,na.rm=TRUE)),]
                dim(rho)            
                barplot(t(head(rho,50)), las=3)
                jj <- match(gene, rownames(rho))
                rname <- rownames(R)
                R <- rho[jj,]
                rownames(R) <- rname
            }
        }
        
        dbg("calcVariableImportance:: 5: drawing tree\n")
        progress$inc(3/10, detail = "drawing tree")
        
        ##------------------------------
        ## create partition tree
        ##------------------------------
        require(rpart)
        R <- R[order(-rowSums(R)),,drop=FALSE]
        sel <- head(rownames(R),100)
        sel <- head(rownames(R),NFEATURES)  ## top50 features
        sel <- intersect(sel,rownames(X))
        tx <- t(X[sel,,drop=FALSE])
        dim(tx)
        
        ## formula wants clean names, so save original names
        colnames(tx) <- gsub("[: +-.,]","_",colnames(tx))
        colnames(tx) <- gsub("[')(]","",colnames(tx))
        colnames(tx) <- gsub("\\[|\\]","",colnames(tx))
        orig.names <- sel
        names(orig.names) <- colnames(tx)
        jj <- names(y)
        ##ny <- length(unique(y))
        ##if(length(jj) < ny*20) jj <- c(jj,jj,jj)
        
        if(do.survival) {
            require(survival)
            time <- abs(y)
            status <- (y>0)    ## dead if positive time
            df <- data.frame( time=time+0.001, status=status, tx)
            ##df <- df[jj,]
            rf <- rpart( Surv(time,status) ~ ., data=df)
        } else {
            df <- data.frame( y=y, tx)
            ##df <- df[jj,]
            rf <- rpart( y ~ ., data=df)
        }
        table(rf$where)
        rf$cptable
        rf$orig.names <- orig.names
        
        rf.nsplit <- rf$cptable[,"nsplit"]
        max(rf.nsplit)
        length(unique(y)) 
        if(grepl("survival",ct)) {
            MAXSPLIT = 4  ## maximum five groups....
        } else {
            MAXSPLIT = 1.5*length(unique(y))  ## maximum N+1 groups
        }
        if( max(rf.nsplit) > MAXSPLIT) {
            cp.idx <- max(which(rf.nsplit <= MAXSPLIT))
            cp0 <- rf$cptable[cp.idx,"CP"]
            ##rf <- prune(rf, cp=0.05)
            rf <- prune(rf, cp=cp0)
        }
        table(rf$where) 
        
        dbg("calcVariableImportance:: *** done!!! ***\n")
        progress$inc(2/10, detail = "done")

        ##y <- y[rownames(ngs$samples)]
        ##tx <- tx[names(y),]
        y <- y[rownames(tx)]
        colnames(tx) <- orig.names[colnames(tx)]
        ##res <- list(P=P, R=R, y=y, X=t(tx), rf=rf)
        res <- list(R=R, y=y, X=t(tx), rf=rf)
        return(res)
    })

    ##================================================================================
    ##==================================== PLOTS =====================================
    ##================================================================================
    
    pdx_importance.RENDER %<a-% reactive({
        
        res <- calcVariableImportance()
        if(is.null(res)) return(NULL)

        dbg("[BiomarkerModule::pdx_importance.RENDER] called\n")
        
        R <- res$R
        R <- R[order(-rowSums(R,na.rm=TRUE)),,drop=FALSE]
        R <- pmax(R,0.05)

        if(FALSE && isolate(input$pdx_level=="geneset")) {
            par(mfrow=c(1,2), oma=c(1,1,1,1)*0.5, mgp=c(2.2,0.8,0))
            par(mar=c(4,8,2,4))
            frame()
            R.top <- head(R,15)
            rownames(R.top) <- tolower(rownames(R.top))
            rownames(R.top) <- substring(rownames(R.top),1,60)
            barplot( t(R.top), las=1, horiz=TRUE,
                    cex.names=0.75, xlab="cumulative importance" )
            klr <- grey.colors(ncol(R))
            legend("topright",
                   legend = colnames(R), fill=klr,
                   cex=0.75, y.intersp=0.8, inset=c(-0.25,0), xpd=NA )

        } else {
            par(mfrow=c(1,1), oma=c(1,1,1,1)*0.2)
            par(mar=c(7,4,0,4))
            R.top <- head(R,40)
            barplot( t(R.top), las=3, horiz=FALSE,
                    cex.names=0.75, ylab="cumulative importance" )
            klr <- grey.colors(ncol(R))
            legend("topright",
                   legend=rev(colnames(R)), fill=rev(klr),
                   cex=0.8, y.intersp=0.8 )
        }
        ## title("Variable importance", cex.main=1.2, line=1.2, adj=0.3)    
    })


    pdx_heatmap.RENDER %<a-% reactive({

        dbg("[BiomarkerModule::pdx_heatmap] reacted\n")
        
        ngs <- inputData()
        alertDataLoaded(session, ngs) 
        req(ngs)

        dbg("[BiomarkerModule::pdx_heatmap] called\n")
        
        res <- calcVariableImportance()
        if(is.null(res)) {
            sendSweetAlert(session=session, title = NULL,
                           text="Please select a target variable to predict, then hit compute.")
            return(NULL)
        }

        cat("<predict:pdx_heatmap> called\n")
        
        gg <- NULL
        if(FALSE && isolate(input$pdx_level=="geneset")) {
            gg <- rownames(res$X)
            gg <- intersect(gg, rownames(ngs$gsetX))
            X <- ngs$gsetX[gg,]
            rownames(X) <- tolower(rownames(X))
        } else {
            gg <- rownames(res$X)
            gg <- intersect(gg, rownames(ngs$X))
            X <- ngs$X[gg,]
        }

        X <- head(X[order(-apply(X,1,sd)),],40)  ## top50

        splitx <- NULL
        ct <- input$pdx_predicted
        do.survival <- grepl("survival",ct,ignore.case=TRUE)
        dbg("calcVariableImportance","do.survival = ",do.survival,"\n")

        splitx <- ngs$Y[colnames(X),ct]
        dbg("<predict:pdx_heatmap> 1:table(splitx)=",table(splitx),"\n")
        if(!is.categorical(splitx) || do.survival) {
            splitx <- NULL
        }
        dbg("<predict:pdx_heatmap> 2:table(splitx)=",table(splitx),"\n")
        
        rownames(X) <- substring(rownames(X),1,40)    
        annot <- ngs$Y[colnames(X),]
        dbg("<predict:pdx_heatmap> dim(annot)=",dim(annot),"\n")
        
        gx.splitmap( X, split=NULL, splitx=splitx, main="  ",
                    show_colnames = FALSE,  ## save space, no sample names
                    show_legend = ifelse(is.null(splitx),TRUE,FALSE),
                    col.annot=annot, annot.ht=2.5, show_rownames=50,
                    lab.len=50, cexRow=0.85, mar=c(2,8))

        dbg("<predict:pdx_heatmap> done!\n")
        
    })    

    pdx_decisiontree.RENDER %<a-% reactive({

        dbg("[BiomarkerModule::pdx_decisiontree] reacted")
        
        res <- calcVariableImportance()
        if(is.null(res)) return(NULL)

        require(rpart)
        require(rpart.plot)
        par(mfrow=c(1,1), mar=c(1,0,2,0))

        is.surv <- grepl("Surv",res$rf$call)[2]
        is.surv
        if(is.surv) {
            require(survival)
            require("partykit")
            (rf <- as.party(res$rf))
            ##(rf <- as.party(prune(res$rf, cp=0.05)))
            plot(rf)
            ##title("Survival tree",cex=1.2,line=0.9,adj=0.35)
            ##table(res$rf$where)
        } else {
            ##visNetworkOutput("pdx_decisiontreeClass")
            ##plotOutput("pdx_decisiontreeClass")
            ##require(visNetwork)
            if(1) {
                rpart.plot(res$rf)
                title("Classification tree",cex=1.2,line=3,adj=0.35)
            } else {
                visTree(res$rf, main="Classification tree",width="100%",legend=FALSE) %>%
                    visInteraction(tooltipStyle='position:fixed;visibility:hidden;padding:5px;white-space:nowrap;font-family:helvetica;font-size:10px;background-color:lightgrey;')
            }
        } 
        dbg("[BiomarkerModule::pdx_decisiontree] done")
        
    })
    
    pdx_boxplots.RENDER %<a-% reactive({
        res <- calcVariableImportance()
        if(is.null(res)) return(NULL)

        dbg("pdx_boxplots:: called\n")    
        vars <- setdiff(res$rf$frame$var,"<leaf>")
        vars <- res$rf$orig.names[vars]
        if(length(vars)==0) return(NULL)
        
        ##vars0 <- setdiff(vars,rownames(res$X))
        ##dbg("pdx_boxplots:: vars0=",vars0,"\n")
        if(FALSE && isolate(input$pdx_level=="geneset")) {
            ##xvars <- res$rf$orig.names[vars]
            vars <- intersect(vars,rownames(res$X))
        } else {
            vars <- intersect(vars,rownames(res$X))
        }

        ## add some other variables
        if(1 && length(vars)<8) {
            jj <- order(-rowSums(res$R,na.rm=TRUE))
            other.vars <- rownames(res$R)[jj]
            vars <- head(unique(c(vars,other.vars)),8)
        }
        
        y <- res$y
        is.surv <- grepl("Surv",res$rf$call)[2]
        is.surv
        if(is.surv) 
        {
            y <- paste0("N",res$rf$where)
            ##y <- as.integer(factor(y, levels=sort(unique(y))))
            names(y) <- colnames(res$X)
            table(y)
        }

        ny <- length(unique(y))    
        par(mfrow=c(2,4), mar=c(3.5,3,2,0.5),
            mgp=c(2,0.8,0), oma=c(0.5,1,1,1)*0)
        if(length(vars)>8) par(mfrow=c(3,4), mar=c(2.8,3,2,0.3),)
        i=1
        for(i in 1:min(12,length(vars))) {
            ##g <- rownames(res$R)[i]
            g <- vars[i]
            gx <- res$X[g,]
            boxplot( gx ~ y, col="grey85", ylim = range(gx),
                    ylab="expression", xlab="", cex.axis=0.001)
            axis(2, cex.axis=0.9)
            ##axis(1, at=1:ny, labels=levels(factor(y)), cex.axis=0.6)
            cex1 <- ifelse( ny >= 8, 0.65, 0.8)
            title(g, cex.main=ifelse(nchar(g)>20,0.85,1))
            nchar(y)
            too.big = (max(nchar(y))>=6 && ny==2) ||
                (max(nchar(y))>=4 && ny>2 || ny>=6 )
            if(too.big) {
                dy <- min(gx) - 0.16*diff(range(gx))
                nn <- max(nchar(y))
                text(1:ny, dy, levels(factor(y)), xpd=NA,
                     cex=cex1, srt=30, adj=1)
            } else {
                mtext(levels(factor(y)), side=1, line=0.7,
                      cex=cex1*0.85, las=1, at=1:ny)
            }

        }

    })


    pdx_importance.opts = tagList()
    callModule(
        plotModule,
        id = "pdx_importance",
        func = pdx_importance.RENDER,
        func2 = pdx_importance.RENDER,
        title = "Variable importance",
        info.text = "A variable importance score for each feature is calculated using multiple machine learning algorithms, including LASSO, elastic nets, random forests, and extreme gradient boosting. By combining several methods, the platform aims to select the best possible biomarkers. The top features are plotted according to cumulative ranking by the algorithms.",
        ##options = pdx_importance.opts,
        label="a",
        pdf.width=10, pdf.height=5,
        height = 270, res=80
        )
    

    pdx_heatmap.opts = tagList()
    callModule(
        plotModule,
        id = "pdx_heatmap",
        func = pdx_heatmap.RENDER,
        func2 = pdx_heatmap.RENDER,
        title = "Heatmap", label="b",
        info.text = "The heatmap shows the expression of genes of the top most important features.",
        ##options = pdx_heatmap.opts,
        pdf.width=10, pdf.height=10,
        height = 370, res=72
    )
    
    pdx_decisiontree.opts = tagList()
    callModule(
        plotModule,
        "pdx_decisiontree",
        func = pdx_decisiontree.RENDER,
        func2 = pdx_decisiontree.RENDER,
        title = "Decision tree", label="c",
        info.text = "The decision tree shows a tree solution for classification based on the top most important features.",
        ##options = pdx_decisiontree.opts,
        pdf.width=10, pdf.height=6,
        height = 315, res=72
    )

    pdx_boxplots.opts = tagList()
    callModule(
        plotModule,
        id = "pdx_boxplots",
        func = pdx_boxplots.RENDER,
        func2 = pdx_boxplots.RENDER,
        title = "Biomarker expression", label="d",
        info.text = "These boxplots shows the expression of genes/samples of the identified features.",
        ##options = pdx_boxplots.opts,
        pdf.width=10, pdf.height=5.5,
        height = 320, res=90
    )

    pdx_biomarker_caption = "<b>Biomarker selection</b>. The expression of certain genes may be used as <i>markers</i> to predict a certain phenotype such as response to a therapy. Finding such <i>biomarkers</i> are of high importance in clinical applications. <b>(a)</b> An importance score for each feature is calculated using multiple machine learning algorithms, including LASSO, elastic nets, random forests, and extreme gradient boosting. The top features are plotted  according to cumulative ranking by the algorithms. <b>(b)</b> The heatmap shows the expression distribution for the top most important features. <b>(c)</b> The decision tree shows (one) tree solution for classification based on the top most important features. <b>(d)</b> Boxplots show the expression of biomarker genes across the groups."

    output$pdx_biomarker_UI <- renderUI({
        fillCol(
            height = fullH,
            flex = c(NA,0.05,1),
            div(HTML(pdx_biomarker_caption),class="caption"),
            br(),
            fillRow(
                flex = c(1,0.1,1),
                height = fullH - 100,                
                fillCol(
                    flex = c(0.7,1),
                    plotWidget(ns("pdx_importance")),
                    plotWidget(ns("pdx_heatmap"))
                ),
                br(), ## spacer
                fillCol(
                    flex = c(1,0.9), 
                    plotWidget(ns("pdx_decisiontree")),
                    plotWidget(ns("pdx_boxplots"))
                )
            )
        )
    })
    outputOptions(output, "pdx_biomarker_UI", suspendWhenHidden=FALSE) ## important!!!


    ## return(NULL)
} ## end-of-Module
