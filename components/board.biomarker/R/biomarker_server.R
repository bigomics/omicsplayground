##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

BiomarkerBoard <- function(id, inputData)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns ## NAMESPACE
    
    fullH = 800  ## full height of panel
    rowH  = 320  ## row height of panel
    imgH  = 260


    pdx_infotext =
        "The <strong>Biomarker Board</strong> performs the biomarker selection that can be used for classification or prediction purposes.

<br><br>To better understand which genes, mutations, or gene sets influence the final phenotype the most, Playground calculates a variable importance score for each feature using state-of-the-art machine learning algorithms, including LASSO, elastic nets, random forests, and extreme gradient boosting, and provides the top 50 features according to cumulative ranking by the algorithms. By combining several methods, the platform aims to select the best possible biomarkers.

<br><br>The phenotype of interest can be multi-categorical classes or patient survival data. Instead of choosing a phenotype, users can also specify a particular contrast from the analysis and perform biomarker selection. The platform also provides a heatmap of samples based on identified top features. 

<br><br>In addition, it generates a classification tree using top features and provides expression boxplots by phenotype classes for features present in the tree. The platform can also provide a survival tree analysis using top features and provides expression boxplots by phenotype classes for features present in the tree."

    ##================================================================================
    ##======================= REACTIVE/OBSERVE FUNCTIONS =============================
    ##================================================================================
    
    shiny::observeEvent( input$pdx_info, {
        dbg("<module-biomarker::observe pdxinfo> reacted")
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Biomarker Board</strong>"),
            shiny::HTML(pdx_infotext),
            easyClose = TRUE, size="l"))
    })

    input_pdx_select <- shiny::reactive({
        dbg("[BiomarkerBoard:<input_pdx_select>]  reacted")
        gg <- input$pdx_select
        if(is.null(gg)) return(NULL)
     
        gg <- strsplit(as.character(gg), split="[, \n\t]")[[1]]
        if(length(gg)==0) return(NULL)
        if(length(gg)==1 && gg[1]!="") gg <- c(gg,gg)  ## hack to allow single gene....
        return(gg)
    }) %>% shiny::debounce(1000)

    shiny::observe({
        ngs <- inputData()
        ##if(is.null(ngs)) return(NULL)
        shiny::req(ngs)
        dbg("[BiomarkerBoard::observe1] reacted")
        dbg("[BiomarkerBoard::observe1] dim(ngs$Y) = ",dim(ngs$Y))
        ct <- colnames(ngs$Y)
        ## ct <- grep("group|sample|patient|donor",ct,value=TRUE,invert=TRUE)
        ## ct <- grep("sample|patient|donor",ct,value=TRUE,invert=TRUE)
        shiny::updateSelectInput(session, "pdx_predicted", choices=ct )
    })

    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        ## input$pdx_runbutton
        dbg("[BiomarkerBoard::observe2] reacted")

        if(FALSE && shiny::isolate(input$pdx_level=="geneset")) {
            ft <- names(COLLECTIONS)
            nn <- sapply(COLLECTIONS, function(x) sum(x %in% rownames(ngs$gsetX)))
            ft <- ft[nn >= 10]
        } else {
            ## gene level
            ##ft <- pgx.getFamilies(ngs,nmin=10,extended=FALSE)
            ft <- names(ngs$families)
        }
        ft <- sort(ft)
        ##if(input$pdx_level == "gene") ft = sort(c("<custom>",ft))
        ft = sort(c("<custom>",ft))
        shiny::updateSelectInput(session, "pdx_filter", choices=ft, selected="<all>")    
    })

    calcVariableImportance <- shiny::eventReactive( input$pdx_runbutton, {
        ## 
        ## This code also features a progress indicator.
        ##
        
        ## input$pdx_runbutton
        dbg("[calcVariableImportance] reacted on runbutton")
        
        ngs <- inputData()
        if(is.null(ngs)) return(NULL)
        shiny::req(ngs, input$pdx_predicted)

        dbg("[calcVariableImportance] 0: ")
        
        ct=2
        ct=12
        colnames(ngs$Y)    
        shiny::isolate(ct <- input$pdx_predicted)        
        do.survival <- grepl("survival",ct,ignore.case=TRUE)
        if(is.null(ct)) return(NULL)

        dbg("[calcVariableImportance] 1: called!")
        
        NFEATURES=50
        NFEATURES=60

        ## Create a Progress object
        progress <- shiny::Progress$new()
        ## Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
        progress$set(message = "Variable importance", value = 0)
        
        if(!(ct %in% colnames(ngs$Y))) return(NULL)
        y0 <- as.character(ngs$Y[,ct])
        names(y0) <- rownames(ngs$Y)
        y <- y0[!is.na(y0)]


        ## augment to 100 samples        
        table(y)
        ##if(length(y)<40) y <- head(rep(y,10),100)
        if(length(y)<100) y <- head(rep(y,100),100)  
        table(y)
        
        ##-------------------------------------------
        ## select features
        ##-------------------------------------------
        ## group prediction
        if(FALSE && shiny::isolate(input$pdx_level)=="geneset") {
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
        ft <- shiny::isolate(input$pdx_filter)
        if(is.null(ft)) return(NULL)
        shiny::isolate(sel <- input_pdx_select())

        is.family <- (ft %in% c(names(ngs$families),names(iGSETS)))
        
        if(ft=='<custom>' && !is.null(sel) && length(sel)>0) {
            ## ------------- filter with user selection
            if(sel[1]!="") {
                dbg("[calcVariableImportance] 2: using custom list of variable ")                
                ##pp <- intersect(rownames(X),sel)
                pp <- rownames(X)[which(toupper(rownames(X)) %in% toupper(sel))]
                X <- X[pp,,drop=FALSE]
            }
        } else if(is.family) {
            pp <- rownames(X)
            if(ft %in% names(ngs$families)) {
                dbg("[calcVariableImportance] 2: using ngs$families")
                gg <- ngs$families[[ft]]
                pp <- filterProbes(ngs$genes, gg)
            } else if(ft %in% names(iGSETS)) {
                dbg("[calcVariableImportance] 2: using genesets")                    
                gg <- unlist(getGSETS(ft))
                pp <- filterProbes(ngs$genes, gg)                
            }
            pp <- intersect(pp,rownames(X))
            X <- X[pp,,drop=FALSE]
        } else {
            dbg("[calcVariableImportance] 2: using all features")            
        }

        dbg("[calcVariableImportance] 3: dim.X = ",dim(X))                
        
        ## ----------- restrict to top 100
        dim(X)
        X <- head(X[order(-apply(X,1,sd)),,drop=FALSE], 10*NFEATURES)  ## top 100
        sdx <- mean(apply(X,1,sd))
        X <- X + 0.25*sdx*matrix(rnorm(length(X)),nrow(X),ncol(X))  ## add some noise
        dim(X)

        dbg("[calcVariableImportance] 4: dim.X = ",dim(X))
        
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
        
        if(FALSE && DEV) {
            is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]",rownames(R)))
            is.multiomics    
            do.multiomics <- (is.multiomics && shiny::isolate(input$pdx_multiomics))
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

        R <- R[order(-rowSums(R)),,drop=FALSE]
        sel <- head(rownames(R),100)
        sel <- intersect(sel,rownames(X))
        sel <- head(rownames(R),NFEATURES)  ## top50 features
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
            time <- abs(y)
            status <- (y>0)    ## dead if positive time
            df <- data.frame( time=time+0.001, status=status, tx)
            ##df <- df[jj,]
            rf <- rpart::rpart( survival::Surv(time,status) ~ ., data=df)
        } else {
            df <- data.frame( y=y, tx)
            ##df <- df[jj,]
            rf <- rpart::rpart( y ~ ., data=df)
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
            ##rf <- rpart::prune(rf, cp=0.05)
            rf <- rpart::prune(rf, cp=cp0)
        }
        table(rf$where) 
        
        dbg("[calcVariableImportance] done!!!\n")
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
    
    pdx_importance.RENDER <- shiny::reactive({
        
        res <- calcVariableImportance()
        if(is.null(res)) return(NULL)

        dbg("[BiomarkerBoard::pdx_importance.RENDER] called\n")
        
        R <- res$R
        R <- R[order(-rowSums(R,na.rm=TRUE)),,drop=FALSE]
        R <- pmax(R,0.05)

        if(FALSE && shiny::isolate(input$pdx_level=="geneset")) {
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
            par(mar=c(5,4,0,4))
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


    pdx_heatmap.RENDER <- shiny::reactive({

        dbg("[BiomarkerBoard::pdx_heatmap] reacted\n")
        
        ngs <- inputData()
        alertDataLoaded(session, ngs) 
        shiny::req(ngs)

        dbg("[BiomarkerBoard::pdx_heatmap] called\n")
        
        res <- calcVariableImportance()

        dbg("[BiomarkerBoard::pdx_heatmap] is.null(res) = ", is.null(res) )
        
        if(is.null(res)) {
            ##shinyWidgets::sendSweetAlert( session=session, title = NULL,
            ##                             text="Please select a target variable to predict, then hit compute.")
            return(NULL)
        }

        cat("<predict:pdx_heatmap> called\n")        
        gg <- NULL
        if(FALSE && shiny::isolate(input$pdx_level=="geneset")) {
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
        ## X <- X + 1e-4 * matrix(rnorm(length(X)),nrow(X),ncol(X)) ## add noise...
        
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
        sdx <- apply(X,1,sd)
        
        dbg("[pdx_heatmap] dim(annot)=",dim(annot),"\n")
        dbg("[pdx_heatmap] sum.NA(X)=",sum(is.na(X)),"\n")
        dbg("[pdx_heatmap] sum.sdx0=",sum(sdx < 1e-3),"\n")
        
        gx.splitmap( X, split=NULL, splitx=splitx, main="  ",
                    dist.method = "euclidean",                    
                    show_colnames = FALSE,  ## save space, no sample names
                    show_legend = ifelse(is.null(splitx),TRUE,FALSE),
                    key.offset = c(0.05,1.03),
                    ## col.annot=annot, annot.ht=2.5,
                    show_rownames = 99,
                    lab.len=50, cexRow=0.88, mar=c(2,8))

        dbg("<predict:pdx_heatmap> done!\n")
        
    })    

    pdx_decisiontree.RENDER <- shiny::reactive({

        dbg("[BiomarkerBoard::pdx_decisiontree] reacted")
        
        res <- calcVariableImportance()
        if(is.null(res)) return(NULL)

        par(mfrow=c(1,1), mar=c(1,0,2,0))
        is.surv <- grepl("Surv",res$rf$call)[2]
        is.surv
        if(is.surv) {

            (rf <- partykit::as.party(res$rf))
            ##(rf <- partykit::as.party(rpart::prune(res$rf, cp=0.05)))
            partykit::plot.party(rf)
            ##title("Survival tree",cex=1.2,line=0.9,adj=0.35)
            ##table(res$rf$where)
        } else {
            ##visNetworkOutput("pdx_decisiontreeClass")
            ##plotOutput("pdx_decisiontreeClass")

            if(1) {
                rpart.plot::rpart.plot(res$rf)
                title("Classification tree",cex=1.2,line=3,adj=0.35)
            } else {
                visNetwork::visTree(res$rf, main="Classification tree",width="100%",legend=FALSE) %>%
                    visNetwork::visInteraction(tooltipStyle='position:fixed;visibility:hidden;padding:5px;white-space:nowrap;font-family:helvetica;font-size:10px;background-color:lightgrey;')
            }
        } 
        dbg("[BiomarkerBoard::pdx_decisiontree] done")
        
    })
    
    pdx_boxplots.RENDER <- shiny::reactive({

        res <- calcVariableImportance()
        if(is.null(res)) return(NULL)

        dbg("[BiomarkerBoard::pdx_boxplots] called\n")    
        vars <- setdiff(res$rf$frame$var,"<leaf>")
        vars <- res$rf$orig.names[vars]
        if(length(vars)==0) return(NULL)
        
        ##vars0 <- setdiff(vars,rownames(res$X))
        ##dbg("pdx_boxplots:: vars0=",vars0,"\n")
        if(FALSE && shiny::isolate(input$pdx_level=="geneset")) {
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
        if(length(vars)>8) par(mfrow=c(3,4), mar=c(2.8,3,2,0.3))
        i=1
        for(i in 1:min(12,length(vars))) {
            ##g <- rownames(res$R)[i]
            g  <- vars[i]
            gx <- res$X[g,]
            boxplot( gx ~ y, col="grey85", ylim = range(gx),
                    ylab="expression", xlab="", cex.axis=0.001)
            axis(2, cex.axis=0.9)
            ##axis(1, at=1:ny, labels=levels(factor(y)), cex.axis=0.6)
            cex1 <- ifelse( ny >= 8, 0.65, 0.8)
            title(g, cex.main=ifelse(nchar(g)>20,0.85,1))
            nchar(y)
            too.big = (max(nchar(y))>=8 && ny==2) ||
                (max(nchar(y))>=5 && ny %in% c(3,4) ||
                 ny>=5 )
            if(too.big) {
                dy <- min(gx) - 0.12*diff(range(gx))
                text(1:ny, dy, levels(factor(y)), xpd=NA,
                     cex=cex1, srt=30, adj=1)
            } else {
                mtext(levels(factor(y)), side=1, line=0.7,
                      cex=cex1*0.85, las=1, at=1:ny)
            }

        }

    })


    pdx_importance.opts = shiny::tagList()
    shiny::callModule(
        plotModule,
        id = "pdx_importance",
        func = pdx_importance.RENDER,
        func2 = pdx_importance.RENDER,
        title = "Variable importance",
        info.text = "<b>Variable importance.</b>. An importance score for each variable is calculated using multiple machine learning algorithms, including LASSO, elastic nets, random forests, and extreme gradient boosting. By combining several methods, the platform aims to select the best possible biomarkers. The top features are plotted according to cumulative ranking by the algorithms.",
        ##options = pdx_importance.opts,
        label="a",
        pdf.width=10, pdf.height=5,
        height = 235, res = 78,
        add.watermark = WATERMARK
    )
    

    pdx_heatmap.opts = shiny::tagList()
    shiny::callModule(
        plotModule,
        id = "pdx_heatmap",
        func = pdx_heatmap.RENDER,
        func2 = pdx_heatmap.RENDER,
        title = "Heatmap", label="b",
        info.text = "<b>Biomarker heatmap.</b> Expression heatmap of top gene features according to their variable importance.",
        ##options = pdx_heatmap.opts,
        pdf.width=10, pdf.height=10,
        height = 435, res=72,
        add.watermark = WATERMARK
    )
    
    pdx_decisiontree.opts = shiny::tagList()
    shiny::callModule(
        plotModule,
        "pdx_decisiontree",
        func = pdx_decisiontree.RENDER,
        func2 = pdx_decisiontree.RENDER,
        title = "Decision tree", label="c",
        info.text = "The decision tree shows a tree solution for classification based on the top most important features.",
        ##options = pdx_decisiontree.opts,
        pdf.width=10, pdf.height=6,
        height = 315, res=72,
        add.watermark = WATERMARK
    )

    pdx_boxplots.opts = shiny::tagList()
    shiny::callModule(
        plotModule,
        id = "pdx_boxplots",
        func = pdx_boxplots.RENDER,
        func2 = pdx_boxplots.RENDER,
        title = "Biomarker expression", label="d",
        info.text = "These boxplots shows the expression of genes/samples of the identified features.",
        ##options = pdx_boxplots.opts,
        pdf.width=10, pdf.height=5.5,
        height = 320, res=90,
        add.watermark = WATERMARK
    )


  })
} ## end-of-Board
