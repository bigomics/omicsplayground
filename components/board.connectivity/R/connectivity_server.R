##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


ConnectivityBoard <- function(id, inputData)
{
  moduleServer(id, function(input, output, session)
  {
      
    ns <- session$ns ## NAMESPACE
    fullH = 750       # row height of panel
    tabH = '70vh'

    cmap_infotext =
        "The <strong>Experiment connectivity</strong> module enables users to compare their data to other datasets. For the selected contrast, this module provides pairwise correlation plots and/or enrichment plots with signatures from other data sets. The <strong>Connectivity map</strong> shows the similarity of the contrasts profiles as a t-SNE plot.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=5' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>

"
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$cmap_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Connectivity Analysis Board</strong>"),
            shiny::HTML(cmap_infotext),
            easyClose = TRUE, size="l" ))
    })

    ## update choices upon change of data set 
    shiny::observe({
        ngs <- inputData()
        ##req(ngs)
        if(is.null(ngs)) return(NULL)
        comparisons <- colnames(ngs$model.parameters$contr.matrix)
        comparisons <- sort(comparisons)
        shiny::updateSelectInput(session, "cmap_contrast", choices=comparisons,
                          selected = head(comparisons,1))

        sigdb <- c("A","B","C")
        sigdb0 <- dir(SIGDB.DIR, pattern="sigdb-.*h5")
        sigdb <- names(ngs$connectivity)  ## only precomputed inside PGX object??
        sigdb <- sort(intersect(sigdb, sigdb0))
        sigdb
        sel <- sigdb[1]
        ##if(any(grepl("sigdb",sigdb))) sel <- grep("sigdb",sigdb,value=TRUE)[1]
        shiny::updateSelectInput(session, "cmap_sigdb", choices=sigdb, selected=sel)
    })
    
    shiny::observe({
        ## reset CMap threshold zero/max
        res <- getConnectivityScores()
        shiny::req(res)
        max <- round(0.999*max(abs(res$score),na.rm=TRUE),digits=1)
        max <- round(0.999*tail(sort(abs(res$score)),10)[1],digits=1)
        shiny::updateSliderInput(session, 'cmap_scorethreshold', value=0, max=max)
        
    })
    
    ## update choices upon change of chosen contrast
    ##observeEvent( input$cmap_level, {
    shiny::observeEvent( input$cmap_contrast, {

        ngs <- inputData()
        shiny::req(ngs)

        ## reset CMap threshold zero/max
        res <- getConnectivityScores()  ## result gets cached
        shiny::req(res)
        max <- round(0.999*max(abs(res$score),na.rm=TRUE),digits=1)
        shiny::updateSliderInput(session, 'cmap_scorethreshold', value=0, max=max)
                               
    })
    
    getCurrentContrast <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs, input$cmap_contrast)        
        ct=1
        ct <- input$cmap_contrast
        fc <- ngs$gx.meta$meta[[ct]]$meta.fx
        names(fc) <- rownames(ngs$gx.meta$meta[[ct]])
        gs <- ngs$gset.meta$meta[[ct]]$meta.fx
        names(gs) <- rownames(ngs$gset.meta$meta[[ct]])
        names(fc) <- toupper(names(fc))  ## de-MOUSE
        list(name=ct, fc=fc, gs=gs)
    })
    
    ##================================================================================
    ##========================= REACTIVE FUNCTIONS ===================================
    ##================================================================================

    getConnectivityFullPath <- function(sigdb) {
        db.exists <- sapply( SIGDB.DIR, function(d) file.exists(file.path(d,sigdb)))
        db.exists
        db.dir <- names(which(db.exists))[1]
        db.dir
        file.path(db.dir, sigdb)                
    }
    
    getConnectivityContrasts <- function(sigdb) {
        dbg("[getConnectivityContrasts] sigdb=",sigdb)
        if(length(sigdb)==0 || is.null(sigdb) || sigdb=="" ) return(NULL)

        db <- getConnectivityFullPath(sigdb)
        dbg("[getConnectivityContrasts] db=",db)
        cn <- NULL
        if(file.exists(db)) {

            cn <- rhdf5::h5read(db, "data/colnames")
        }
        cn
    }
    
    getConnectivityMatrix <- function(sigdb, select=NULL, genes=NULL)
    {

        if(0) {
            dbg("[getConnectivityMatrix] reacted")
            sigdb = "sigdb-archs4.h5"
            sigdb = "sigdb-creeds.h5"
            sigdb <- input$cmap_sigdb
            shiny::req(sigdb)
        }
        if(sigdb=="" || is.null(sigdb)) {
            dbg("[getConnectivityMatrix] ***WARNING*** sigdb=",sigdb)
            return(NULL)
        }
                
        db.exists <- sapply( SIGDB.DIR, function(d) file.exists(file.path(d,sigdb)))
        db.exists
        X <- NULL
        if(any(db.exists)) {
            db.dir <- names(which(db.exists))[1]
            db.dir
            if(grepl("csv$",sigdb)) {
                X <- read.csv(file.path(db.dir, sigdb), row.names=1, check.names=FALSE)
                X <- as.matrix(X)
                X <- X[,colMeans(is.na(X)) < 0.99,drop=FALSE]  ## omit empty columns
                if(!is.null(genes)) X <- X[intersect(genes,rownames(X)),,drop=FALSE]
                if(!is.null(select)) X <- X[, intersect(select,colnames(X))]
            }
            if(grepl("h5$",sigdb)) {
                h5.file <- file.path(db.dir, sigdb)                
                cn <- rhdf5::h5read(h5.file, "data/colnames")
                rn <- rhdf5::h5read(h5.file, "data/rownames")
                rowidx <- 1:length(rn)
                colidx <- 1:length(cn)
                if(!is.null(genes)) rowidx <- match(intersect(genes,rn),rn)                
                if(!is.null(select)) colidx <- match(intersect(select,cn),cn)

                nr <- length(rowidx)
                nc <- length(colidx)
                dbg("*** WARNING *** reading large H5 file:",nr,"x",nc,"")

                X  <- rhdf5::h5read(h5.file, "data/matrix", index = list(rowidx,colidx) )
                rownames(X) <- rn[rowidx]
                colnames(X) <- cn[colidx]
            }
        } else {
            cat("[getConnectivityMatrix] WARNING: could not retrieve matrix\n")
            ## X <- as.matrix(PROFILES$FC)
            ## X <- X[,colMeans(is.na(X)) < 0.99,drop=FALSE]  ## omit empty columns
            ## if(!is.null(genes)) X <- X[intersect(genes,rownames(X)),,drop=FALSE]
            ## if(!is.null(select)) X <- X[, intersect(select,colnames(X))]
        }
        class(X)        
        return(X)
    }

    getEnrichmentMatrix <- function(sigdb, select=NULL, nc=-1)
    {

        if(sigdb=="" || is.null(sigdb)) {
            dbg("[getEnrichmentMatrix] ***WARNING*** sigdb=",sigdb)
            return(NULL)
        }
        if(!grepl("h5$",sigdb)) {
            stop("getEnrichmentMatrix:: only for H5 database files")
            return(NULL)
        }

        h5exists <- function(h5.file, obj) {
            xobjs <- apply(rhdf5::h5ls(h5.file)[,1:2],1,paste,collapse="/")
            obj %in% gsub("^/|^//","",xobjs)
        }
                
        db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d,sigdb)))
        db.exists
        Y <- NULL
        if(any(db.exists)) {            
            db.dir <- names(which(db.exists))[1]
            db.dir
            h5.file <- file.path(db.dir, sigdb)                
            cn <- rhdf5::h5read(h5.file, "data/colnames")

            has.gs   <- h5exists(h5.file, "enrichment/genesets")
            has.gsea <- h5exists(h5.file, "enrichment/GSEA") 
            if(!has.gs && has.gsea) {
                dbg("[getEnrichmentMatrix] WARNING: PGX object has no enrichment results")
                return(NULL)
            }

            rn <- rhdf5::h5read(h5.file, "enrichment/genesets")
            rowidx <- 1:length(rn)
            colidx <- 1:length(cn)
            if(!is.null(select)) colidx <- match(intersect(select,cn),cn)
            Y  <- rhdf5::h5read(h5.file, "enrichment/GSEA", index = list(rowidx,colidx) )
            rownames(Y) <- rn[rowidx]
            colnames(Y) <- cn[colidx]
            dim(Y)            
            sdy <- apply(Y,1,sd)
            Y <- Y[order(-sdy),]
        }

        ## cluster genesets into larger groups
        if(nc>0) {
            hc <- hclust(dist(Y[,]))
            idx <- paste0("h",cutree(hc, nc))
            Y2 <- tapply( 1:nrow(Y), idx, function(i) colMeans(Y[i,,drop=FALSE]))
            Y2 <- do.call(rbind, Y2)
            idx.names <- tapply(rownames(Y),idx,paste,collapse=",")
            idx.names <- gsub("H:HALLMARK_","",idx.names)
            idx.names <- gsub("C2:KEGG_","",idx.names)
            rownames(Y2) <- as.character(idx.names[rownames(Y2)])
            Y <- Y2
        }

        if(nrow(Y)==0) {
            return(NULL)
        }

        class(Y)        
        return(Y)
    }
    
    getSignatureMatrix <- function(sigdb) {

        if(sigdb=="" || is.null(sigdb)) {
            dbg("[getSignatureMatrix] ***WARNING*** sigdb=",sigdb)
            return(NULL)
        }
        ##if(!is.null(select)) dbg("[getEnrichmentMatrix] length(select)=",length(select))

        if(!grepl("h5$",sigdb)) {
            stop("getEnrichmentMatrix:: only for H5 database files")
        }
        
        db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d,sigdb)))
        db.exists
        up=dn=NULL
        if(any(db.exists)) {

            db.dir <- names(which(db.exists))[1]
            db.dir
            h5.file <- file.path(db.dir, sigdb)
            rhdf5::h5ls(h5.file)
            cn <- rhdf5::h5read(h5.file, "data/colnames")
            dn <- rhdf5::h5read(h5.file, "signature/sig100.dn")
            up <- rhdf5::h5read(h5.file, "signature/sig100.up")
            colnames(dn) <- cn
            colnames(up) <- cn
        }
        list(up=up, dn=dn)
    }
    
    getConnectivityPositions <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)
        
        ## get the foldchanges of selected comparison and neighbourhood
        method="umap";dims=2
        method <- input$cmap_layout
        ##if(!method %in% c("pca","tsne","umap","volcano")) return(NULL)

        if(method=="volcano") {
            res <- getConnectivityScores()
            pos <- data.frame( rho=res$rho, NES=res$NES )
            rownames(pos) <- res$pathway
            colnames(pos) <- c("rho","NES")
            return(pos)
        }
        
        do3d = FALSE
        do3d = TRUE
        ##do3d = input$cmap_3d
        do3d = "3D" %in% input$cmap_plotoptions
        dims <- 2 + 1*do3d

        sigdb = "sigdb-archs4.h5"
        sigdb <- input$cmap_sigdb
        shiny::req(sigdb)
        h5.ref <- grepl("h5$",sigdb)
        h5.file <- NULL
        pos <- NULL                
        
        if(!h5.ref) {
            stop("sigdb must be H5 format!")
        }
            
        h5.file = "/home/kwee/Playground/omicsplayground/libx/sigdb-archs4.h5"
        h5.file
        db.exists <- sapply(SIGDB.DIR, function(d) file.exists(file.path(d,sigdb)))
        db.exists
        if(!any(db.exists)) {
            cat("*** WARNING *** cannot locate signature matrix file")
            return(NULL)
        }
        db.dir <- names(which(db.exists))[1]
        h5.file <- file.path(db.dir, sigdb)
        h5.file


        rhdf5::h5closeAll()
        rhdf5::h5ls(h5.file)
        
        pos <- NULL
        if(method == "pca" && dims==2) try(pos  <- rhdf5::h5read(h5.file, "clustering/pca2d"))
        if(method == "pca" && dims==3) try(pos  <- rhdf5::h5read(h5.file, "clustering/pca3d"))
        if(method == "tsne" && dims==2) try(pos <- rhdf5::h5read(h5.file, "clustering/tsne2d"))
        if(method == "tsne" && dims==3) try(pos <- rhdf5::h5read(h5.file, "clustering/tsne3d"))
        if(method == "umap" && dims==2) try(pos <- rhdf5::h5read(h5.file, "clustering/umap2d"))
        if(method == "umap" && dims==3) try(pos <- rhdf5::h5read(h5.file, "clustering/umap3d"))
        dim(pos)

        ## normalize
        pos <- scale(pos)

        ## set row/colnames
        cn <- rhdf5::h5read(h5.file, "data/colnames")
        rownames(pos) <- cn
        xyz <- c("x","y","z")
        colnames(pos) <- paste0(toupper(method),"-",xyz[1:ncol(pos)])
        dbg("[getConnectivityPositions] retrieved",nrow(pos),"cluster positions from H5 file")

        return(pos)
    })
    
    getConnectivityScores <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs, input$cmap_contrast)
        shiny::validate(shiny::need("connectivity" %in% names(ngs), "no 'connectivity' in object."))
        
        ntop = 1000
        sigdb = "sigdb-gtex.h5"
        sigdb = "sigdb-lincs.h5"
        sigdb = "sigdb-lincs-cp.h5"
        sigdb = "sigdb-creeds.h5"
        sigdb = "sigdb-archs4.h5"
        sigdb <- input$cmap_sigdb
        shiny::req(sigdb)

        all.scores <- NULL
        if(sigdb %in% names(ngs$connectivity)) {
            dbg("[getConnectivityScores] extracting precomputed scores")
            all.scores <- ngs$connectivity[[sigdb]]
        } else {
            dbg("[getConnectivityScores] ERROR : could not get scores")
            return(NULL)
        }
        
        ct = names(all.scores)[1]
        ct <- input$cmap_contrast
        if(!ct %in% names(all.scores)) {
            dbg("[getConnectivityScores] ERROR : contrast not in connectivity scores")
            return(NULL)
        }
        
        names(all.scores)
        scores <- as.data.frame(all.scores[[ct]])
        if(input$cmap_abs_score==FALSE) {
            ## put sign back!!!
            scores$score <- scores$score * sign(scores$rho) 
        }
        scores <- scores[order(-abs(scores$score)),]
        scores <- scores[!duplicated(scores$pathway),]
        rownames(scores) <- scores$pathway
        dim(scores)

        if(nrow(scores)==0 || ncol(scores)==0 ) {
            dbg("[getConnectivityScores] ERROR : scores has zero dimensions")
            return(NULL)
        }
        
        if(input$cmap_hideclustcontrasts) {
            sel <- grep("cluster[:]",scores$pathway,invert=TRUE)
            scores <- scores[sel,,drop=FALSE]
        }
        
        ## only those in existing database
        cts <- getConnectivityContrasts(sigdb)
        scores <- scores[which(rownames(scores) %in% cts),,drop=FALSE]
        dim(scores)
        
        ## filter on significance
        qsig = 0.05
        qsig <- input$connectivityScoreTable_qsig
        scores <- scores[which(scores$padj <= qsig),,drop=FALSE]        
        scores <- scores[order(-scores$score),,drop=FALSE]
        
        no.le <- !("leadingEdge" %in% colnames(scores))
        no.le
        abs_score = TRUE
        abs_score <- input$cmap_abs_score
        ntop = 100
        ##if(DEV) ntop <- input$cmap_le_ntop
        
        if(no.le && abs_score==TRUE) {            
            
            dbg("[getConnectivityScores] recomputing leading edge")
            ## recreate "leadingEdge" list
            sig <- getSignatureMatrix(sigdb)
            fc  <- getCurrentContrast()$fc
            fc  <- fc[order(-abs(fc))]

            fc.up <- head(names(fc[fc>0]),ntop)
            fc.dn <- head(names(fc[fc<0]),ntop)
            ff <- c(fc.up, fc.dn)            
            e1 <- apply(sig$up, 2, function(g) intersect(ff,g))
            e2 <- apply(sig$dn, 2, function(g) intersect(ff,g))
            ee <- mapply(c, e1, e2)
            ee <- ee[match(scores$pathway,names(ee))]
            scores$leadingEdge <- ee
        }
        if(no.le && abs_score==FALSE) {            
            dbg("[getConnectivityScores] recomputing leading edge")
            ## recreate "leadingEdge" list
            sig <- getSignatureMatrix(sigdb)
            fc  <- getCurrentContrast()$fc
            fc  <- fc[order(-abs(fc))]
            
            fc.up <- head(names(fc[fc>0]), ntop)
            fc.dn <- head(names(fc[fc<0]), ntop)
            
            p1  <- apply(sig$up, 2, function(g) intersect(fc.up,g))
            p2  <- apply(sig$dn, 2, function(g) intersect(fc.dn,g))
            pp <- mapply(c, p1, p2)

            n1  <- apply(sig$up, 2, function(g) intersect(fc.dn,g))
            n2  <- apply(sig$dn, 2, function(g) intersect(fc.up,g))
            nn <- mapply(c, n1, n2)

            ee <- vector("list",nrow(scores))
            pos.rho <- which(scores$rho >= 0)
            neg.rho <- which(scores$rho < 0)
            ee[pos.rho] <- pp[match(scores$pathway[pos.rho],names(pp))]
            ee[neg.rho] <- nn[match(scores$pathway[neg.rho],names(nn))]
            scores$leadingEdge <- ee
        }

        ## bail out
        if(nrow(scores)==0) {
            return(NULL)
        }

        return(scores)
    })


    getThresholdedConnectivityScores <- shiny::reactive({

        dbg("[getThresholdedConnectivityScores] reacted")
        
        res <- getConnectivityScores()
        if(is.null(res)) return(NULL)

        ## add all points (not within score table)
        pos0 <- getConnectivityPositions()  ## all points
        pp <- rownames(pos0)
        res <- res[match(pp,rownames(res)),,drop=FALSE]
        rownames(res) <- rownames(pos0)
        res$pathway <- rownames(pos0)
        res$score[is.na(res$score)] <- 0
        res$rho[is.na(res$rho)] <- 0

        ## select on minimum score
        minscore <- input$cmap_scorethreshold
        shiny::req(input$cmap_scorethreshold)
        minscore <- min(minscore, 0.999*max(abs(res$score),na.rm=TRUE))
        sel <- which(abs(res$score) >= minscore)
        if(length(sel)==0) return(NULL)
        res <- res[sel,,drop=FALSE]
        
        ## rownames(res) <- res$pathway
        cmap_grouped <- "grouped" %in% input$cmap_plotoptions
        if(cmap_grouped) {
            res <- res[order(-res$score),]
            dset <- sub("].*","]",res$pathway)
            names(dset) <- res$pathway
            idx <- which(!duplicated(dset))
            ##idx <- tapply(1:nrow(res), res$pathway, function(i) i)
            resx <- res[idx,,drop=FALSE]
            tt <- table(dset)[dset[resx$pathway]]
            resx$pathway <- paste0(resx$pathway," (+",tt," contrasts)")
            resx$datasets <- sub("].*","]",resx$pathway)
            rownames(resx) <- resx$datasets
            res <- resx
        }

        dbg("[getThresholdedConnectivityScores] dim(res)=",dim(res))
        
        return(res)
    })


    ##================================================================================
    ## Correlation score table
    ##================================================================================

    PERTINFO <- NULL
    pert_info.file = file.path(FILESX,"GSE92742_Broad_LINCS_pert_info.txt")
    if(file.exists(pert_info.file)) {
        PERTINFO <- read.csv(pert_info.file,sep="\t",row.names=1)
    }
    
    connectivityScoreTable.RENDER <- shiny::reactive({
        
        df <- getConnectivityScores()
        if(is.null(df)) return(NULL)

        kk <- c("pathway","score","rho","NES","padj","size","leadingEdge")
        kk <- c("pathway","score","rho","NES","padj","leadingEdge")
        kk <- intersect(kk,colnames(df))
        ##df <- df[,kk,with=FALSE]
        df <- df[,kk]
        df <- df[abs(df$score)>0,,drop=FALSE]
        
        ##--------- temporarily add LINCS descriptive name !!!!!!!!!!!!!! -----------------
        if(DEV && input$cmap_sigdb=="sigdb-lincs.h5" && !is.null(PERTINFO)) {
            dd <- sub("\\|.*","", df$pathway)
            pert_iname = PERTINFO[match(dd,rownames(PERTINFO)),"pert_iname"]
            df$pathway <- paste0(df$pathway," (",pert_iname,")")
        }
        ##---------- temporarily add LINCS descriptive name !!!!!!!!!!!!!! -----------------
        
        ##colnames(df) <- sub("padj","NES.q",colnames(df))
        df$leadingEdge <- shortstring(sapply(df$leadingEdge,paste,collapse=","),40)
        df$pathway <- shortstring(df$pathway,100)
        df$leadingEdge <- NULL
        
        colnames(df) <- sub("pathway","dataset/contrast",colnames(df))
        score.col <- which(colnames(df) == "score")
        numcols <- c("score","pval","padj","NES.q","ES","NES","rho","R2")
        numcols <- intersect(numcols, colnames(df))
        
        DT::datatable( df,
                      rownames = FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip',
                          pageLength = 99999,
                          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = tabH,
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      )  %>%
            DT::formatSignif(numcols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle( "score",
                                ##background = DT::styleColorBar( score.col, 'lightblue'),
                                background = color_from_middle(
                                    df[,"score"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    })

    connectivityScoreTable.RENDER2 <- shiny::reactive({
        connectivityScoreTable.RENDER() %>%
            DT::formatStyle(0, target='row', fontSize='13px', lineHeight='90%')
    })

    
    connectivityScoreTable_info = "<b>Similarity scores.</b> Normalized enrichment scores (NES) and Pearson correlation (rho) of reference profiles with respect to the currently selected contrast. The top 100 up/down genes are considered for the calculation of rho or NES. The score is calculated as rho^2*NES. "

    connectivityScoreTable_opts = shiny::tagList(
        shiny::selectInput(ns("connectivityScoreTable_qsig"),"threshold (padj)",
                    c(0.01,0.05,0.2,1),selected=1)
    )
    
    connectivityScoreTable <- shiny::callModule(
        tableModule,
        id = "connectivityScoreTable", label = "b",
        func = connectivityScoreTable.RENDER,
        info.text = connectivityScoreTable_info,
        options = connectivityScoreTable_opts,
        info.width = "300px",
        title = "Similarity scores",
        height = c(260,720), width = c('auto',1280)
    )
    
    ##============================================================================
    ## FC correlation/scatter plots
    ##============================================================================
    
    mfplots=c(4,5)
    cmap_FCFCscatter <- function(fc, F, mfplots, ylab)
    {
        ## get the foldchanges of selected comparison and neighbourhood
        F0 <- F
        F[is.na(F)] <- 0  ## really??
        names(fc) <- toupper(names(fc))
        gg <- intersect(names(fc),rownames(F)) ## uppercase for MOUSE
        fc <- fc[gg]
        
        nplots <- mfplots[1]*mfplots[2]
        F <- F[gg,1:min(nplots,ncol(F)),drop=FALSE]        
        F0 <- F0[gg,colnames(F),drop=FALSE]        
        i=1
        par(mfrow=mfplots, mar=c(5.1,1.6,0.2,0.5),
            mgp=c(2.6,0.7,0), oma=c(0,3,0,0))
        i=1
        for(i in 1:ncol(F)) {
            ct1 <- colnames(F)[i]
            ct1x <- sub("\\]","]\n",ct1)
            ##ct1x <- substring(sub("[:]",":\n",ct1),1,50)
            nna <- (is.na(fc) | is.na(F0[,ct1]))
            col <- c("grey15","grey70")[1 + nna]            
            base::plot( F[,ct1], fc, pch=20, cex=0.5,
                 cex.lab=0.9, cex.axis=0.9,
                 xlab = ct1x, ylab = "", col=col)
            abline(v=0, h=0, lty=2, lwd=0.5)
            abline(lm( fc ~ F0[,ct1]), col="red")
            if(i%%mfplots[2]==1) {
                mtext( ylab, 2, line=3, cex=0.60 )
            }            
        }        
    }

    mfplots=c(4,5)
    cmap_FCFCenplot <- function(fc, F, mfplots, ylab, res)
    {
        names(fc) <- toupper(names(fc))        
        nplots <- mfplots[1]*mfplots[2]        
        i=1
        par(mfrow=mfplots, mar=c(0.1,4,2.6,1))        
        for(i in 1:min(ncol(F),nplots)) {

            j1 <- head(order(F[,i]),100)
            j2 <- head(order(-F[,i]),100)
            gset.dn <- rownames(F)[j1]
            gset.up <- rownames(F)[j2]
            gset.both <- c(gset.dn,gset.up)
            rnk <- fc
            pw <- colnames(F)[i]
            if(0) {
                gsea.enplot.UPDN(
                    rnk, gset.up, gset.dn, sum.trace=FALSE,
                    xlab="", main=pw, cex.main=0.8, len.main=32)
            } else {
                gsea.enplot(abs(rnk), gset.both, xlab="",
                            main=pw, cex.main=0.8, len.main=32)
            }
            R <- res[match(pw,res$pathway),,drop=FALSE]
            legend("topright", cex = 0.75, y.intersp=0.85, bty='n',
                   c( paste("NES=", round(R$NES[1],3)),
                     paste("padj=", round(R$padj[1],4)))
                   )
        }
    }
    
    getTopProfiles <- shiny::reactive({
        ## Get profiles of top-enriched contrasts (not all genes...)
        ##
        ##        
        df <- getConnectivityScores()
        ##if(is.null(df)) return(NULL)
        ##pw <- head(df$pathway,28)

        ii=1:100;sigdb="sigdb-archs4.h5"        

        ii <- connectivityScoreTable$rows_all()
        shiny::req(ii,input$cmap_sigdb)        
        ii <- head(ii,50)  ## 50??
        pw <- df$pathway[ii]
        
        sigdb <- input$cmap_sigdb
        shiny::req(sigdb)
        
        fc <- getCurrentContrast()$fc
        ngenes = 1000
        ngenes = 500
        var.genes <- head(names(sort(-abs(fc))),ngenes)
        var.genes <- unique(c(var.genes, sample(names(fc),ngenes)))  ## add some random        
        F <- getConnectivityMatrix(sigdb, select=pw, genes=var.genes)
        dbg("[getTopProfiles] *** READING H5 *** dim(F)=",dim(F))
        ## F <- t(t(F) * sign(score))        
        pw <- intersect(pw, colnames(F))
        F <- F[,pw,drop=FALSE]
        return(F)        
    })
    
    cmap_FCFCplots.RENDER <- shiny::reactive({        

        ngs <- inputData()
        alertDataLoaded(session,ngs)

        shiny::req(ngs, input$cmap_contrast)        
        ##ct <- input$cmap_contrast
        ##fc <- ngs$gx.meta$meta[[ct]]$meta.fx
        ##names(fc) <- rownames(ngs$gx.meta$meta[[ct]])
        res1 <- getCurrentContrast()
        fc <- res1$fc
        ct <- res1$name
        F <- getTopProfiles()
        if(NCOL(F)==0) return(NULL)
        F <- F[,1:min(ncol(F),10),drop=FALSE]
        
        if(input$fcfc_plottype=="scatter") {
            mfplots <- c(2,5)
            cmap_FCFCscatter(fc, F, mfplots, ylab=ct)
        } else {
            mfplots <- c(3,4)
            df <- getConnectivityScores()            
            cmap_FCFCenplot(fc, F, mfplots, ylab, df)
        }
        
    })
    
    cmap_FCFCplots.opts <- shiny::tagList(        
        shiny::radioButtons(ns("fcfc_plottype"),"Plot type:",c("scatter","enrichment"),
                     inline=TRUE)
    )
    
    cmap_FCFCplots_info = "<b>FC scatter plots.</b> Scatter plots of gene expression foldchange values between two contrasts. Foldchanges that are similar show high correlation, i.e. are close to the diagonal. You can switch to enrichment type plots in the plot settings."

    cmap_FCFCplots_caption = "<b>FC scatter plots.</b> Scatter plots of gene expression foldchange values between two contrasts. Foldchanges that are similar show high correlation, i.e. are close to the diagonal."
    
    shiny::callModule(
        plotModule,
        "cmap_FCFCplots", label = "a",
        func = cmap_FCFCplots.RENDER,
        func2 = cmap_FCFCplots.RENDER,
        options = cmap_FCFCplots.opts,
        title = "FC scatter plots",
        info.text = cmap_FCFCplots_info,
        pdf.height=4.5, pdf.width=10, 
        height = c(360,600), width=c("auto",1280),
        res = c(90,110),
        add.watermark = WATERMARK
    )

      
    ##================================================================================
    ## Cumulative FC barplot
    ##================================================================================

    cumulativeFCtable <- shiny::reactive({
        
        F <- getTopProfiles()
        F[is.na(F)] <- 0

        ## maximum 10??
        MAXF = 20
        ## F <- F[,1:min(MAXF,ncol(F)),drop=FALSE]  

        ## multiply with sign of rho
        df <- getConnectivityScores()
        rho1 <- df$rho[match(colnames(F),df$pathway)]
        F <- t(t(F) * sign(rho1))
        
        ## add current contrast
        cc <- getCurrentContrast()
        shiny::req(cc)
        fc <- cc$fc[rownames(F)]
        fc[is.na(fc)] <- 0
        F <- cbind(fc[rownames(F)], F)
        colnames(F)[1] <- "thisFC"
        colnames(F)[1] <- cc$name

        if(input$cumFCplot_absfc) {
            F <- abs(F)  
        }
        F <- F[order(-rowMeans(F**2)),,drop=FALSE]
        F
    })
    
    cumFCplot.RENDER <- shiny::reactive({
        
        F <- cumulativeFCtable()
        shiny::req(F)

        MAXF=10
        NGENES=64
        
        F <- F[,1:min(MAXF,ncol(F)),drop=FALSE]
        if(input$cumFCplot_order=="FC") {
            F <- F[order(-abs(F[,1])),]
            F1 <- head(F,NGENES)
            F1 <- F1[order(F1[,1]),,drop=FALSE]
        } else {
            F1 <- head(F,NGENES)
            F1 <- F1[order(rowMeans(F1)),,drop=FALSE]
        }
        
        par(mfrow=c(1,1), mar=c(7,3.5,0,0), mgp=c(2.4,1,0))
        maxfc <- max(abs(rowSums(F1,na.rm=TRUE)))
        ylim <- c(-1*(min(F1,na.rm=TRUE)<0),1.2)*maxfc

        col1 = grey.colors(ncol(F1),start=0.15)
        ##barplot(t(F1), las=3, cex.names=0.85, col=col1,
        ##        ylim=ylim, ylab="cumulative logFC")
        pgx.stackedBarplot(F1, cex.names=0.85, col=col1,
                           ## ylim=ylim,
                           ylab="cumulative logFC")        
        F1
    })

    cumFCplot.RENDER2 <- shiny::reactive({
        F1 <- cumFCplot.RENDER()
        col1 = grey.colors(ncol(F1),start=0.15)
        legend("topleft", legend = rev(colnames(F1)),
               fill = rev(col1), cex=0.72, y.intersp=0.85)
        
    })

    cumFCplot.opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns('cumFCplot_absfc'),'Absolute foldchange',FALSE),
               "Take the absolute foldchange for calculating the cumulative sum.",
               placement="right", options = list(container = "body")),
        withTooltip( shiny::radioButtons(ns('cumFCplot_order'),'Order:',
                             choiceValues = c("FC","cumFC"),
                             choiceNames = c("this FC","cumFC"),
                             selected="cumFC", inline=TRUE),
               "How to order the cumulative barplot.",
               placement="right", options = list(container = "body"))
        ## withTooltip( shiny::radioButtons(ns('cumFCplot_ntop'),'N-best:', c(5,10,20),
        ##                      selected=10, inline=TRUE),
        ##        "How many closest profiles to consider for calculating the cumulative sum.",
        ##        placement="right", options = list(container = "body"))
    )

    cumFCplot.caption = "<b>Meta-foldchange.</b> The barplot visualizes the cumulative foldchange between the top-10 most similar profiles."

    cumFCplot.info = "<b>Meta-foldchange.</b> The barplot visualizes the cumulative foldchange between the top-10 most similar profiles. Genes that are frequently shared with high foldchange will show a higher cumulative score. You can choose between signed or absolute foldchange in the options."
    
    shiny::callModule(
        plotModule,
        id = "cumFCplot",
        title="Cumulative foldchange", label="a",
        func = cumFCplot.RENDER,
        func2 = cumFCplot.RENDER2,
        csvFunc = cumulativeFCtable,
        info.text =  cumFCplot.info,
        ## caption =  cumFCplot.caption,
        options = cumFCplot.opts,
        pdf.height = 6, pdf.width = 9, 
        height = c(300, 600), width = c('auto',1300),
        res = c(72,90),
        add.watermark = WATERMARK
    )    
        
    ##================================================================================
    ## Cumulative enrichment barplot
    ##================================================================================
    
    cumEnrichmentTable <- shiny::reactive({

        shiny::req(input$cmap_sigdb)        
        if(!grepl(".h5$",input$cmap_sigdb)) return(NULL)
        
        df <- getConnectivityScores()
        if(is.null(df)) return(NULL)        
        ii <- connectivityScoreTable$rows_all()
        shiny::req(ii)
        
        sel <- head(df$pathway[ii],10)
        sigdb <- input$cmap_sigdb
        F <- getEnrichmentMatrix(sigdb, select=sel)
        if(is.null(F)) return(NULL)

        ## multiply with sign of enrichment
        rho1 <- df$rho[match(colnames(F),df$pathway)]
        F <- t(t(F) * sign(rho1))

        if(input$cumgsea_absfc) {
            F <- abs(F)  
        }
        F <- F[order(-rowMeans(F**2)),,drop=FALSE]

        ## add current contrast
        if(1) {
            ct <- getCurrentContrast()
            gx <- ct$gs[match(rownames(F),names(ct$gs))]
            names(gx) <- rownames(F)
            gx[is.na(gx)] <- 0
            F <- cbind(gx, F)
            colnames(F)[1] <- ct$name
        }
        
        F
    })
    
    cumEnrichmentPlot.RENDER <- shiny::reactive({
        
        ##
        F <- cumEnrichmentTable()
        if(is.null(F)) {
            frame()
            return(NULL)
        }

        NSETS=20
        if(input$cumgsea_order == "FC") {
            F <- F[order(-abs(F[,1])),]
            F <- head(F,NSETS)
            F <- F[order(F[,1]),,drop=FALSE]
        } else {
            F <- F[order(-rowMeans(F**2)),]
            F <- head(F,NSETS)
            F <- F[order(rowMeans(F)),,drop=FALSE]
        }
        
        rownames(F) <- gsub("H:HALLMARK_","",rownames(F))
        rownames(F) <- gsub("C2:KEGG_","",rownames(F))
        rownames(F) <- shortstring(rownames(F),72)
        maxfc <- max(abs(rowSums(F,na.rm=TRUE)))
        xlim <- c(-1*(min(F,na.rm=TRUE)<0),1.2)*maxfc
        
        par(mfrow=c(1,2), mar=c(4,1,0,0.5), mgp=c(2.4,1,0))
        frame()
        col1 = grey.colors(ncol(F),start=0.15)
        pgx.stackedBarplot(
            F, las=1, cex.names=0.8, col=col1, hz = TRUE,
            ## las=1, xlim=xlim,
            xlab="cumulative enrichment")

        fname <- sub("\\].*","]",colnames(F))
        
    })

    cumEnrichmentPlot.RENDER2 <- shiny::reactive({
        
        ##
        F <- cumEnrichmentTable()
        if(is.null(F)) {
            frame()
            return(NULL)
        }

        NSETS=40
        if(input$cumgsea_order == "FC") {
            F <- F[order(-abs(F[,1])),]
            F <- head(F,NSETS)
            F <- F[order(F[,1]),,drop=FALSE]
        } else {
            F <- F[order(-rowMeans(F**2)),]
            F <- head(F,NSETS)
            F <- F[order(rowMeans(F)),,drop=FALSE]
        }
        
        rownames(F) <- gsub("H:HALLMARK_","",rownames(F))
        rownames(F) <- gsub("C2:KEGG_","",rownames(F))
        rownames(F) <- shortstring(rownames(F),72)
        maxfc <- max(abs(rowSums(F,na.rm=TRUE)))
        xlim <- c(-1*(min(F,na.rm=TRUE)<0),1.2)*maxfc
        
        par(mfrow=c(1,2), mar=c(4.5,1,0.4,1), mgp=c(2.4,1,0))
        frame()
        col1 = grey.colors(ncol(F),start=0.15)
        pgx.stackedBarplot(
            F, las=1, cex.names=0.78, col=col1, hz = TRUE,
            ## las=1, xlim=xlim,
            xlab="cumulative enrichment")

        fname <- sub("\\].*","]",colnames(F))
        
    })
    
    
    cumEnrichmentPlot.opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns('cumgsea_absfc'),'Absolute foldchange',FALSE),
               "Take the absolute foldchange for calculating the cumulative sum.",
               placement="right", options = list(container = "body")),
        withTooltip( shiny::radioButtons(ns('cumgsea_order'),'Order:',
                             choiceValues = c("FC","cumFC"),
                             choiceNames = c("this FC","cumFC"),
                             selected="cumFC", inline=TRUE),
               "How to order the cumulative barplot.",
               placement="right", options = list(container = "body"))        
        ## withTooltip( shiny::radioButtons(ns('cumFCplot_ntop'),'N-best:', c(5,10,20),
        ##                      selected=10, inline=TRUE),
        ##        "How many closest profiles to consider for calculating the cumulative sum.",
        ##        placement="right", options = list(container = "body"))
    )
    
    cumEnrichmentPlot.caption = "<b>Meta-enrichment.</b> The barplot visualizes the cumulative enrichment of the top-10 most similar profiles."

    cumEnrichmentPlot.info = "<b>Meta-enrichment.</b> The barplot visualizes the cumulative enrichment of the top-10 most similar profiles. Gene sets that are frequently shared with high enrichment will show a higher cumulative scores. You can choose between signed or absolute enrichment in the options."
    
    shiny::callModule(
        plotModule,
        id = "cumEnrichmentPlot",
        title="Cumulative enrichment", label="b",
        func = cumEnrichmentPlot.RENDER,
        func2 = cumEnrichmentPlot.RENDER2,
        csvFunc = cumEnrichmentTable,
        info.text =  cumEnrichmentPlot.info,
        ## caption =  cumFCplot.caption,
        options = cumEnrichmentPlot.opts,
        pdf.height = 8, pdf.width = 12, 
        height = c(300,720), width = c('auto',1000),
        res = c(72,90),
        add.watermark = WATERMARK
    )    


    ##=============================================================================
    ## CONNECTIVITY MAP
    ##=============================================================================
           
    connectivityMap.RENDER <- shiny::reactive({
        
        ngs <- inputData()
        shiny::req(ngs, input$cmap_contrast)
       
        ## get positions
        res0 <- getThresholdedConnectivityScores()
        if(is.null(res0)) return(NULL)        
        res0 <- res0[,c("pathway","score","rho")]
        pos0 <- getConnectivityPositions()
        if(is.null(pos0)) return(NULL)
        if(is.null(res0)) return(NULL)
        
        cmap_grouped <- "grouped" %in% input$cmap_plotoptions
        if(cmap_grouped) {
            ##res0 <- collapseByDataset(res0)
            pset <- sub("].*","]",rownames(pos0))            
            pos0 <- apply(pos0,2,function(x) tapply(x,pset,mean)) 
        }
        if(is.null(res0) || nrow(res0)==0) return(NULL)
        if(is.null(pos0) || nrow(pos0)==0) return(NULL)
        
        ##res <- res[match(rownames(pos0),res$pathway),,drop=FALSE]
        pp <- intersect(rownames(pos0), rownames(res0))
        ##pp <- head(sample(pp),4000)
        length(pp)
        res <- res0[match(pp,rownames(res0)),,drop=FALSE]
        pos <- pos0[match(pp,rownames(pos0)),,drop=FALSE]

        ## draw high score last
        jj <- order(res$score)
        res <- res[jj,,drop=FALSE]
        pos <- pos[jj,,drop=FALSE]
        
        ## parse groups/dataset from name
        jj <- grep("\\]",rownames(pos))
        dset <- rep("[this_data]",nrow(pos))
        pw <- res$pathway  ## full name (may be modified by collapsing)
        dset[jj] <- sub("].*","]",pw[jj])
        ct.name <- sub(".*\\]","",pw)
        table(dset)
                
        ## make dataframe
        score=rho=NULL
        df <- data.frame(pos, dataset=dset, contrast=ct.name, name=res$pathway,
                         score=res$score, rho=res$rho, check.names=FALSE)
        dim(df)
        rownames(df) <- rownames(pos)
        
        do3d <- ncol(pos)==3
        mode = "markers"
        ann.text = rep(NA,nrow(df))
        ##showlabels <- input$cmap_showlabel
        showlabels <- "label" %in% input$cmap_plotoptions
        if(showlabels) ann.text = df$name


        tt.info <- paste(
            'Contrast:', df$contrast,
            '</br>Dataset:', df$dataset,
            '</br>Similarity score:', round(df$score,3),
            '</br>Correlation:', round(df$rho,3)
        )
        cex1 = c(1.0,0.8,0.6,0.4)[1 + 1*(nrow(pos)>30) + 1*(nrow(pos)>200) +
                                  1*(nrow(pos)>500)]
        this.data <- 1*grepl("this_data",dset)
        shapevar <- 1 + 1*this.data
        symbols = c("circle","x")
        
        bluered64 <- colorRampPalette(
            colors = c("navyblue","royalblue4","grey90","indianred3","firebrick4"))(64)
        greyred64 <- colorRampPalette(colors=c("grey85","grey70","indianred3","firebrick4"))(64)
        if("dark" %in% input$cmap_plotoptions) {
            dbg("*** DARK MODE ***\n")       
            greyred64 <- colorRampPalette(colors=c("grey15","grey30","indianred3","firebrick4"))(64)
            ##greyred64 <- colorRampPalette(colors=c("grey15","grey15","grey30","green4"))(64)
        }

        ## add transparency
        greyred64 <- sapply(1:64,function(i) paste0(greyred64[i],sprintf("%02X",0+4*i-1)))
        bluered64 <- sapply(1:64,function(i) paste0(bluered64[i],sprintf("%02X",abs(-259+8*i-1)) ))        
        colorby <- "score"
        colorby <- input$cmap_cmapcolorby
        sizevar <- (df$score)**2
        colorpal <- NULL
        marker.col <- NULL
        
        if(colorby == "dataset") {
            ## dataset
            colorvar <- dset  ## [GSE1234-xxx]
            marker.col <- list()
            colorpal <- rep(RColorBrewer::brewer.pal(8,"Set2"),99)

        } else if(colorby == "hallmark") {
            Y <- t(getEnrichmentMatrix(input$cmap_sigdb,nc=15))
            Y <- Y[match(rownames(df),rownames(Y)),]
            colorvar <- colnames(Y)[max.col(Y)]
            colorvar <- shortstring(colorvar,40) ## too long!!!
            marker.col <- list()
            colorpal <- rep(RColorBrewer::brewer.pal(8,"Set2"),99)

        } else if(colorby == "score") {
            dbg("*** SCORE ***\n")
            gamma1 <- input$cmap_scoregamma
            score1 <- sign(df$score) * (abs(df$score)/max(abs(df$score),na.rm=TRUE))**gamma1
            ##score1 <- abs(score1) ## absolute??
            colorvar <- score1
            colorpal <- greyred64
            cmin = 0
            if(min(score1)<0) {
                colorpal <- bluered64
                cmin = -1
            }
            marker.col <- list( colorscale = colorpal,
                               cmin = cmin, cmax=1, colorbar = list(title="score", len=0.33), 
                               opacity=0.966, reversescale = FALSE)
        } else {
            stop("[connectivityMap.RENDER] FATAL:: invalid colorby") ## should not come here
        }
        
        if(do3d ) {
            ## 3D plot
            sizeref <- 0.06 * max(1,nrow(df)/1000)**0.33
            sizeref
            if("large" %in% input$cmap_plotoptions) sizeref <- 0.25*sizeref
            
            plt <- plotly::plot_ly(
                df, x = df[,1], y = df[,2], z = df[,3], 
                mode = 'markers', type="scatter3d",
                ##mode = 'markers', type = "pointcloud",
                color = colorvar, colors = colorpal, 
                size = sizevar, sizes = c(5,35)*cex1,                
                marker = c(marker.col,
                           list( sizeref = sizeref,
                                line = list(color="grey50", width=0, opacity=0.5)
                                )),
                ## symbol = shapevar, symbols = symbols,
                text = tt.info
            )                        
            if(FALSE && showlabels && nrow(pos)<100) {
                ## slow...
                plt <- plt %>%
                    plotly::add_annotations(
                        x = pos[,1], y = pos[,2], z = pos[,3],
                        text = ann.text,
                        ##xref = "x", yref = "y",
                        showarrow = FALSE)
            }
            plt <- plt %>%
                plotly::layout(
                    xaxis = list(title=colnames(pos)[1]),
                    yaxis = list(title=colnames(pos)[2]),
                    zaxis = list(title=colnames(pos)[3]) ) 
            
        } else {
            ## 2D plot
            ##
            sizeref <- 0.08 * max(1,nrow(df)/1000)**0.33
            sizeref
            if("large" %in% input$cmap_plotoptions) sizeref <- 0.85*sizeref

            plt <- plotly::plot_ly(
                df, x = df[,1], y = df[,2],
                source="cmap2d", key=rownames(df),
                mode = 'markers', type = "scattergl",
                ## type = "pointcloud",
                color = colorvar, colors = colorpal, 
                marker = c( marker.col,
                           ##sizemin = 3, sizemax = 100, ## for pointcloud
                           list( sizeref = sizeref,
                                line = list(color="grey50", width=0, opacity=0.5)
                                )),
                ## symbol = shapevar, symbols=symbols,
                size = sizevar, sizes = c(5,30)*cex1,
                text = tt.info ) 
            
            if(showlabels && nrow(pos)<100) {
                ## this is slow if not careful. too many annotation labels slows down
                ##
                plt <- plt  %>%
                    plotly::add_annotations(
                        x = pos[,1], y = pos[,2],
                        text = ann.text,
                        ##xref = "x", yref = "y",
                        yanchor = "bottom",
                        yshift = 2,
                        textposition = 'top',
                        showarrow = FALSE)
            }
            plt <- plt %>%
                plotly::layout(
                    xaxis = list(title=colnames(pos)[1]),
                    yaxis = list(title=colnames(pos)[2]) ) 
            
        }

        ## scale range to cover 95%
        xrng <- quantile( pos[,1], probs=c(0.025,0.975))
        yrng <- quantile( pos[,2], probs=c(0.025,0.975))
        xrng <- xrng + 0.15*c(-1,1)*diff(xrng)
        yrng <- yrng + 0.15*c(-1,1)*diff(yrng)
        if(0 && ncol(pos)==2) {
            plt <- plt %>%
                plotly::layout(
                    xaxis = list(range = xrng),
                    yaxis = list(range = yrng) ) 
        }
        if(0 && ncol(pos)==3) {
            zrng <- quantile( pos[,3], probs=c(0.025,0.975))
            zrng <- zrng + 0.15*c(-1,1)*diff(zrng)
            plt <- plt %>%
                plotly::layout(
                    xaxis = list(range = xrng),
                    yaxis = list(range = yrng),
                    zaxis = list(range = zrng) ) 
        }
        
        plt <- plt %>%
            plotly::config( toImageButtonOptions =
                        list(format='svg', height=800, width=800, scale=1.1)) %>%
            plotly::event_register('plotly_selected')

        if("dark" %in% input$cmap_plotoptions) {
            plt <- darkmode(plt, dim=ncol(pos) )
        }
        
        return(plt)
    })

    connectivityMap.RENDER1 <- shiny::reactive({
        ## Hide colorbar
        ##
        plt <- connectivityMap.RENDER() 
        if(is.null(plt)) return(NULL)
        plt <- plt %>% plotly::layout(showlegend = FALSE) %>% plotly::hide_colorbar()
        return(plt)
    })
    
    connectivityMap.opts = shiny::tagList(
        withTooltip(shiny::radioButtons(
            ##"Choose the plot layout: t-SNE, PCA or UMAP",
            ns('cmap_layout'),"Layout:",c("pca","tsne","volcano"), selected="tsne", inline=TRUE),
            "Choose the plot layout: t-SNE, PCA, or volcano-type",
            placement="right", options = list(container = "body")),
        withTooltip(shiny::sliderInput(ns('cmap_scorethreshold'),"Score threshold:", 0, 1, 0, step=0.01),
               "Threshold the points by minimum score",
               placement="right", options = list(container = "body")),
        withTooltip(shiny::radioButtons(
            ns('cmap_cmapcolorby'),"Color by:", c("score","dataset","hallmark"),
            inline=TRUE), "Color the points by score, dataset or hallmark",
            placement="right", options = list(container = "body")),
        withTooltip(shiny::sliderInput(ns('cmap_scoregamma'),"Color gamma:", 0.1, 2, 0.5, step=0.1),
               "Gamma for color adjustments",
               placement="right", options = list(container = "body")),
        withTooltip( shiny::checkboxGroupInput(
            ns('cmap_plotoptions'),"Other options:",
            choiceValues = c("label","grouped","3D","dark","large"),
            choiceNames = c("show label","group by dataset","3D plot", "dark mode",
                            "larger points"),
            selected = c("label","3D") ),
            "Show labels, group by dataset, show 3D plot, dark mode.",
            placement="top", options = list(container = "body"))
    )

    connectivityMap_info = "<b>The Connectivity Map</b> shows the similarity of the contrasts profiles as a t-SNE plot. Contrasts that are similar will be clustered close together, contrasts that are different are placed farther away."

    
    shiny::callModule(
        plotModule,
        "connectivityMap", label = "a",
        func = connectivityMap.RENDER1, plotlib="plotly",
        func2 = connectivityMap.RENDER, 
        options = connectivityMap.opts,
        download.fmt = c("pdf","png","html"),  ## PDF/PNG does not work???
        title = "Connectivity map",
        info.text = connectivityMap_info,
        pdf.width=8, pdf.height=8,
        height = c(fullH-100,750), width = c('auto',1000),
        res=90,
        add.watermark = WATERMARK
    )

    ##-------------------------------------------------------------------------------
    ## Leading-edge graph 
    ##-------------------------------------------------------------------------------    
    
    getLeadingEdgeGraph <- shiny::reactive({
            
        dbg("[getLeadingEdgeGraph] reacted")

        ##df <- (ngs$connectivity[[2]][[1]])
        df <- getConnectivityScores()
        ##df <- getThresholdedConnectivityScores()        
        if(is.null(df)) return(NULL)        
        df$score[is.na(df$score)] <- 0
        df <- df[which(df$score>0),]
        dim(df)
        
        ## always only top-N from selection
        ntop <- as.integer(input$LEgraph_ntop)
        ##ii <- order(-df$score)
        ii <- connectivityScoreTable$rows_all()
        shiny::req(ii)
        ii <- head(order(-abs(df$score)),25)
        ii <- head(ii,ntop)
        df <- df[ii,,drop=FALSE]

        le.genes <- sort(unique(unlist(df$leadingEdge)))        
        A <- 1*sapply( df$leadingEdge, function(g) le.genes %in% g)
        rownames(A) <- le.genes
        A <- head( A[order(-rowMeans(A)),,drop=FALSE], 100)        
        head(A)
        ##adjM <- pmax( cor(t(A)), 0)  ## only positive correlation???
        adjM <- (A %*% t(A))
        adjM <- adjM / max(adjM,na.rm=TRUE)

        gr <- igraph::graph_from_adjacency_matrix(
            adjM, mode="undirected", weighted=TRUE, diag=FALSE)

        ## set graph threshold to some sensible value [0,1]
        wt0 <- tail(sort(abs(igraph::E(gr)$weight)),150)[1] ## about 150 edges
        ##dbg("[getLeadingEdgeGraph] setting edge threshold to wt=",wt0)
        ##updateSliderInput(session, "LEgraph_threshold", value = 0.99*wt0)
        shiny::updateSliderInput(session, "LEgraph_threshold", value = 0)
        
        return(gr)
    })

    leadingEdgeGraph.VISNETWORK <- shiny::reactive({

        dbg("[leadingEdgeGraph.VISNETWORK] reacted")
        
        gr <- getLeadingEdgeGraph()
        if(is.null(gr)) return(NULL)
        
        max(abs(igraph::E(gr)$weight))
        minwt <- 0.5
        minwt <- input$LEgraph_threshold
        minwt <- min(c(minwt, 0.99*max(abs(igraph::E(gr)$weight),na.rm=TRUE)))
        gr <- igraph::subgraph.edges(gr, which(abs(igraph::E(gr)$weight) >= minwt) )
        if(length(igraph::V(gr))==0) return(NULL)

        fc=cumFC=NULL
        fc <- getCurrentContrast()$fc        
        fc <- fc[igraph::V(gr)$name]
        cumFC <- cumulativeFCtable()
       
        cumFC <- cumFC[igraph::V(gr)$name,]
        fontsize = 22
        fc <- fc / max(abs(fc))
        
        sizevar <- input$LEgraph_sizevar
        vsize = 15
        if(sizevar=="centrality") {
            vsize <- log(1 + igraph::betweenness(gr))
        } else if(sizevar=="cumFC") {
            fc1 <- rowMeans(cumFC)
            vsize <- abs(fc1[match(igraph::V(gr)$name,names(fc1))])**2
        } else if(sizevar=="FC") {
            vsize <- abs(fc[match(igraph::V(gr)$name,names(fc))])**2
        } else {
            vsize <- 1
        }
        vsize <- 3 + 12 * (abs(vsize) / max(abs(vsize),na.rm=TRUE))**0.5

        bluered.pal <- colorRampPalette(
            colors = c("royalblue4","royalblue2","grey90","indianred3","firebrick4")
        )       
        vcolor <- bluered.pal(65)[ 33 + round(32*fc) ]
        vcolor <- paste0(vcolor,"AA") ## add transparency
        
        ## defaults graph parameters
        gene <- igraph::V(gr)$name
        igraph::V(gr)$label <- igraph::V(gr)$name
        igraph::V(gr)$title <- paste0("<b>",gene,"</b><br>",GENE.TITLE[toupper(gene)])
        igraph::V(gr)$size  <- vsize      ## rather small
        igraph::V(gr)$color  <- vcolor

        ew = abs(igraph::E(gr)$weight)
        igraph::E(gr)$width <- 1.5 * (0.2 + 10*(ew/max(ew,na.rm=TRUE))**2)
        igraph::E(gr)$color <- "#DDD"  ## lightgrey
        
        visdata <- visNetwork::toVisNetworkData(gr, idToLabel=FALSE)
        
        ## ------------------ plot using visNetwork (zoomable) -----------------
        graph <- visNetwork::visNetwork(
            nodes = visdata$nodes,
            edges = visdata$edges) %>%
            visNetwork::visNodes(font = list(size = fontsize))  %>%
            visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            visNetwork::visIgraphLayout(layout="layout_nicely")
        graph
        
    })
    
    leadingEdgeGraph.opts <- shiny::tagList(
        withTooltip( shiny::sliderInput(ns("LEgraph_threshold"),"edge threshold:",0, 1, 0, 0.01),
               "Threshold value for edges."),
        withTooltip( shiny::radioButtons(ns("LEgraph_ntop"),"N-neighbours:",c(5,10,25,100),selected=10,inline=TRUE),
               "Number of simlar experiments to consider."),
        withTooltip( shiny::radioButtons(ns("LEgraph_sizevar"),"Size:",c("FC","cumFC","centrality"),
                             selected="cumFC", inline=TRUE),
               "Parameter for node size.")        
    )
    
    leadingEdgeGraph_info = "<b>Leading-edge graph.</b> Network of shared leading-edge genes between top-N most similar signatures. The edge width corresponds to the number of signatures that share that pair of genes in their top differentially expressed genes. In the plot options you can set the threshold of the edges."
    
    shiny::callModule(
        plotModule,
        "leadingEdgeGraph", label = "a",
        func = leadingEdgeGraph.VISNETWORK,
        plotlib="visnetwork",
        options = leadingEdgeGraph.opts,
        title = "Leading-edge graph",
        info.text = leadingEdgeGraph_info,
        pdf.height=8, pdf.width=8, 
        height = c(720,720), width=c("auto",1300),
        res = c(90,100),
        add.watermark = WATERMARK
    )

    ##-------------------------------------------------------------------------------
    ## Enrichment graph 
    ##-------------------------------------------------------------------------------    
    
    getEnrichmentGraph <- shiny::reactive({
        
        dbg("[getEnrichmentGraph] reacted")
                
        ## get enrichment scores
        F <- cumEnrichmentTable()

        if(input$enrichGraph_oddweighting) {
            gr2 <- getLeadingEdgeGraph()
            le.genes <- igraph::V(gr2)$name
            ##gsets <- GSETS[rownames(F)]
            gsets <- getGSETS(rownames(F))
            ##gsets <- lapply(gsets, function(x) intersect(x,le.genes))
            gsets <- gsets[sapply(gsets,length)>=5]
            bg <- unique(unlist(gsets))
            ft <- gset.fisher(le.genes, gsets, fdr=1.0,
                              min.genes=3, max.genes=99999,
                              background=bg)
            ft <- ft[match(rownames(F),rownames(ft)),]            
            F <- F * log(1+ft$odd.ratio)            
        }

        ## get top average-enriched
        ntop <- as.integer(input$enrichGraph_ntop)
        le.sets <- apply(F, 2, function(x) head(order(-abs(x)), ntop))
        le.sets <- apply(le.sets, 2, function(i) rownames(F)[i])
        A <- 1*apply( le.sets, 2, function(g) rownames(F) %in% g)
        rownames(A) <- rownames(F)

        ## create graph
        A <- head( A[order(-rowMeans(A)),,drop=FALSE], 100)        
        head(A)
        ##adjM <- pmax( cor(t(A)), 0)  ## only positive correlation???
        adjM <- (A %*% t(A))
        adjM <- adjM / max(adjM,na.rm=TRUE)

        gr <- igraph::graph_from_adjacency_matrix(
            adjM, mode="undirected", weighted=TRUE, diag=FALSE)

        ## set graph threshold to some sensible value [0,1]
        wt0 <- tail(sort(abs(igraph::E(gr)$weight)),100)[1] ## about 80 edges
        ##dbg("[getEnrichmentGraph] setting edge threshold to wt=",wt0)
        ##updateSliderInput(session, "enrichGraph_threshold", value = 0.99*wt0)
        shiny::updateSliderInput(session, "enrichGraph_threshold", value = 0)
        
        return(gr)
    })
    
    enrichmentGraph.VISNETWORK <- shiny::reactive({




        ## return(NULL)

        gr <- getEnrichmentGraph()
        if(is.null(gr)) return(NULL)
                
        gr2 <- getLeadingEdgeGraph()
        pw <- igraph::V(gr)$name
        le.genes <- igraph::V(gr2)$name

        ##gsets <- GSETS[pw]
        gsets <- getGSETS(pw)
        pw.genes <- sapply( gsets, function(gs) intersect(gs,le.genes))
        pw.genes <- sapply(pw.genes, paste, collapse=" ")
        
        max(abs(igraph::E(gr)$weight))
        minwt <- 0.5
        minwt <- input$enrichGraph_threshold
        minwt <- min(c(minwt, 0.99*max(abs(igraph::E(gr)$weight),na.rm=TRUE)))
        gr <- igraph::subgraph.edges(gr, which(abs(igraph::E(gr)$weight) >= minwt) )
        if(length(igraph::V(gr))==0) return(NULL)

        fc=cumFC=NULL
        cumFC <- cumEnrichmentTable()
        cumFC <- cumFC[igraph::V(gr)$name,]
        fc <- cumFC[,1]
        fontsize = 18
        fc <- fc / max(abs(fc))
        
        sizevar <- input$enrichGraph_sizevar
        vsize = 15
        if(sizevar=="centrality") {
            vsize <- log(1 + igraph::betweenness(gr))
        } else if(sizevar=="cumFC") {
            fc1 <- rowMeans(cumFC)
            vsize <- abs(fc1[match(igraph::V(gr)$name,names(fc1))])**2
        } else if(sizevar=="FC") {
            vsize <- abs(fc[match(igraph::V(gr)$name,names(fc))])**2
        } else {
            vsize <- 1
        }
        vsize <- 3 + 12 * (abs(vsize) / max(abs(vsize),na.rm=TRUE))**0.5

        bluered.pal <- colorRampPalette(
            ##colors = c("navyblue","royalblue4","grey90","indianred3","firebrick4")
            colors = c("royalblue4","royalblue2","grey90","indianred3","firebrick4")
        )       
        vcolor <- bluered.pal(65)[ 33 + round(32*fc) ]
        vcolor <- paste0(vcolor,"AA") ## add transparency
        
        ## defaults graph parameters
        vname <- sub("H:HALLMARK_|C2:KEGG_","",igraph::V(gr)$name)
        igraph::V(gr)$label <- vname
        igraph::V(gr)$title <- paste0("<b>",vname,"</b><br>",pw.genes)
        igraph::V(gr)$size  <- vsize      ## rather small
        ##V(gr)$color  <- c("skyblue","salmon")[1 + 1*(sign(fc)>0)]       ## rather small
        igraph::V(gr)$color  <- vcolor

        ew = abs(igraph::E(gr)$weight)
        igraph::E(gr)$width <- 1.5 * (0.2 + 10*(ew/max(ew,na.rm=TRUE))**2)
        igraph::E(gr)$color <- "#DDD"  ## lightgrey

        visdata <- visNetwork::toVisNetworkData(gr, idToLabel=FALSE)
        
        ## ------------------ plot using visNetwork (zoomable) -----------------
        graph <- visNetwork::visNetwork(
            nodes = visdata$nodes,
            edges = visdata$edges) %>%
            visNetwork::visNodes(font = list(size = fontsize))  %>%
            ## visNetwork::visEdges(hidden=FALSE, width=2, color=list(opacity=0.9))  %>%
            visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE)) %>%
            ## visNetwork::visPhysics(enabled=TRUE)  %>%
            ## visNetwork::visInteraction(hideEdgesOnDrag = TRUE) %>%
            ## visNetwork::visIgraphLayout(layout="layout.norm", layoutMatrix=pos)
            visNetwork::visIgraphLayout(layout="layout_nicely")
        graph
        
    })
    
    enrichmentGraph.opts <- shiny::tagList(
        withTooltip( shiny::sliderInput(ns("enrichGraph_threshold"),"edge threshold:",0, 1, 0, 0.01),
               "Threshold value for edges."),
        withTooltip( shiny::radioButtons(ns("enrichGraph_ntop"),"N-neighbours:",c(5,10,25,100),
                             selected=10, inline=TRUE),
               "Number of simlar experiments to consider."),
        withTooltip( shiny::checkboxInput(ns("enrichGraph_oddweighting"),"Odd ratio weighting",FALSE),
               "Odds ratio weighting."),
        withTooltip( shiny::radioButtons(ns("enrichGraph_sizevar"),"Size:",c("FC","cumFC","centrality"),
                             selected="cumFC", inline=TRUE),
               "Parameter for node size.")        
    )
    
    enrichmentGraph_info = "<b>Enrichment graph.</b> Network of shared enriched genesets between top-N most similar signatures. The edge width corresponds to the number of signatures that share that pair of genesets in their top enriched genesets. In the plot options you can set the threshold the edges."
    
    shiny::callModule(
        plotModule,
        "enrichmentGraph", label = "b",
        func = enrichmentGraph.VISNETWORK, plotlib="visnetwork",
        options = enrichmentGraph.opts,
        title = "Enrichment graph",
        info.text = enrichmentGraph_info,
        pdf.height=8, pdf.width=8, 
        height = c(720,720), width=c("auto",1200),
        res = c(90,100),
        add.watermark = WATERMARK
    )
    

    ##======================================================================    
    ## Pairs
    ##======================================================================

    ##----------------------------------------------------------------------
    ## Scatterplot matrix in plotly
    ##
    ## From: https://plot.ly/r/splom/
    ##----------------------------------------------------------------------

    cmapPairsPlot.PLOT <- shiny::reactive({



        
        ##res = pgx.getMetaFoldChangeMatrix(ngs, what="meta")
        ##fc1 <- getCurrentContrast()$fc

        ngs <- inputData()
        shiny::req(ngs, input$cmap_contrast)        
        ##ct <- input$cmap_contrast
        sigdb = "sigdb-archs4.h5"
        sigdb <- input$cmap_sigdb

        all.ct <- getConnectivityContrasts(sigdb)
        ct1 <- all.ct[1]
        fc1 <- ngs$gx.meta$meta[[1]]$meta.fx
        names(fc1) <- rownames(ngs$gx.meta$meta[[1]])
        
        fc1 <- getCurrentContrast()$fc
        ct1 <- getCurrentContrast()$name

        sigdb <- input$cmap_sigdb
        shiny::req(sigdb)        
        ct2 = all.ct[1]
        ##sel.row <- connectivityScoreTable2$rows_selected()
        sel.row <- connectivityScoreTable$rows_selected()
        shiny::req(sel.row)        
        df <- getConnectivityScores()
        df <- df[abs(df$score)>0,,drop=FALSE]
        ct2 <- rownames(df)[sel.row]
        ##ct2 <- input$cmap_splom_comparewith
        fc2 <- getConnectivityMatrix(sigdb, select=ct2)[,1]
        
        ## match with selection filter        
        gg <- unique(c(names(fc1),names(fc2))) ## union or intersection??
        fc1 <- fc1[match(gg,names(fc1))]
        fc2 <- fc2[match(gg,names(fc2))]
        df <- data.frame(fc2, fc1)
        rownames(df) <- gg
        colnames(df) <- c(ct2, ct1)
        df[is.na(df)] <- 0  ## missing as zero???
        df <- df[order(-rowMeans(abs(df**2),na.rm=TRUE)),]
        dim(df)

        ## Number of selected genes
        sel.genes = grep("^CD",rownames(df),value=TRUE)
        ##sel.genes = head(rownames(df),ntop)  ## high-light top100
        logfc <- as.numeric(input$cmap_logFC)
        sel.genes <- rownames(df)[ rowSums(abs(df) > logfc)>=1 ]  ## minimum FC
        head(sel.genes)

        not.sel <- setdiff(rownames(df),sel.genes)
        if(length(not.sel)>0) {
            nr <- 1000
            ii <- c(sel.genes, sample(not.sel,nr,replace=TRUE))
            df <- df[unique(ii),,drop=FALSE]
        }
        
        ## remove NA genes from selection
        if(1) {
            na.fc <- rownames(df)[rowSums(is.na(df))>0] ## probably was missing
            na.zero <- rownames(df)[rowSums(df==0)>0] ## probably was missing
            sel.genes <- setdiff(sel.genes, c(na.fc,na.zero))
        }
        
        ## resort selection so that selected genes are drawn last to avoid
        ## covering them up.
        is.sel = (rownames(df) %in% sel.genes)
        highlight = TRUE
        ##highlight = input$cmap_splom_highlight
        if(highlight) {
            df.color = c("#00000033","#0066FF")[1 + is.sel]
            df.color = c("#AAAAAA55","#1e60BB88")[1 + is.sel]
            ##df.color = sample(c("#AAAAAA66","#1e60BBCC"),nrow(df),replace=TRUE)
        } else {            
            df.color = rep("#77777788",nrow(df))
            df.color = rep("#1e60BB88",nrow(df))
        }
        
        ## Labels for top 50 
        label.text0 = head(rownames(df)[which(is.sel)],50)
        label.text <- shortstring(label.text0,30)
        if(sum(is.na(label.text))) label.text[is.na(label.text)] <- ""
        
        ## reorder so the selected genes don't get overlapped
        jj <- order(is.sel)
        df <- df[jj,]
        df.color <- df.color[jj]
        sel1 <- match(label.text0, rownames(df))  ## index for labeled
        
        ## Tooltip text for all 
        gg <- rownames(df)
        tt <- paste0("<b>",gg,"</b> ", ngs$genes[gg,"gene_title"])
        tt <- gsub("_"," ",tt)
        tt <- sapply(tt,breakstring2,50,brk="<br>")

        ## plotly
        ##   
        axis = list( showline=TRUE,
                    zeroline=TRUE,
                    gridcolor='#dddf',
                    ticklen=4
                    )
        
        if(ncol(df)>=3) {
            dimensions = lapply(colnames(df), function(a) list(label=a, values = df[,a]))

            ## compute correlations 
            rho = cor(df, use="pairwise")
            rho.text = paste("r=",as.vector(round(rho,digits=3)))
            n = ncol(df)

            ## annotation positions (approximated by eye...)
            xann = 1.02*(as.vector(mapply(rep,seq(0,0.98,1/n),n)) + 0.05*1/n)
            ##xann = as.vector(mapply(rep,seq(0,1,1/(n-1)),n))
            yann = 1.08*(as.vector(rep(seq(1,0.02,-1/n),n)) - 0.15*1/n - 0.04)
            ##yann = as.vector(rep(seq(1,0.0,-1/(n-1)),n))
            
            p <- plotly::plot_ly(df, source="cmapSPLOM", key=rownames(df) ) %>%
                plotly::add_trace(
                    type = 'splom',
                    dimensions = dimensions,
                    text = tt,
                    hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
                    marker = list(
                        color = df.color,
                        ## colorscale = pl_colorscale,
                        size = 5,
                        line = list(
                            width = 0.3,
                            ##color = 'rgb(230,230,230)'
                            color = 'rgb(0,0,0)'
                        )
                    )
                ) %>%
                plotly::add_annotations(
                    x = xann,
                    y = yann,
                    text = rho.text,
                    font = list(size=11),
                    xanchor = "left",
                    align = 'left',
                    showarrow = FALSE,
                    xref = 'paper',
                    yref = 'paper',
                    borderpad = 3, 
                    bordercolor = 'black',
                    borderwidth = 0.6) %>%
                plotly::layout(
                    ##title= 'Scatterplot matrix',
                    hovermode='closest',
                    dragmode= 'select',
                    ##annotations = annot,
                    ## plot_bgcolor='rgba(240,240,240, 0.95)',
                    ## template = "plotly_dark",
                    xaxis = c( domain=NULL, axis),
                    yaxis = c( domain=NULL, axis),
                    xaxis2=axis, xaxis3=axis, xaxis4=axis, xaxis5=axis, xaxis6=axis, xaxis7=axis,
                    yaxis2=axis, yaxis3=axis, yaxis4=axis, yaxis5=axis, yaxis6=axis, yaxis7=axis
                )
            ## %>% plotly::style(diagonal = list(visible = F))

        } else {

            rho = cor(df[,1], df[,2], use="pairwise")
            rho
            annot.rho <- list(
                text = paste("r=",round(rho,4)),
                font = list(size=14),
                align = 'left',
                showarrow = FALSE,
                xref = 'paper',
                yref = 'paper',
                x = 0.03,
                y = 0.97,
                borderpad = 8, 
                bordercolor = 'black',
                borderwidth = 0.6)

            p <- plotly::plot_ly( data=df[,1:2], x = df[,1], y= df[,2],
                         type = 'scattergl', mode = 'markers',
                         source='cmapSPLOM', key = rownames(df),
                         text = tt,
                         hovertemplate = paste0("<br>%{text}<br>x: %{x}<br>y: %{y}<extra></extra>"),
                         marker = list(
                             color = df.color,
                             size = 8,
                             line = list(
                                 width = 0.3,
                                 ##color = 'rgb(230,230,230)'
                                 color = 'rgb(0,0,0)'
                             ))
                         )
            if(length(sel1)>0) {
                p <- p %>%
                    plotly::add_annotations(
                        x = df[sel1,1],
                        y = df[sel1,2],
                        text = as.character(label.text),
                        xanchor = 'center',
                        yanchor = 'top',
                        font = list(size=14),
                        xref = "x", yref = "y",
                        showarrow = FALSE,
                        ax = 20, ay = -40
                    )
            }
            
            p <- p %>%
                plotly::layout(
                    ## title= 'Scatterplot',
                    annotations = annot.rho,
                    hovermode = 'closest',
                    ## dragmode= 'select',  ## default dragging mode
                    ## plot_bgcolor='rgba(240,240,240, 0.95)',
                    ## template = "plotly_dark",
                    xaxis = c(title = paste(colnames(df)[1],"   (logFC)"), axis),
                    yaxis = c(title = paste(colnames(df)[2],"   (logFC)"), axis)
                ) 
        }

        p <- p %>%
            ## plotly::config(displayModeBar = FALSE) %>% ## disable buttons
            plotly::config( toImageButtonOptions = list(format='svg', height=800, width=800, scale=1.1)) %>%
            plotly::event_register('plotly_selected')

        dbg("cmapPairsPlot:: done\n")
        p    
    })
    
    cmapPairsPlot.opts = shiny::tagList(
        withTooltip( shiny::selectInput(ns("cmap_logFC"),"logFC threshold:",c(0,0.5,1,2,3,4),selected=1), "Threshold for (log) foldchange to highlight in plot.",  placement="right", options = list(container = "body"))
    )

    cmapPairsPlot_info = "The <strong>Pairs</strong> panel provides pairwise scatterplots of differential expression profiles for the selected contrasts. The main purpose of this panel is to identify similarity or dissimilarity between selected contrasts. The pairs plot is interactive and shows information of each gene with a mouse hover-over."

    cmapPairsPlot_caption = "<b>Pairs plot.</b> Pairwise scatterplots of two differential expression profiles for selected contrasts. Similar profiles will show high correlation with points close to the diagonal."

    shiny::callModule(
        plotModule,
        id = "cmapPairsPlot", label="c",
        func = cmapPairsPlot.PLOT,
        plotlib="plotly",
        title = "Scatterplot matrix (pairs)",
        options = cmapPairsPlot.opts,
        download.fmt = c("html"),
        pdf.width=8, pdf.height=8,
        height = c(fullH-80,700), width=c('auto',1000), res = 95,
        info.text = cmapPairsPlot_info,
        add.watermark = WATERMARK
    )

    connectivityScoreTable2.RENDER <- shiny::reactive({
        
        df <- getConnectivityScores()
        if(is.null(df)) return(NULL)

        kk <- c("pathway","score","rho","NES","padj","size","leadingEdge")
        kk <- c("score","pathway","rho","NES","padj")
        kk <- intersect(kk,colnames(df))
        ##df <- df[,kk,with=FALSE]
        df <- df[,kk]
        df <- df[abs(df$score)>0,,drop=FALSE]
        
        ##colnames(df) <- sub("padj","NES.q",colnames(df))
        ##df$leadingEdge <- shortstring(sapply(df$leadingEdge,paste,collapse=","),40)
        df$pathway <- shortstring(df$pathway,110)
        df$leadingEdge <- NULL
        
        dbg("[connectivityScoreTable.RENDER] dim(df)=",dim(df))        
        
        colnames(df) <- sub("pathway","dataset/contrast",colnames(df))
        score.col <- which(colnames(df) == "score")
        numcols <- c("score","pval","padj","NES.q","ES","NES","rho","R2")
        numcols <- intersect(numcols, colnames(df))
        
        DT::datatable( df,
                      rownames = FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip',
                          pageLength = 99999,
                          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scrollX = TRUE,
                          scrollY = "100vh",
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      )  %>%
            DT::formatSignif(numcols,3) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')  %>%
                DT::formatStyle( "score",
                                ##background = DT::styleColorBar( score.col, 'lightblue'),
                                background = color_from_middle(
                                    df[,"score"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    })
    
    connectivityScoreTable2_opts = shiny::tagList(
        shiny::selectInput(ns("connectivityScoreTable2_qsig"),"threshold (padj)",
                    c(0.01,0.05,0.2,1),selected=1)
    )
    
    connectivityScoreTable2 <- shiny::callModule(
        tableModule,
        id = "connectivityScoreTable2", label = "b",
        func = connectivityScoreTable2.RENDER, 
        info.text = connectivityScoreTable_info,
        options = connectivityScoreTable2_opts,
        info.width = "150px",
        title = "Similarity scores",
        height = c(660,700), width = c('auto',1280)
    )
    
    ##=============================================================================
    ## CONNECTIVITY HEATMAP
    ##=============================================================================

    connectivityHeatmap.RENDER <- shiny::reactive({
        ##
        F <- cumulativeFCtable()
        shiny::req(F)
        F <- F[,1:min(NCOL(F),25),drop=FALSE]
        if(input$cumFCplot_order=="FC") {
            F <- F[order(-abs(F[,1])),]
        }
        F1 <- head(F,80)
        par(mfrow=c(1,1), mar=c(0,0,0,0))
        ##gx.heatmap(t(F1), keysize=0.85, mar=c(6,30))
        gx.splitmap(t(F1), split=1,
                    ##cluster_columns = FALSE,
                    cluster_columns = TRUE,
                    cluster_rows = TRUE,
                    rowlab.maxlen = 80,                    
                    ## zsym = TRUE,
                    symm.scale = TRUE,
                    mar = c(15,0,0,60), 
                    key.offset = c(0.90, 0.2),
                    cexRow=0.9, cexCol=0.75)

    })
    
    connectivityHeatmap.opts = shiny::tagList(
    )
    
    connectivityHeatmap_info = "<b>The Connectivity Heatmap</b> shows the most similar profiles as a heatmap. Contrasts that are similar will be clustered close together."

    connectivityHeatmap_caption = "<b>Connectivity Heatmap.</b> Similarity of the contrasts profiles as a heatmap. Contrasts that are similar will be clustered close together."
    
    shiny::callModule(
        plotModule,
        "connectivityHeatmap", label = "c",
        func = connectivityHeatmap.RENDER,
        func2 = connectivityHeatmap.RENDER,
        options = connectivityHeatmap.opts,
        title = "Connectivity Heatmap",
        info.text = connectivityHeatmap_info,
        pdf.width=14, pdf.height=5.5,
        height = c(480,550), width = c('auto',1400),
        res = c(90,90),
        add.watermark = WATERMARK
     )
  })      
} ## end-of-Board 
