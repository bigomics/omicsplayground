##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

DrugConnectivityBoard <- function(id, inputData)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns ## NAMESPACE
    
    fullH = 750
    rowH = 660  ## row height of panel
    tabH = 200  ## row height of panel
    tabH = '60vh'  ## row height of panel  
    
    dsea_infotext = paste("<b>This module performs drug enrichment analysis</b> to see if certain drug activity or drug
        sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate
        biological functions through mechanism-of-action (MOA) and known drug molecular targets.<br><br>
        In the <a href='https://portals.broadinstitute.org/cmap/'>Drug Connectivity Map</a> panel,
        you can correlate your signature with known drug profiles from the L1000 database.
        An activation-heatmap compares drug activation profiles across multiple contrasts.
        This facilitates to quickly see and detect the similarities between contrasts for certain drugs.<br><br><br><br>
        <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6'
        frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture'
        allowfullscreen></iframe></center>"
    )

    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        ct <- names(ngs$drugs)
        shiny::updateSelectInput(session, "dsea_method", choices=ct)
    })
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$dsea_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Drug Connectivity Analysis Board</strong>"),
            shiny::HTML(dsea_infotext),
            easyClose = TRUE, size="l" )
        )
    })

    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        ct <- colnames(ngs$model.parameters$contr.matrix)
        shiny::updateSelectInput(session, "dsea_contrast", choices=sort(ct) )
    })
    
    ##================================================================================
    ## Reactive functions
    ##================================================================================

    getActiveDSEA <- shiny::reactive({
        
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        shiny::req(ngs)        
        shiny::req(input$dsea_contrast, input$dsea_method)
        
        dbg("[getActiveDSEA] reacted")
        
        names(ngs$drugs)
        dmethod="CTRPv2/sensitivity"
        dmethod="L1000/activityXL"
        dmethod="L1000/gene"
        colnames(ngs$contrasts)
        contr="treatment:Gefitinib_vs_CT" 
        
        contr = input$dsea_contrast
        if(is.null(contr)) return(NULL)
        
        dmethod <- input$dsea_method
        if(is.null(dmethod)) return(NULL)

        dr    <- ngs$drugs[[dmethod]]

        dbg("[getActiveDSEA] dim(dr$X) = ",dim(dr$X))
        dbg("[getActiveDSEA] dim(dr$P) = ",dim(dr$P))
        dbg("[getActiveDSEA] dim(dr$Q) = ",dim(dr$Q))
        dbg("[getActiveDSEA] dim(dr$stats) = ",dim(dr$stats))                
        
        nes   <- round(dr$X[,contr],4)
        pv    <- round(dr$P[,contr],4)
        qv    <- round(dr$Q[,contr],4)
        drug  <- rownames(dr$X)
        stats <- dr$stats[,contr]
        annot <- dr$annot
        nes[is.na(nes)] <- 0
        qv[is.na(qv)] <- 1
        pv[is.na(pv)] <- 1
        
        ## !!!SHOULD MAYBE BE DONE IN PREPROCESSING???
        if(is.null(annot)) {
            message("[getActiveDSEA] WARNING:: missing drug annotation in PGX file!")
            annot <- read.csv(file.path(FILESX,"sig/L1000_repurposing_drugs.txt"),
                              sep="\t", comment.char="#")
            rownames(annot) <- annot$pert_iname
        }

        ## compile results matrix
        jj <- match(toupper(drug), toupper(rownames(annot)))
        annot <- annot[jj,c("moa","target")]        
        dt <- data.frame( drug=drug, NES=nes, pval=pv, padj=qv, annot)
        ##dt <- dt[order(-abs(dt$NES)),]
        ##dt <- dt[order(dt$pval,-abs(dt$NES)),]
        dt <- dt[order(-dt$NES),]
        
        if(length(input$dseatable_filter)==0) {
            ## sometimes UI is not ready
            return(NULL)
        }
        
        if(input$dseatable_filter) {
            sel <- which(dt$moa!='' | dt$target!='')
            dt <- dt[sel,,drop=FALSE]
        }

        dbg("[getActiveDSEA] done!")
        dsea <- list(table=dt, clust=dr$clust, stats=stats)        

        return(dsea)
    })

    getMOA.target <- shiny::reactive({
        ## meta-GSEA on molecular targets
        dsea <- getActiveDSEA()
        dt <- dsea$table
        shiny::req(dt)        
        targets.list <- lapply(as.character(dt$target),
                               function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]) )
        names(targets.list) <- rownames(dt)
        targets <- setdiff(unique(unlist(targets.list)),c(NA,""," "))
        gmt <- lapply(targets, function(g)
            names(which(sapply(targets.list,function(t) (g %in% t)))))
        names(gmt) <- targets
        
        rnk <- dt$NES
        names(rnk) <- rownames(dt)
        suppressWarnings(
            moa.target <- fgsea::fgsea( gmt, rnk, nperm=20000)
        )
        moa.target <- moa.target[order(-abs(moa.target$NES)),]
        ##head(moa.target)
        return(moa.target)
    })
    
    getMOA.class <- shiny::reactive({
        ## meta-GSEA on MOA terms
        dsea <- getActiveDSEA()
        dt <- dsea$table
        shiny::req(dt)
        moa.list <- lapply(as.character(dt$moa),
                           function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]))
        names(moa.list) <- rownames(dt)
        moa <- setdiff( unlist(moa.list), c("",NA," "))
        gmt <- lapply(moa, function(g) names(which(sapply(moa.list,function(t) (g %in% t)))))
        names(gmt) <- moa
        rnk <- dt$NES
        names(rnk) <- rownames(dt)
        suppressWarnings(
            moa.class <- fgsea::fgsea( gmt, rnk, nperm=20000)
        )
        moa.class <- moa.class[order(-abs(moa.class$NES)),]        
        return(moa.class)
    })

    getMOA <- shiny::reactive({
        moatype <- input$dsea_moatype
        res <- NULL
        if(moatype=='target gene') res <- getMOA.target()
        if(moatype=='drug class')  res <- getMOA.class()
        res
    })
    
    ##================================================================================
    ## PLOTTING
    ##================================================================================
    
    dsea_enplots.RENDER <- shiny::reactive({

        ngs <- inputData()
        if(is.null(ngs$drugs)) return(NULL)        
        shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))        
        shiny::req(input$dsea_contrast, input$dsea_method)

        dbg("[dsea_enplots.RENDER] called!")
        
        dsea <- getActiveDSEA()
        dt <- dsea$table        
        
        ## filter with table selection/search
        ii  <- dsea_table$rows_selected()
        jj  <- dsea_table$rows_all()
        shiny::req(jj)  ## must have non-empty table
        
        if(length(ii)>0) {
            dt <- dt[ii,,drop=FALSE]
        }
        if(length(ii)==0 && length(jj)>0) {
            dt <- dt[jj,,drop=FALSE]
        }

        if(nrow(dt)==0) return(NULL)
        
        ## rank vector for enrichment plots
        dmethod <- input$dsea_method        
        rnk <- dsea$stats
        if(length(rnk)==0) return(NULL)
        
        ## ENPLOT TYPE
        if(nrow(dt)==1) {
            par(oma=c(1,1,1,1))
            par(mfrow=c(1,1), mar=c(4,4,1.1,2), mgp=c(2.3,0.9,0))
            lab.cex = 1            
            xlab="Rank in ordered dataset"
            ylab="Rank metric"
            nc=1
        } else {
            dt <- head(dt, 16)
            lab.cex = 0.75
            xlab=''
            ylab=""
            nc = ceiling(sqrt(nrow(dt)))           
            par(oma=c(0,1.6,0,0))
            par(mfrow=c(nc,nc), mar=c(0.3,1.0,1.3,0), mgp=c(1.9,0.6,0))
        }
        
        i=1
        for(i in 1:nrow(dt)) {
            dx <- rownames(dt)[i]
            dx
            gmtdx <- grep(dx,names(rnk),fixed=TRUE,value=TRUE)  ## L1000 naming
            length(gmtdx)
            ##if(length(gmtdx) < 3) { frame(); next }
            dx1 <- substring(dx,1,26)
            par(cex.axis=0.001)
            if(i%%nc==1) par(cex.axis=0.98)
            suppressWarnings(
                gsea.enplot( rnk, gmtdx, main=dx1, cex.main=1.2,
                            xlab=xlab, ylab=ylab)
            )
            nes <- round(dt$NES[i],2)
            qv  <- round(dt$padj[i],3)
            tt <- c( paste("NES=",nes), paste("q=",qv) )
            legend("topright", legend=tt, cex=0.8, y.intersp=0.85, bty='n')
            if(i%%nc==1 && nrow(dt)>1) {
                mtext('rank metric', side=2, line=1.8, cex=lab.cex)
            }
        }
        
    })    
    
    ##dsea_moaplot.RENDER <- shiny::reactive({
    dsea_moaplot.RENDER <- shiny::reactive({    

        ngs <- inputData()
        shiny::req(ngs)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))    

        dbg("[dsea_moaplot.RENDER] reacted")
        
        res <- getMOA()
        ntop = 16
        jj <- unique(c(head(order(-res$NES),ntop),tail(order(-res$NES),ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        par(mfrow=c(2,1), mar=c(4,3.5,0.1,0), mgp=c(1.7,0.65,0))
        barplot(moa.top, horiz=FALSE, las=3,
                ylab="enrichment  (NES)",
                cex.names=0.96 )
        
    })    

    ##dsea_moaplot.RENDER2 <- shiny::reactive({
    dsea_moaplot.RENDER2 <- shiny::reactive({    

        ngs <- inputData()
        shiny::req(ngs, input$dsea_contrast, input$dsea_method)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))    

        res <- getMOA()
        ntop = 32
        jj <- unique(c(head(order(-res$NES),ntop),tail(order(-res$NES),ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        par(mfrow=c(2,1), mar=c(4,3.5,0.1,0), mgp=c(1.7,0.65,0))
        barplot(moa.top, horiz=FALSE, las=3,
                ylab="enrichment  (NES)",
                cex.names=1 )
        
    })    
        
    dsea_table.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        if(is.null(ngs$drugs)) return(NULL)
        
        dsea <- getActiveDSEA()
        shiny::req(dsea)
        res <- dsea$table
        res$moa <- shortstring(res$moa,60)
        res$target <- shortstring(res$target,30)
        res$drug   <- shortstring(res$drug,60)

        colnames(res) <- sub("moa","MOA",colnames(res))
        DT::datatable( res, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=NULL),
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'Blfrtip', buttons = c('copy','csv','pdf'),
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = '70vh',
                          scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "NES",
                                background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })

    dseaPlotActmap <- function(ngs, dmethod, contr, nterms, nfc) {

        if(is.null(ngs$drugs)) return(NULL)
        ##dmethod="activity/L1000";contr=1        
        nes <- ngs$drugs[[dmethod]]$X
        qv  <- ngs$drugs[[dmethod]]$Q
        score <- nes * (1 - qv)**2
        score[is.na(score)] <- 0
        ## if(NCOL(score)==1) score <- cbind(score,score)  ## UGLY....
        
        ## reduce score matrix
        score = score[order(-score[,contr]**2),,drop=FALSE] ## sort by score
        dsea <- getActiveDSEA()
        res <- dsea$table
        ## filter with table selection/search
        ii  <- dsea_table$rows_all()
        shiny::req(ii)
        if(length(ii)>0) {
            res <- res[ii,,drop=FALSE]
            dd <- intersect(res$drug,rownames(score))
            score = score[dd,,drop=FALSE]
        }
        if(nrow(score) <= 1) return(NULL)
        
        score = head(score,nterms) ## max number of terms    
        score = score[,head(order(-colSums(score**2)),nfc),drop=FALSE] ## max contrs/FC        
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        
        if(input$dsea_normalize) score <- t(t(score) / (1e-8+sqrt(colMeans(score**2))))
        score <- sign(score) * abs(score)**3   ## fudging
        score <- score / (1e-8 + max(abs(score),na.rm=TRUE))
        
        if(NCOL(score)>1) {
            d1 <- as.dist(1-cor(t(score),use="pairwise"))
            d2 <- as.dist(1-cor(score,use="pairwise"))

            d1[is.na(d1)] <- 1
            d2[is.na(d2)] <- 1
            jj=1;ii=1:nrow(score)
            ii <- hclust(d1)$order
            jj <- hclust(d2)$order
            score <- score[ii,jj,drop=FALSE]
        } else {
            score <- score[order(-score[,1]),,drop=FALSE]
        }
        
        cex2=1
        colnames(score) = substring(colnames(score),1,30)
        rownames(score) = substring(rownames(score),1,50)
        cex2=0.85        
        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,1,0,0))

        corrplot::corrplot( score, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex = 0.9*cex2, tl.col = "grey20", tl.srt = 90)

    }      
        
    dsea_actmap.RENDER <- shiny::reactive({


        ngs <- inputData()
        shiny::req(ngs, input$dsea_contrast, input$dsea_method)

        shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        if(is.null(ngs$drugs)) return(NULL)
        
        dmethod="activity/L1000";contr=1
        dmethod <- input$dsea_method        
        contr = input$dsea_contrast
        if(is.null(contr)) return(NULL)

        dseaPlotActmap(ngs, dmethod, contr, nterms=50, nfc=20)

    })    

    dsea_actmap.RENDER2 <- shiny::reactive({


        ngs <- inputData()
        shiny::req(ngs, input$dsea_contrast, input$dsea_method)

        shiny::validate(shiny::need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        if(is.null(ngs$drugs)) return(NULL)
        
        dmethod="activity/L1000";contr=1
        dmethod <- input$dsea_method        
        contr = input$dsea_contrast
        if(is.null(contr)) return(NULL)

        dseaPlotActmap(ngs, dmethod, contr, nterms=50, nfc=100)
        
    })    

        
    ##--------- DSEA enplot plotting module
    dsea_enplots.opts = shiny::tagList()
    
    shiny::callModule(
        plotModule,
        id = "dsea_enplots",
        func = dsea_enplots.RENDER,
        func2 = dsea_enplots.RENDER,         
        title = "Drug connectivity", label="a",
        info.text = "<strong>Drug connectivity</strong> correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space.",
        options = dsea_enplots.opts,
        pdf.height = 6.5, pdf.width = 12.8, 
        height = c(0.54*rowH,750), width=c('auto',1280),
        res = c(78,110),
        add.watermark = WATERMARK        
    )
    
    ##---------- DSEA Activation map plotting module
    dsea_moaplot.opts = shiny::tagList(
        withTooltip( shiny::radioButtons(ns('dsea_moatype'),'Plot type:',c("drug class","target gene"),inline=TRUE),
               "Select plot type of MOA analysis: by class description or by target gene.")
    )
    shiny::callModule(
        plotModule,
        id = "dsea_moaplot",
        func = dsea_moaplot.RENDER,
        func2 = dsea_moaplot.RENDER2, 
        ## csvFunc = getMOA,
        title = "Mechanism of action", label="c",
        info.text = "This plot visualizes the <strong>mechanism of action</strong> (MOA) across the enriched drug profiles. On the vertical axis, the GSEA normalized enrichment score of the MOA class or gene target is plotted. You can switch to visualize between MOA class or target gene.",
        options = dsea_moaplot.opts,
        pdf.width=6, pdf.height=6,
        height = c(0.54*rowH,700), width=c('auto',1400),
        res=c(70,110),
        add.watermark = WATERMARK
    )

    ##-------- Activation map plotting module
    dsea_actmap.opts = shiny::tagList(
        withTooltip(shiny::checkboxInput(ns('dsea_normalize'),'normalize activation matrix',FALSE), "Normalize columns of the activation matrix.")
    )
    shiny::callModule(
        plotModule,
        id = "dsea_actmap",
        func = dsea_actmap.RENDER,
        func2 = dsea_actmap.RENDER2, 
        title = "Activation matrix", label="d",
        info.text = "The <strong>Activation Matrix</strong> visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = dsea_actmap.opts,
        pdf.width=6, pdf.height=9,
        height = c(fullH,750), width=c("100%",1400),        
        res = 72,
        add.watermark = WATERMARK
    )

    ##--------buttons for table
    dsea_table.opts = shiny::tagList(
        withTooltip(shiny::checkboxInput(ns('dseatable_filter'),'only annotated drugs',TRUE),
               "Show only annotated drugs.")
    )  
    dsea_table <- shiny::callModule(
        tableModule,
        id = "dsea_table", label="b",
        func = dsea_table.RENDER, 
        options = dsea_table.opts,
        info.text="<b>Enrichment table.</b> Enrichment is calculated by correlating your signature with known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug.", 
        title = "Enrichment table",
        height = c(360,700)
    )
       
    ##=======================================================================================
    ## CONNECTIVITY MAP
    ##=======================================================================================
    
    ##---------------------------------------------------------------------
    ## Enrichment plot
    ##---------------------------------------------------------------------
    
    cmap_enplot.RENDER <- shiny::reactive({
        
        pgx <- inputData()
                
        ## get selected drug (from table)
        dsea <- getActiveDSEA()
        dt <- dsea$table
        ii <- 1:nrow(dt)
        ii <- 1
        ii <- cmap_table$rows_selected()
        if(length(ii)==0) {
            ii <- cmap_table$rows_all()
        }
        if(length(ii)==0) return(NULL)

        ## draw enrichemtn plot
        d <- dt$drug[ii[1]]        
        rnk <- dsea$stats
        dtype <- sub("[@_].*$","",names(rnk))
        gmt = names(rnk)[dtype==d]
        ## gsea.enplot(rnk, gmt, main=d)
        p1 <- gsea.enplotly(rnk, gmt, main=d)        
        return(p1)
    })
        
    cmap_enplot.info = "<strong>Connectivity map.</strong> correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space."
    cmap_enplot.opts = shiny::tagList()
    
    shiny::callModule(
        plotModule,
        id = "cmap_enplot",
        func = cmap_enplot.RENDER,
        func2 = cmap_enplot.RENDER,
        plotlib = "plotly",
        ## plotlib2 = "plotly",
        title = "ENRICHMENT PLOT", label="a",
        info.text = cmap_enplot.info,
        options = cmap_enplot.opts,
        pdf.height = 6, pdf.width = 10, 
        height = c(305,600), width=c('auto',1000),
        res = c(80,105),
        add.watermark = WATERMARK        
    )

    ##---------------------------------------------------------------------
    ## Enrichment UMAP
    ##---------------------------------------------------------------------
    
    plotCMAP <- function(pgx, db, contr, moa.target, moa.class,
                         labtype = 'drugs', showlabel=TRUE,  
                         npoints=100, nlabel=10,
                         lab.wt=TRUE, lab.gamma=1, lab.cex=1,
                         opacity = 0.15, softmax=1, 
                         title=NULL, plotlib='base')
    {     
        if(0) {
            db="L1000/activity"
            res <- pgx$drugs[[db]]
            contr=colnames(res$stats)[1]
            top = array(1, dimnames=list(rownames(res$X)[1]))
            showlabel=lab.wt=TRUE;opacity=0.9;
            lab.wt=lab.gamma=1;lab.cex=1
            npoints=100;nlabel=10
            labtype = 'drugs';plotlib='base'
        }

        if(!"drugs" %in% names(pgx)) {
            frame()
            text(0.5,0.5,"Error: PGX object does not have CMAP results",col='red3')
            return(NULL)
        }

        res <- pgx$drugs[[db]]
        if(!"clust" %in% names(res)) {
            frame()
            text(0.5,0.5,"Error: PGX object does not have CMAP cluster positions",col='red3')
            return(NULL)
        }
        
        pos <- res$clust
        var <- res$stats[,contr]
        ## highlight genes from table
        nes1  <- res$X[,contr]  ## NES for drugs
                
        ## compute median position of drugs
        xdrugs <- gsub("[_@].*","",rownames(pos))
        pos1 <- apply(pos,2,function(x) tapply(x,xdrugs,median))
        pos1 <- pos1[names(nes1),]
        var1 <- nes1
        
        ## compute median position of moa class
        xmoa <- res$annot[xdrugs,"moa"]
        xmoa <- strsplit(xmoa, split="[\\|;,]")
        nmoa <- sapply(xmoa, length)
        ii   <- unlist(mapply(rep, 1:nrow(pos), nmoa))
        pos2 <- apply(pos[ii,],2,function(x) tapply(x,unlist(xmoa),median))
        moa2 <- moa.class
        nes2 <- array(moa2$NES, dimnames=list(moa2$pathway))               
        var2 <- nes2[rownames(pos2)]
        
        ## compute median position of targets
        xtarget <- res$annot[xdrugs,"target"]
        xtarget <- strsplit(xtarget, split="[\\|;,]")
        ntarget <- sapply(xtarget, length)
        ii   <- unlist(mapply(rep, 1:nrow(pos), ntarget))
        pos3 <- apply(pos[ii,],2,function(x) tapply(x,unlist(xtarget),median))
        ##var3 <- tapply(var[ii], unlist(xtarget), median, na.rm=TRUE)
        moa3 <- moa.target
        nes3 <- array(moa3$NES, dimnames=list(moa3$pathway))               
        var3 <- nes3[rownames(pos3)]
        
        ## create extended positions and variable (drugs, moa, target)
        var  <- var / max(abs(var),na.rm=TRUE)     ## replicate level
        var1 <- var1 / max(abs(var1),na.rm=TRUE)   ## drug level
        var2 <- var2 / max(abs(var2),na.rm=TRUE)   ## MOA class level
        var3 <- var3 / max(abs(var3),na.rm=TRUE)   ## target level                     
        
        rownames(pos)  <- paste0("0:",rownames(pos))
        rownames(pos1) <- paste0("1:",rownames(pos1))
        rownames(pos2) <- paste0("2:",rownames(pos2))
        rownames(pos3) <- paste0("3:",rownames(pos3))        
        names(var)  <- paste0("0:",names(var))
        names(var1) <- paste0("1:",names(var1))
        names(var2) <- paste0("2:",names(var2))
        names(var3) <- paste0("3:",names(var3))        
        
        xpos <- rbind(pos, pos1, pos2, pos3)
        xvar <- c(var, var1, var2, var3)
        sum(duplicated(names(xvar)))
          
        h1=h2=NULL
        ## limit number of labels/points
        labtype        
        if(labtype=='replicate') {
            xx <- res$stats[,contr]
            names(xx) <- paste0("0:",names(xx))
        } else if(labtype=='drugs') {
            xx <- res$X[,contr]
            names(xx) <- paste0("1:",names(xx))            
        } else if(toupper(labtype)=='MOA') {
            xx <- moa.class$NES
            names(xx) <- moa.class$pathway
            names(xx) <- paste0("2:",names(xx))            
        } else if(labtype=='target') {
            xx <- moa.target$NES
            names(xx) <- moa.target$pathway
            names(xx) <- paste0("3:",names(xx))            
        } else {
            message("ERROR:: labtype switch error")
            return(NULL)
        }
        xx <- xx[order(-abs(xx))]
        h1 <- head(names(xx),npoints)
        h2 <- head(h1,nlabel)
        if(!showlabel) {
            h2 <- NULL
        }
                
        if(is.null(title)) title <- contr
        
        plt <- NULL
        if(plotlib=='base') {
            
            wcex = lab.cex 
            if(lab.wt) {
                wcex <- 1.2*(abs(xvar)/max(abs(xvar),na.rm=TRUE))**lab.gamma
                wcex = lab.cex * wcex
                names(wcex) <- names(xvar)
            }
            wcex[is.na(wcex)] <- 1
            
            pgx.scatterPlotXY.BASE(
                xpos, var=xvar, title = title,
                xlab = "UMAP-x", ylab = "UMAP-y",
                labels = sub(".*:","",rownames(xpos)),
                hilight = h1, hilight2 = h2, hilight.cex=1.1,
                cex = 1, cex.lab = wcex, cex.title = 0.95,
                legend = TRUE, zsym = TRUE,
                rstep = 0.2, dlim = 0.01,
                softmax=softmax, opacity = opacity 
            )

        } else {
            plt <- pgx.scatterPlotXY(
                xpos, var=xvar, plotlib=plotlib, title = title,
                xlab = "UMAP-x", ylab = "UMAP-y",
                hilight = h1, hilight2 = h2, hilight.cex=1.1,
                cex = 1, cex.lab = cex.lab, cex.title = 1.0,
                legend = TRUE, zsym = TRUE,
                softmax=1, opacity = opacity)
        }        
        plt
    }
    
    dsea_cmap.RENDER <- shiny::reactive({

        dbg("[dsea_cmap.RENDER] reacted!")
        
        pgx <- inputData()
        shiny::req(pgx)
        db="L1000/gene";contr="treatment:Gefitinib_vs_CT"

        db <- input$dsea_method
        contr <- input$dsea_contrast
        shiny::req(db)
        shiny::req(contr)

        dsea <- getActiveDSEA()
        
        ## get reactive values
        rows_selected=1;rows_all=1:nrow(dsea$table)
        rows_selected <- cmap_table$rows_selected()
        rows_all <- cmap_table$rows_all()                
        
        if(is.null(rows_all) || length(rows_all)==0) {
            return(NULL)
        }

        drugs_all <- rownames(dsea$table)[rows_all]
       
        moa.class  <- getMOA.class()
        moa.target <- getMOA.target()    
               
        labtype <- input$cmap_labeltype
        nlabel  <- as.integer(input$cmap_nlabel)
        showlab <- ("show" %in% input$cmap_labeloptions)
        lab.wt <- !("fixed" %in% input$cmap_labeloptions)
                
        ##---------------  plot -------------------
        all.contr <- colnames(pgx$contrasts)
        all.contr
        contr <- all.contr[1]
        
        nr <- ceiling(sqrt(length(all.contr)))
        par(mfrow=c(nr,nr))            
        for(contr in all.contr) {                
            tt <- paste0(contr," (",toupper(labtype),")")
            ##tt <- toupper(labtype)                
            plotCMAP(pgx, db, contr,
                     moa.target, moa.class,
                     labtype=labtype, showlabel=showlab,
                     lab.wt=lab.wt, lab.gamma=1, lab.cex=1.6,
                     opacity=0.15, softmax=0, 
                     npoints=nlabel, nlabel=nlabel,
                     title=tt, plotlib='base')                           
        }
            
        dbg("[dsea_cmap.RENDER] done!")
        
    })
    
    dsea_cmap.info = "<strong>Connectivity map.</strong> correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space."

    dsea_cmap.opts = shiny::tagList(
        tipifyL(shiny::radioButtons(
            ns('cmap_labeltype'),'label type:',
            c("drugs","MOA","target"),inline=TRUE),
            "Label point with drugs, MOA terms or targets (if no drug selected)."),
        tipifyL(shiny::radioButtons(ns('cmap_nlabel'),'number of labels:',c(3,10,20,100),
                            selected=10, inline=TRUE),
               "Number of labels to show."),
        tipifyL(shiny::checkboxGroupInput(
            ns('cmap_labeloptions'),'label options:',  choices=c("show","fixed"),
            selected=c("show"), inline=TRUE ), "Other labels options.")
    )
    
    shiny::callModule(
        plotModule,
        id = "dsea_cmap",
        func = dsea_cmap.RENDER,
        func2 = dsea_cmap.RENDER,
        title = "CONNECTIVITY MAP", label="b",
        info.text = dsea_cmap.info,
        options = dsea_cmap.opts,
        pdf.height = 8, pdf.width = 8, 
        height = c(750,750), width=c('auto',900),
        res = c(80,105),
        add.watermark = WATERMARK        
    )

    ##---------------------------------------------------------------------
    ## Enrichment table
    ##---------------------------------------------------------------------
    
    cmap_table.RENDER <- shiny::reactive({
        pgx <- inputData()
        shiny::req(pgx)
        if(is.null(pgx$drugs)) return(NULL)
        
        dsea <- getActiveDSEA()
        shiny::req(dsea)
        res <- dsea$table
        res$moa    <- shortstring(res$moa,30)
        res$target <- shortstring(res$target,80)
        res$drug   <- shortstring(res$drug,60)
        res$pval   <- NULL
        res$padj   <- NULL                
        
        colnames(res) <- sub("moa","MOA",colnames(res))
        DT::datatable( res, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=NULL),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = '70vh', scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
            DT::formatStyle( "NES",
                            background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                            backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center') 
    })
    
    cmap_table <- shiny::callModule(
        tableModule,
        id = "cmap_table", label="c",
        func = cmap_table.RENDER,
        info.text="<b>Enrichment table.</b> Enrichment is calculated by correlating your signature with known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug.", 
        title = "CONNECTIVITY TABLE",
        height = c(380,740)
    )

  })
}
