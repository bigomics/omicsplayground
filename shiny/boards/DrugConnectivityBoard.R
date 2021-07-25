##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing DrugConnectivityBoard")

DrugConnectivityInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

DrugConnectivityUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        flex = c(1),
        height = 750,
        tabsetPanel(
            id = ns("tabs"),
            tabPanel("Drug enrichment",uiOutput(ns("DSEA_enrichment_UI"))),
            tabPanel("Connectivity map",uiOutput(ns("DSEA_cmap_UI")))            
            ## tabPanel("Fire plot (dev)",uiOutput(ns("fireplot_UI")))            
        )
    )
}

DrugConnectivityBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    inputData <- env[["load"]][["inputData"]]
    fullH = 750
    rowH = 660  ## row height of panel
    tabH = 200  ## row height of panel
    tabH = '60vh'  ## row height of panel    
    description = "<h3>Drug Connectivity</h3> Perform drug connectivity analysis
to see if certain drug activity or drug sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate biological functions through mechanism-of-action (MOA) and known drug molecular targets."
    output$description <- renderUI(HTML(description))
    
    dsea_infotext = paste("<b>This module performs drug enrichment analysis</b> to see if certain drug activity or drug sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate biological functions through mechanism-of-action (MOA) and known drug molecular targets.

<br><br> In the <a href='https://portals.broadinstitute.org/cmap/'>Drug Connectivity Map</a> panel, you can correlate your signature with known drug profiles from the L1000 database. An activation-heatmap compares drug activation profiles across multiple contrasts. This facilitates to quickly see and detect the similarities between contrasts for certain drugs.


<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
")

    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("dsea_info"), "Youtube", icon = icon("youtube") ),
                   "Show more information about this module."),
            hr(), br(),             
            tipify( selectInput(ns("dsea_contrast"),"Contrast:", choices=NULL),
                   "Select the contrast corresponding to the comparison of interest.",
                   placement="top"),
            tipify( selectInput(ns('dsea_method'),"Analysis type:", choices = ""),
                   "Select type of drug enrichment analysis: activity or sensitivity (if available).",
                   placement="top"),
            ##tipify( actionLink(ns("dsea_options"), "Options", icon=icon("cog", lib = "glyphicon")),
            ##       "Show/hide advanced options", placement="top"),
            ## br(),
            ## conditionalPanel(
            ##     "input.dsea_options % 2 == 1", ns=ns,
            ##     tagList()
            ## )
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!

    observe({
        ngs <- inputData()
        req(ngs)
        ct <- names(ngs$drugs)
        updateSelectInput(session, "dsea_method", choices=ct)
    })
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$dsea_info, {
        showModal(modalDialog(
            title = HTML("<strong>Drug Connectivity Analysis Board</strong>"),
            HTML(dsea_infotext),
            easyClose = TRUE, size="l" ))
    })

    observe({
        ngs <- inputData()
        req(ngs)
        ct <- colnames(ngs$model.parameters$contr.matrix)
        ##ct <- c(ct,"<sd>")
        updateSelectInput(session, "dsea_contrast", choices=sort(ct) )
    })
    
    ##================================================================================
    ## Reactive functions
    ##================================================================================

    getDseaTable <- reactive({
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        req(ngs)        
        req(input$dsea_contrast, input$dsea_method)
        
        dbg("[getDseaTable] reacted")
        
        dmethod="CTRPv2/sensitivity";comparison=1
        dmethod="L1000/activity";comparison=1        
        dmethod="L1000/gene";comparison=1        

        names(ngs$gx.meta$meta)
        comparison = input$dsea_contrast
        if(is.null(comparison)) return(NULL)

        dbg("[getDseaTable] 1: ")
        
        names(ngs$drugs)
        dmethod <- input$dsea_method
        if(is.null(dmethod)) return(NULL)

        dbg("[getDseaTable] 2: ")
        
        dr <- ngs$drugs[[dmethod]]
        nes <- round(dr$X[,comparison],4)
        pv  <- round(dr$P[,comparison],4)
        qv  <- round(dr$Q[,comparison],4)
        drug <- rownames(dr$X)
        stats <- dr$stats
        annot <- dr$annot
        nes[is.na(nes)] <- 0
        qv[is.na(qv)] <- 1
        pv[is.na(pv)] <- 1

        dbg("[getDseaTable] 3: dim.annot = ",dim(annot))        
        
        ## !!!SHOULD MAYBE BE DONE IN PREPROCESSING???
        if(is.null(annot)) {
            message("[getDseaTable] WARNING:: missing drug annotation in PGX file!")
            annot <- read.csv(file.path(FILESX,"sig/L1000_repurposing_drugs.txt"),
                              sep="\t", comment.char="#")
            rownames(annot) <- annot$pert_iname
        }

        dbg("[getDseaTable] 4: ")

        jj <- match(toupper(drug), toupper(rownames(annot)))
        annot <- annot[jj,c("moa","target")]        
        dt <- data.frame( drug=drug, NES=nes, pval=pv, padj=qv, annot)
        ##dt <- dt[order(-abs(dt$NES)),]
        ##dt <- dt[order(dt$pval,-abs(dt$NES)),]
        dt <- dt[order(-dt$NES),]

        dbg("[getDseaTable] 5: ")
        
        req(input$dseatable_filter)
        if(input$dseatable_filter) {
            sel <- which(dt$moa!='' | dt$target!='')
            dt <- dt[sel,,drop=FALSE]
        }

        dbg("[getDseaTable] done!")
        
        return(dt)
    })

    getMOA.target <- reactive({
        ## meta-GSEA on molecular targets
        dt <- getDseaTable()
        require(dt)        
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
            res <- fgsea( gmt, rnk, nperm=20000)
        )
        res <- res[order(res$pval),]
        return(res)
    })
    
    getMOA.class <- reactive({
        ## meta-GSEA on MOA terms
        dt <- getDseaTable()
        require(dt)
        moa.list <- lapply(as.character(dt$moa),
                           function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]))
        names(moa.list) <- rownames(dt)
        moa <- setdiff( unlist(moa.list), c("",NA," "))
        gmt <- lapply(moa, function(g) names(which(sapply(moa.list,function(t) (g %in% t)))))
        names(gmt) <- moa
        rnk <- dt$NES
        names(rnk) <- rownames(dt)
        suppressWarnings(
            res <- fgsea( gmt, rnk, nperm=20000)
        )
        res <- res[order(res$pval),]
        return(res)
    })

    getMOA <- reactive({
        moatype <- input$dsea_moatype
        if(moatype=='target gene') res <- getMOA.target()
        if(moatype=='drug class')  res <- getMOA.class()
        res
    })

    ##================================================================================
    ## PLOTTING
    ##================================================================================

    
    dsea_enplots.RENDER %<a-% reactive({

        ngs <- inputData()
        if(is.null(ngs$drugs)) return(NULL)        
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))        
        req(input$dsea_contrast, input$dsea_method)

        dbg("[dsea_enplots.RENDER] called!")
        
        comparison=1
        comparison = input$dsea_contrast
        if(is.null(comparison)) return(NULL)

        dbg("[dsea_enplots.RENDER] 1:")
        
        dt <- getDseaTable()
        ## filter with table selection/search
        ii  <- dsea_table$rows_selected()
        jj  <- dsea_table$rows_all()
        req(jj)
        if(length(jj)>0) {
            dt <- dt[jj,,drop=FALSE]
        }

        dbg("[dsea_enplots.RENDER] dim.dt", dim(dt))
        if(nrow(dt)==0) return(NULL)
        
        ## rank vector for enrichment plots
        dmethod <- input$dsea_method        
        dbg("[dsea_enplots.RENDER] dmethod = ",dmethod)
        dbg("[dsea_enplots.RENDER] comparison = ",comparison)

        stats <- ngs$drugs[[dmethod]]$stats
        dbg("[dsea_enplots.RENDER] dim(stats) = ",dim(stats))
        dbg("[dsea_enplots.RENDER] colnames(stats) = ",colnames(stats))        

        rnk <- stats[,comparison]
        dctype <- sub("_.*$","",names(rnk))
        ##table(rownames(dt) %in% dctype)
        ##table(sapply(rownames(dt), function(g) sum(grepl(g,names(rnk),fixed=TRUE))))

        dbg("[dsea_enplots.RENDER] len.rnk = ", length(rnk))
        if(length(rnk)==0) return(NULL)
        
        ## ENPLOT TYPE
        ##itop <- c( head(order(-dt$NES),10), tail(order(-dt$NES),10))
        itop <- 1:min(16,nrow(dt))
        par(oma=c(0,1.6,0,0))
        par(mfrow=c(4,4), mar=c(0.3,1.0,1.3,0), mgp=c(1.9,0.6,0))
        i=1
        for(i in itop) {

            dx <- rownames(dt)[i]
            dx
            gmtdx <- grep(dx,names(rnk),fixed=TRUE,value=TRUE)  ## L1000 naming
            length(gmtdx)
            ##if(length(gmtdx) < 3) { frame(); next }
            dx1 <- substring(dx,1,26)
            par(cex.axis=0.001)
            if(i%%4==1) par(cex.axis=0.98)
            suppressWarnings(
                gsea.enplot( rnk, gmtdx, main=dx1, cex.main=1.2, xlab="", ylab="")
            )
            nes <- round(dt$NES[i],2)
            qv  <- round(dt$padj[i],3)
            tt <- c( paste("NES=",nes), paste("q=",qv) )
            legend("topright", legend=tt, cex=0.8, y.intersp=0.85, bty='n')
            if(i%%4==1) {
                mtext('rank metric', side=2, line=1.8, cex=0.75)
            }
        }
        
    })    
    

    ##dsea_moaplot.RENDER %<a-% reactive({
    dsea_moaplot.RENDER <- reactive({    

        ngs <- inputData()
        req(ngs)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    

        dbg("[dsea_moaplot.RENDER] reacted")
        
        res <- getMOA()
        dbg("[dsea_moaplot.RENDER] dim(res)=",dim(res))
        ntop = 16
        jj <- unique(c(head(order(-res$NES),ntop),tail(order(-res$NES),ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        
        par(mfrow=c(2,1), mar=c(4,3.5,0.1,0), mgp=c(1.7,0.65,0))
        barplot(moa.top, horiz=FALSE, las=3,
                ylab="enrichment  (NES)", cex.names=cex )
        
    })    

    ##dsea_moaplot.RENDER2 %<a-% reactive({
    dsea_moaplot.RENDER2 <- reactive({    

        ngs <- inputData()
        req(ngs, input$dsea_contrast, input$dsea_method)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    

        dbg("[dsea_moaplot.RENDER2] reacted")
        
        res <- getMOA()
        dbg("[dsea_moaplot.RENDER2] dim(res)=",dim(res))
        ntop = 32
        jj <- unique(c(head(order(-res$NES),ntop),tail(order(-res$NES),ntop)))
        moa.top <- res$NES[jj]
        names(moa.top) <- res$pathway[jj]
        
        par(mfrow=c(2,1), mar=c(4,3.5,0.1,0), mgp=c(1.7,0.65,0))
        barplot(moa.top, horiz=FALSE, las=3,
                ylab="enrichment  (NES)", cex.names=cex )
        
    })    
        
    dsea_table.RENDER <- reactive({
        ngs <- inputData()
        req(ngs)
        if(is.null(ngs$drugs)) return(NULL)
        
        res <- getDseaTable()
        req(res)
        res$moa <- shortstring(res$moa,60)
        res$target <- shortstring(res$target,30)
        res$drug   <- shortstring(res$drug,60)

        colnames(res) <- sub("moa","MOA",colnames(res))
        DT::datatable( res, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
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

    dseaPlotActmap <- function(ngs, dmethod, comparison, nterms, nfc) {

        if(is.null(ngs$drugs)) return(NULL)
        ##dmethod="activity/L1000";comparison=1        
        nes <- ngs$drugs[[dmethod]]$X
        qv  <- ngs$drugs[[dmethod]]$Q
        score <- nes * (1 - qv)**2
        score[is.na(score)] <- 0
        ## if(NCOL(score)==1) score <- cbind(score,score)  ## UGLY....
        
        ## reduce score matrix
        score = score[order(-score[,comparison]**2),,drop=FALSE] ## sort by score
        if(1) {
            res <- getDseaTable()
            ## filter with table selection/search
            ii  <- dsea_table$rows_all()
            req(ii)
            if(length(ii)>0) {
                res <- res[ii,,drop=FALSE]
                dd <- intersect(res$drug,rownames(score))
                score = score[dd,,drop=FALSE]
            }
        }
        if(nrow(score) <= 1) return(NULL)
        
        score = head(score,nterms) ## max number of terms    
        score = score[,head(order(-colSums(score**2)),nfc),drop=FALSE] ## max comparisons/FC        
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
        require(corrplot)
        corrplot( score, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex = 0.9*cex2, tl.col = "grey20", tl.srt = 90)

    }      
        
    dsea_actmap.RENDER <- reactive({

        require(igraph)
        ngs <- inputData()
        req(ngs, input$dsea_contrast, input$dsea_method)

        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        if(is.null(ngs$drugs)) return(NULL)
        
        dmethod="activity/L1000";comparison=1
        dmethod <- input$dsea_method        
        comparison = input$dsea_contrast
        if(is.null(comparison)) return(NULL)

        dseaPlotActmap(ngs, dmethod, comparison, nterms=50, nfc=20)

    })    

    dsea_actmap.RENDER2 <- reactive({

        require(igraph)
        ngs <- inputData()
        req(ngs, input$dsea_contrast, input$dsea_method)

        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        if(is.null(ngs$drugs)) return(NULL)
        
        dmethod="activity/L1000";comparison=1
        dmethod <- input$dsea_method        
        comparison = input$dsea_contrast
        if(is.null(comparison)) return(NULL)

        dseaPlotActmap(ngs, dmethod, comparison, nterms=50, nfc=100)
        
    })    

        
    ##--------- DSEA enplot plotting module
    dsea_enplots.opts = tagList()
    
    callModule(
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
    ##outputOptions(output, "dsea_enplots", suspendWhenHidden=FALSE) ## important!!!
    
    ##---------- DSEA Activation map plotting module
    dsea_moaplot.opts = tagList(
        tipify( radioButtons(ns('dsea_moatype'),'Plot type:',c("drug class","target gene"),inline=TRUE),
               "Select plot type of MOA analysis: by class description or by target gene.")
        ##tipify( radioButtons(ns('dsea_moamethod'),'Method:',c("count","enrichment"),inline=TRUE),
        ##       "Select method of MOA analysis: count significant terms or enrichment test.")
    )
    callModule(
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
    dsea_actmap.opts = tagList(
        tipify(checkboxInput(ns('dsea_normalize'),'normalize activation matrix',FALSE), "Normalize columns of the activation matrix.")
    )
    callModule(
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
    dsea_table.opts = tagList(
        tipify(checkboxInput(ns('dseatable_filter'),'only annotated drugs',TRUE),
               "Show only annotated drugs.")
    )  
    dsea_table <- callModule(
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
    
    cmap_enplot.RENDER <- reactive({
        
        pgx <- inputData()
        
        method="L1000/activity"
        method <- input$dsea_method
        res <- pgx$drugs[[method]]
        names(res)

        ## get selected drug (from table)
        dt <- getDseaTable()
        ii <- 1:nrow(dt)
        ii <- 1
        ii <- cmap_table$rows_selected()
        if(length(ii)==0) {
            ii <- cmap_table$rows_all()
        }
        if(length(ii)==0) return(NULL)

        ## draw enrichemtn plot
        d <- dt$drug[ii[1]]        
        contr=colnames(res$stats)[1]
        contr <- input$dsea_contrast        
        rnk <- res$stats[,contr]
        dtype <- sub("_.*$","",names(rnk))
        gmt = names(rnk)[dtype==d]
        ## gsea.enplot(rnk, gmt, main=d)
        p1 <- gsea.enplotly(rnk, gmt, main=d)        
        return(p1)
    })
        
    cmap_enplot.info = "<strong>Connectivity map.</strong> correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space."
    cmap_enplot.opts = tagList()
    
    callModule(
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
    
    getCmapHilights <- reactive({

        dt <- getDseaTable()
        dt.nes <- array(dt$NES, dimnames=list(dt$drug))       

        ii <- 1:nrow(dt)
        ii <- 1
        ii <- cmap_table$rows_selected()
        is.sel = (length(ii) >0)
        ltype <- input$cmap_labeltype
        if(!is.sel && ltype=='drugs') {
            ## show the top 20 drugs (drug level)
            ii <- cmap_table$rows_all()
            nes <- dt.nes[ii]
            nes <- nes[order(-abs(nes))]            
            d1 = d2 = nes
        } else if(!is.sel && ltype=='MOA') {
            ## show the top 20 drugs (drug level)
            ii <- cmap_table$rows_all()
            mm <- unlist(strsplit(dt$moa[ii],split='\\|'))
            moa <- getMOA.class()
            nes <- array(moa$NES, dimnames=list(moa$pathway))
            nes <- nes[intersect(mm,names(nes))]
            nes <- nes[order(-abs(nes))]
            d1 = d2 = nes              
        } else if(!is.sel && ltype=='target') {
            ## show the top 20 drugs (drug level)
            ii <- cmap_table$rows_all()
            mm <- unlist(strsplit(dt$target[ii],split='\\|'))
            moa <- getMOA.target()
            nes <- array(moa$NES, dimnames=list(moa$pathway))
            nes <- nes[intersect(mm,names(nes))]
            nes <- nes[order(-abs(nes))]            
            d1 = d2 = nes
        } else {
            ## show the replicates of selected drug (replicate level)            
            ii <- ii[1]
            d1 = d2 = dt.nes[ii]
        }
        ht <- list(d1=d1, d2=d2, is.sel=is.sel)
        ht
    })
    
    plotCMAP <- function(pgx, db, contr, ht, label, wt.font, plotlib)
    {     
        res <- pgx$drugs[[db]]
        names(res)
        if(0) {
            db="L1000/activity"
            contr=colnames(res$stats)[1]
            g = array(1, dimnames=list(rownames(res$X)[1]))
            ht=list(d1=g, d2=g, is.sel=FALSE)
            label=wt.font=TRUE
            plotlib='base'
        }
        
        if(!"clust" %in% names(res)) {
            frame()
            text(0.5,0.5,"Error: PGX object does not have CMAP cluster positions",col='red3')
            return(NULL)
        }
        
        pos <- res$clust
        var <- res$stats[,contr]
        ## highlight genes from table
        x1 <- res$X[,contr]
        xdrugs <- gsub("_.*","",rownames(pos))

        h1=h2=head(rownames(res$X))        
        h1 <- names(ht$d1)
        h2 <- names(ht$d2)
        if(ht$is.sel) {
            j1 <- which(xdrugs %in% h1)
            h1 <- rownames(pos)[j1]
            j2 <- which(xdrugs %in% h2)
            h2 <- rownames(pos)[j2]
        }

        ## limit number of labels/points
        h1 <- head(h1,100)
        h2 <- head(h2,20)
        
        if(!label) {
            h2 <- NULL
        }
        
        dbg("[plotCMAP] dim.pos = ", dim(pos) )
        dbg("[plotCMAP] len.xdrugs = ", length(xdrugs) )
        dbg("[plotCMAP] len.var = ", length(var) )                
        var <- var[match(rownames(pos),names(var))]
        names(var) <- rownames(pos)
        dbg("[plotCMAP] 2: len.var = ", length(var) )
        
        ## compute median position of drugs
        xdrugs <- sub("_.*","",rownames(pos))
        pos1 <- apply(pos,2,function(x) tapply(x,xdrugs,median))
        var1 <- tapply(var, xdrugs, median, na.rm=TRUE)

        ## compute median position of moa class
        xmoa <- res$annot[xdrugs,"moa"]
        xmoa <- strsplit(xmoa, split="\\|")
        nmoa <- sapply(xmoa, length)
        ii   <- unlist(mapply(rep, 1:nrow(pos), nmoa))
        pos2 <- apply(pos[ii,],2,function(x) tapply(x,unlist(xmoa),median))
        ##var2 <- tapply(var[ii], unlist(xmoa), median, na.rm=TRUE)
        moa2 <- getMOA.class()
        nes2 <- array(moa2$NES, dimnames=list(moa2$pathway))               
        var2 <- nes2[rownames(pos2)]
        
        ## compute median position of moa targets
        xtarget <- res$annot[xdrugs,"target"]
        xtarget <- strsplit(xtarget, split="\\|")
        ntarget <- sapply(xtarget, length)
        ii   <- unlist(mapply(rep, 1:nrow(pos), ntarget))
        pos3 <- apply(pos[ii,],2,function(x) tapply(x,unlist(xtarget),median))
        ##var3 <- tapply(var[ii], unlist(xtarget), median, na.rm=TRUE)
        moa3 <- getMOA.target()
        nes3 <- array(moa3$NES, dimnames=list(moa3$pathway))               
        var3 <- nes3[rownames(pos3)]
        dbg("[plotCMAP] head.var3 = ", head(var3))
        
        ## create extended positions and variable (drugs, moa, target)
        xpos <- rbind(pos, pos1, pos2, pos3)
        var <- var / max(abs(var),na.rm=TRUE)
        var1 <- var1 / max(abs(var1),na.rm=TRUE)
        var2 <- var2 / max(abs(var2),na.rm=TRUE)
        var3 <- var3 / max(abs(var3),na.rm=TRUE)                        
        xvar <- c(var, var1, var2, var3)

        dbg("[plotCMAP] len.xvar = ", length(xvar))
        dbg("[plotCMAP] dim.xpos = ", dim(xpos))
        dbg("[plotCMAP] names.xvar = ", head(names(xvar)))
        dbg("[plotCMAP] names.var1 = ", head(names(var1)))
        dbg("[plotCMAP] names.var2 = ", head(names(var2)))        
        dbg("[plotCMAP] names.var3 = ", head(names(var3)))
        dbg("[plotCMAP] xvar.h2 = ", head(xvar[h2]))            

        dbg("[plotCMAP] h1 = ", h1)
        dbg("[plotCMAP] h2 = ", h2)
        dbg("[plotCMAP] length.h1 = ", length(h1))
        dbg("[plotCMAP] length.h2 = ", length(h2))
        dbg("[plotCMAP] h2.in.pos1 = ", table(h2 %in% rownames(xpos) ))
        
        cex.lab <- ifelse(length(h2) > 10, 1.05, 1.2)
        cex.lab <- ifelse(length(h2) > 20, 0.95, cex.lab)
        cex.lab
        
        plt <- NULL
        if(plotlib=='base') {
            
            wcex.lab = cex.lab 
            if(wt.font) {
                wcex <- 1.2*(abs(xvar)/max(abs(xvar),na.rm=TRUE))**0.3
                if(!ht$is.sel) wcex <- 1.3*wcex
                wcex.lab = wcex.lab * wcex
                names(wcex.lab) <- names(xvar)
            }
            wcex.lab[is.na(wcex.lab)] <- 1

            dbg("[plotCMAP] len.wcexlab = ", length(wcex.lab))
            dbg("[plotCMAP] min.wcexlab = ", min(wcex.lab,na.rm=TRUE))
            dbg("[plotCMAP] max.wcexlab = ", max(wcex.lab,na.rm=TRUE))                        
            dbg("[plotCMAP] head.wcexlab = ", head(wcex.lab))
            dbg("[plotCMAP] head.wcexlab.h2 = ", head(wcex.lab[h2]))                        
            dbg("[plotCMAP] sum.isna.wcexlab = ", sum(is.na(wcex.lab)))
            
            pgx.scatterPlotXY.BASE(
                xpos, var=xvar, title = contr,
                xlab = "UMAP-x", ylab = "UMAP-y",
                hilight = h1, hilight2 = h2, hilight.cex=1.1,
                cex = 1, cex.lab = wcex.lab, cex.title = 1.0,
                legend = TRUE, zsym = TRUE,
                rstep=0.2, dlim=0.01,
                ## col=c("grey80","grey80"),
                softmax=0, opacity = 0.2)
            ## title(contr, cex.main=1, line=0.3)

        } else {
            plt <- pgx.scatterPlotXY(
                xpos, var=xvar, plotlib=plotlib,
                title = contr,
                xlab = "UMAP-x", ylab = "UMAP-y",
                hilight = h1, hilight2 = h2, hilight.cex=1.1,
                cex = 1, cex.lab = cex.lab, cex.title = 1.0,
                legend = TRUE, zsym = TRUE,
                ##rstep=0.2, dlim=0.01
                softmax=0, opacity = 0.2)
        }        
        plt
    }
    
    dsea_cmap.RENDER <- reactive({
        pgx <- inputData()
        ht <- getCmapHilights()
        req(pgx,ht)
        db <- input$dsea_method
        contr <- input$dsea_contrast
        lab <- input$cmap_showlabels
        wt.font <- !input$cmap_fixedlabel
        plotCMAP(pgx, db, contr, ht, lab, wt.font, plotlib='base')        
    })
    
    dsea_cmap.RENDER2 <- reactive({
        pgx <- inputData()
        ht <- getCmapHilights()
        req(pgx,ht)
        db <- input$dsea_method
        contr <- input$dsea_contrast
        lab <- input$cmap_showlabels
        wt.font <- FALSE
        ## wt.font <- input$cmap_wtfont
        plotCMAP(pgx, db, contr, ht, lab, wt.font, plotlib='plotly')        
    })
    
    dsea_cmap.info = "<strong>Connectivity map.</strong> correlates your signature with known drug profiles from the L1000 database, and shows similar and opposite profiles by running the GSEA algorithm on the drug profile correlation space."

    dsea_cmap.opts = tagList(
        tipify(radioButtons(ns('cmap_labeltype'),'label type:',
                            c("drugs","MOA","target"),inline=TRUE),
               "Label point with drugs, MOA terms or targets (if no drug selected)."),                
        tipify(checkboxInput(ns('cmap_showlabels'),'show labels',TRUE),
               "Show/hide labels."),
        tipify(checkboxInput(ns('cmap_fixedlabel'),'fixed label size',FALSE),
               "Fix size of labels (not proportional to enrichment score).")
    )
    
    callModule(
        plotModule,
        id = "dsea_cmap",
        func = dsea_cmap.RENDER,
        func2 = dsea_cmap.RENDER2,
        ##plotlib = "plotly",
        plotlib2 = "plotly",
        title = "CONNECTIVITY MAP", label="b",
        info.text = dsea_cmap.info,
        options = dsea_cmap.opts,
        pdf.height = 8, pdf.width = 8, 
        height = c(700,750), width=c('auto',900),
        res = c(80,105),
        add.watermark = WATERMARK        
    )

    ##---------------------------------------------------------------------
    ## Enrichment table
    ##---------------------------------------------------------------------
    
    cmap_table.RENDER <- reactive({
        pgx <- inputData()
        req(pgx)
        if(is.null(pgx$drugs)) return(NULL)
        
        res <- getDseaTable()
        req(res)
        res$moa    <- shortstring(res$moa,30)
        res$target <- shortstring(res$target,80)
        res$drug   <- shortstring(res$drug,60)
        ## res$target <- NULL
        res$pval   <- NULL
        res$padj   <- NULL                
        
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
                          scrollY = '70vh', scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
            DT::formatStyle( "NES",
                            background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                            backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center') 
    })
    
    cmap_table <- callModule(
        tableModule,
        id = "cmap_table", label="c",
        func = cmap_table.RENDER, 
        ## options = cmap_table.opts,
        info.text="<b>Enrichment table.</b> Enrichment is calculated by correlating your signature with known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug.", 
        title = "CONNECTIVITY TABLE",
        height = c(350,740)
    )
    
    ##=======================================================================================
    ## PAGE LAYOUT
    ##=======================================================================================
    
    dsea_enrichment_caption = "<b>(a)</b> <b>Drug connectivity</b> correlates your signature with known drug perturbation profiles from the L1000 database. The figures show the most similar (or opposite) profiles by running the GSEA algorithm on the profile correlation space. <b>(b)</b> <b>Enrichment table</b> summarizing the statistical results of the drug enrichment analysis. <b>(c)</b> <b>Mechanism-of-action</b> plot showing the top most frequent drug class (or target genes) having similar or opposite enrichment compared to the query signature. <b>(d)</b> <b>Activation matrix</b> visualizing enrichment levels of drug signatures across multiple contrast profiles." 

    output$DSEA_enrichment_UI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,            
            div(HTML(dsea_enrichment_caption),class="caption"),
            br(),
            fillRow(
                height = rowH,
                flex = c(2.8,1), 
                fillCol(
                    flex = c(1.4,0.15,1),
                    height = rowH,
                    fillRow(
                        flex=c(1.1,1),
                        plotWidget(ns("dsea_enplots")),
                        plotWidget(ns("dsea_moaplot"))
                    ),
                    br(),  ## vertical space
                    tableWidget(ns("dsea_table"))        
                ),
                plotWidget(ns("dsea_actmap"))
            )
        )
    })
    outputOptions(output, "DSEA_enrichment_UI", suspendWhenHidden=FALSE) ## important!!!

    dsea_cmap_caption = "<b>(a)</b> <b>Enrichment plot.</b> Enrichment of the selected drug perturbation profile with your signature. <b>(b)</b> <b>Enrichment table</b> summarizing the statistical results of the drug enrichment analysis. <b>(c)</b> <b>Connectivity map.</b> Plot showing the top signatures as UMAP. Each point is one L1000 experiment. The color corresponds to the rank correlation between the drug signatures and your selected contrast." 

    output$DSEA_cmap_UI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,            
            div(HTML(dsea_cmap_caption),class="caption"),
            br(),
            fillRow(
                height = rowH,
                flex = c(1,0.06,1.5),
                fillCol(
                    flex = c(1.1,0.05,1),                    
                    plotWidget(ns("cmap_enplot")),
                    br(),
                    tableWidget(ns("cmap_table"))                    
                ),
                br(),
                plotWidget(ns("dsea_cmap"))
            )
        )
    })
    outputOptions(output, "DSEA_cmap_UI", suspendWhenHidden=FALSE) ## important!!!
    

}
