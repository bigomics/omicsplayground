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
            tabPanel("Drug CMap",uiOutput(ns("DSEA_analysis_UI")))
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
    description = "<b>Drug Connectivity Analysis</b>. <br> Perform drug connectivity analysis
to see if certain drug activity or drug sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate biological functions through mechanism-of-action (MOA) and known drug molecular targets. "
    output$description <- renderUI(HTML(description))
    
    dr_infotext = paste("<b>This module performs drug enrichment analysis</b> to see if certain drug activity or drug sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate biological functions through mechanism-of-action (MOA) and known drug molecular targets.

<br><br> In the <a href='https://portals.broadinstitute.org/cmap/'>Drug Connectivity Map</a> panel, you can correlate your signature with known drug profiles from the L1000 database. An activation-heatmap compares drug activation profiles across multiple contrasts. This facilitates to quickly see and detect the similarities between contrasts for certain drugs.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
")

    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("dr_info"), "Youtube", icon = icon("youtube") ),
                   "Show more information about this module."),
            hr(), br(),             
            tipify( selectInput(ns("dr_contrast"),"Contrast:", choices=NULL),
                   "Select the contrast corresponding to the comparison of interest.",
                   placement="top"),
            tipify( selectInput(ns('dsea_method'),"Analysis type:", choices = ""),
                   "Select type of drug enrichment analysis: activity or sensitivity (if available).",
                   placement="top"),
                tipify( actionLink(ns("dr_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Show/hide advanced options", placement="top"),
            br(),
            conditionalPanel(
                "input.dr_options % 2 == 1", ns=ns,
                tagList(
                    tipify(checkboxInput(ns('dr_normalize'),'normalize activation matrix',TRUE),
                           "Click to fine-tune the coloring of an activation matrices.")
                )
            )
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
    
    observeEvent( input$dr_info, {
        showModal(modalDialog(
            title = HTML("<strong>Drug Connectivity Analysis Board</strong>"),
            HTML(dr_infotext),
            easyClose = TRUE, size="l" ))
    })

    observe({
        ngs <- inputData()
        req(ngs)
        ct <- colnames(ngs$model.parameters$contr.matrix)
        ##ct <- c(ct,"<sd>")
        updateSelectInput(session, "dr_contrast", choices=sort(ct) )
    })
    
    ##================================================================================
    ## Drug signature enrichment analysis L1000
    ##================================================================================

    getDseaTable <- reactive({
        ngs <- inputData()
        alertDataLoaded(session,ngs)
        req(ngs)        
        req(input$dr_contrast, input$dsea_method)
        
        comparison=1
        names(ngs$gx.meta$meta)
        comparison = input$dr_contrast
        if(is.null(comparison)) return(NULL)
        
        names(ngs$drugs)
        dmethod = "activity-combo/L1000"
        dmethod = "activity/L1000"
        dmethod <- input$dsea_method
        if(is.null(dmethod)) return(NULL)
        
        fc <- ngs$gx.meta$meta[[comparison]]$meta.fx
        names(fc) <- rownames(ngs$gx.meta$meta[[1]])

        nes <- round(ngs$drugs[[dmethod]]$X[,comparison],4)
        pv  <- round(ngs$drugs[[dmethod]]$P[,comparison],4)
        qv  <- round(ngs$drugs[[dmethod]]$Q[,comparison],4)
        drug <- rownames(ngs$drugs[[dmethod]]$X)
        stats <- ngs$drugs[[dmethod]]$stats
        annot <- ngs$drugs[[dmethod]]$annot
        nes[is.na(nes)] <- 0
        qv[is.na(qv)] <- 1
        pv[is.na(pv)] <- 1
        
        ## SHOULD MAYBE BE DONE IN PREPROCESSING???
        if(is.null(annot)) {
            annot <- read.csv(file.path(FILES,"L1000_repurposing_drugs.txt"),
                              sep="\t", comment.char="#")
        }
        
        jj <- match( toupper(drug), toupper(rownames(annot)) )
        annot <- annot[jj,c("moa","target")]        
        res <- data.frame( drug=drug, NES=nes, pval=pv, padj=qv, annot)
        ##res <- res[order(-abs(res$NES)),]
        ##res <- res[order(res$pval,-abs(res$NES)),]
        res <- res[order(-res$NES),]
        
        return(res)
    })

    dsea_enplots.RENDER %<a-% reactive({

        ngs <- inputData()
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))        
        req(input$dr_contrast, input$dsea_method)

        comparison=1
        comparison = input$dr_contrast
        if(is.null(comparison)) return(NULL)

        res <- getDseaTable()
        ## filter with table selection/search
        ii  <- dsea_table$rows_all()
        req(ii)
        if(length(ii)>0) {
            res <- res[ii,,drop=FALSE]
        }
        
        ## rank vector for enrichment plots
        dmethod <- input$dsea_method
        rnk <- ngs$drugs[[dmethod]]$stats[,comparison]
        dctype <- sub("_.*$","",names(rnk))
        ##table(rownames(res) %in% dctype)
        ##table(sapply(rownames(res), function(g) sum(grepl(g,names(rnk),fixed=TRUE))))
        
        ## ENPLOT TYPE
        ##itop <- c( head(order(-res$NES),10), tail(order(-res$NES),10))
        itop <- 1:min(20,nrow(res))
        par(oma=c(0,1,0,0))
        par(mfrow=c(4,5), mar=c(1,1.5,1.8,1))
        i=1
        for(i in itop) {
            dx <- rownames(res)[i]
            dx
            gmtdx <- grep(dx,names(rnk),fixed=TRUE,value=TRUE)  ## L1000 naming
            length(gmtdx)
            ##if(length(gmtdx) < 3) { frame(); next }
            dx1 <- substring(dx,1,26)
            gsea.enplot( rnk, gmtdx, main=dx1, cex.main=1.1, xlab="")
            nes <- round(res$NES[i],2)
            qv  <- round(res$padj[i],3)
            tt <- c( paste("NES=",nes), paste("q=",qv) )
            legend("topright", legend=tt, cex=0.85, y.intersp=0.85)
        }
        
    })    

    dsea_moaplot.RENDER %<a-% reactive({

        ngs <- inputData()
        req(ngs, input$dr_contrast, input$dsea_method)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    

        dbg("[dsea_moaplot.RENDER] reacted")
        
        comparison=1
        comparison = input$dr_contrast
        if(is.null(comparison)) return(NULL)
        
        res <- getDseaTable()
        dbg("[dsea_moaplot.RENDER] dim(res)=",dim(res))
        
        dmethod="combo"
        moatype="target gene"
        names(ngs$drugs)
        dmethod <- input$dsea_method
        moatype <- input$dsea_moatype
        NTOP=12
        dtg.top = moa.top = NULL
        
        if(moatype=="target gene") {
            
            ## GSEA on molecular targets
            targets.list <- lapply(as.character(res$target),
                                   function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]) )
            names(targets.list) <- rownames(res)
            targets <- setdiff(unique(unlist(targets.list)),c(NA,""," "))
            gmt <- lapply(targets, function(g) names(which(sapply(targets.list,function(t) (g %in% t)))))
            names(gmt) <- targets
            
            rnk <- res$NES
            names(rnk) <- rownames(res)
            suppressWarnings( res1 <- fgsea( gmt, rnk, nperm=20000) )
            res1 <- res1[order(res1$pval),]
            head(res1)            
            
            jj <- unique(c(head(order(-res1$NES),NTOP),tail(order(-res1$NES),NTOP)))
            dtg.top <- res1$NES[jj]
            names(dtg.top) <- res1$pathway[jj]
            
        }
        if(moatype=="drug class") {
            
            ## GSEA on moa terms
            moa.list <- lapply(as.character(res$moa),
                               function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]))
            names(moa.list) <- rownames(res)
            moa <- setdiff( unlist(moa.list), c("",NA," "))
            gmt <- lapply(moa, function(g) names(which(sapply(moa.list,function(t) (g %in% t)))))
            names(gmt) <- moa
            
            rnk <- res$NES
            names(rnk) <- rownames(res)
            suppressWarnings( res1 <- fgsea( gmt, rnk, nperm=20000) )
            res1 <- res1[order(res1$pval),]
            head(res1)            
            
            jj <- unique(c(head(order(-res1$NES),NTOP),tail(order(-res1$NES),NTOP)))
            moa.top <- res1$NES[jj]
            names(moa.top) <- res1$pathway[jj]                                    
        }
        
        if(1) {

            ##layout(matrix(1:2,nrow=1),widths=c(1.4,1))
            ##par(mfrow=c(2,1))
            par(mar=c(4,15,5,0.5), mgp=c(2,0.7,0))
            par(mfrow=c(1,1))
            ylab = "enrichment  (NES)"

            if(moatype=="drug class") {
                par(mfrow=c(2,1), mar=c(4,3.8,1,0.2), mgp=c(1.9,0.7,0))
                barplot(moa.top, horiz=FALSE, las=3,
                        ylab=ylab, cex.names = 0.8 )
                ##title(main="MOA", line=1 )
            } else {
                par(mfrow=c(2,1), mar=c(0,3.8,1,0.2), mgp=c(1.9,0.7,0))
                barplot(dtg.top, horiz=FALSE, las=3, ## ylab="drugs (n)",
                        ylab=ylab, cex.names = 0.8 )
                ##title(main="target gene", line=1 )
            }
        }
        
    })    

    ##
    ## TEMPORARY: please merge with dsea_moaplot.RENDER
    ##
    dsea_moaplot.RENDER2 %<a-% reactive({

        ngs <- inputData()
        req(ngs, input$dr_contrast, input$dsea_method)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    

        dbg("[dsea_moaplot.RENDER] reacted")
        
        comparison=1
        comparison = input$dr_contrast
        if(is.null(comparison)) return(NULL)
        
        res <- getDseaTable()
        dbg("[dsea_moaplot.RENDER] dim(res)=",dim(res))
        
        dmethod="combo"
        moatype="target gene"
        names(ngs$drugs)
        dmethod <- input$dsea_method
        moatype <- input$dsea_moatype
        NTOP = 24
        dtg.top = moa.top = NULL
        
        if(moatype=="target gene") {
            
            ## GSEA on molecular targets
            targets.list <- lapply(as.character(res$target),
                                   function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]) )
            names(targets.list) <- rownames(res)
            targets <- setdiff(unique(unlist(targets.list)),c(NA,""," "))
            gmt <- lapply(targets, function(g) names(which(sapply(targets.list,function(t) (g %in% t)))))
            names(gmt) <- targets
            
            rnk <- res$NES
            names(rnk) <- rownames(res)
            suppressWarnings( res1 <- fgsea( gmt, rnk, nperm=20000) )
            res1 <- res1[order(res1$pval),]
            head(res1)            
            
            jj <- unique(c(head(order(-res1$NES),NTOP),tail(order(-res1$NES),NTOP)))
            dtg.top <- res1$NES[jj]
            names(dtg.top) <- res1$pathway[jj]
            
        }
        if(moatype=="drug class") {
            
            ## GSEA on moa terms
            moa.list <- lapply(as.character(res$moa),
                               function(s) trimws(strsplit(s,split="[\\|;,]")[[1]]))
            names(moa.list) <- rownames(res)
            moa <- setdiff( unlist(moa.list), c("",NA," "))
            gmt <- lapply(moa, function(g) names(which(sapply(moa.list,function(t) (g %in% t)))))
            names(gmt) <- moa
            
            rnk <- res$NES
            names(rnk) <- rownames(res)
            suppressWarnings( res1 <- fgsea( gmt, rnk, nperm=20000) )
            res1 <- res1[order(res1$pval),]
            head(res1)            
            
            jj <- unique(c(head(order(-res1$NES),NTOP),tail(order(-res1$NES),NTOP)))
            moa.top <- res1$NES[jj]
            names(moa.top) <- res1$pathway[jj]                                    
        }
        
        if(1) {

            ##layout(matrix(1:2,nrow=1),widths=c(1.4,1))
            ##par(mfrow=c(2,1))
            par(mar=c(4,15,5,0.5), mgp=c(2,0.7,0))
            par(mfrow=c(1,1))
            ylab = "enrichment  (NES)"

            if(moatype=="drug class") {
                par(mfrow=c(2,1), mar=c(4,4,1,0.5), mgp=c(2,0.7,0))
                barplot(moa.top, horiz=FALSE, las=3,
                        ylab=ylab, cex.names = 0.9 )
                ##title(main="MOA", line=1 )
            } else {
                par(mfrow=c(2,1), mar=c(0,4,1,0.5), mgp=c(2,0.7,0))
                barplot(dtg.top, horiz=FALSE, las=3, ## ylab="drugs (n)",
                        ylab=ylab, cex.names = 0.9 )
                ##title(main="target gene", line=1 )
            }
        }
        
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
                          scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "NES",
                                background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })

    dsea_actmap.plotdata <- reactive({
        require(igraph)
        ngs <- inputData()
        req(ngs, input$dr_contrast, input$dsea_method)

        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        if(is.null(ngs$drugs)) return(NULL)
        
        dmethod="activity/L1000";comparison=1
        dmethod <- input$dsea_method        
        comparison = input$dr_contrast
        if(is.null(comparison)) return(NULL)
        
        nes <- ngs$drugs[[dmethod]]$X
        qv  <- ngs$drugs[[dmethod]]$Q
        score <- nes * (1 - qv)**2
        score[is.na(score)] <- 0
        if(NCOL(score)==1) score <- cbind(score,score)
        
        ## reduce score matrix
        ##score = head(score[order(-rowSums(abs(score))),],40)
        ##score = score[head(order(-rowSums(score**2)),50),] ## max number of terms
        score = score[head(order(-score[,comparison]**2),50),,drop=FALSE] ## max number of terms    
        score = score[,head(order(-colSums(score**2)),25),drop=FALSE] ## max comparisons/FC

        cat("dsea_actmap:: dim(score)=",dim(score),"\n")
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        jj=1;ii=1:nrow(score)
        ii <- hclust(d1)$order
        jj <- hclust(d2)$order
        score <- score[ii,jj,drop=FALSE]
        
        cex2=1
        colnames(score) = substring(colnames(score),1,25)
        rownames(score) = substring(rownames(score),1,42)
        if(ncol(score)>15) {
            rownames(score) = substring(rownames(score),1,34)
            cex2=0.85
        }
        if(ncol(score)>25) {
            rownames(score) = substring(rownames(score),1,25)
            colnames(score) <- rep("",ncol(score))
            cex2=0.7
        }

        score2 <- score
        if(input$dr_normalize) score2 <- t( t(score2) / apply(abs(score2),2,max)) 
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**3   ## fudging

        list( score=score2, cex=cex2)
    })        
        
    dsea_actmap.RENDER <- reactive({

        plt <- dsea_actmap.plotdata()

        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,1,0,0))
        require(corrplot)
        corrplot( plt$score, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex = 0.9*plt$cex, tl.col = "grey20", tl.srt = 90)
    })    

    dsea_actmap.RENDER2 <- reactive({

        plt <- dsea_actmap.plotdata()
        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,1,0,0))
        require(corrplot)
        corrplot( t(plt$score), is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex = 0.9*plt$cex, tl.col = "grey20", tl.srt = 90)
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
        pdf.width=11, pdf.height=7,
        height = c(0.54*rowH,650), width=c('auto',1280),
        res = c(72,90)
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
        title = "Mechanism of action", label="c",
        info.text = "This plot visualizes the <strong>mechanism of action</strong> (MOA) across the enriched drug profiles. On the vertical axis, the number of drugs with the same MOA are plotted. You can switch to visualize between MOA or target gene.",
        options = dsea_moaplot.opts,
        pdf.width=4, pdf.height=6,
        height = c(0.54*rowH, 600), width=c('auto',1000),
        res=c(72,100)
    )

    ##-------- Activation map plotting module
    dsea_actmap.opts = tagList()
    callModule(
        plotModule,
        id = "dsea_actmap",
        func = dsea_actmap.RENDER,
        func2 = dsea_actmap.RENDER2, 
        title = "Activation matrix", label="d",
        info.text = "The <strong>Activation Matrix</strong> visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = dsea_actmap.opts,
        pdf.width=6, pdf.height=10,
        height = c(fullH-120,600), width = c('auto',1200),
        res=72
    )

    ##--------buttons for table
    dsea_table <- callModule(
        tableModule,
        id = "dsea_table", label="b",
        func = dsea_table.RENDER, 
        info.text="<b>Enrichment table.</b> Enrichment is calculated by correlating your signature with known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug.", 
        title = "Enrichment table",
        height = c(240,700)
    )
       
    ##-----------------------------------------
    ## Page layout
    ##-----------------------------------------

    dsea_analysis_caption = "<b>(a)</b> <b>Drug connectivity</b> correlates your signature with known drug perturbation profiles from the L1000 database. The figures show the most similar (or opposite) profiles by running the GSEA algorithm on the profile correlation space. <b>(b)</b> <b>Enrichment table</b> summarizing the statistical results of the drug enrichment analysis. <b>(c)</b> <b>Mechanism-of-action</b> plot showing the top most frequent drug class (or target genes) having similar or opposite enrichment compared to the query signature. <b>(d)</b> <b>Activation matrix</b> visualizing enrichment levels of drug signatures across multiple contrast profiles." 

    output$DSEA_analysis_UI <- renderUI({
        fillCol(
            flex = c(NA,0.035,1),
            height = fullH,            
            div(HTML(dsea_analysis_caption),class="caption"),
            br(),
            fillRow(
                height = rowH,
                flex = c(2.8,1), 
                fillCol(
                    flex = c(1.4,0.15,1),
                    height = rowH,
                    fillRow(
                        flex=c(2.0,1),
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
    outputOptions(output, "DSEA_analysis_UI", suspendWhenHidden=FALSE) ## important!!!


}
