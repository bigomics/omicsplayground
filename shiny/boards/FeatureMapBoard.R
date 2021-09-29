##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

message(">>> sourcing FeatureMapBoard")

FeatureMapInputs <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    shiny::tagList(
        shiny::uiOutput(ns("description")),
        shiny::uiOutput(ns("inputsUI"))
    )
}

FeatureMapUI <- function(id) {
    ns <- shiny::NS(id)  ## namespace
    ui <- shiny::fillCol(
        height = 750,
        shiny::tabsetPanel(
            id = ns("tabs"),
            shiny::tabPanel("Gene", shiny::uiOutput(ns("geneUMAP_UI"))),
            shiny::tabPanel("Geneset", shiny::uiOutput(ns("gsetUMAP_UI")))            
        )
    )
    ui
}

FeatureMapBoard <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE

    inputData <- env[["load"]][["inputData"]]

    fullH = 800  ## full height of page
    rowH1 = 220  ## row 1 height
    rowH2 = 460  ## row 2 height
    
    description = "<h3>Feature Maps</h3> Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). It allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers."
    output$description <- shiny::renderUI(shiny::HTML(description))
    
    infotext ="Visually explore and compare expression signatures on UMAP plots. Feature-level clustering is based on pairwise co-expression between genes (or genesets). This is in contrast to sample-level clustering which clusters samples by similarity of their expression profile. Feature-level clustering allows one to detect gene modules, explore gene neighbourhoods, and identify potential drivers, to study the relationships between features.
<br><br>The tabs present Gene Maps and Geneset Maps and are computed for gene and geneset features, respectively. The clustering of features is computed using UMAP from either the normalized log-expression matrix (logCPM) or the log-foldchange matrix (logFC), with the covariance as distance metric. The UMAP from the logCPM is the default, but in cases of strong batch/tissue effects the UMAP from the logFC matrix is a better choice. We prefer the covariance distance metric instead of the correlation because it takes the size of the foldchange into account. Doing so, genes that are close together in corners in the outer rim are those with high pairwise covariance, i.e. have high correlation and high FC.
<br><br>The maps can be colored according to the foldchange signature of the group contrasts (i.e. comparisons), or colored by the average relative log-expression according to some phenotype condition. Multiple signatures can then be easily compared by visually inspection of the colors.
"
    
    ##================================================================================
    ##========================= OUTPUT UI ============================================
    ##================================================================================
    
    umap_caption = "<h4>Gene maps</h4>"
    output$geneUMAP_UI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.02,1,0.3),
            height = 1.1*fullH,
            shiny::div(shiny::HTML(umap_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.03,1.2),
                height = 0.85*fullH,                
                plotWidget(ns('geneUMAP')),
                shiny::br(),                
                plotWidget(ns('geneTopPlots')) ## %>% shinycssloaders::withSpinner()
            ),
            tableWidget(ns('geneTable'))
        )
    })
    ## important for plot options to be updated correctly...
    shiny::outputOptions(output, "geneUMAP_UI", suspendWhenHidden=FALSE) 


    ##----------------------------------------------
    gset_caption = "<h4>Geneset maps</h4>"    
    output$gsetUMAP_UI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.02,1,0.3),
            height = 1.1*fullH,
            shiny::div(shiny::HTML(gset_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.03,1.2),
                plotWidget(ns('gsetUMAP')),
                shiny::br(),
                plotWidget(ns('gsetTopPlots')) ##  %>% shinycssloaders::withSpinner()
            ),
            tableWidget(ns('gsetTable'))
        )

    })
    shiny::outputOptions(output, "gsetUMAP_UI", suspendWhenHidden=FALSE)     
    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    
    output$inputsUI <- shiny::renderUI({
        ui <- shiny::tagList(
            shiny::actionLink(ns("info"), "Info", icon=icon("info-circle")),
            shiny::hr(), shiny::br(),             
            
            ## data set parameters
            shiny::selectInput(ns('sigvar'),'Show phenotype:', choices=NULL, multiple=FALSE),
            shiny::br(),
            shiny::br(),
            shiny::actionLink(ns("options"), "Options", icon=icon("cog", lib = "glyphicon")),
            shiny::br(),br(),
            shiny::conditionalPanel(
                "input.options % 2 == 1", ns=ns,
                shiny::tagList(
                    tipifyR(shiny::radioButtons(ns('umap_type'),'UMAP datatype:',
                                         choices=c('logCPM','logFC'), inline=TRUE),
                            "The UMAP can be computed from the normalized log-expression (logCPM), or from the log-foldchange matrix (logFC). Clustering based on logCPM is the default, but when batch/tissue effects are present the logFC might be better."),
                    tipifyR(shiny::selectInput(ns('filter_genes'),'Show genes:',
                                        choices=NULL, multiple=FALSE),
                           "Filter the genes to highlight on the map."),
                    tipifyR(shiny::selectInput(ns('filter_gsets'),'Show genesets:',
                                        choices=NULL, multiple=FALSE),
                            "Filter the genesets to highlight on the map.")
                )
            )
        )
        
    })
    shiny::outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!
    
    ##================================================================================
    ##======================= PRECOMPUTE FUNCTION ====================================
    ##================================================================================

    if(0) {
        load("~/Playground/omicsplayground/data/geiger2016-arginine-test.pgx")
        GSETxGENE <- readRDS("~/Playground/omicsplayground/lib/gset-sparseG-XL.rds")
    }
            
    
    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================

    shiny::observeEvent( input$info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Feature Map Analysis</strong>"),
            shiny::HTML(infotext),
            easyClose = TRUE ))
    })

    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        
        families <- names(FAMILIES)
        shiny::updateSelectInput(session, "filter_genes", choices=families,
                          selected = '<all>')
        
        gsetcats <- sort(unique(gsub(":.*","",rownames(ngs$gsetX))))
        shiny::updateSelectInput(session, "filter_gsets", choices=gsetcats,
                          selected = 'H')        

        cvar <- pgx.getCategoricalPhenotypes(ngs$samples, max.ncat=99)
        xpheno <- c("<foldchange>", cvar)
        shiny::updateSelectInput(session, "sigvar", choices=xpheno,
                          selected = "<foldchange>")        
        
    })
    
    ##================================================================================
    ##============================= FUNCTIONS ========================================
    ##================================================================================
    
    ##hilight=hilight2=NULL;source="";plotlib='base';cex=0.9
    plotUMAP <- function(pos, var, hilight=NULL, nlabel=20,  title="",
                         cex=0.9, source="", plotlib='base')
    {        

        dbg("[plotUMAP] called")

        if(!is.null(hilight)) {
            hilight <- intersect(hilight,rownames(pos))
            hilight <- intersect(hilight,names(var))            
            hilight <- hilight[order(-var[hilight])]        
            if(min(var,na.rm=TRUE)<0) {
                hilight2 <- c(head(hilight,nlabel/2),tail(hilight,nlabel/2))
                hilight2 <- unique(hilight2)
            } else {
                hilight2 <- head(hilight,nlabel)
            }
        }
        if(length(hilight) > 0.33*length(var)) hilight <- hilight2

        cexlab  = ifelse(length(hilight2) <= 20, 1, 0.85)
        cexlab  = ifelse(length(hilight2) <= 8, 1.15, cexlab)
        opacity = ifelse(length(hilight2)>0, 0.4, 0.90)
        ##cex = 0.9
        ## opacity = ifelse(length(hilight)>0, 0.15, 1)
        if(plotlib=='plotly') opacity <- sqrt(opacity) ## less opacity..
        
        p <- pgx.scatterPlotXY(
            pos,
            var = var,
            plotlib = plotlib,            
            softmax = TRUE,
            cex.lab = 1.3*cexlab,
            opacity = opacity,
            cex = cex,
            zsym = (min(var,na.rm=TRUE)<0),
            hilight.cex = cex,
            ##hilight.col = 'red',
            hilight.col = NULL,
            hilight.lwd = 0.8,
            hilight = hilight,
            hilight2 = hilight2,
            title = title
            ##legend.pos = 'bottomright',
            ## source = source,
            ## key = rownames(pos)            
        )

        if(0) {
            p <- p %>%
                plotly::event_register('plotly_selected') %>%
                ## plotly::config(displayModeBar = TRUE) %>%
                plotly::layout(dragmode= 'select')
        }
        p 
    }
    
    plotFeaturesPanel <- function(pos, F, ntop, nr, nc, sel, progress) {        

        dbg("[plotFeatures] called")

        par(mar=c(1.6,1.5,0.5,0), oma=c(1,1,0,0)*2)
        par(mar=c(1.1,1.0,0.5,0), oma=c(1,1,0,0)*2)
        par(mgp=c(1.35,0.5,0), las=0, cex.axis=0.85, cex.lab=0.9, xpd=TRUE)
        ##ntop <- 10        
        cex = ifelse(nc>3, 0.5, 0.7)
        jj <- 1:nrow(F)
        if(ncol(F)>4 && nrow(F)>8000) jj <- sample(1:nrow(F),8000)  ## subsample for speed
        if(ncol(F)>9 && nrow(F)>4000) jj <- sample(1:nrow(F),4000)  ## subsample for speed
        if(ncol(F)>16 && nrow(F)>2000) jj <- sample(1:nrow(F),2000)  ## subsample for speed

        i=1
        for(i in 1:ncol(F)) {                            
            if(!interactive()) progress$inc(1/ncol(F))            
            var <- F[,i]
            var <- var[match(rownames(pos),names(var))]                
            zsym <- ifelse(min(var,na.rm=TRUE)>=0, FALSE, TRUE)
            hmarks <- NULL
            if(!is.null(sel)) {
                hmarks <- intersect(sel,names(var))
                hmarks <- head(hmarks[order(var[hmarks])],ntop)
            }
            opacity = ifelse(is.null(hmarks),0.9,0.4)
            xlab = "UMAP-x"
            ylab = "UMAP-y"
            xaxs = yaxs = TRUE
            if(i%%nc != 1) {
                ylab = ''
                yaxs = FALSE
            }
            if((i-1)%/%nc < (nr-1)) {
                xlab = ''
                xaxs = FALSE                
            }

            pgx.scatterPlotXY.BASE(
                pos[jj,], var=var[jj], 
                zsym=zsym, set.par=FALSE, softmax=1,
                cex=cex, cex.legend = 0.9, cex.lab=1.2, bty='n',
                col='grey70', dlim=c(0.05,0.05),
                hilight=hmarks, hilight2=NULL,
                hilight.col=NULL, opacity=opacity,
                ##xlab = xlab, ylab = ylab,
                xlab = '', ylab = '',
                xaxs = xaxs, yaxs = yaxs,
                hilight.lwd=0.5, hilight.cex=1.3)

            cex1 <- ifelse(ncol(F) <= 16, 1.2, 1)
            title(colnames(F)[i], cex.main=cex1, line= -0.75)                                
            mtext(xlab,1,line=1.5,cex=0.6)
            mtext(ylab,2,line=1.6,cex=0.6)            

        }               
    }

    getGeneUMAP_FC <- shiny::reactive({
        ## buffered reactive
        ngs <- inputData()
        shiny::withProgress({            
            F <- pgx.getMetaMatrix(ngs, level='gene')$fc
            F <- scale(F, center=FALSE)
            pos <- pgx.clusterBigMatrix(t(F), methods='umap', dims=2)[[1]]
            pos <- pos.compact(pos)
        }, message="computing foldchange UMAP", value=0.5 )            
        pos
    })
       
    getGeneUMAP <- shiny::reactive({
        ngs <- inputData()
        if(input$umap_type=='logFC') {
            message("[getGeneUMAP] computing foldchange UMAP")
            pos <- getGeneUMAP_FC()
        } else {
            pos <- ngs$cluster.genes$pos[['umap2d']]
        }
        pos
    })
    
    getGsetUMAP_FC <- shiny::reactive({
        ## buffered reactive
        ngs <- inputData()
        shiny::withProgress({
            F <- pgx.getMetaMatrix(ngs, level='geneset')$fc
            F <- scale(F, center=FALSE)
            dbg("getGsetUMAP_FC:: dim.F = ",dim(F))
            pos <- pgx.clusterBigMatrix(t(F), methods='umap', dims=2)[[1]]
            pos <- pos.compact(pos)
        }, message="computing foldchange UMAP (genesets)", value=0.5 )
        pos
    })

    getGsetUMAP <- shiny::reactive({
        ngs <- inputData()
        if(input$umap_type=='logFC') {
            message("[getGsetUMAP] computing foldchange UMAP (genesets)")            
            pos <- getGsetUMAP_FC()
        } else {
            pos <- ngs$cluster.gsets$pos[['umap2d']]
        }        
        pos
    })
        
    ## ================================================================================
    ## ========================= PLOTTING MODULES =====================================
    ## ================================================================================
    
    ##geneUMAP.RENDER %<a-% shiny::reactive({
    geneUMAP.RENDER <- shiny::reactive({            

        dbg("[geneUMAP.RENDER] reacted")
        ngs <- inputData()
        shiny::req(ngs)
        
        pos <- getGeneUMAP()
        hilight <- NULL
        colgamma <- as.numeric(input$umap_gamma)
        
        ## select on table filter
        F  <- pgx.getMetaMatrix(ngs)$fc
        F <- scale(F,center=FALSE)
        colorby <- input$umap_colorby
        if(colorby=='sd.FC') {
            fc <- (rowMeans(F**2))**0.2
        } else if(colorby=='mean.FC') {
            fc <- rowMeans(F)
        } else {
            ## sdX
            cX <- ngs$X - rowMeans(ngs$X,na.rm=TRUE)
            fc <- sqrt(rowMeans(cX**2))
        }
        fc <- sign(fc)*abs(fc/max(abs(fc)))**colgamma
        hilight <- names(sort(-fc))

        ## filter on table
        hilight <- selGenes()
        nlabel <- as.integer(input$umap_nlabel)
        
        p <- plotUMAP(pos, fc, hilight, nlabel=nlabel, title=colorby,
                      cex=0.9, source="", plotlib='base')        
        p
    })

    
    geneUMAP.RENDER2 <- shiny::reactive({        

        dbg("[geneUMAP.RENDER] reacted")
        ngs <- inputData()
        shiny::req(ngs)
        
        ##pos <- ngs$cluster.genes$pos[['umap2d']]
        pos <- getGeneUMAP()        
        hilight <- NULL
        colgamma <- as.numeric(input$umap_gamma)
        
        ## select on table filter
        F <- pgx.getMetaMatrix(ngs)$fc
        F <- scale(F,center=FALSE)
        colorby <- input$umap_colorby
        if(colorby=='var.FC') {
            fc <- (rowMeans(F**2))**0.2
        } else if(colorby=='mean.FC') {
            fc <- rowMeans(F)
        } else {
            cX <- ngs$X - rowMeans(ngs$X,na.rm=TRUE)
            fc <- sqrt(rowMeans(cX**2))
        }
        fc <- sign(fc)*abs(fc/max(abs(fc)))**colgamma

        ## filter on table
        hilight <- selGenes()
        nlabel <- as.integer(input$umap_nlabel)
        
        p <- plotUMAP(pos, fc, hilight, nlabel=nlabel, title=colorby,
                      cex=1.2, source="", plotlib='plotly')        
        p
    })

    
    geneUMAP_info = "<b>Gene mapping.</b> Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis)."

    geneUMAP.opts <- shiny::tagList(
        shiny::selectInput(ns('umap_nlabel'),'nr labels:',
                    c(0,10,20,50,100,1000), selected=50),
        shiny::sliderInput(ns('umap_gamma'),'color gamma:',
                    min=0.2, max=2, value=1, step=0.2),
        shiny::radioButtons(ns('umap_colorby'),'color by:',
                     choices = c("var.FC","mean.FC","sd.X"),
                     selected = "mean.FC", inline=TRUE )
    )

    shiny::callModule(
        plotModule, 
        id = "geneUMAP", ##ns=ns,
        ##plotlib = 'plotly',
        ##plotlib = 'ggplot',
        title = "GENE MAP", label="a",
        func = geneUMAP.RENDER,
        ##outputFunc = "function(x,...) shiny::plotOutput(x,brush='geneUMAP_brush',...)",        
        outputFunc = sub("XXX",ns("geneUMAP_brush"),"function(x,...)plotOutput(x,brush='XXX',...)"),
        func2 = geneUMAP.RENDER2, 
        plotlib2 = 'plotly',
        download.fmt = c("png","pdf"),
        options = geneUMAP.opts,
        info.text = geneUMAP_info,        
        height = c(600, 750), width = c('auto',1200),
        pdf.width=10, pdf.height=8, res=c(72,100),
        add.watermark = WATERMARK
    )
    
    ## getSelectedGenes <- shiny::reactive({
    ##     ev <- plotly::event_data("plotly_selected", source='geneUMAP')
    ##     sel <- ev$key
    ##     ##message('[getSelectedGenes] sel = ', paste(sel,collapse=' '))
    ##     sel
    ## })
    
    selGenes <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)        
        sel <- input$filter_genes
        selgenes <- FAMILIES[[sel]]
        selgenes
    }) ## %>% shiny::debounce(1000)
    

    geneTopPlots.RENDER %<a-% shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)

        ##pos <- ngs$cluster.genes$pos[['umap2d']]
        pos <- getGeneUMAP()

        pheno='tissue'
        pheno <- input$sigvar
        if(pheno %in% colnames(ngs$samples)) {
            X <- ngs$X - rowMeans(ngs$X)
            y <- ngs$samples[,pheno]
            F <- do.call(cbind,tapply(1:ncol(X),y,function(i)
                rowMeans(X[,i,drop=FALSE])))
        } else {
            F <- pgx.getMetaMatrix(ngs, level='gene')$fc
        }
        
        dbg("[geneTopPlots.RENDER] dim.F = ",dim(F))
        if(nrow(F)==0) return(NULL)
        
        ## ntop=15
        nc = ceiling(sqrt(ncol(F)))
        nr = ceiling(ncol(F)/nc)
        nr
        nc
        nr2 <- ifelse(nr <= 3, nc, nr)
        par(mfrow=c(nr2,nc), mar=c(2,1,1,0), mgp=c(1.6,0.55,0), las=0)
        progress <- NULL
        if(!interactive())  {
            progress <- shiny::Progress$new()
            on.exit(progress$close())    
            progress$set(message = "Computing feature plots...", value = 0)
        }
        plotFeaturesPanel(pos, F, ntop=ntop, nr, nc, sel=NULL, progress)                
        
    })

    
    geneTopPlots_info = "<b>Gene signature maps.</b> Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis)."
    
    geneTopPlots.opts <- shiny::tagList(
        ##radioButtons(ns('gene_plottype'),'plot type:',c("umap","bar"),inline=TRUE)        
    )

    shiny::callModule(
        plotModule, 
        id = "geneTopPlots", ##ns=ns,
        ##plotlib = 'plotly',
        ##plotlib = 'ggplot',
        title = "GENE SIGNATURES", label="b",
        func = geneTopPlots.RENDER,
        func2 = geneTopPlots.RENDER, 
        download.fmt = c("png","pdf"),
        options = geneTopPlots.opts,
        info.text = geneTopPlots_info,        
        height = c(600, 750), width = c('auto',1200),
        pdf.width=11, pdf.height=9,
        res=c(80,90),
        add.watermark = WATERMARK
    )
        
    ##-----------------------------------------------------------------
    ##----------------------  Geneset UMAP ----------------------------
    ##-----------------------------------------------------------------
        
    gsetUMAP.RENDER %<a-% shiny::reactive({        

        dbg("[gsetUMAP.RENDER] reacted")
        ngs <- inputData()

        ##pos <- ngs$cluster.gsets$pos[['umap2d']]
        pos <- getGsetUMAP()
        hilight <- NULL
        colgamma <- as.numeric(input$gsmap_gamma)
        
        F <- pgx.getMetaMatrix(ngs, level='geneset')$fc
        F <- scale(F,center=FALSE)
        colorby <- input$gsmap_colorby
        if(colorby=='sd.FC') {
            fc <- (rowMeans(F**2))**0.2
        } else if(colorby=='mean.FC') {
            fc <- rowMeans(F)
        } else {
            cX <- ngs$gsetX - rowMeans(ngs$gsetX,na.rm=TRUE)
            fc <- sqrt(rowMeans(cX**2))
        }
        fc <- sign(fc)*abs(fc/max(abs(fc)))**colgamma

        ## filter on table
        hilight <- selGsets()
        nlabel <- as.integer(input$gsmap_nlabel)
        
        par(mfrow=c(1,1))
        p <- plotUMAP(pos, fc, hilight, nlabel=nlabel, title=colorby,
                      cex=0.9, source="", plotlib='base')        
        p

    })

    gsetUMAP.RENDER2 %<a-% shiny::reactive({        

        dbg("[gsetUMAP.RENDER] reacted")
        ngs <- inputData()

        ##pos <- ngs$cluster.gsets$pos[['umap2d']]
        pos <- getGsetUMAP()        
        hilight <- NULL
        colgamma <- as.numeric(input$gsmap_gamma)
        
        F <- pgx.getMetaMatrix(ngs, level='geneset')$fc
        F <- scale(F,center=FALSE)
        colorby <- input$gsmap_colorby
        if(colorby=='var.FC') {
            fc <- (rowMeans(F**2))**0.2
        } else if(colorby=='mean.FC') {
            fc <- rowMeans(F)
        } else {
            cX <- ngs$gsetX - rowMeans(ngs$gsetX,na.rm=TRUE)
            fc <- sqrt(rowMeans(cX**2))
        }
        fc <- sign(fc)*abs(fc/max(abs(fc)))**colgamma

        ## filter on table
        hilight <- selGsets()
        nlabel <- as.integer(input$gsmap_nlabel)
        
        par(mfrow=c(1,1))
        p <- plotUMAP(pos, fc, hilight, nlabel=nlabel, title=colorby,
                      cex=1.2, source="", plotlib='plotly')        
        p

    })
        
    ## getSelectedGsets <- shiny::reactive({
    ##     ev <- plotly::event_data("plotly_selected", source='gsetUMAP')
    ##     sel <- ev$key
    ##     ##message('[getSelectedGenes] sel = ', paste(sel,collapse=' '))
    ##     sel
    ## })
    
    gsetUMAP_info = "<b>Geneset maps.</b> Analysis of network topology for various soft-thresholding powers. The left panel shows the scale-free fit index (y-axis) as a function of the soft-thresholding power (x-axis). The right panel displays the mean connectivity (degree, y-axis) as a function of the soft-thresholding power (x-axis)."

    gsetUMAP.opts <- shiny::tagList(
        shiny::selectInput(ns('gsmap_nlabel'),'nr labels:',
                    choices=c(0,10,20,50,100,1000),selected=20),        
        shiny::sliderInput(ns('gsmap_gamma'),'color gamma:',
                    min=0.2, max=2, value=1, step=0.2),
        shiny::radioButtons(ns('gsmap_colorby'),'color by:',
                     choices = c("sd.FC","mean.FC","sd.X"),
                     selected = "sd.FC", inline=TRUE )
    )

    shiny::callModule(
        plotModule, 
        id = "gsetUMAP", ##ns=ns,
        ##plotlib = 'plotly',
        title = "GENESET UMAP", label="a",
        func = gsetUMAP.RENDER,
        func2 = gsetUMAP.RENDER2,
        outputFunc = sub("XXX",ns("gsetUMAP_brush"),"function(x,...)plotOutput(x,brush='XXX',...)"),
        plotlib2 = 'plotly',
        download.fmt = c("png","pdf"),
        options = gsetUMAP.opts,
        info.text = gsetUMAP_info,        
        height = c(600, 750), width = c('auto',1200),
        pdf.width=8, pdf.height=8,
        res = c(72,80),
        add.watermark = WATERMARK
    )
    
    ##-----------------------------------------------------------------
    ##----------------------  Enrichment plots ------------------------
    ##-----------------------------------------------------------------       

    selGsets <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        db <- input$filter_gsets
        gsets <- rownames(ngs$gsetX)
        gsets <- grep(paste0("^",db,":"),gsets,value=TRUE)
        gsets
    }) ## %>% shiny::debounce(1000)
    
    
    gsetTopPlots.RENDER %<a-% shiny::reactive({        

        dbg("[gsetTopPlots.RENDER] reacted")
        ngs <- inputData()
        shiny::req(ngs)
        
        ## pos <- ngs$cluster.gsets$pos[['umap2d']]
        pos <- getGsetUMAP()        
        hilight <- NULL

        pheno='tissue'
        pheno <- input$sigvar
        if(pheno %in% colnames(ngs$samples)) {
            X <- ngs$gsetX - rowMeans(ngs$gsetX)
            y <- ngs$samples[,pheno]
            F <- do.call(cbind,tapply(1:ncol(X),y,function(i)
                rowMeans(X[,i,drop=FALSE])))
        } else {
            F <- pgx.getMetaMatrix(ngs, level='geneset')$fc
        } 
        if(nrow(F)==0) return(NULL)        

        ntop = 15
        nc = ceiling(sqrt(ncol(F)))
        nr = ceiling(ncol(F)/nc)
        nr2 <- ifelse(nr <= 3, nc, nr)        
        par(mfrow=c(nr2,nc), mar=c(3,2,1.5,2), mgp=c(1.6,0.6,0))
        par(mfrow=c(nr2,nc), mar=c(3,1,1,0.5), mgp=c(1.6,0.55,0))        
        progress <- NULL
        if(!interactive())  {
            progress <- shiny::Progress$new()
            on.exit(progress$close())    
            progress$set(message = "Computing feature plots...", value = 0)
        }
        plotFeaturesPanel(pos, F, ntop, nr, nc, sel=NULL, progress)                
        
    })
    
    gsetTopPlots.info = "<b>Module enrichment.</b> "
    gsetTopPlots.opts = shiny::tagList(
        ## shiny::radioButtons(ns('gset_plottype'),'plot type:',c("umap","bar"),inline=TRUE)        
    )
    
    shiny::callModule(
        plotModule, 
        id = "gsetTopPlots", ##ns=ns,
        title = "GENESET SIGNATURES", label="b",
        func  = gsetTopPlots.RENDER,
        func2 = gsetTopPlots.RENDER, 
        download.fmt = c("png","pdf"),
        options = gsetTopPlots.opts,
        info.text = gsetTopPlots.info,        
        height = c(600, 750),
        width = c('auto',1200),
        pdf.width = 11, pdf.height = 9, 
        res = c(80,90),
        add.watermark = WATERMARK
    )    


    ##-----------------------------------------------------------------
    ##-------------------------  Tables -------------------------------
    ##-----------------------------------------------------------------       

    geneTable.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        if(is.null(ngs$drugs)) return(NULL)
        
        pos <- getGeneUMAP()

        ## detect brush
        sel.genes <- NULL
        ##b <- input$ftmap-geneUMAP_brush  ## ugly??
        dbg("[geneTable.RENDER] names.input = ",names(input))        
        b <- input[['geneUMAP_brush']]  ## ugly??        
        dbg("[geneTable.RENDER] is.null(b) = ",is.null(b))
        dbg("[geneTable.RENDER] length(b) = ",length(b))
        dbg("[geneTable.RENDER] names(b) = ",names(b))
        
        if(!is.null(b) & length(b)) {
            sel <- which( pos[,1] > b$xmin & pos[,1] < b$xmax &
                          pos[,2] > b$ymin & pos[,2] < b$ymax )
            sel.genes <- rownames(pos)[sel]
        }
        dbg("[geneTable.RENDER] sel.genes = ",sel.genes)
        
        pheno='tissue'
        pheno <- input$sigvar
        is.fc=FALSE
        if(pheno %in% colnames(ngs$samples)) {
            X <- ngs$X - rowMeans(ngs$X)
            y <- ngs$samples[,pheno]
            F <- do.call(cbind,tapply(1:ncol(X),y,function(i)
                rowMeans(X[,i,drop=FALSE])))
            is.fc=FALSE            
        } else {
            F <- pgx.getMetaMatrix(ngs, level='gene')$fc
            is.fc=TRUE
        }
        
        if(!is.null(sel.genes)) {
            sel.genes <- intersect(sel.genes,rownames(F))
            F <- F[sel.genes,]
        }
        F <- F[order(-rowMeans(F**2)),]
        
        tt <- shortstring(ngs$genes[rownames(F),'gene_title'],60)
        F <- cbind(sd.X = sqrt(rowMeans(F**2)), F)
        if(is.fc) colnames(F)[1] = "sd.FC"
        F  <- round(F, digits=3)
        df <- data.frame(gene=rownames(F), title=tt, F, check.names=FALSE)
        
        DT::datatable( df, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=NULL),
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'Blfrtip', buttons = c('copy','csv','pdf'),
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = '70vh',
                          scroller = TRUE,
                          deferRender = TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') ## %>% 
                ## DT::formatStyle( "NES",
                ##                 background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                ##                 backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                ##                 backgroundPosition = 'center') 
    })

    
    ##--------buttons for table
    geneTable.opts = shiny::tagList()  
    geneTable <- shiny::callModule(
        tableModule,
        id = "geneTable", label="c",
        func = geneTable.RENDER, 
        options = geneTable.opts,
        info.text="<b>Gene table.</b>", 
        title = "Gene table",
        height = c(280,750), width=c('auto',1400)
    )

    ##--------------------------------------------------------------------------------
    ##---------------------------- GENESET -------------------------------------------
    ##--------------------------------------------------------------------------------    
    
    gsetTable.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        if(is.null(ngs$drugs)) return(NULL)
        
        pos <- getGsetUMAP()

        ## detect brush
        sel.gsets <- NULL
        ##b <- input$ftmap-geneUMAP_brush  ## ugly??
        dbg("[gsetTable.RENDER] names.input = ",names(input))        
        b <- input[['gsetUMAP_brush']]  ## ugly??        
        dbg("[gsetTable.RENDER] is.null(b) = ",is.null(b))
        dbg("[gsetTable.RENDER] length(b) = ",length(b))
        dbg("[gsetTable.RENDER] names(b) = ",names(b))
        
        if(!is.null(b) & length(b)) {
            sel <- which( pos[,1] > b$xmin & pos[,1] < b$xmax &
                          pos[,2] > b$ymin & pos[,2] < b$ymax )
            sel.gsets <- rownames(pos)[sel]
        }
        dbg("[gsetTable.RENDER] sel.gsets = ",sel.gsets)
        
        pheno='tissue'
        pheno <- input$sigvar
        is.fc=FALSE
        if(pheno %in% colnames(ngs$samples)) {
            X <- ngs$gsetX - rowMeans(ngs$gsetX)
            y <- ngs$samples[,pheno]
            F <- do.call(cbind,tapply(1:ncol(X),y,function(i)
                rowMeans(X[,i,drop=FALSE])))
            is.fc=FALSE            
        } else {
            F <- pgx.getMetaMatrix(ngs, level='geneset')$fc
            is.fc=TRUE
        }
        
        if(!is.null(sel.gsets)) {
            sel.gsets <- intersect(sel.gsets,rownames(F))
            F <- F[sel.gsets,]
        }
        F <- F[order(-rowMeans(F**2)),]
        
        ## tt <- shortstring(ngs$genes[rownames(F),'gene_title'],60)
        F <- cbind(sd.X = sqrt(rowMeans(F**2)), F)
        if(is.fc) colnames(F)[1] = "sd.FC"
        F  <- round(F, digits=3)
        gs <- substring(rownames(F),1,100)
        df <- data.frame(geneset=gs, F, check.names=FALSE)
        
        DT::datatable( df, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=NULL),
                      fillContainer = TRUE,
                      options=list(
                          ##dom = 'Blfrtip', buttons = c('copy','csv','pdf'),
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = '70vh',
                          scroller = TRUE,
                          deferRender = TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') ## %>% 
                ## DT::formatStyle( "NES",
                ##                 background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                ##                 backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                ##                 backgroundPosition = 'center') 
    })

    ##--------buttons for table
    gsetTable.opts = shiny::tagList(
    )  
    gsetTable <- shiny::callModule(
        tableModule,
        id = "gsetTable", label="c",
        func = gsetTable.RENDER, 
        options = gsetTable.opts,
        info.text="<b>Gene table.</b>", 
        title = "Geneset table",
        height = c(280,750), width=c('auto',1400)
    )
    

    
    ## ========================================================================
    ## ========================================================================
    ## ========================================================================

} ## end of Board



if(0) {

    x <- 1:20
    y <- runif(20)
    KEGGgraph::plot(x, y, axes=FALSE, frame.plot=TRUE)
    Axis(side=1, labels=FALSE)
    Axis(side=2, labels=TRUE)
    
}
    
    
