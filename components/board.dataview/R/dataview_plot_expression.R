##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_expression_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    info_text = paste0('Expression barplot of grouped samples (or cells) for the gene selected in the <code>Search gene</code> Samples (or cells) in the barplot can be ungrouped by setting the <code>grouped</code> under the main <i>Options</i>.')

    opts <- shiny::tagList(
        shiny::radioButtons(ns('geneplot_type'),'plot type (grouped)', c('bar','violin','box'),
                            inline=TRUE)
    )
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Abundance/expression",
        label = label,
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,        
        info.text = info_text,
        caption = NULL,
        caption2 = NULL,        
        options = opts,
        download.fmt=c("png","pdf","csv"),         
        width = c("auto","1200"),
        height = height
    )
    
}

dataview_plot_expression_server <- function(id, pgxdata, parent.input, watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        dbg("[dataview_expressionplot_server] created!")
        
        plot_data <- shiny::reactive({
           
             dbg("[dataview_expressionplot_server:plot_data] reacted! ")

            pgx <- pgxdata()
            shiny::req(pgx)
            shiny::req(parent.input$data_groupby,
                       parent.input$search_gene,
                       parent.input$data_type)
            
            dbg("[dataview_expressionplot_server:plot_data] calling... ")
            
            search_gene <- parent.input$search_gene
            samples = colnames(pgx$X)
            if(!is.null(parent.input$data_samplefilter)) {
                samples <- selectSamplesFromSelectedLevels(pgx$Y, parent.input$data_samplefilter)
            }
            nsamples = length(samples)
            
            grpvar=1
            grp <- rep(NA,length(samples))
            grpvar <- parent.input$data_groupby
            if(grpvar != "<ungrouped>") {
                grp  = factor(as.character(pgx$Y[samples,grpvar]))
            }

            pp <- rownames(pgx$genes)[match(search_gene, pgx$genes$gene_name)]            
            gx = NULL
            ylab = NULL
            if(parent.input$data_type=="counts") {
                gx = pgx$counts[pp,samples]
                ylab="expression (counts)"
            } else if(parent.input$data_type=="CPM") {
                gx = 2**pgx$X[pp,samples]
                ylab="expression (CPM)"
            } else if(parent.input$data_type=="logCPM") {
                gx = pgx$X[pp,samples]
                ylab="expression (log2CPM)"
            }
            
            pd <- list(
                df = data.frame(
                    x = gx,
                    group = grp,
                    samples = samples
                ),
                geneplot_type = input$geneplot_type,
                data_groupby = parent.input$data_groupby,
                ylab = ylab,
                gene = search_gene
            )
            return(pd)
        })


        plot.RENDER <- function() {
            
            pd  <- plot_data()
            df <- pd[['df']]
            
            par(mar=c(7,3.5,2,1), mgp=c(2.1,0.8,0))
            
            BLUE = rgb(0.2,0.5,0.8,0.8)
            bee.cex = ifelse(length(df$x)>500,0.1,0.2)
            bee.cex = c(0.3,0.1,0.05)[cut(length(df$x),c(0,100,500,99999))]
            
            ##if(parent.input$data_grouped) {
            if(pd$data_groupby != "<ungrouped>") {
                nnchar = nchar(paste(unique(df$group),collapse=''))
                srt = ifelse(nnchar < 20, 0, 35)
                ngrp <- length(unique(df$group))
                cx1 = ifelse( ngrp < 10, 1, 0.8)
                cx1 = ifelse( ngrp > 20, 0.6, cx1)
                if(pd$geneplot_type == 'bar') {
                    gx.b3plot(
                        df$x,
                        df$group,
                        las=3,
                        main=pd$gene,
                        ylab=pd$ylab,
                        cex.main=1,
                        col.main="#7f7f7f",
                        bar=TRUE,
                        border=NA,
                        ## bee = ifelse(length(df$x) < 500,TRUE,FALSE),
                        bee.cex = bee.cex,
                        ## sig.stars=TRUE,
                        ## max.stars=5,
                        xlab="",
                        names.cex=cx1,
                        srt=srt,
                        col = rgb(0.4,0.6,0.85,0.85)
                    )
                } else if(pd$geneplot_type == 'violin') {
                    pgx.violinPlot(
                        df$x,
                        df$group,
                        main = pd$gene,
                        cex.main=1,
                        xlab = '',
                        ylab = ylab,
                        ##vcol = rgb(0.2,0.5,0.8,0.8),
                        vcol = rgb(0.4,0.6,0.85,0.85),
                        srt = srt
                    )
                    
                } else {
                    boxplot(
                        df$x ~ df$group,
                        main = pd$gene,
                        cex.main = 1.0,
                        ylab = pd$ylab,
                        xlab = '',
                        xaxt = 'n',
                        col = rgb(0.4,0.6,0.85,0.85)
                    )
                    yy <- sort(unique(df$group))
                    text(x = 1:length(yy),
                         y = par("usr")[3] - 0.03*diff(range(df$x)),
                         labels = yy,
                         xpd = NA,
                         srt = srt,
                         adj = ifelse(srt==0, 0.5, 0.965),
                         cex = cx1)
                }
            }  else {

                ## plot as bars
                barplot(
                    df$x,
                    col = BLUE,
                    las = 3,
                    cex.names = 0.8,
                    ylab = pd$ylab,
                    xlab = "",
                    main = pd$gene,
                    cex.main = 1,
                    col.main = "#7f7f7f",
                    border = NA,
                    names.arg = rep(NA,length(df$x))
                )

                ## add labels if needed
                nx <-  length(df$x)
                if(nx < 100) {
                    cx1 = ifelse(nx > 20, 0.8, 0.9)
                    cx1 = ifelse(nx > 40, 0.6, cx1)
                    cx1 = ifelse(nx < 10, 1, cx1)
                    text(
                        x = (1:nx-0.5)*1.2,
                        y = -0.04*max(df$x),
                        labels = names(df$x),
                        las = 3,
                        cex = cx1,
                        pos = 2,
                        adj = 0,
                        offset = 0,
                        srt = 45,
                        xpd = TRUE
                    )
                }
            }

        }
        
        modal_plot.RENDER <- function() {
            plot.RENDER()
        }
        
        PlotModuleServer(
            "pltmod",
            plotlib = "base",
            plotlib2 = "base",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            renderFunc = shiny::renderPlot,
            renderFunc2 = shiny::renderPlot,        
            ##renderFunc = shiny::renderCachedPlot,
            ##renderFunc2 = shiny::renderCachedPlot,        
            res = c(96,120)*1,                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


