##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_histogram_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    info_text = paste0('Histogram of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the <code>grouped</code> setting uder the main <i>Options</i>.')
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Counts histogram",
        label = label,
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,        
        info.text = info_text,
        options = NULL,
        download.fmt = c("png","pdf","csv"),         
        width = c("auto","100%"),
        height = height
    )
    
}

dataview_plot_histogram_server <- function(id, getCountsTable, watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        gx.hist <- function(gx, n=1000, main="",ylim=NULL, plot=TRUE) {
            jj <- 1:nrow(gx)
            if(length(jj)>n) jj <- sample(jj,n,replace=TRUE)
            h0 <- hist(as.vector(c(gx[jj],min(gx),max(gx))),
                       breaks = 120,
                       plot = plot,
                       main = main,
                       border = FALSE,
                       col = "grey",
                       freq = FALSE, ## ylim=ylim,
                       xlim = c(min(gx),max(gx)),
                       xlab = "expression (log2)",
                       cex.lab = 1)
            i = 1
            H <- c()
            for(i in 1:ncol(gx)) {
                h1 <- hist(gx[jj,i], breaks=h0$breaks,plot=FALSE)
                lines( h0$mids, h1$density, col="black", lwd=0.5 )
                H <- cbind(H, h1$density)
            }
            colnames(H) <- colnames(gx)
            data.frame(mids=h0$mids, density=h0$density, H)
        }

        ## extract data from pgx object
        plot_data  <- shiny::reactive({
            dbg("[dataview_counts_histplot_server:plot_data] reacted!")
            res <- getCountsTable()
            shiny::req(res)
            hh <- gx.hist(gx=res$log2counts, n=2000, plot=FALSE)
            list(
                histogram = hh,
                log2counts = res$log2counts
            )
        })
        
        plot.RENDER <- function() {
            res <- plot_data()
            shiny::req(res)
            par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
            hh <- gx.hist(gx=res$log2counts, n=2000) #, main="histogram")
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
            res = c(90,170)*1,                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


