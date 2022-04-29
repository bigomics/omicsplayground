##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_averagecounts_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    menu_grouped = '<code>grouped</code>'
    info_text = paste0('Barplot showing the average count levels of major gene types such as CD molecules, kinanses or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Average count by gene type",
        label = label,
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,        
        info.text = info_text,
        options = NULL,
        download.fmt=c("png","pdf","csv"),         
        width = c("auto","100%"),
        height = height
    )
    
}

dataview_plot_averagecounts_server <- function(id, 
                                               getCountsTable,
                                               watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        ## extract data from pgx object
        plot_data  <- shiny::reactive({
            dbg("[dataview_counts_histplot_server:plot_data] reacted!")
            res <- getCountsTable()
            shiny::req(res)
            list(
                avg.counts = res$avg.counts
            )
        })        

        plot.RENDER <- function() {
            
            res <- plot_data()
            shiny::req(res)
            
            par(mar=c(8,3.5,2,0.5), mgp=c(2.2,0.8,0))

            klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)),
                                    alpha = TRUE)(nrow(res$avg.counts))
            
            names.arg = colnames(res$avg.counts)
            if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)) }
            cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
            ## ---- xlab ------ ###
            ymax = max(colSums(res$avg.counts, na.rm=TRUE))
            barplot(res$avg.counts, las=3, #main="average counts by gene type", cex.main=1.6,
                                        #cex.names=res$cx1+0.04,
                    border=NA, cex.lab=1.0,
                    names.arg=names.arg, cex.names=cex.names,
                    ylim=c(0,ymax)*1.6, ylab="average count", col=klr)
            leg <- legend("topleft", legend=rev(rownames(res$avg.counts)),
                          fill=rev(klr), cex=1, y.intersp=0.75, bty="n", plot = FALSE)
            leftx <- leg$rect$left*0.9
            rightx <- leg$rect$right*0.9
            topy <- leg$rect$top
            bottomy <- leg$rect$bottom
            legend(x = c(leftx, rightx), y = c(topy, bottomy),
                   legend=rev(rownames(res$avg.counts)),
                   fill=rev(klr), cex=0.9, y.intersp=0.75, bty="n")
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
            res = c(90,170),                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


