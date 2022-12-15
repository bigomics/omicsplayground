##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_abundance_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    menu_grouped = '<code>grouped</code>'
    info_text = paste0('Barplot showing the percentage of counts in terms of major gene types such as ribosomal protein genes, kinases or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Abundance of major gene types",
        label = label,
        plotlib = "plotly",
        info.text = info_text,
        options = NULL,
        download.fmt=c("png","pdf","csv"),         
        width = c("auto","100%"),
        height = height
    )
    
}

dataview_plot_abundance_server <- function(id, 
                                           getCountsTable,
                                           watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {


        ## extract data from pgx object
        plot_data  <- shiny::reactive({
            dbg("[dataview_counts_histplot_server:plot_data] reacted!")
            res <- getCountsTable()
            shiny::req(res)
            res <- list(
                prop.counts = res$prop.counts
            )
            res
        })

        plot.RENDER <- function() {
            
            res <- plot_data()
            shiny::req(res)            
            
            klr <- colorRampPalette(
                c(rgb(0.2,0.5,0.8,0.8),
                  rgb(0.2,0.5,0.8,0.1)), alpha = TRUE)(nrow(res$prop.counts))
        
            ymax = max(colSums(res$prop.counts, na.rm=TRUE))
            names.arg = colnames(res$prop.counts)
            if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)) }
            cex.names <- ifelse(length(names.arg)>10,0.8,0.9)

            par(mar=c(8,3.5,2,0), mgp=c(2.2,0.8,0))
            barplot(res$prop.counts, las=3,
                    cex.lab=1.0, border = NA,
                    ylim = c(0,ymax)*1.6, ylab = "abundance (%)",
                    names.arg = names.arg, cex.names = cex.names,
                    col = klr)
            leg <- legend("topleft", legend=rev(rownames(res$prop.counts)),
                          fill=rev(klr),cex=1, y.intersp=0.75, bty="n", plot = FALSE)
            leftx  <- leg$rect$left*0.9
            rightx <- leg$rect$right*0.9
            topy <- leg$rect$top
            bottomy <- leg$rect$bottom
            legend(x = c(leftx, rightx), y = c(topy, bottomy),
                   legend=rev(rownames(res$prop.counts)),
                   fill=rev(klr), bty="n", cex=0.9, y.intersp=0.75)
        }

        modal_plot.RENDER <- function() {
            plot.RENDER()
        }

        plotly.RENDER <- function() {

            res <- plot_data()
            shiny::req(res)            
            
            long.data <- reshape2::melt( head(res$prop.counts,5) )
            colnames(long.data) <- c("gene","sample","value")

            ## stacked barchart
            fig <- 
              plotly::plot_ly(
                data = long.data,
                x = ~sample,
                y = ~value,
                type = 'bar',
                color = ~gene,
                colors = omics_pal_d("muted")(length(unique(long.data$gene))), 
                hovertemplate = ~paste0(
                  "Sample: <b>", sample, "</b><br>",
                  "Gene: <b>", gene, "</b><br>",
                  "Cum. proportion: <b>", sprintf("%2.1f", value), "%</b>",
                  "<extra></extra>"
                )
              ) %>% 
              plotly::layout(
                barmode = 'stack',
                xaxis = list(title = FALSE),
                yaxis = list(title = "Cumulative proportion", ticksuffix = "%"),
                font = list(family = "Lato"),
                margin = list(l = 10, r = 10, b = 10, t = 10)   
              ) %>% 
              plotly_default1()
            fig
            
        }
        
        modal_plotly.RENDER <- function() {
            fig <- plotly.RENDER() %>%
                plotly::layout(
                    showlegend = TRUE, ## TODO: really TRUE here?
                    font = list(
                        size = 18
                    )
                )
            ## fig <- plotly::style(fig, marker.size = 14)
            fig
        }

        PlotModuleServer(
            "pltmod",
            plotlib = "plotly",
            ##plotlib2 = "base",
            func = plotly.RENDER,
            func2 = modal_plotly.RENDER,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            res = c(90,170),                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


