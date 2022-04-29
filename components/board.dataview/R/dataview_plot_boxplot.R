##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_boxplot_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    menu_grouped='<code>grouped</code>'
    info_text = paste0('Boxplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Counts distribution",
        label = label,
        info.text = info_text,
        height = height
    )
    
}

dataview_plot_boxplot_server <- function(id, parent.input, getCountsTable, watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        ## extract data from pgx object
        plot_data  <- shiny::reactive({
            res = getCountsTable()
            req(res)
            list(
                log2counts = res$log2counts
            )
        })
            
        plot.RENDER <- function() {
            res <- plot_data()
            shiny::req(res)
            
            par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
            ## ---- xlab ------ ###
            xaxt="l"
            names.arg = colnames(res$log2counts)
            if( length(names.arg) > 20){
                names.arg = rep("",length(names.arg))
                xaxt = "n"
            }
            
            cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
            boxplot(
                res$log2counts,
                col=rgb(0.2,0.5,0.8,0.4),
                names = names.arg,
                cex.axis = cex.names,
                border = rgb(0.824,0.824,0.824,0.9),
                xaxt = xaxt,
                las = 3,
                cex.lab = 1,
                ylab = "counts (log2)",
                outline = FALSE,
                varwidth = FALSE
            )
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


