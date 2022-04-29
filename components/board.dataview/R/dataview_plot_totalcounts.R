##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_totalcounts_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)

    menu_grouped='<code>grouped</code>'
    info_text = paste('Barplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the',menu_grouped, 'setting uder the main <i>Options</i>.')       
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Total counts",
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

dataview_plot_totalcounts_server <- function(id,
                                             getCountsTable,
                                             r.data_groupby = reactive(""),
                                             watermark=FALSE)
{
  moduleServer( id, function(input, output, session) {
    
    plot_data  <- shiny::reactive({
      
      data_groupby <- r.data_groupby()
      
      tbl = getCountsTable()
      req(tbl)
      
      ylab = "counts (million)"
      if(data_groupby != "<ungrouped>") {
        ylab = "average counts (million)"
      }
      
      res <- list(
        total.counts = tbl$total.counts,
        ylab = ylab
      )
      return(res)
    })
    

    plot.RENDER <- function() {

      res <- plot_data()
      shiny::req(res)
      
      ## ---- xlab ------ ###
      names.arg = names(res$total.counts)
      if( length(names.arg) > 20){ names.arg = "" }
      cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
      
      par(mar=c(8,4,2,0.5), mgp=c(2.2,0.8,0))
      barplot(
        res$total.counts/1e6, las=3, border = NA,
        col=rgb(0.2,0.5,0.8,0.8), 
        cex.names = cex.names,
        cex.lab = 1,
        ylab = res$ylab,
        ylim = c(0,max(res$total.counts)/1e6)*1.1,
        names.arg = names.arg
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
      res = c(90,170)*1,                ## resolution of plots
      pdf.width = 6, pdf.height = 6,
      add.watermark = watermark
    )

  })  ## end of moduleServer
}


