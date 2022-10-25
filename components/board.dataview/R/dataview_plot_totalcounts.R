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
        plotlib = "plotly",
        #outputFunc = plotOutput,
        #outputFunc2 = plotOutput,        
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
          
          ylab = "Total counts"
          if(data_groupby != "<ungrouped>") {
              ylab = "Mean total counts"
          }
          
          res <- list(
              df = data.frame(
                  sample = names(tbl$total.counts),
                  counts = tbl$total.counts
              ),
              ylab = ylab
          )
          return(res)
      })
      
      
      plot.RENDER <- function() {
          
          res <- plot_data()
          shiny::req(res)
          df <- res[[1]]
          
          ## ---- xlab ------ ###
          names.arg = df$sample
          if( length(names.arg) > 20){ names.arg = "" }
          cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
          
          par(mar=c(8,4,2,0.5), mgp=c(2.2,0.8,0))
          barplot(
              df$counts/1e6, las=3, border = NA,
              col=rgb(0.2,0.5,0.8,0.8), 
              cex.names = cex.names,
              cex.lab = 1,
              ylab = paste(res$ylab,"(M)"),
              ylim = c(0,max(df$counts)/1e6)*1.1,
              names.arg = names.arg
          )
      }
      
      modal_plot.RENDER <- function() {
          plot.RENDER()
      }

      plotly.RENDER <- function() {
          
          res <- plot_data()
          shiny::req(res)
          df <- res[[1]]

          fig <- 
            plotly::plot_ly(
              data = df,
              x = ~sample,
              y = ~counts,
              type = 'bar',
              marker = list(
                color = omics_colors("mid_blue")
              ), 
              hovertemplate = ~paste0(
                "Sample: <b>", sample,
                "</b><br>", res$ylab, ": <b>", sprintf("%8.0f", counts),
                "</b><extra></extra>"
              )
            ) %>% 
            plotly::layout(
              xaxis = list(title = FALSE),
              yaxis = list(title = res$ylab),
              font = list(family = "Lato"),
              margin = list(l = 10, r = 10, b = 10, t = 10)   
            ) %>% 
            plotly_default1()
            
          fig
      }
      
      modal_plotly.RENDER <- function() {
          plotly.RENDER()
      }
      
    PlotModuleServer(
        "pltmod",
        plotlib = "plotly",
        func = plotly.RENDER,
        func2 = modal_plotly.RENDER,
        csvFunc = plot_data,   ##  *** downloadable data as CSV
        res = c(90,170)*1,                ## resolution of plots
        pdf.width = 6, pdf.height = 6,
        add.watermark = watermark
    )
      
  })  ## end of moduleServer
}


