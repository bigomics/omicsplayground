##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

intersection_plot_venn_diagram_ui <- function(id, label='', height=c(600,800)) {
  
  ns <- shiny::NS(id)
  
  info_text = "The Venn diagram visualizes the number of intersecting genes between the profiles. The list of intersecting genes with further details is also reported in an interactive table below, where users can select and remove a particular contrasts from the intersection analysis."
  
  venndiagram.opts = shiny::tagList(
    shiny::radioButtons(ns('include'),'Counting:', choices=c("both","up/down"), inline=TRUE)
  )
  
  PlotModuleUI(
    ns("vennplot"),
    title = "Venn Diagram",
    label = "b",
    info.text = info_text,
    options = venndiagram.opts,
    download.fmt = c("png","pdf","csv"),
    height = height
  )
  
}


intersection_plot_venn_diagram_server <- function(id, 
                                                  getSignificantTable,
                                                  watermark=FALSE){
  moduleServer( id, function(input, output, session) {
    
    dbg("[intersection_plot_venn_diagram_server] created!")
    
    venndiagram.RENDER <- shiny::reactive({
      dt <- getSignificantTable()
      shiny::req(dt)
      
      if(is.null(dt) || nrow(dt)==0) return(NULL)
      
      dt1 = dt[,2:ncol(dt),drop=FALSE]
      label = LETTERS[1:ncol(dt1)]
      colnames(dt1) = label
      include = "both"
      if(input$include=="up/down") {
        include = c("up/down")
      }
      
      par(mfrow=c(1,1), mar=c(1,1,3,1)*0, bty="n")
      par(oma=c(0.0,0,0,0))
      if(ncol(dt1)==1) {
        frame()
        text(0.5, 0.5, "Error: Venn diagram needs at least two groups")
        return(NULL)            
      } 
      if(ncol(dt1)>7) {
        frame()
        text(0.5, 0.5, "Error: too many groups for Venn diagram")
        return(NULL)
      } 
      
      p <- NULL
      if(0) {
        ## dt1 = dt1[,1:min(5,ncol(dt1))]
        limma::vennDiagram(
          dt1,  main=NULL, cex.main=0.2, cex=1.2, mar=c(0,0,2,0),
          include=include, bty="n", fg=grey(0.7),
          circle.col=c("turquoise", "salmon","lightgreen","orange") )
        tt = paste(label,"=",colnames(dt)[-1])
        legend("topleft", legend=tt, bty='n', cex=0.9, y.intersp=0.95,
               inset=c(0.04,-0.01), xpd=TRUE)       
        
      } else {
        
        ##colnames(dt1) <- colnames(dt)[-1]
        
        x <- apply(dt1, 2, function(x) rownames(dt1)[which(x!=0)])
        
        xlen <- sapply(x,length) 
        dbg("[venndiagram.RENDER] list.len = ",xlen)
        
        if(length(x)==0 || all(xlen==0)) {
          frame()
          text(0.5, 0.5, "Error: no valid genes. Please adjust thresholds.",
               col='red')
          return(NULL)
        }
        
        if(include=="up/down"){
          x_up <- apply(dt1, 2, function(x) paste(rownames(dt1)[which(x>0)], x[which(x>0)]))
          x_down <- apply(dt1, 2, function(x) paste(rownames(dt1)[which(x<0)], x[which(x<0)]))
          
          count_up <- ggVennDiagram::process_region_data(ggVennDiagram::Venn(x_up))$count
          count_down <- ggVennDiagram::process_region_data(ggVennDiagram::Venn(x_down))$count
          count_both <- paste0(count_up, "\n", count_down)
          
          p <- ggVennDiagram::ggVennDiagram(
            x,
            label = "both",
            edge_size = 0.4
          ) +
            ggplot2::scale_fill_gradient(low="grey90",high = "red") +
            ggplot2::theme( legend.position = "none",
                            plot.margin = ggplot2::unit(c(1,1,1,1)*0.3, "cm"))
          
          p$layers[[4]]$data$both <- count_both
        } else{
          p <- ggVennDiagram::ggVennDiagram(
            x,
            label = "count",
            edge_size = 0.4
          ) +
            ggplot2::scale_fill_gradient(low="grey90",high = "red") +
            ggplot2::theme( legend.position = "none", 
                            plot.margin = ggplot2::unit(c(1,1,1,1)*0.3, "cm"))
        }
        
        ## legend
        tt = paste(label,"=",colnames(dt)[-1])
        n1  <- ceiling(length(tt)/2)
        tt1 <- tt[1:n1]
        tt2 <- tt[(n1+1):length(tt)]
        if(length(tt2) < length(tt1)) tt2 <- c(tt2,"   ")
        tt1 <- paste(tt1, collapse='\n')
        tt2 <- paste(tt2, collapse='\n')            
        
        xlim <- ggplot2::ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
        ylim <- ggplot2::ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
        x1 = xlim[1] - 0.1*diff(xlim)
        x2 = xlim[1] + 0.6*diff(xlim)
        y1 = ylim[2] + 0.12*diff(xlim)
        
        p <- p +
          ggplot2::annotate("text", x = x1, y = y1, hjust="left",
                            label = tt1, size=4, lineheight=0.83) +
          ggplot2::annotate("text", x = x2, y = y1, hjust="left",
                            label = tt2, size=4, lineheight=0.83) +
          ggplot2::coord_sf(clip="off")
      }        
      p
    })
    
    PlotModuleServer(
      "vennplot",
      func = venndiagram.RENDER,
      res = c(60,120),                ## resolution of plots
      pdf.width = 8, pdf.height = 5,
      add.watermark = watermark
    )
    
  }
  )
}