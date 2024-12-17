##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_enrichment_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::selectInput(ns("labeltype"), "Label type",
      choices = c("feature","symbol","gene_title"))
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_enrichment_server <- function(id,
                                        pgx,
                                        gsea,
                                        input_k = reactive(1),
                                        req.selection = FALSE,
                                        select = reactive(NULL),
                                        ntop = 15, 
                                        watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function(ntop=20, par=TRUE) {
      gsea <- gsea()      
      validate(need(!is.null(gsea), "missing MOFA data."))      

      k <- input_k()
      shiny::req(k)
      shiny::req(k %in% names(gsea))
      sel <- select()
      
      nx <- grep("[:]",rownames(pgx$X),value=TRUE)      
      dtypes <- unique(sub(":.*","",nx))
      dtypes
      dbg("[mofa_plot_enrichment_server] dtypes = ",dtypes)
      if(length(sel)==1) {
        if(length(dtypes) > 1) {
          par(mfrow = c(1,length(dtypes)))
          par(mar=c(4,4,0.5,1))
          if(length(dtypes)>3) par(mfrow=c(2,ceiling(length(dtypes)/2)))
          for(dt in dtypes) {
            tt <- paste0(k, " (",toupper(dt),")")
            playbase::pgx.Volcano(pgx, contrast=k,
              hilight=sel, label=sel, labeltype=input$labeltype,
              ntop=10, plotlib="base",
              datatype=dt, cex=0.8, fc=0.5, title=tt)
          }
        } else {
          playbase::pgx.Volcano(pgx, contrast=k, hilight=sel,
            labeltype=input$labeltype, datatype=NULL, cex=0.8)
        }
      } else {
        if(par) {
          par(mfrow=c(1,2), mar=c(4,4,0.5,1))
          plot.new()
        }
        playbase::mofa.plot_enrichment(
          gsea[[k]],
          ntop = ntop,
          select = select(),
          strip.names = TRUE,
          par = FALSE,
          title = "")
      }
    }

    plot.RENDER2 <- function() {
      par(mfrow=c(1,2), mar=c(4,14,2,2))
      plot.new()
      plot.RENDER(ntop=25, par=FALSE) 
    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      func2 = plot.RENDER2,
      pdf.width = 9, pdf.height = 5,
      res = c(80, 130),
      add.watermark = watermark
    )

    
  })
}
