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

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_enrichment_server <- function(id,
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
      ## shiny::req(gsea$gsea)
      ##if(req.selection) shiny::req(select())
      
      shiny::req(k %in% names(gsea))
      
      if(FALSE && length(select()==1)) {
        ## DOES NOT WORK...
        gs <- select()
        gset <- names(which(PATHBANK[gs,]!=0))
        datatypes <- colnames(gsea[[k]]$score)
        dt=datatypes[1]
        for(dt in datatypes) {
          rnk <- gsea$rnk
          gset1 <- intersect(gset, names(rnk))
          gsea.enplot( rnk, gset1, main=dt)
        }
        
      } else {
        if(par) {
          par(mfrow=c(1,2), mar=c(4,4,2,2))
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
      res = c(72, 110),
      add.watermark = watermark
    )

    
  })
}
