##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_mgsea_ui <- function(
    id,
    title = "",
    info.text = "",
    info.methods,
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::radioButtons(ns("size.par"), "size by", c("p-value"="p", "q-value"="q"))
  )
  
  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    options = options,
    info.text = info.text,
    info.methods = info.methods,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}


mofa_plot_mgsea_server <- function(id,
                                   mgsea,
                                   input_k = reactive(1),
                                   select = reactive(NULL),
                                   watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {
      mgsea <- mgsea()
      validate(need(!is.null(mgsea), "missing GSEA data."))        
      k <- input_k()
      shiny::req(k)

      cn <- colnames(mgsea[[k]])
      types <- sub("^pval.","",grep("^pval",cn,value=TRUE))
      types <- intersect(c("mx","px","gx",types), types)
      if(length(types)==1) types <- rep(types,2)
      hilight <- NULL
      selected <- select()
      if(!is.null(selected) && length(selected)<100) {
        hilight <- selected
      }
      
      par(mar=c(4.5,4.5,1,0.5))
      saveRDS(list(mgsea,types,input$size.par,hilight), "~/Desktop/MNT/plot.RDS")
      playbase::mgsea.plot_scatter(
        mgsea[[k]], type1=types[1], type2=types[2],
        size.par = input$size.par,
        hilight = hilight)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 8,
      res = c(72, 110),
      add.watermark = watermark
    )

    
  })
}
