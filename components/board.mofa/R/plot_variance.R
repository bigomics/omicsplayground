##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_variance_ui <- function(
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

mofa_plot_variance_server <- function(id,
                                      mofa,
                                      type = c("factorxview","factor","view")[1],
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      res <- mofa()
      shiny::req(res$model)
      if( type == "factorxview") {
        plt <- NULL
        if( inherits(res$model, "MOFA")) {
          plt <- MOFA2::plot_variance_explained(res$model, x="view", y="factor")
        } else {
          par(mar=c(3,1,2,5))                  
          purples <- RColorBrewer::brewer.pal(9,"Purples")
          ##col <- colorRampPalette(c("grey98", purples))(64)
          playbase::gx.imagemap(t(res$V), clust=FALSE, col=purples)
        }
        return(plt)
      }
      if( type == "view") {
        y <- rowSums( res$V )
        par(mar=c(3,4,2,0))        
        barplot(y, names.arg = names(y), ylab="Var. (%)")
      }
      if( type == "factor") {
        y <- colSums( res$V )
        par(mar=c(3,4,2,0))
        barplot(y, names.arg = names(y), ylab="Var. (%)")
      }
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 5, pdf.height = 5,
      res = c(80, 110),
      add.watermark = watermark
    )

    
  })
}
