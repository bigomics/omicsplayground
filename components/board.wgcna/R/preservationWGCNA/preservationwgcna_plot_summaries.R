##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_plot_summaries_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("showiqr"),
      label = "Show IQR",
      value = FALSE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

preservationWGCNA_plot_summaries_server <- function(id,
                                                    rwgcna
                                                    ) {
  moduleServer(id, function(input, output, session) {
    
    plot.RENDER <- function() {

      res <- rwgcna()
      shiny::req(res)

      nsets <- ncol(res$Zsummary)
      nr <- ceiling(sqrt(nsets*3))
      par(mfrow=c(nr,nr), mar=c(5,5,4,1)) 
      playbase::wgcna.plotPreservationSummaries(
        res, setpar=FALSE)

    }
    
    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8,
      pdf.height = 10,
      res = c(80, 100),
      add.watermark = FALSE
    )
    
  })
}



