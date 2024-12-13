##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_boxplots_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width  = 400) {

  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::selectizeInput(
      ns("selected_pheno"), "Select phenotype",
      choices=NULL, multiple=TRUE), 
    shiny::checkboxInput(
      inputId = ns("collapse"),
      label = "Collapse covariates"
    )
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

mofa_plot_boxplots_server <- function(id,
                                      mofa,
                                      input_factor = reactive(1),
                                      watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    observeEvent( mofa(),{
      data <- mofa()
      phenos <- colnames(data$samples)
      updateSelectInput(session, "selected_pheno",
        choices = phenos,
        selected = head(phenos,12))
    })

    plot.RENDER <- function() {
      res <- mofa()
      shiny::req(res)
      k <- input_factor()
      if(!is.null(k)) shiny::req(k %in% colnames(res$F))

      pheno <- input$selected_pheno
      if(!is.null(pheno)) shiny::req(all(pheno %in% colnames(res$samples)))
      nph <- length(pheno)
      shiny::validate( shiny::need(nph>0,"Must select at least one phenotype"))

      samples <- data.frame( res$samples )
      nr <- ceiling(sqrt(nph))
      nc <- ceiling(nph / nr)
      par(mfrow=c(nr,nc), mar=c(4,4,2.8,0.5))      
      for(ph in pheno) {
        y <- samples[,ph]
        if(all(is.na(y))) next
        f1 <- res$F[,k]
        if(class(y) %in% c("numeric","integer") &&
             all(y %in% c(NA,0:9))) y <- factor(y)
        isfactor <- (class(y) %in% c("character","factor","logical"))
        ylab <- paste(k,"score")
        if(isfactor) {
          y <- factor(y)
          boxplot( f1 ~ y, main="", ylab=ylab, xlab=ph)
        } else {
          plot( y, f1, main="", ylab=ylab, xlab=ph)
        }
        title(ph, cex.main=1.2)
      }
        
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 10, pdf.height = 6,
      res = c(85, 120),
      add.watermark = watermark
    )
  })
}
