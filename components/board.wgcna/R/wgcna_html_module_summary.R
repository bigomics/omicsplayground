##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_html_module_summary_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height,
    width) {
  ns <- shiny::NS(id)

  opts <- shiny::tagList(
    shiny::checkboxInput(ns("show_cov"), "covariance", FALSE)
  )

  PlotModuleUI(
    ns("text"),
    outputFunc = htmlOutput,
    title = title,
    label = label,
    info.text = info.text,
    #options = opts,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_html_module_summary_server <- function(id,
                                             wgcna,
                                             selected_module,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {


    info_text <- function() {
      wgcna <- wgcna()
      module <- selected_module()
      ## trait <- selected_trait()      
      
      top <- playbase::wgcna.getTopGenesAndSets(wgcna, module=module)      
      topgenes <- paste( top$genes[[module]], collapse=";")
      res1 <- paste("<b>Key genes:</b>", topgenes)

      res2 <- "Sorry. Description not available."
      if("summary" %in% names(wgcna) && module %in% names(wgcna$summary)) {
        res2 <- paste("<b>Summary:</b>", wgcna$summary[[module]])
      }
      res2 <- gsub("\n","<p>",res2)
      res2 <- gsub(" [*]{2}","<b>",res2)
      res2 <- gsub("[*]{2} ","</b>",res2)

      res3 <- "<i>[AI generated]</i>"
      
      res <- paste(res1, res2, res3, sep="<br><br><p>")
      return(res)
    }

    info.RENDER <- function() {
      res <- info_text()
      div(class = "gene-info", shiny::HTML(res))
    }

    info.RENDER2 <- function() {
      res <- info_text()
      div(shiny::HTML(res))
    }
    
    PlotModuleServer(
      "text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = info.RENDER,
      func2 = info.RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
