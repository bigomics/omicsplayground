##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
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

  PlotModuleUI(
    ns("text"),
    outputFunc = htmlOutput,
    title = title,
    label = label,
    info.text = info.text,
    #options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_html_module_summary_server <- function(id,
                                             wgcna,
                                             multi = FALSE,
                                             r_module,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {


    info_text <- function() {
      wgcna <- wgcna()
      module <- r_module()
      shiny::req(wgcna)      
      res <- "Summary not available"
      if("summary" %in% names(wgcna) && module %in% names(wgcna$summary)) {
        res <- wgcna$summary[[module]]
        res <- gsub("\n","<p>",res)
        res <- gsub(" [*]{2}","<b>",res)
        res <- gsub("[*]{2} ","</b>",res)
        res <- paste0("<b>",module,":</b> ", res)
      }
      return(res)
    }

    info.RENDER <- function() {
      res <- info_text()
      shiny::div(class="gene-info", shiny::HTML(res))
    }

    info.RENDER2 <- function() {
      res <- info_text()
      shiny::div( shiny::HTML(res), style="font-size:22px;" )
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
