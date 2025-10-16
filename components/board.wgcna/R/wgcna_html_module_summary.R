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

  opts <- shiny::tagList(
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
                                             pgx,
                                             multi = FALSE,
                                             r_module,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {


    info_text <- function() {
      wgcna <- wgcna()
      module <- r_module()
      ##pgx <- pgx()
      shiny::req(wgcna)
      
      shiny::validate(shiny::need("gsea" %in% names(wgcna),
        "Error: object has not enrichment results (gsea)")
      )

      shiny::req(pgx)
      annot <- pgx$genes
      
      if(multi) {
        top <- playbase::wgcna.getConsensusTopGenesAndSets(
          wgcna, annot=annot, module=module, ntop=25)
      } else {
        top <- playbase::wgcna.getTopGenesAndSets(
          wgcna, annot=annot, module=module, ntop=25)
      }
      topgenes <- paste(top$genes[[module]], collapse=";")
      res1 <- paste("<b>Key genes:</b>", topgenes)      
      res2 <- "<b>Summary:</b> not available"
      if("summary" %in% names(wgcna) && module %in% names(wgcna$summary)) {
        res2 <- wgcna$summary[[module]]
        res2 <- gsub("\n","<p>",res2)
        res2 <- gsub(" [*]{2}","<b>",res2)
        res2 <- gsub("[*]{2} ","</b>",res2)
        res2 <- paste("<b>Summary:</b>", res2)
      }
      
      res <- paste(res1, res2, sep="<br><br><p>")
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
