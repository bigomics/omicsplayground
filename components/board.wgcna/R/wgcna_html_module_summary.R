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
  width
) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("text"),
    outputFunc = htmlOutput,
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf", "svg")
  )
}

wgcna_html_module_summary_server <- function(id,
                                             wgcna,
                                             multi = FALSE,
                                             r_annot = reactive(NULL),
                                             r_module,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    contents_text <- shiny::reactive({
      module <- isolate(r_module())
      paste0(
        "<b>", module, " module</b><br><br>",
        "<div class='alert alert-info'>AI module summary refactor ongoing.</div>"
      )
    })

    text.RENDER <- function() {
      res <- contents_text()
      shiny::div(class = "gene-info", shiny::HTML(res))
    }

    text.RENDER2 <- function() {
      res <- contents_text()
      shiny::div(shiny::HTML(res), style = "font-size:22px;")
    }

    PlotModuleServer(
      "text",
      plotlib = "generic",
      plotlib2 = "generic",
      func = text.RENDER,
      func2 = text.RENDER2,
      renderFunc = shiny::renderUI,
      renderFunc2 = shiny::renderUI,
      pdf.width = 8, pdf.height = 5,
      res = c(75, 100),
      add.watermark = watermark
    )
  })
}
