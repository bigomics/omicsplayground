##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

wgcna_table_genes_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    shiny::checkboxInput(ns("showpvalues"),"Show p-values", FALSE)
  )
  
  TableModuleUI(
    ns("datasets"),
    options = options,
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

wgcna_table_genes_server <- function(id,
                                     wgcna,
                                     selected_module,
                                     selected_trait
                                     ) {
  moduleServer(id, function(input, output, session) {

    RENDER <- function() {
      res <- wgcna()
      module <- selected_module()
      trait <- selected_trait()      
      
      df <- playbase::wgcna.getGeneStats(res, module=module, trait=trait, plot=FALSE) 
      if(input$showpvalues==FALSE) {
        df <- df[, grep("Pvalue", colnames(df), invert=TRUE), drop=FALSE]
      }
      df <- df[ which(df$module == module), , drop=FALSE ]
      
      numeric.cols <- grep("^module$", colnames(df), invert=TRUE)

      DT::datatable(
        df,
        rownames = TRUE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        plugins = "scrollResize",
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", #
          scrollX = TRUE, #
          scrollY = "70vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE
        ) 
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "10px", lineHeight = "70%")
    }

    RENDER_modal <- function() {
      dt <- RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    }

    geneTable_module <- TableModuleServer(
      "datasets",
      func = RENDER,
      func2 = RENDER_modal,
      selector = "none"
    )

    return(geneTable_module)
  })
}
