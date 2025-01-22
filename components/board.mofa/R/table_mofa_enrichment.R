##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_mofa_enrichment_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

mofa_table_mofa_enrichment_server <- function(id,
                                              gsea,
                                              input_k = reactive(1),
                                              watermark = TRUE
                                              )
{
  moduleServer(id, function(input, output, session) {

    table.RENDER <- function(full=FALSE) {
      gsea <- gsea()

      validate(need(!is.null(gsea), "missing GSEA data."))            
      k=1
      k <- input_k()  ## which factor/phenotype
      shiny::req(k)

      df <- gsea[[k]]
      numeric.cols <- grep("score|pval|p$|q$|rho",colnames(df))
      
      DT::datatable(
        df,
        rownames = FALSE, #
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        plugins = "scrollResize",
        options = list(
          dom = "lfrtip", #
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, #
          scrollY = "70vh",
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }
    
    table <- TableModuleServer(
      "table",
      func = table.RENDER,
      selector = "single"
    )

    return(table)
  })
}
