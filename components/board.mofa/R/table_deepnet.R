##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

table_deepnet_gradients_ui <- function(
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

table_deepnet_gradients_server <- function(id,
                                           net,
                                           pgx,
                                           phenoFC,
                                           conditions,
                                           datatypes
                                           )
{
  moduleServer(id, function(input, output, session) {

    table.RENDER <- function(full=TRUE) {
      net <- net()
      grad <- net$get_gradients()[[1]]
      all.fc <- phenoFC()

      dbg("[table_deepnet_gradients_server] reacted!")
      
      cond <- conditions()
      dtypes <- datatypes()
      shiny::validate( shiny::need( length(cond)>0, "please select a condition"))
      shiny::validate( shiny::need( length(dtypes)>0, "please select a datatype"))      
      
      ## check if we are too early after a change
      net <- net()
      shiny::req( dtypes %in% names(net$X))
      shiny::req( cond %in% net$labels[[1]])

      grad <- grad[dtypes]
      grad <- lapply(grad, function(g) g[,cond,drop=FALSE] )
      
      df <- playbase::deep.plotGradientVSFoldchange(
        grad, fc = all.fc, data=TRUE, par=FALSE)
      
      numeric.cols <- colnames(df)
      
      ## add gene title for readable names (e.g. for metabolomics)
      feature <- stringr::str_trunc(rownames(df),20)
      symbol  <- stringr::str_trunc(pgx$genes[feature,"symbol"],20)
      title   <- stringr::str_trunc(pgx$genes[feature,"gene_title"],60)
      if(full) {
        df <- cbind( feature=feature, symbol=symbol, title=title, df)
      } else {
        df <- cbind( feature=feature, symbol=symbol, df)
      }

      DT::datatable(
        df,
        rownames = FALSE, 
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        plugins = "scrollResize",
        options = list(
          dom = "lfrtip", #
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

    table.RENDER2 <- function() {
      table.RENDER(full=TRUE)
    }
    
    table <- TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER2,      
      selector = "single"
    )

    return(table)
  })
}
