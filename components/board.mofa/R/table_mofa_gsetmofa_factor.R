##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_gsetmofa_factor_ui <- function(
    id,
    label = "",
    title = "",
    info.text = "",
    caption = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("tablemodule"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

mofa_table_gsetmofa_factor_server <- function(id,
                                              mofa,
                                              selected_factor = reactive(NULL),
                                              selected_trait = reactive(NULL)
                                              )
{
  moduleServer(id, function(input, output, session) {

    table.RENDER <- function() {
      mofa <- mofa()
      validate(need(!is.null(mofa), "missing MOFA data."))            
      
      m=1
      k <- selected_factor()  ## which factor
      ph <- selected_trait()  ## which phenotype            
      shiny::req(k,ph)

      w <- mofa$W[,k]
      g <- rownames(mofa$W)
      Y <- mofa$Y[,ph]
      rho <- cor(t(mofa$X), Y)[,1]
      
      df <- data.frame( factor=k, gene = g, weight=w, rho=rho, check.names=FALSE)
      df <- df[order(-df$w),]
      numeric.cols <- grep("score|weight|rho",colnames(df))
      
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
      "tablemodule",
      func = table.RENDER,
      selector = "single"
    )

    return(table)
  })
}
