##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_gsetmofa_ui <- function(
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

mofa_table_gsetmofa_server <- function(id,
                                       mofa,
                                       selected_module = reactive(NULL),
                                       selected_trait = reactive(NULL)
                                       )
{
  moduleServer(id, function(input, output, session) {

    table.RENDER <- function() {
      mofa <- mofa()
      validate(need(!is.null(mofa), "missing MOFA data."))            

      dbg("[table_gsetmofa] 1:")
      m=1
      m  <- selected_module()  ## which factor/phenotype
      ph <- selected_trait()  ## which factor/phenotype            
      shiny::req(m,ph)

      gset.mofa <- mofa$gset.mofa      
      w <- gset.mofa$W[,m]
      gs <- rownames(gset.mofa$W)
      Y <- gset.mofa$Y[,ph]
      rho <- cor(t(gset.mofa$X), Y)[,1]
      
      df <- data.frame( module=m, geneset = gs, weight=w, rho=rho, check.names=FALSE)
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
          deferRender = TRUE,
          columnDefs = list(
            list(
              targets = "geneset", ## with no rownames column 1 is column 2
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 80 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 80) + '...</span>' : data;",
                "}"
              )
            )
          )
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
