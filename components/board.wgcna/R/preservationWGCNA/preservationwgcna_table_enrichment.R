##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

preservationWGCNA_table_enrichment_ui <- function(
    id,
    label = "a",
    title = "Title",
    info.text = "Info",
    caption = "Caption",
    height = 400,
    width = 400 ) {
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

preservationWGCNA_table_enrichment_server <- function(id,
                                               rwgcna,
                                               rmodule = reactive(NULL)
                                               )
{
  moduleServer(id, function(input, output, session) {

    table_df <- function() {
      res <- rwgcna()
      module <- rmodule()
      shiny::req(module)      
      df <- res$gsea[[module]]
      return(df)
    }
    
    render_table <- function(full=TRUE) {
      df <- table_df()
      shiny::req(df)
      ##df$module <- factor(df$module)
      df$module <- NULL ## don't show
      
      if(!full) {
        sel <- c("geneset","module","score","q.value","overlap")
        sel <- intersect(sel, colnames(df))
        df <- df[,sel]
      }
      
      ## set correct types for filter
      numeric.cols <- which(sapply(df, class) == "numeric")
      
      dt <- DT::datatable(
        df,
        rownames = FALSE, #
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        plugins = c("scrollResize","ellipsis"),
        #filter = 'top',
        options = list(
          dom = "lfrtip", #
          scrollX = TRUE, #
          scrollY = "70vh",
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE,
          autoWidth = TRUE,          
          columnDefs = list(list(
            targets = c(0), ## without rownames column 1 is target 0
            render = DT::JS("$.fn.dataTable.render.ellipsis( 60, false )")
          ))          
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
      
      if(1) {
        dt <- dt %>%
          DT::formatStyle(
            "score",
            background = color_from_middle(df$score, "lightblue", "#f5aeae"),
            backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
          )
      }

      return(dt)
    }

    table.RENDER <- function() {
      render_table(full=FALSE)
    }

    table.RENDER2 <- function() {
      render_table(full=TRUE)
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
