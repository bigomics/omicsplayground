##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

consensusWGCNA_table_modulegenes_ui <- function(
    id,
    label = "a",
    title = "Title",
    info.text = "Info",
    caption = "Caption",
    height = 400,
    width = 400 ) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("fulltable"),
      label = "Show full table",
      value = TRUE
    ),
    shiny::checkboxInput(
      inputId = ns("showallmodules"),
      label = "Show all modules",
      value = FALSE
    )
  )

  
  TableModuleUI(
    ns("table"),
    info.text = info.text,
    options = options,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

consensusWGCNA_table_modulegenes_server <- function(id,
                                                    mwgcna,
                                                    r_annot,
                                                    r_trait = reactive(NULL),
                                                    r_module = reactive(NULL)
                                                    )
{
  moduleServer(id, function(input, output, session) {

    table_df <- function() {

      cons <- mwgcna()
      
      trait <- r_trait()
      module <- r_module()
      annot <- r_annot()
      
      shiny::req(cons)
      shiny::req(trait)
      shiny::req(module)
      shiny::req(annot)      

      if(input$showallmodules) module <- NULL
      
      stats <- playbase::wgcna.getConsensusGeneStats(
        cons,
        stats = cons$stats,
        trait = trait,
        module = module
      )
      
      which_table <- ifelse(input$fulltable, "full", "consensus")
      df <- stats[[which_table]]
      
      if(!is.null(annot)) {
        df$title <- playbase::probe2symbol(df$feature, annot, query="gene_title")
        symbol <- playbase::probe2symbol(df$feature, annot, query="symbol")
        if(mean(df$feature == symbol) < 0.2) df$symbol <- symbol        
      }
      
      return(df)
    }
    
    render_table <- function(full=TRUE) {

      df <- table_df()      
      
      ## set correct types for filter
      df$module <- factor(df$module)
      
      score.cols <- grepl("^score", colnames(df)) & !grepl("Pvalue", colnames(df))
      if(!full) {
        cols <- c("module","feature","symbol","title")
        cols <- c(cols, colnames(df)[which(score.cols)], "consensus")
        cols <- intersect(cols, colnames(df))
        df <- df[,cols]
      }

      if(!input$showallmodules) df$module <- NULL
      
      ## order name first
      cols <- unique(c("module","feature","symbol","title",colnames(df)))
      cols <- intersect(cols, colnames(df))
      df <- df[,cols]      
      
      ## rename
      colnames(df) <- sub("moduleMembership","MM",colnames(df))
      colnames(df) <- sub("traitSignificance","TS",colnames(df))
      
      numeric.cols <- which(sapply(df, class) == "numeric")
      score.cols <- grepl("^score", colnames(df)) & !grepl("Pvalue", colnames(df))
      score.vals <- df[,score.cols,drop=FALSE]
      
      DT::datatable(
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
          columnDefs = list(
            list(
              targets = c(1), ## without rownames column 2 is target 1
              render = DT::JS("$.fn.dataTable.render.ellipsis( 20, false )")
            ),
            list(
              targets = c(2), ## without rownames column 2 is target 1
              render = DT::JS("$.fn.dataTable.render.ellipsis( 40, false )")
            )
          )                    
        ) ## end of options
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          score.cols,
          background = color_from_middle(score.vals, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) 
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
