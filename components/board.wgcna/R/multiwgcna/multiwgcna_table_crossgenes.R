##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_table_crossgenes_ui <- function(
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
      inputId = ns("showallmodules"),
      label = "Show all modules",
      value = FALSE
    )
  )
  
  TableModuleUI(
    ns("table"),
    info.text = info.text,
    #options = options,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

multiwgcna_table_crossgenes_server <- function(id,
                                               mwgcna,
                                               r_annot,
                                               r_module
                                               )
{
  moduleServer(id, function(input, output, session) {

    table_df <- function() {

      wgcna <- mwgcna()
      module <- r_module()
      annot <- r_annot()

      shiny::req(wgcna)
      shiny::req(module)
      shiny::req(annot)      
      
      df <- playbase::wgcna.getModuleCrossGenes(
        wgcna,
        ref = NULL,
        multi = TRUE,
        ngenes = 100,
        modules = module
      )[[1]]
            
      if(!is.null(annot)) {
        ## if we have titles, add them
        title2 <- playbase::probe2symbol(df$gene, annot, query="gene_title")
        if(!all(title2 %in% c(NA,"","NA"))) df$title <- title2
        ## if we have symbols (and they differ from features), add them
        symbol <- playbase::probe2symbol(df$gene, annot, query="symbol")        
        df$symbol <- NULL
        nqq <- mean(symbol != df$gene & !is.na(symbol), na.rm=TRUE)
        if(nqq > 0.5) df$symbol <- symbol        
      }
      
      ## df <- df[df$module %in% module,]  ## select??
      df <- df[order(-df$rho),]
      
      return(df)
    }
    
    render_table <- function(full=TRUE) {

      df <- table_df()      

      ## set correct types for filter
      df$module <- factor(df$module)

      if(!full) {
        sel <- c("module","gene","symbol","title","rho")
        sel <- intersect(sel, colnames(df))
        if(all(c("symbol","gene") %in% sel)) sel <- setdiff(sel, "gene")
        df <- df[,sel]
      }
      
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
          columnDefs = list(
            list(
              targets = c("title"), ## without rownames column 3 is target 2
              render = DT::JS("$.fn.dataTable.render.ellipsis( 60, false )")
            )
          )                    
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          "rho",
          background = color_from_middle(df$rho, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )

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
      csvFunc = table_df
      ##selector = "single"
    )

    return(table)
  })
}
