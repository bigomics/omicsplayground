##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_table_modulegenes_ui <- function(
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
    options = options,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = label
  )
}

multiwgcna_table_modulegenes_server <- function(id,
                                                mwgcna,
                                                r_annot,
                                                r_phenotype = reactive(NULL),
                                                r_module = reactive(NULL)
                                                )
{
  moduleServer(id, function(input, output, session) {
    
    table_df <- function() {

      wgcna <- mwgcna()

      pheno <- r_phenotype()
      module <- r_module()
      annot <- r_annot()

      shiny::req(wgcna)
      shiny::req(pheno)
      shiny::req(module)
      shiny::req(annot)      
      
      if(input$showallmodules) module <- NULL
      
      df <- list()
      for(dt in names(wgcna)) {
        stats <- playbase::wgcna.getGeneStats(
          wgcna = wgcna[[dt]],
          trait = pheno,
          plot = FALSE,
          module = module,
          col = NULL,
          main = NULL)
        
        if(nrow(stats)>0) {
          ## make feature first column
          sel <- unique(c("feature",colnames(stats)))
          sel <- intersect(sel, colnames(stats))
          df[[dt]] <- stats[,sel]
        }
      }
      
      cols <- Reduce(intersect, lapply(df,colnames))
      df <- lapply(df, function(a) a[,cols] ) 
      for(i in 1:length(df)) rownames(df[[i]]) <- paste0(names(df)[i],":",rownames(df[[i]]))
      names(df) <- NULL
      df <- do.call(rbind, df)

      if(!is.null(annot)) {
        ## if we have titles, add them
        title2 <- playbase::probe2symbol(df$feature, annot, query="gene_title")
        if(!all(title2 %in% c(NA,"","NA"))) df$title <- title2

        ## if we have symbols (and they differ from features), add them
        symbol <- playbase::probe2symbol(df$feature, annot, query="symbol")        
        df$symbol <- NULL
        nqq <- mean(symbol != df$feature & !is.na(symbol), na.rm=TRUE)
        if( nqq > 0.5) df$symbol <- symbol        
      }
      
      ## df <- df[df$module %in% module,]  ## select??
      df <- df[order(-df$score),]
      
      return(df)
    }
    
    render_table <- function(full=TRUE) {
      
      df <- table_df()      
      
      ## set correct types for filter
      df$module <- factor(df$module)
      
      if(!full) {
        sel <- c("module","feature","symbol","title","score","traitSignificance","moduleMembership")
      } else {
        sel <- c("module","feature","symbol","title","score")
        sel <- unique(c(sel, colnames(df)))
      }
      
      sel <- intersect(sel, colnames(df))
      if(!full) {
        if(all(c("symbol","feature") %in% sel)) sel <- setdiff(sel, "feature")
      }
      df <- df[,sel]
      
      ## rename
      colnames(df) <- sub("moduleMembership","MM",colnames(df))
      colnames(df) <- sub("traitSignificance","TS",colnames(df))
      
      numeric.cols <- which(sapply(df, class) == "numeric")

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
              #targets = c(1), ## without rownames column 2 is target 1
              targets = c("feature"), ## without rownames column 2 is target 1
              render = DT::JS("$.fn.dataTable.render.ellipsis( 80, false )")
            )
          )                    
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          "score",
          background = color_from_middle(df$score, "lightblue", "#f5aeae"),
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
      csvFunc = table_df,
      selector = "single"
    )

    return(table)
  })
}
