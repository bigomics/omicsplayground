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
      
      df <- list()
      for(dt in names(wgcna)) {
        stats <- playbase::wgcna.getGeneStats(
          wgcna = wgcna[[dt]],
          trait = pheno,
          plot = FALSE,
          showallmodules = input$showallmodules,
          module = module,
          col = NULL,
          main = NULL)
        if(nrow(stats)>0) {
          stats$feature <- NULL
          df[[dt]] <- cbind(feature=rownames(stats), stats)
        }
      }
      
      cols <- Reduce(intersect, lapply(df,colnames))
      df <- lapply(df, function(a) a[,cols] ) 
      for(i in 1:length(df)) rownames(df[[i]]) <- paste0(names(df)[i],":",rownames(df[[i]]))
      names(df) <- NULL
      df <- do.call(rbind, df)
      
      if(!is.null(annot)) {
        tt <- playbase::probe2symbol(rownames(df), annot, query="gene_title")
        df$title <- tt
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
        sel <- c("module","feature","title","score","traitSignificance","moduleMembership")
      } else {
        sel <- c("module","feature","title","score","traitSignificance","moduleMembership")
        sel <- unique(c(sel, colnames(df)))
      }
      
      sel <- intersect(sel, colnames(df))
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
        plugins = "scrollResize",
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
      selector = "single"
    )

    return(table)
  })
}
