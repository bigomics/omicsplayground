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
                                              pgx,
                                              selected_factor = reactive(NULL),
                                              selected_trait = reactive(NULL),
                                              selected_pathway = reactive(NULL)
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
      
      sel <- selected_pathway()      
      ##if(length(sel)==0) return(NULL)
      validate(need(length(sel)>0, "Please select a geneset"))
      
      sel <- sub("^gset.[pmg]x:","",sel)
      gs.genes <- names(which(mofa$GMT[,sel[1]] != 0))
      pp <- playbase::map_probes(pgx$genes, gs.genes)
      pp <- intersect(pp, rownames(mofa$W))
      
      w <- mofa$W[pp,k]
      s <- playbase::probe2symbol(pp, pgx$genes, "symbol")
      f <- playbase::probe2symbol(pp, pgx$genes, "feature")
      tt <- playbase::probe2symbol(pp, pgx$genes, "gene_title")            
      tt <- stringr::str_trunc(tt,40)
      Y <- mofa$Y[,ph]
      rho <- cor(t(mofa$X[pp,,drop=FALSE]), Y)[,1]
      
      df <- data.frame( factor=k, symbol=s, title=tt, weight=w,
                       rho=rho, check.names=FALSE)
      df <- df[order(-df$w),,drop=FALSE]
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
          scrollX = TRUE, #
          scrollY = "70vh",
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE
          ## columnDefs = list(
          ##   list(
          ##     targets = "title", 
          ##     render = DT::JS(
          ##       "function(data, type, row, meta) {",
          ##       "return type === 'display' && data.length > 35 ?",
          ##       "'<span title=\"' + data + '\">' + data.substr(0, 35) + '...</span>' : data;",
          ##       "}"
          ##     )
          ##   )
          ## )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          "weight",
          background = color_from_middle(df$weight, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) 
  
    }
    
    table <- TableModuleServer(
      "tablemodule",
      func = table.RENDER,
      selector = "single"
    )

    return(table)
  })
}
