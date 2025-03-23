##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_enrichmentgenes_ui <- function(
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

mofa_table_enrichmentgenes_server <- function(id,
                                              pgx,
                                              selected_factor = reactive(NULL),
                                              selected_pathway = reactive(NULL)
                                              )
{
  moduleServer(id, function(input, output, session) {

    table_df <- function(full) {
      validate(need("mofa" %in% names(pgx), "No MOFA slot"))      

      k <- selected_factor()  
      shiny::req(k)
      gset <- selected_pathway()  ## which factor/phenotype      
      shiny::validate(shiny::need(length(gset)==1, "Please select a gene set"))      
      shiny::req(gset %in% colnames(pgx$GMT))
      genes <- names(which(pgx$mofa$GMT[,gset]!=0))
      pp <- playbase::map_probes(pgx$genes, genes, "symbol")
      W <- pgx$mofa$W
      pp <- intersect(pp, rownames(W))
      ww <- sort( W[pp,k], decreasing = TRUE )
      
      if(full) {
        aa <- pgx$genes[pp,c("feature","symbol","gene_title")]
      } else {
        aa <- pgx$genes[pp,c("feature","symbol")]
      }
      if(all(aa$feature == aa$symbol)) aa$symbol <- NULL      
      df <- data.frame(aa, weight=ww, check.names=FALSE )
      return(df)
    }
    
    table.RENDER <- function(full=FALSE) {

      df <- table_df(full)
      numeric.cols <- grep("score|weight|centrality",colnames(df))
      
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
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          "weight",
          background = color_from_middle(df$weight, "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ) 

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
