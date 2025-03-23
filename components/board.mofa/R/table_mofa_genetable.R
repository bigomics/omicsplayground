##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_genetable_ui <- function(
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

mofa_table_genetable_server <- function(id,
                                        mofa,
                                        selected_factor = reactive(NULL),
                                        annot 
                                        )
{
  moduleServer(id, function(input, output, session) {

    table.RENDER <- function(full=FALSE) {
      mofa <- mofa()
      validate(need(!is.null(mofa), "missing MOFA data."))            
      k=1
      k <- selected_factor()  ## which factor/phenotype
      shiny::req(k)

      wct <- playbase::mofa.plot_centrality(mofa, k, justdata=TRUE)
      wct <- wct[order(-wct$weight),]
      if(full) {
        aa <- annot()[wct$feature,c("feature","symbol","gene_title")]
      } else {
        aa <- annot()[wct$feature,c("feature","symbol")]        
      }
      if(all(aa$feature == aa$symbol)) aa$symbol <- NULL      
      wct$feature <- NULL
      df <- data.frame( factor=k, aa, wct )

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
        ) %>%
        DT::formatStyle(
          "centrality",
          background = color_from_middle(df$centrality, "white", "#fec34d"),
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
