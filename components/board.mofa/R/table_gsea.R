##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_table_mgsea_ui <- function(
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

mofa_table_mgsea_server <- function(id,
                                   gsea,
                                   datatypes = NULL,
                                   input_k = reactive(1),
                                   top = 20
                                   )
{
  moduleServer(id, function(input, output, session) {

    table.RENDER <- function(full=FALSE) {
      gsea <- gsea()

      validate(need(!is.null(gsea), "missing GSEA data."))            
      k=1
      k <- input_k()  ## which factor/phenotype
      shiny::req(k)
      
      if(is.null(datatypes)) datatypes <- head(colnames(gsea[[1]]$pval),2)
      if(length(datatypes)==1) datatypes <- rep(datatypes,2)
      
      shiny::req(gsea[[k]])
      
      S <- as.matrix(unclass(gsea[[k]]$rho[,datatypes,drop=FALSE]))
      P <- as.matrix(unclass(gsea[[k]]$pval[,datatypes,drop=FALSE]))
      multi.score <- gsea[[k]]$multi$score
      multi <- data.frame(unclass(gsea[[k]]$multi))
      colnames(S) <- paste0(colnames(S),".rho")
      colnames(P) <- paste0(colnames(P),".pval")
      colnames(multi) <- paste0("multi.",colnames(multi))
      snames <- stringr::str_trunc(rownames(S),72)      
      sign <- paste0( c("-","+")[1+(sign(S[,1])==1)],
                     c("-","+")[1+(sign(S[,2])==1)] )
      if(full) {
        df <- data.frame( pathway = snames, multi, S, P)
        rownames(df) <- rownames(S)
        df <- df[order(-multi$multi.score),]
      } else {
        df <- data.frame( pathway = snames, multi.score, S)
        rownames(df) <- rownames(S)
        df <- df[order(-multi.score),]
      }
      numeric.cols <- grep("score|pval|p$|q$|rho",colnames(df))
      
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

    table.RENDER2 <- function(full=FALSE) {
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
