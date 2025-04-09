##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
pcsf_table_centrality_ui <- function(
    id,
    title,
    info.text,
    caption,
    width,
    height) {
  ns <- shiny::NS(id)

  table_opts <- shiny::tagList()

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    ##    options = table_opts,
    title = title,
  )
}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
pcsf_table_centrality_server <- function(id,
                                         pgx,
                                         r_contrast,
                                         r_pcsf) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    table_data <- function() {
      
      df <- playbase::pgx.getPCSFcentrality(
        pgx,
        contrast = r_contrast(),
        pcsf = r_pcsf(),
        n = 100,
        plot = FALSE
      )

      ## warning. sometimes NaN
      df$centrality[is.nan(df$centrality)] <- 0
      
      return(df)
    }

    
    table.RENDER <- function(full=FALSE) {

      df <- table_data()
      
      if(full==FALSE) {
        cols <- c("symbol","gene_title","centrality","logFC")
        cols <- intersect(cols, colnames(df))
        df <- df[,cols, drop=FALSE]
      }      
      #num.cols <- match(c("centrality", "logFC"), colnames(df))
      num.cols <- c("centrality", "logFC")

      dt <- ui.DataTable(
        df,
        rownames = FALSE,
        num.cols = num.cols,
        color.cols = num.cols,
        substr.cols = c("gene_title"),
        substr.len = 50
      )
            
      return(dt)
    }

    table.RENDER_modal <- function() {
      dt <- table.RENDER(full=TRUE)
      dt$x$options$scrollY <- SCROLLY_MODAL
      return(dt)
    }

    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "none"
    )
  }) # end module server
} # end server
