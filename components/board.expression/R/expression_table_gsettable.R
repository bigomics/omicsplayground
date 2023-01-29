##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
expression_table_gsettable_ui <- function(id,
                                          label='') {

  ns <- shiny::NS(id)

  tableWidget(ns("table"))

}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
expression_table_gsettable_server <- function(id,
                                              gx_related_genesets,
                                              height,
                                              width,
                                              watermark=FALSE){
  moduleServer( id, function(input, output, session) {

    gsettable.RENDER <- shiny::reactive({
      df <- gx_related_genesets()
      if (is.null(df)) {
        return(NULL)
      }

      df$geneset <- wrapHyperLink(df$geneset, rownames(df))

      DT::datatable(df,
                    ## class = 'compact cell-border stripe',
                    class = "compact",
                    rownames = FALSE, escape = c(-1, -2),
                    extensions = c("Scroller"),
                    fillContainer = TRUE,
                    options = list(
                      ## dom = 'lfrtip',
                      dom = "frtip",
                      paging = TRUE,
                      pageLength = 16, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                      scrollX = TRUE,
                      ## scrollY = tabV,
                      scrollY = FALSE,
                      scroller = FALSE,
                      deferRender = TRUE,
                      search = list(
                        regex = TRUE,
                        caseInsensitive = TRUE
                        ## search = 'GOBP:'
                      )
                    ), ## end of options.list
                    selection = list(mode = "single", target = "row", selected = NULL)
      ) %>%
        ## formatSignif(1:ncol(df),4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("fx", background = color_from_middle(df$fx, "lightblue", "#f5aeae"))
      # }, server=FALSE)
    })

    gsettable_text <- "By clicking on a gene in the Table <code>I</code>, it is possible to see which genesets contain that gene in this table, and check the differential expression status in other comparisons from the <code>Gene in contrasts</code> plot under the <code>Plots</code> tab."

    gsettable <- shiny::callModule(
      tableModule,
      id = "gsettable",
      func = gsettable.RENDER,
      info.text = gsettable_text, label = "II",
      title = "Gene sets with gene",
      height = height, width = width
    )
    return(gsettable)
  }) #end module server
} #end server