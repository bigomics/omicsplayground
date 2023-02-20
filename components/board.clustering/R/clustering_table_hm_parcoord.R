##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' UI code for table code: clustering board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
clustering_table_hm_parcoord_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("tablemod"))
}

#' Server side table code: clustering board
#'
#' @param id
#' @param watermark
#'
#' @export
clustering_table_hm_parcoord_server <- function(id = "hm_parcoord_table",
                                                hm_parcoord.selected,
                                                watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    hm_parcoord_table.RENDER <- shiny::reactive({
      hm_parcoord.selected <- hm_parcoord.selected()

      mat = hm_parcoord.selected$mat
      clust = hm_parcoord.selected$clust
      df <- data.frame(cluster=clust, mat, check.names=FALSE)
      numeric.cols <- 2:ncol(df)
      DT::datatable(
        df, rownames=TRUE, ## escape = c(-1,-2),
        extensions = c('Buttons','Scroller'),
        selection=list(mode='single', target='row', selected=NULL),
        class = 'compact hover',
        fillContainer = TRUE,
        options=list(
          dom = 'lfrtip', ##buttons = c('copy','csv','pdf'),
          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ##scrollY = TRUE,
          ##scrollY = 170,
          scrollY = '70vh',
          scroller=TRUE, deferRender=TRUE
        )  ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols,3) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    hm_parcoord_table_info = "In this table, users can check mean expression values of features across the conditions for the selected genes."

    shiny::callModule(
      tableModule,
      id = "tablemod",
      func = hm_parcoord_table.RENDER,
      info.text = hm_parcoord_table_info,
      ##options = clustannot_table_opts,
      title="Selected genes",
      label="b",
      height = c(240,700),
      width=c('auto',1000),
      ##caption = clustannot_caption
    )
  }) # end module server
} # end server
