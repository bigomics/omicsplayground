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
clustering_table_clustannot_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  clustannot_table_info_text = "In this table, users can check mean correlation values of features in the clusters with respect to the annotation references database selected in the settings."

  TableModuleUI(
    ns("datasets"),
    width = width,
    height = height,
    info.text = clustannot_table_info_text,
    title = "Annotation scores",
    label = "b"
  )

}

#' Server side table code: expression board
#'
#' @param id
#' @param watermark
#'
#' @export
clustering_table_clustannot_server <- function(id,
                                               getClustAnnotCorrelation,
                                               xann_level,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    clustannot_table.RENDER <- shiny::reactive({

      rho = getClustAnnotCorrelation()
      xann_level <- xann_level()
      if(is.null(rho)) return(NULL)

      ##rownames(rho) = shortstring(rownames(rho),50)
      rho.name = shortstring(sub(".*:","",rownames(rho)),60)
      ##rho = data.frame(cbind( name=rho.name, rho))
      df = data.frame( feature=rho.name, round(as.matrix(rho),digits=3))
      rownames(df) = rownames(rho)
      if(xann_level=="geneset") {
        df$feature <- wrapHyperLink(df$feature, rownames(df))
      }

      DT::datatable(
        df, rownames=FALSE, escape = c(-1,-2),
        extensions = c('Buttons','Scroller'),
        selection=list(mode='single', target='row', selected=c(1)),
        class = 'compact hover',
        fillContainer = TRUE,
        options=list(
          dom = 'lfrtip', buttons = c('copy','csv','pdf'),
          ##pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ##scrollY = TRUE,
          ##scrollY = 170,
          scrollY = '23vh',
          scroller=TRUE,
          deferRender=TRUE
        )  ## end of options.list
      ) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    clustannot_table.RENDER_modal <- shiny::reactive({
      dt <- clustannot_table.RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = clustannot_table.RENDER,
      func2 = clustannot_table.RENDER_modal,
      selector = "none"
    )

  }) # end module server
} # end server
