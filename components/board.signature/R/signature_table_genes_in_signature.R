##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

signature_table_genes_in_signature_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("table"))
}

signature_table_genes_in_signature_server <- function(id,
                                                      getEnrichmentGeneTable,
                                                      tabH) {
  moduleServer(id, function(input, output, session) {

    enrichmentGeneTable.RENDER <- shiny::reactive({

      df <- getEnrichmentGeneTable()
      shiny::req(df)

      color_fx = as.numeric(df[,3:ncol(df)])
      color_fx[is.na(color_fx)] <- 0  ## yikes...

      numeric.cols <- colnames(df)[3:ncol(df)]
      numeric.cols

      DT::datatable(df, class='compact cell-border stripe',
                    rownames=FALSE,
                    extensions = c('Scroller'),
                    ## selection='none',
                    selection = list(mode='single', target='row', selected=NULL),
                    fillContainer=TRUE,
                    options=list(
                      dom = 'lrftip',
                      ## pageLength = 40, ##lengthMenu = c(20, 30, 40, 60, 100, 250),
                      scrollX = TRUE, scrollY = tabH, scroller=TRUE,
                      deferRender=FALSE
                    )) %>%  ## end of options.list
        DT::formatSignif(numeric.cols,4) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
        DT::formatStyle(
          numeric.cols,
          background = color_from_middle(color_fx,'lightblue','#f5aeae'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center')
    })

    info.text2 = "<b>Gene table.</b> Genes of the current signature corresponding to the selected contrast. Genes are sorted by decreasing (absolute) fold-change."
    enrichmentGeneTable <- shiny::callModule(
      tableModule,
      id = "table",
      func = enrichmentGeneTable.RENDER,
      info.text = info.text2,
      caption2 = info.text2,
      title = tags$div(
        HTML('<span class="module-label">(b)</span>Genes in signature')
      ),
      height = c(360,700)
    )
    return(enrichmentGeneTable)
  })
}
