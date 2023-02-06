##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wgcna_table_enrichment_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("enrichTable"))
}

wgcna_table_enrichment_server <- function(id,
                                          enrich_table) {
  moduleServer(id, function(input, output, session) {

    enrichTable.RENDER <- shiny::reactive({

      df <- enrich_table()
      numeric.cols <- grep("score|value|ratio",colnames(df))

      DT::datatable(
        df, rownames=FALSE,
        extensions = c('Buttons','Scroller'),
        selection = list(mode='single', target='row', selected=NULL),
        class = 'compact cell-border stripe hover',
        fillContainer = TRUE,
        options=list(
          dom = 'lfrtip',
          scrollX = TRUE,
          scrollY = '70vh',
          scroller=TRUE, deferRender=TRUE
        )  ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols,3) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    enrichTable_info = "In this table, users can check mean expression values of features across the conditions for the selected genes."

    enrichTable_module <- shiny::callModule(
      tableModule, id = "enrichTable",
      func  = enrichTable.RENDER,
      info.text = enrichTable_info,
      title = tags$div(
        HTML('<span class="module-label">(e)</span>Module enrichment')
      ),
      height = c(250,650)
    )

    return(enrichTable_module)
  })
}
