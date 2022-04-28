##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_table_samples_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("tbl"))  
}


dataview_table_samples_server <- function(id,
                                          pgx,
                                          r.samples = reactive("")
                                          )
{
  moduleServer(id, function(input, output, session) {
    
    table_data <- shiny::reactive({
      shiny::req(pgx$Y,pgx$samples,r.samples())
      samples <- r.samples()
      dt <- pgx$samples[samples,,drop=FALSE]
      dt
    }) 
   
    table.RENDER <- function() {
      dt <- table_data()
      req(dt)      
      DT::datatable( dt,
                    class = 'compact cell-border stripe hover',
                    rownames = TRUE,
                    extensions = c('Buttons','Scroller'),
                    selection = list(mode='single', target='row', selected=1),
                    options=list(
                      dom = 'lfrtip',
                      scroller=TRUE, scrollX = TRUE, scrollY = 150,
                      deferRender=TRUE
                    )) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')         
    } 

    modal_table.RENDER <- function() {
      dt <- table_data()
      req(dt)      
      DT::datatable( dt,
                    class = 'compact cell-border stripe hover',
                    rownames = TRUE,
                    extensions = c('Buttons','Scroller'),
                    selection = list(mode='single', target='row', selected=1),
                    options=list(
                      dom = 'lfrtip',
                      scroller=TRUE, scrollX = TRUE, scrollY = 600,
                      deferRender=TRUE
                    )) %>%
        DT::formatStyle(0, target='row', fontSize='20px', lineHeight='70%')         
    }
      
    info_text = "<b>Sample information table.</b> Phenotype information about the samples. Phenotype variables
                 starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."

    shiny::callModule(
      tableModule, "tbl", label="",
      func = table.RENDER,
      func2 = modal_table.RENDER,
      title = "Sample information",
      filename = "samples.csv",
      info.text = info_text,
      caption2 = info_text
      ##height = c(280,750),
      ##width=c('auto','100%')
    )
      
  })  ## end of moduleServer
} ## end of server

