##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_table_rescources_ui <- function(id)
{
  ns <- shiny::NS(id)

  shiny::fillCol(
    flex = c(NA,0.02,1),
    height = 750,
    tags$div(
      HTML("<b>Resource information.</b> Details about the execution times of the methods,
                     dimensions and memory sizes of objects.")
    ),
    shiny::br(),
    shiny::fillRow(
      flex = c(5,1, 2,1, 1.5, 2),
      tableWidget(ns("timings")),
      shiny::br(),
      tableWidget(ns("objectdims")),
      shiny::br(),
      tableWidget(ns("objectsizes")),
      shiny::br()
    )
  )
}


dataview_table_resources_server <- function(id, pgx)
{
  moduleServer(id, function(input, output, session) {
    
    ##================================================================================
    ## Timings
    ##================================================================================

    datatable_timings.RENDER <- shiny::reactive({

      shiny::req(pgx$timings)

      dbg("[datatable_timings.RENDER] reacted")

      ##if(is.null(pgx$timings)) return(NULL)
      D <- data.frame()
      if(!is.null(pgx$timings)) {
        D <- round(pgx$timings[,1:3],digits=3)
        D <- apply(D, 2, function(x) tapply(x, rownames(D), sum))
        catg <- gsub("^\\[|\\].*","",rownames(D))
        metd <- gsub("^.*\\]","",rownames(D))
        D <- data.frame(category=catg, method=metd, D)
      }
      D
    })

    datatable_timings.RENDER <- function(){
      D <- timings_data()
      req(D)
      DT::datatable( D, rownames=FALSE,
                    options = list(dom='tp', pageLength = 100),
                    class = 'compact cell-border stripe hover' ) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    }
    
    timings_text = 'The <b>timings</b> table reports more detailed information about the object dimensions, object sizes and execution times of the methods.'

    datatable_timings <- shiny::callModule(
      tableModule, "timings",
      func = datatable_timings.RENDER,
      info.text = timings_text,
      options = NULL, title='Timings'
    )

    ##================================================================================
    ## Object dimensions
    ##================================================================================

    
    objectdims_data <- reactive({
      shiny::req(pgx$X)
      dims1 <- lapply( pgx, dim)
      lens <- sapply( pgx, length)
      dims2 <- t(sapply( pgx[which(!sapply(dims1,is.null)) ], dim))
      kk <- which(sapply(dims1,is.null))
      dims2 <- rbind(dims2, cbind(lens[kk],0))
        colnames(dims2) = c("nrows","ncols")
      D = data.frame( object=rownames(dims2), dims2, check.names=FALSE)
    })

    objectdims.RENDER <- function() {      
      D <- objectdims_data()
      req(D)
      DT::datatable( D, rownames=FALSE,
                    options = list(dom='t', pageLength = 50),
                    class = 'compact cell-border stripe hover') %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    }

    objectdims_text = 'This table provides details about the data dimensions of objects.'

    datatable_objectdims <- shiny::callModule(
      tableModule, "objectdims",
      func = objectdims.RENDER,
      info.text = objectdims_text,
      options = NULL, title='Object dimensions'
    )

    ##================================================================================
    ## Object sizes
    ##================================================================================

    objectsize_data <- shiny::reactive({
      shiny::req(pgx$name)
      objsize <- sapply(pgx,object.size)
      objsize <- round( objsize/1e6, digits=2)
      data.frame( object=names(pgx), "size.Mb"=objsize, check.names=FALSE)
    })
    
    objectsize.RENDER <- function() {
      D <- objectsize_data()
      req(D)
      DT::datatable( D, rownames=FALSE,
                    options = list(dom='t', pageLength = 50),
                    class = 'compact cell-border stripe hover') %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    }

    objectsize_text = "This table provides information about  about the memory sizes of objects"

    datatable_objectsize <- shiny::callModule(
      tableModule, "objectsize",
      func = objectsize.RENDER,
      options = NULL, title='Object sizes',
      info.text = objectsize_text
    )
      
  })  ## end of moduleServer
} ## end of server

