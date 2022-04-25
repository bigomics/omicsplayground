##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_table_contrasts_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("tbl"))  
}


dataview_table_contrasts_server <- function(id,
                                            pgx,
                                            r.samples = reactive("")
                                            )
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns
    
    contrasts_data <- shiny::reactive({
      
      dbg("[DataView:contrasts:contrast_data] reacted")

      shiny::req(pgx$Y,pgx$model.parameters)
      shiny::req(r.samples(),!is.null(input$ctbygroup))

      dbg("[DataView:contrasts:contrast_data] called!")

      ## dereference reactives
      samples   <- r.samples()
      ctbygroup <- input$ctbygroup
      
      dt <- NULL
      if(ctbygroup=="group") {
        ct <- pgx$model.parameters$contr.matrix
        kk <- which(rownames(ct) %in% pgx$model.parameters$group[samples])
        dt <- ct[kk,,drop=FALSE]
      } else {
        dt <- pgx$model.parameters$exp.matrix[samples,,drop=FALSE]
      }
      dt <- sign(dt)
      dt[dt==0] <- NA
      dt
    }) ## %>% bindCache(pgx$Y, r.samples(), ctbygroup)

    
    table.RENDER <- function() {
      
      dt <- contrasts_data()
      req(dt)

      dbg("[DataView:contrasts:table.RENDER] reacted!")
      
      tabH = 600  ## height of tables
      colnames(dt) <- sub("[_. ]vs[_. ]","\nvs ",colnames(dt))
      
      DT::datatable( dt,
                    class = 'compact cell-border stripe hover',
                    rownames = TRUE,
                    extensions = c('Buttons','Scroller'),
                    selection = list(mode='single', target='row', selected=1),
                    options=list(
                      dom = 'lfrtip',
                      scroller=TRUE, scrollX = TRUE, scrollY = tabH,
                      deferRender=TRUE,
                      autoWidth = TRUE
                    )) %>%
        DT::formatStyle(0, target='row', fontSize='12px', lineHeight='70%') %>%
        DT::formatStyle(colnames(dt),
                        background = color_from_middle(c(-1,1), 'lightblue', '#f5aeae'),
                        backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',
                        backgroundPosition = 'center')

    }

    info_text = "<b>Contrast table.</b> Table summarizing the contrasts of all comparisons. Here, you can check which samples belong to which groups for the different comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison."

    opts = shiny::tagList(
      withTooltip( shiny::radioButtons(
        ns('ctbygroup'),
        "Show by:", choices=c("sample","group")),
        "Show contrasts by group or by samples.",
        placement="right", options = list(container = "body")
        )
    )
    
    contrastTable <- shiny::callModule(
      tableModule, "tbl",
      func = table.RENDER,
      csvFunc = contrasts_data,
      options = opts,
      title = "Contrast table",
      filename = "contrasts.csv",
      info.text = info_text
    )


  })  ## end of moduleServer
} ## end of server

