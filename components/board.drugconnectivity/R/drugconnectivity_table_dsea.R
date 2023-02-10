##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


drugconnectivity_table_dsea_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("dsea_table"))
}


drugconnectivity_table_dsea_server <- function(id,
                                               getActiveDSEA)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns

    table_data <- shiny::reactive({
      dsea <- getActiveDSEA()
      shiny::req(dsea)

      dt <- dsea$table
      return(dt)
    })

    table.RENDER <- function() {
      res <- table_data()
      res$moa <- shortstring(res$moa, 60)
      res$target <- shortstring(res$target, 30)
      res$drug <- shortstring(res$drug, 60)

      colnames(res) <- sub("moa", "MOA", colnames(res))
      DT::datatable(res,
                    rownames = FALSE,
                    class = "compact cell-border stripe hover",
                    extensions = c("Scroller"),
                    selection = list(mode = "single",
                                     target = "row",
                                     selected = NULL),
                    fillContainer = TRUE,
                    options = list(
                      dom = "lfrtip",
                      scroller = TRUE, scrollX = TRUE,
                      scrollY = "70vh",
                      deferRender = TRUE
                    )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px",
                        lineHeight = "70%") %>%
        DT::formatStyle("NES",
                        background = color_from_middle(res[, "NES"],
                                                       "lightblue",
                                                       "#f5aeae"),
                        backgroundSize = "98% 88%",
                        backgroundRepeat = "no-repeat",
                        backgroundPosition = "center"
        )
    }

    info_text <- strwrap("<b>Enrichment table.</b> Enrichment is calculated by
                         correlating your signature with known drug profiles
                         from the L1000 database. Because the L1000 has multiple
                         perturbation experiment for a single drug, drugs are
                         scored by running the GSEA algorithm on the
                         contrast-drug profile correlation space. In this way,
                         we obtain a single score for multiple profiles of a
                         single drug.")

    table.opts <- shiny::tagList()
    dsea_table <- shiny::callModule(
      tableModule,
      id = "dsea_table",
      func = table.RENDER,
      options = table.opts,
      info.text = info_text,
      selector = "single",
      title = "Enrichment table",
      height = c(360, 700)
    )

    return(dsea_table)

  })  ## end of moduleServer
} ## end of server
