##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


connectivity_table_connectivityScoreTable_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("table"))
}


connectivity_table_connectivityScoreTable_server <- function(id,
                                               inputData,
                                               getFilteredKeggTable,
                                               fa_contrast,
                                               tabH)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns

    table_data <- shiny::reactive({
      res <- list(
        pgx = inputData(),
        df = getFilteredKeggTable(),
        fa_contrast = fa_contrast
      )
      return(res)
    })

    table_RENDER <- function() {
      res <- table_data()
      pgx <- res$pgx
      df <- res$df
      comparison <- res$fa_contrast

      if (is.null(pgx$meta.go)) return(NULL)
      if (is.null(comparison)) return(NULL)
      if (is.null(df)) return(NULL)
      if (nrow(df) == 0) return(NULL)

      ## add hyperlink
      url <- paste0("https://www.genome.jp/kegg-bin/show_pathway?map=hsa",
                    df$kegg.id,
                    "&show_description=show")
      df$kegg.id <- paste0("<a href='", url, "' target='_blank'>",
                           df$kegg.id, "</a>")

      numeric.cols <- colnames(df)[which(sapply(df, is.numeric))]

      DT::datatable(df,
                    rownames = FALSE, escape = c(-1, -2),
                    class = "compact cell-border stripe hover",
                    extensions = c("Scroller"),
                    selection = list(mode = "single", target = "row",
                                     selected = 1),
                    fillContainer = TRUE,
                    options = list(
                      dom = "lfrtip",
                      scrollX = TRUE,
                      scrollY = tabH, scroller = TRUE, deferRender = TRUE
                    ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px",
                        lineHeight = "70%") %>%
        DT::formatStyle("logFC",
                        background = color_from_middle(df[, "logFC"],
                                                       "lightblue",
                                                       "#f5aeae"),
                                                       backgroundSize = "98% 88%",
                        backgroundRepeat = "no-repeat",
                        backgroundPosition = "center"
        )
    }

    info_text <- strwrap("<strong>Enrichment table.</strong> The table is
                         interactive; enabling user to sort on different
                         variables and select a pathway by clicking on the row
                         in the table. The scoring is performed by considering
                         the total number of genes in the pathway (n), the
                         number of genes in the pathway supported by the contrast
                         profile (k), the ratio of k/n, and the ratio of
                         |upregulated or downregulated genes|/k. Additionally,
                         the table contains the list of the upregulated and
                         downregulated genes for each pathway and a q value from
                         the Fisher’s test for the overlap.")

    table_opts <- shiny::tagList()

    my_table <- shiny::callModule(
      tableModule,
      id = "table",
      label = "",
      func = table_RENDER,
      options = table_opts,
      info.text = info_text,
      info.width = '350px',
      title = "Enrichment table",
      height = c(270, 700)
    )

    return(my_table)

  })  ## end of moduleServer
} ## end of server
