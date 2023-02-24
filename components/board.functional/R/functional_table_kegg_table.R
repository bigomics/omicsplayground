##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


functional_table_kegg_table_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

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

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "Enrichment table",
    label = "b"
  )

}


functional_table_kegg_table_server <- function(id,
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
                      scrollY = "15vh", scroller = TRUE, deferRender = TRUE
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

    table_RENDER_modal <- shiny::reactive({
      dt <- table_RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    my_table <- TableModuleServer(
      "datasets",
      func = table_RENDER,
      func2 = table_RENDER_modal,
      selector = "none"
    )

    return(my_table)

  })  ## end of moduleServer
} ## end of server
