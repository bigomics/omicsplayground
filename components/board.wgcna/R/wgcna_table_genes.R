##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

wgcna_table_genes_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("geneTable"))
}

wgcna_table_genes_server <- function(id,
                                     wgcna.compute,
                                     selected_module) {
  moduleServer(id, function(input, output, session) {
    geneTable.RENDER <- shiny::reactive({
      out <- wgcna.compute()

      k <- selected_module()
      genes <- out$me.genes[[k]]
      tt <- GENE.TITLE[toupper(genes)]
      rho <- cor(out$datExpr[, genes], out$net$MEs[, k])[, 1]

      df <- data.frame(module = k, gene = genes, me.rho = rho, title = tt)
      numeric.cols <- grep("score|value|ratio|rho", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE, ## escape = c(-1,-2),
        ## filter = 'top',
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact cell-border stripe hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip", ## buttons = c('copy','csv','pdf'),
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE, ## scrollY = TRUE,
          ## scrollY = 170,
          scrollY = "70vh",
          scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })

    geneTable_info <- "Genes in the selected WGCNA module."

    geneTable_module <- shiny::callModule(
      tableModule,
      id = "geneTable",
      func = geneTable.RENDER, ## ns=ns,
      info.text = geneTable_info,
      title = tags$div(
        HTML('<span class="module-label">(d)</span>Module genes')
      ),
      height = c(250, 650)
    )

    return(geneTable_module)
  })
}
