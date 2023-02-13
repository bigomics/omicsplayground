##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

enrichment_table_enrichment_analysis_ui <- function(id) {
  ns <- shiny::NS(id)

  tableWidget(ns("gseatable"))
}

enrichment_table_enrichment_analysis_server <- function(id,
                                                        getFilteredGeneSetTable) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    gseatable.RENDER <- shiny::reactive({
      rpt <- getFilteredGeneSetTable()
      if (is.null(rpt)) {
        return(NULL)
      }
      if (nrow(rpt) == 0) {
        return(NULL)
      }

      if (!("GS" %in% colnames(rpt))) rpt <- cbind(GS = rownames(rpt), rpt)
      if ("GS" %in% colnames(rpt)) rpt$GS <- shortstring(rpt$GS, 72)
      if ("size" %in% colnames(rpt)) rpt$size <- as.integer(rpt$size)

      fx <- NULL
      fx.col <- grep("score|fx|fc|sign|NES|logFC", colnames(rpt))[1]
      if (length(fx.col) > 0) fx <- rpt[, fx.col]

      jj <- which(sapply(rpt, is.numeric))
      if (length(jj) > 0) rpt[, jj] <- round(rpt[, jj], digits = 4)
      jj <- which(sapply(rpt, is.character) | sapply(rpt, is.factor))
      if (length(jj) > 0) rpt[, jj] <- apply(rpt[, jj, drop = FALSE], 2, shortstring, 100)

      if (!input$gs_showqvalues) {
        rpt <- rpt[, grep("^q[.]|^q$", colnames(rpt), invert = TRUE)]
      }

      ## wrap genesets names with known links.
      rpt$GS <- wrapHyperLink(rpt$GS, rownames(rpt))
      selectmode <- "single"

      is.numcol <- sapply(rpt, is.numeric)
      numcols <- which(is.numcol & !colnames(rpt) %in% c("size"))
      numcols <- colnames(rpt)[numcols]

      colnames(rpt) <- sub("GS", "geneset", colnames(rpt))

      DT::datatable(rpt,
        class = "compact cell-border stripe hover",
        rownames = FALSE,
        escape = c(-1, -5),
        extensions = c("Scroller"),
        fillContainer = TRUE,
        selection = list(mode = selectmode, target = "row", selected = NULL),
        options = list(
          dom = "frtip",
          paging = TRUE,
          pageLength = 15, ##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          scrollX = TRUE,
          scrollY = FALSE,
          scroller = FALSE,
          deferRender = TRUE,
          search = list(
            regex = TRUE,
            caseInsensitive = TRUE
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numcols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(fx.col,
          background = color_from_middle(fx, "lightblue", "#f5aeae")
        )
    })

    gseatable_text <- paste("Similar to the differential gene expression analysis, users can perform differential expression analysis on a geneset level that is referred as gene set enrichment analysis. To ensure statistical reliability, the platform performs the gene set enrichment analysis using multiple methods, including", a_Spearman, ", ", a_GSVA, ", ", a_ssGSEA, ", ", a_Fisher, ", ", a_GSEA, ", ", a_camera, " and ", a_fry, ".<br><br>The combined result from the methods is displayed in this table, where for each geneset the <code>meta.q</code> corresponds to the highest <code>q</code> value provided by the methods and the number of <code>stars</code> indicate how many methods identified the geneset as significant (<code>q < 0.05</code>). The table is interactive; users can sort it by <code>logFC</code>, <code>meta.q</code> and <code>starts</code>. Additionally, the list of genes in that geneset are displayed in the second table on the right. Users can filter top N = {10} differently enriched gene sets in the table by clicking the <code>top 10 gene sets</code> from the table <i>Settings</i>.")

    gseatable_opts <- shiny::tagList(
      withTooltip(shiny::checkboxInput(ns("gs_showqvalues"), "show indivivual q-values", FALSE),
        "Show all q-values of each individual statistical method in the table.",
        placement = "top", options = list(container = "body")
      )
    )

    gseatable <- shiny::callModule(
      tableModule,
      id = "gseatable",
      func = gseatable.RENDER,
      info.text = gseatable_text,
      options = gseatable_opts,
      title = tags$div(
        HTML('<span class="module-label">(I)</span>Enrichment analysis')
      ),
      info.width = "500px",
      height = c(285, 700),
      selector = "single"
    )

    return(gseatable)
  })
}
