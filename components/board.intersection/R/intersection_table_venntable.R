##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


intersection_table_venntable_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("venntable"))
}


intersection_table_venntable_server <- function(id,
                                                getSignificanceCalls,
                                                inputData,
                                                level,
                                                getFoldChangeMatrix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    getSignificantFoldChangeMatrix <- shiny::reactive({
      ##
      ## Filters FC matrix with significance and user-defined
      ## intersection region.
      dt <- getSignificanceCalls()
      shiny::req(dt)

      isect <- input$intersection
      fc0 <- getFoldChangeMatrix()$fc
      if (length(isect) == 0) {
        fc1 <- fc0
      } else {
        ## only genes at least significant in one group
        jj <- which(rowSums(dt[, 2:ncol(dt), drop = FALSE] != 0) > 0)
        if (length(jj) == 0) {
          return(NULL)
        }
        dt <- dt[jj, , drop = FALSE]

        ## check same sign
        if (input$include == "up/down") {
          kk <- 1 + match(c("B", "C"), LETTERS[1:10])
          kk <- 1 + match(isect, LETTERS[1:10])
          kk <- intersect(kk, 1:ncol(dt))

          dt1 <- dt[, kk, drop = FALSE]
          jj <- which(rowMeans(sign(dt1) == +1) == 1 |
            (rowMeans(sign(dt1) == -1) == 1))
          dt <- dt[jj, , drop = FALSE]
          remove(dt1)
        }

        ## only genes in the selected intersection
        intersection <- "ABC"
        intersection <- paste0(input$intersection, collapse = "")
        dt <- dt[which(dt$intersection == intersection), , drop = FALSE]
      }

      ## filtered by family/collection
      fc1 <- fc0[intersect(rownames(dt), rownames(fc0)), , drop = FALSE]
      if (nrow(dt) == 1) {
        fc1 <- matrix(fc1, nrow = 1)
        rownames(fc1) <- rownames(dt)
        colnames(fc1) <- colnames(fc0)
      }

      ## filtered by SPLOM selection
      splom.sel <- plotly::event_data("plotly_selected", source = "splom")
      sel.keys <- as.character(splom.sel$key)
      if (1 && length(sel.keys) > 0) {
        sel <- intersect(sel.keys, rownames(fc1))
        fc1 <- fc1[sel, , drop = FALSE]
      }

      ## only active/selected comparisons
      sel <- colnames(dt)[-1]
      kk <- match(sel, gsub(" \\(-\\)", "", colnames(fc1)))
      fc1 <- fc1[, kk, drop = FALSE]

      ## order
      fc1 <- fc1[order(-rowMeans(fc1)), , drop = FALSE]
      fc1 <- round(fc1, digits = 3)
      colnames(fc1) <- LETTERS[1:ncol(fc1)]
      ## fc0 = data.frame(fc0)

      ## add intersection code
      sel <- match(rownames(fc1), rownames(dt))
      fc1 <- data.frame(intersection = dt$intersection[sel], fc = fc1)

      ## filter on user selection
      vv <- input$venntable_intersection
      if (vv != "<all>") {
        sel <- which(fc1$intersection == vv)
        fc1 <- fc1[sel, , drop = FALSE]
      }
      return(fc1)
    })

    venntable.RENDER <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs)

      ## get foldchanges
      fc0 <- getSignificantFoldChangeMatrix() ## isolate??
      if (is.null(fc0) || nrow(fc0) == 0) {
        return(NULL)
      }

      ## add gene name/title
      if (level == "gene") {
        gene <- as.character(ngs$genes[rownames(fc0), "gene_name"])
        gene.tt <- substring(GENE.TITLE[gene], 1, 50)
        gene.tt <- as.character(gene.tt)
        ## fc0 = data.frame( name=name, title=gene.tt, fc0)
        fc0 <- data.frame(name = gene, fc0, check.names = FALSE)
      } else {
        name <- substring(rownames(fc0), 1, 50)
        name[is.na(name)] <- "NA"
        fc0 <- data.frame(name = name, fc0, check.names = FALSE)
      }

      df <- data.frame(fc0, check.names = FALSE)
      nsc <- setdiff(1:ncol(df), 2)
      ## dt <- dt[rownames(fc0),]
      ## D <- cbind(intersection=dt$intersection, D)
      DT::datatable(df,
        class = "compact cell-border stripe",
        rownames = FALSE,
        extensions = c("Scroller"), selection = "none",
        fillContainer = TRUE,
        options = list(
          ## dom = 'lfrtip',
          dom = "tip",
          ## buttons = c('copy','csv','pdf'),
          ## pageLength = 20,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
          ## columnDefs = list(list(targets=nsc, searchable = FALSE)),
          scrollX = TRUE,
          ## scrollY = 150,
          scrollY = "70vh",
          scroller = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    })


    info_text <- "Table of genes in selected intersection."

    venntable_opts <- shiny::tagList(
      shiny::selectInput(ns("venntable_intersection"), "Filter intersection:", choices = NULL)
    )

    shiny::callModule(
      tableModule,
      id = "venntable",
      func = venntable.RENDER,
      ## caption = venntable_buttons,
      options = venntable_opts,
      title = "INTERSECTION",
      label = "c",
      info.text = info_text,
      ## info.width = "400px",
      height = c(260, 750),
      width = c("auto", 1200)
    )
  }) ## end of moduleServer
} ## end of server
