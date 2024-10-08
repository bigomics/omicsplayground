##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_table_rawdata_ui <- function(
    id,
    width,
    height,
    title,
    caption,
    info.text) {
  ns <- shiny::NS(id)


  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_full_table"),
        "show full annotation",
        FALSE
      ),
      "Show full table. Show all feature annotation columns."
    )
  )

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    caption = caption,
    options = options,
    width = width,
    height = height,
    title = title
  )
}

dataview_table_rawdata_server <- function(id,
                                          pgx,
                                          r.gene = reactive(""),
                                          r.data_type = reactive("counts"),
                                          r.samples = reactive(""),
                                          r.groupby = reactive(""),
                                          scrollY = "auto") {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      ## get current view of raw_counts

      ## dereference reactives
      gene <- r.gene()
      data_type <- r.data_type()
      samples <- r.samples()
      groupby <- r.groupby()

      parse_sample <- function(data) {
        if (samples[1] == "") samples <- colnames(data)
        samples <- intersect(colnames(data), samples)
        parsed_data <- data[, samples, drop = FALSE]
      }

      if (is.null(gene) || gene == "" || is.na(gene)) {
        gene <- rownames(pgx$X)[1]
      }

      logx <- parse_sample(pgx$X)
      if (data_type == "counts") {
        x <- parse_sample(pgx$counts)
      } else {
        x <- logx
      }
      x0 <- x

      # Handle to avoid errors on dataset change
      shiny::req(any(rownames(x) == gene))

      k <- which(rownames(x) == gene)
      rho <- cor(t(logx), logx[k, ], use = "pairwise")[, 1]
      rho <- rho[match(rownames(x), names(rho))]
      rho <- round(rho, digits = 3)
      sdx <- round(matrixStats::rowSds(x, na.rm = TRUE), digits = 3)
      avg <- round(rowMeans(x, na.rm = TRUE), digits = 3)

      group <- NULL
      if (groupby %in% colnames(pgx$Y)) {
        group <- pgx$Y[colnames(x), groupby]
      }
      if (length(samples) > 500 && groupby == "<ungrouped>") {
        group <- pgx$model.parameters$group
      }
      do.grouped <- (groupby != "<ungrouped>")
      if (do.grouped && !is.null(group)) {
        allgroups <- sort(unique(group))
        newx <- c()
        for (gr in allgroups) {
          mx <- rowMeans(x[, which(group == gr), drop = FALSE], na.rm = TRUE)
          newx <- cbind(newx, mx)
        }
        rownames(newx) <- rownames(x)
        colnames(newx) <- paste0("avg.", allgroups, "")
        x <- newx
      }

      x <- round(as.matrix(x), digits = 3)
      x95 <- quantile(as.vector(x0[which(x0 > 0)]), probs = 0.95)
      x99 <- quantile(as.vector(x0[which(x0 > 0)]), probs = 0.99)

      if (NCOL(x) == 0 || nrow(x) == 0) {
        return(NULL)
      }

      ## create final dataframe
      pp <- rownames(x)
      annot <- pgx$genes[pp, ]
      if (!input$show_full_table) {
        annot <- annot[, c("feature", "symbol", "gene_title")]
      }
      annot$gene_title <- substring(annot$gene_title, 1, 50)
      if (mean(head(annot$feature, 1000) == head(annot$symbol, 1000), na.rm = TRUE) > 0.8) {
        annot$symbol <- NULL
      }

      df <- data.frame(
        annot,
        rho = rho,
        SD = sdx,
        AVG = avg,
        as.matrix(x),
        check.names = FALSE
      )

      ## if symbol and feature as same, drop symbol column
      df <- df[order(-df$rho, -df$SD), , drop = FALSE]

      list(
        df = df,
        x95 = x95,
        x99 = x99
      )
    })

    rawdataTable.RENDER <- function() {
      dt <- table_data()
      req(dt, dt$df)

      numcols <- grep("gene|title", colnames(dt$df), value = TRUE, invert = TRUE)
      tabH <- 700 ## height of table

      DT::datatable(
        dt$df,
        rownames = FALSE,
        fillContainer = TRUE, #
        class = "compact hover",
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = 1),
        options = list(
          dom = "frtip",
          pageLength = 100,
          lengthMenu = c(25, 40, 100, 250, 1000),
          scroller = TRUE,
          scrollY = scrollY,
          scrollResize = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(numcols,
          background = DT::styleColorBar(data = c(0, dt$x99), color = unname(omics_colors("light_blue"))),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    rawdataTable.RENDER_modal <- shiny::reactive({
      dt <- rawdataTable.RENDER()
      dt$df$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = rawdataTable.RENDER,
      func2 = rawdataTable.RENDER_modal,
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
