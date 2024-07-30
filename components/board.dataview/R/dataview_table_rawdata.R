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

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    caption = caption,
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

      if (data_type %in% c("counts", "abundance")) {
        x <- parse_sample(pgx$counts)
        x_cpm <- parse_sample(pgx$X)
      } else if (data_type == "CPM") {
        x <- edgeR::cpm(pgx$counts, log = FALSE)
      } else {
        ## log2CPM
        x <- parse_sample(pgx$X)
      }

      x0 <- x

      ## ------------------ select samples


      ## Quickly (?) calculated correlation to selected gene

      ## compute statistics
      rho <- sdx <- avg <- NULL

      if (data_type == "counts") {
        xgenes <- pgx$genes[rownames(x), "gene_name"]
        k <- which(xgenes == gene)

        xgenes_cpm <- pgx$genes[rownames(x_cpm), "gene_name"]
        k_cpm <- which(xgenes_cpm == gene)
      } else {
        ## log2CPM
        xgenes <- pgx$genes[rownames(x), "gene_name"]
        k <- which(xgenes == gene)
      }

      if (data_type == "counts") {
        # compute the geometric mean, exp(mean(log(x+1)))
        logx <- log(x[rownames(x), ] + 1)

        # correlation should be equal between counts and logCPM, use logCPM
        rho <- cor(t(x_cpm[, colnames(x)]), x_cpm[k_cpm, colnames(x)], use = "pairwise")[, 1]
        rho <- round(rho[rownames(x_cpm)], digits = 3)
        rho <- rho[match(rownames(x), names(rho))]
        names(rho) <- rownames(x)

        # geometric std deviation
        sdx <- round(exp(apply(logx[, samples], 1, sd, na.rm = TRUE)), digits = 3)
        # geometric mean
        avg <- round(exp(rowMeans(logx, na.rm = TRUE)), digits = 3)
      } else {
        # compute the geometric mean, mean(x)
        logx <- x
        rho <- cor(t(logx[, samples]), logx[k, samples], use = "pairwise")[, 1]
        rho <- round(rho[rownames(logx)], digits = 3)
        sdx <- round(apply(logx[, samples], 1, sd, na.rm = TRUE), digits = 3)
        avg <- round(rowMeans(logx, na.rm = TRUE), digits = 3)
      }

      group <- NULL
      if (groupby %in% colnames(pgx$Y)) {
        group <- pgx$Y[colnames(logx), groupby]
      }
      if (length(samples) > 500 && groupby == "<ungrouped>") {
        group <- pgx$model.parameters$group
      }
      do.grouped <- (groupby != "<ungrouped>")
      if (do.grouped && !is.null(group)) {
        allgroups <- sort(unique(group))
        newx <- c()
        for (gr in allgroups) {
          if (data_type == "counts") {
            mx <- exp(rowMeans(logx[, which(group == gr), drop = FALSE], na.rm = TRUE))
          } else {
            mx <- rowMeans(logx[, which(group == gr), drop = FALSE], na.rm = TRUE)
          }

          newx <- cbind(newx, mx)
        }
        rownames(newx) <- rownames(logx)
        colnames(newx) <- paste0("avg.", allgroups, "")
        x <- newx
      }

      x <- round(as.matrix(x), digits = 3)
      x95 <- quantile(as.vector(x0[which(x0 > 0)]), probs = 0.95)
      x99 <- quantile(as.vector(x0[which(x0 > 0)]), probs = 0.99)

      if (NCOL(x) == 0 || nrow(x) == 0) {
        return(NULL)
      }

      dbg("[dataview_rawdata:table_data] create dataframe")
      pp <- rownames(x)
      feature <- pgx$genes[pp, "feature"]
      symbol <- pgx$genes[pp, "symbol"]      
      gene.title <- pgx$genes[pp, "gene_title"]
      gene.title <- substring(gene.title, 1, 50)
      
      df <- data.frame(
        gene = feature,
        symbol = symbol,          
        title = gene.title,
        rho = rho,
        SD = sdx,
        AVG = avg,
        as.matrix(x),
        check.names = FALSE
      )
      
      ## if symbol and feature as same, drop symbol column
      if(mean(head(df$gene,1000)==head(df$symbol,1000)) > 0.8) {
        df$symbol <- NULL
      }

      if (DATATYPEPGX == "proteomics") {
        colnames(df)[1] <- "protein"
      }
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
