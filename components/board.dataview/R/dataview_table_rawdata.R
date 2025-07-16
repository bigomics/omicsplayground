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
      gene <- r.gene()
      data_scale <- r.data_type()
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
        # So old datasets work (they can be missaligned)
        jj <- which(rownames(pgx$X) %in% rownames(pgx$counts))
        dt <- pgx$counts[rownames(pgx$X)[jj], ]
        ##
        x <- parse_sample(dt)
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

      if (NCOL(x) == 0 || nrow(x) == 0) return(NULL)

      # create final dataframe
      annot <- pgx$genes[rownames(x), ]
      cl <- c("feature", "symbol", "gene_title")
      if (!input$show_full_table) annot <- annot[, cl]

      if ("human_ortholog" %in% colnames(annot))
        colnames(annot)[colnames(annot) == "human_ortholog"] <- "ortholog"

      # hide symbol column if symbol is feature
      if (mean(head(annot$feature, 1000) == head(annot$symbol, 1000), na.rm = TRUE) > 0.8) {
        annot$symbol <- NULL
      }

      pct.na <- round(rowMeans(is.na(pgx$counts[rownames(annot), ])) * 100, 1)

      df <- data.frame(
        annot,
        rho = rho,
        pct.missingness = pct.na,
        SD = sdx,
        AVG = avg,
        as.matrix(x),
        check.names = FALSE
      )
      
      # if symbol and feature as same, drop symbol column
      df <- df[order(-df$rho, -df$SD), , drop = FALSE]

      list(df = df, groupby = groupby, samples = samples, data_scale = data_scale, x95 = x95, x99 = x99)

    })

    rawdataTable.RENDER <- function() {
      dt <- table_data()      
      req(dt, dt$df)
      DF <- dt$df
      groupby <- dt$groupby
      samples <- dt$samples
      data_scale <- dt$data_scale

      numcols <- grep("gene|title", colnames(DF), value = TRUE, invert = TRUE)
      tabH <- 700

      na.map=NULL; rm.cols=NULL
      is.imp <- sum(is.na(pgx$counts))>0 && sum(is.na(pgx$X))==0
      if (is.imp && data_scale == "log2") {
        if(groupby != "<ungrouped>" && groupby %in% colnames(pgx$samples)) {
          group <- pgx$samples[samples, groupby]
          if (!is.null(group)) {
            allgroups <- sort(unique(group))
            i=1; na.map=list()
            for(i in 1:length(allgroups)) {
              jj <- which(pgx$samples[samples, groupby] == allgroups[i])
              samples1 <- intersect(samples, rownames(pgx$samples)[jj]) 
              counts <- pgx$counts[rownames(DF), samples1, drop = FALSE]
              nas <- apply(counts, 2, function(x) unname(which(is.na(x))))
              nas <- nas[sapply(nas, length) > 0]
              na.map[[i]] <- unique(unlist(unname(nas)))
              names(na.map)[i] <- paste0("avg.", allgroups[i], "")
              rm(nas)
            }
          }
        } else {
          counts <- pgx$counts[rownames(DF), , drop = FALSE]
          na.map <- apply(counts, 2, function(x) unname(which(is.na(x))))
          na.map <- na.map[sapply(na.map, length) > 0]
        }
        i=1
        for(i in 1:length(na.map)) {
          imp.info <- ifelse(seq_len(nrow(DF)) %in% na.map[[i]], "yes", "no")
          DF <- cbind(DF, imp.info)
          colnames(DF)[ncol(DF)] <- paste0(names(na.map)[i], ".impinfo")
        }
        rm.cols <- grep(".impinfo", colnames(DF)) - 1 # zero-based index
      }
      
      DTable <- DT::datatable(
        DF,
        rownames = FALSE,
        fillContainer = TRUE,
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
          deferRender = TRUE,
          columnDefs = list(list(targets = rm.cols, visible = FALSE))
        )
      );
      if (!is.null(na.map)) {
        i=1
        for(i in 1:length(na.map)) {
          ss <- paste0(names(na.map)[i], ".impinfo")
          DTable <- DTable %>% DT::formatStyle(names(na.map)[i], valueColumns = ss,
            backgroundColor = DT::styleEqual("yes","#FFE4E1"))
        }
      }
      DTable <- DTable %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle(
          numcols,
          background = DT::styleColorBar(data = c(0, dt$x99), color = unname(omics_colors("light_blue"))),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
      DTable
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
