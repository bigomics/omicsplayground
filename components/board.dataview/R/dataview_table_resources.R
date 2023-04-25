##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


dataview_table_rescources_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_column_wrap(
    width = 1,
    height = "100%",
    heights_equal = "row",          
    style = htmltools::css(grid_template_columns = "4fr 4fr 4fr"),
    tableWidget(ns("timings")),
    tableWidget(ns("pgxobject")),
    tableWidget(ns("objects"))    
  )
}


dataview_table_resources_server <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ## ================================================================================
    ## Timings
    ## ================================================================================

    timings_data <- shiny::reactive({
      shiny::req(pgx$timings)

      ## if(is.null(pgx$timings)) return(NULL)
      D <- data.frame()
      if (!is.null(pgx$timings)) {
        D <- round(pgx$timings[, 1:3], digits = 3)
        D <- apply(D, 2, function(x) tapply(x, rownames(D), sum))
        catg <- gsub("^\\[|\\].*", "", rownames(D))
        metd <- gsub("^.*\\]", "", rownames(D))
        D <- data.frame(category = catg, method = metd, D)
      }
      D
    })

    datatable_timings.RENDER <- function() {
      D <- timings_data()
      req(D)
      DT::datatable(D,
        rownames = FALSE,
        options = list(dom = "t", pageLength = 50),
        class = "compact hover"
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    timings_text <- "The <b>timings</b> table reports more detailed information about the object dimensions, object sizes and execution times of the methods."

    datatable_timings <- shiny::callModule(
      tableModule, "timings",
      func = datatable_timings.RENDER,
      info.text = timings_text,
      options = NULL, title = "Timings"
    )

    ## ================================================================================
    ## PGX Object dimensions
    ## ================================================================================

    pgx_data <- reactive({
      shiny::req(pgx$X)
      dims1 <- lapply(pgx, dim)
      lens  <- sapply(pgx, length)
      sel.matrix <- names(pgx)[which(!sapply(dims1, is.null))]
      dims2 <- do.call(rbind, dims1[sel.matrix])
      kk <- which(sapply(dims1, is.null))
      dims2 <- rbind(dims2, cbind(lens[kk], 0))
      colnames(dims2) <- c("nrows", "ncols")

      objsize <- sapply(pgx, object.size)
      objsize <- round(objsize / 1e6, digits = 2)
      objsize <- objsize[rownames(dims2)]

      data.frame(object = rownames(dims2), dims2, "size.Mb" = objsize, check.names = FALSE)
    })

    pgx.RENDER <- function() {
      D <- pgx_data()
      D <- D[order(-D$size.Mb),]
      req(D)
      DT::datatable(D,
        rownames = FALSE,
        options = list(dom = "t", pageLength = 50),
        class = "compact hover"
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    pgx_text <- "This table provides details about the pgx object."

    datatable_pgxdims <- shiny::callModule(
      tableModule, "pgxobject",
      func = pgx.RENDER,
      info.text = pgx_text,
      options = NULL,
      title = "PGX slot sizes"
    )

    ## ================================================================================
    ## PGX Object dimensions
    ## ================================================================================

    object_data <- reactive({
      shiny::req(pgx$X)
      obj <- ls(envir=.GlobalEnv)
      sizes <- sapply(obj,function(s) object.size(get(s)))
      sizes.Mb <- round(as.numeric(sizes) / 1024**2, digits=2)
      names(sizes) <- names(sizes.Mb) <- obj
      sizes.Mb <- sort(sizes.Mb, decreasing=TRUE)
      data.frame(object = names(sizes.Mb), size.Mb = sizes.Mb)
    })

    object.RENDER <- function() {
      D <- object_data()
      req(D)
      D <- head(D, 35)
      DT::datatable(D,
        rownames = FALSE,
        options = list(dom = "t", pageLength = 50),
        class = "compact hover"
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    object_text <- "This table provides size details about R objects."

    shiny::callModule(
      tableModule, "objects",
      func = object.RENDER,
      info.text = object_text,
      options = NULL,
      title = "R object sizes"
    )
    
  }) ## end of moduleServer
} ## end of server
