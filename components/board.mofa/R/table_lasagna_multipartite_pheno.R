lasagna_multipartite_pheno_table_ui <- function(id,
                                                label = "",
                                                title = "",
                                                info.text = "",
                                                caption = "",
                                                height,
                                                width) {

  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(ns("show_all_ctx"), "Show all contrasts", FALSE),
      "Show all available contrasts"
    )
  )

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    width = width,    
    height = height,
    title = title,
    options = options,
    caption = caption,
    label = label
  )

}

lasagna_multipartite_pheno_table_server <- function(id,
                                                    data,
                                                    pgx,
                                                    input_contrast,
                                                    scrollY = "auto") {

  moduleServer(id, function(input, output, session) {

    table_data <- shiny::reactive({
      res <- data()
      ctx <- input_contrast()
      shiny::req(res, ctx)
      Y=NULL
      cc <- (!is.null(res[["Y"]]) & !is.null(res[["X"]]))
      if (cc) { 
        Y <- as.data.frame(res[["Y"]])
        kk <- intersect(colnames(res[["X"]]), Y)
        if (any(kk)) Y <- Y[kk, , drop = FALSE]
      }
      msg <- "Missing data from LASAGNA multipartite graph!"
      shiny::validate(shiny::need(!is.null(Y), msg))
      return(list(Y = Y, ctx = ctx))
    })
    
    table.RENDER <- function(full = TRUE) {

      dt <- table_data()
      Y <- dt[["Y"]]      
      ctx <- dt[["ctx"]]
      shiny::req(Y, ctx)

      Y <- cbind(Sample = rownames(Y), Y)
      kk <- intersect(c("Sample", ctx), colnames(Y))
      if (!input$show_all_ctx) Y <- Y[, kk, drop = FALSE]
      rownames(Y) <- 1:nrow(Y)

      DTable <- DT::datatable(
        Y,
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
          deferRender = TRUE
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
      
      return(DTable)

    }
    

    table.RENDER2 <- shiny::reactive({ return(table.RENDER()) })
    
    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER2,
      selector = "none"
    )

  })

}
