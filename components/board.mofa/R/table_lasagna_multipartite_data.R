lasagna_multipartite_data_table_ui <- function(id,
                                               label = "",
                                               title = "",
                                               info.text = "",
                                               caption = "",
                                               height = 400,
                                               width = 400) {

  ns <- shiny::NS(id)

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    width = width,    
    height = height,
    title = title,
    caption = caption,
    label = label
  )

}

lasagna_multipartite_data_table_server <- function(id, data, scrollY = "auto") {

  moduleServer(id, function(input, output, session) {

    table_data <- shiny::reactive({
      res <- data()
      shiny::req(res)
      X <- as.data.frame(res[["X"]])
      Y <- as.data.frame(res[["Y"]])
      kk <- intersect(colnames(X), rownames(Y))
      X <- round(X[, kk, drop = FALSE], digits = 2)      
      hh <- grep("PHENO:|SINK|SOURCE", rownames(X))
      if (any(hh)) X <- X[-hh, , drop = FALSE]
      cc <- (!is.null(X) & !is.null(Y) & length(kk) == ncol(X))
      shiny::validate(
        shiny::need(cc, "Missing data from LASAGNA multipartite graph!")
      )
      rm(Y)
      return(X)
    })
    
    table.RENDER <- function(full = TRUE) {

      dt <- table_data()
      shiny::req(dt)

      DTable <- DT::datatable(
        dt,
        rownames = TRUE,
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
