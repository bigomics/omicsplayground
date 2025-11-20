lasagna_multipartite_graph_table_ui <- function(id,
                                                label = "",
                                                title = "",
                                                info.text = "",
                                                caption = "",
                                                height,
                                                width) {

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

lasagna_multipartite_graph_table_server <- function(id,
                                                    data,
                                                    scrollY = "auto") {

  moduleServer(id, function(input, output, session) {

    table_data <- shiny::reactive({
      res <- data()
      shiny::req(res)
      is.igraph <- ("graph" %in% names(res)) & (any(class(res[["graph"]]) == "igraph"))
      ## igraph::graph_attr(res$graph)
      G <- NULL
      if (is.igraph) G <- igraph::as_data_frame(res[["graph"]])
      shiny::validate(
        shiny::need(!is.null(G), "Missing igraph data from LASAGNA multipartite")
      )
      return(G)
    })
    
    table.RENDER <- function(full = TRUE) {

      G <- table_data()
      shiny::req(G)
      num_cols <- sapply(G, function(x) all(grepl("^\\s*-?\\d*\\.?\\d+\\s*$", x)))      
      if (any(num_cols)) {
        num_cols <- names(which(num_cols))
        kk <- match(num_cols, colnames(G))
        for(i in 1:length(kk)) G[,kk[i]] <- round(as.numeric(G[,kk[i]]), digits=2)
      }
      
      DTable <- DT::datatable(
        G,
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
