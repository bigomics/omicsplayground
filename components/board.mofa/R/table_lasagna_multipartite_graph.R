lasagna_multipartite_edges_table_ui <- function(id,
                                                label = "",
                                                title = "",
                                                info.text = "",
                                                caption = "",
                                                height = 400,
                                                width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    withTooltip(
      shiny::checkboxInput(
        ns("show_intracor"),
        "Show intra-layer correlations", FALSE
      ),
      "Show intra-layer edges"
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

lasagna_multipartite_edges_table_server <- function(id, data, scrollY = "auto") {
  moduleServer(id, function(input, output, session) {
    table_data <- shiny::reactive({
      res <- data()
      shiny::req(res)
      G <- visNetwork::toVisNetworkData(res$graph)
      shiny::validate(shiny::need(
        "edges" %in% names(G),
        "Missing edges from LASAGNA multipartite data!"
      ))
      E <- G$edges
      E$rho <- round(E$rho, digits = 2)
      E$weight <- round(E$weight, digits = 2)
      rm(res)
      gc()
      return(E)
    })

    table.RENDER <- function(full = TRUE) {
      dt <- table_data()
      shiny::req(dt)

      jj <- grep("->", as.character(dt$connection_type))
      if (length(jj) > 0) dt1 <- dt[jj, , drop = FALSE]
      if (input$show_intracor & length(jj) > 0) {
        dt1 <- dt
      }

      DTable <- DT::datatable(
        dt1,
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

    table.RENDER2 <- shiny::reactive({
      return(table.RENDER())
    })

    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER2,
      selector = "none"
    )
  })
}
