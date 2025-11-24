lasagna_multipartite_nodes_table_ui <- function(id,
                                               label = "",
                                               title = "",
                                               info.text = "",
                                               caption = "",
                                               height = 400,
                                               width = 400) {

  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::radioButtons(ns("labeltype"), "Label type:",
      c("feature", "symbol", "title"), selected = "feature", inline = TRUE)
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

lasagna_multipartite_nodes_table_server <- function(id,
                                                    data,
                                                    pgx,
                                                    scrollY = "auto") {

  moduleServer(id, function(input, output, session) {

    table_data <- shiny::reactive({      
      
      res <- data()
      shiny::req(res)
      G <- visNetwork::toVisNetworkData(res$graph)
      
      shiny::validate(shiny::need("nodes" %in% names(G),
        "Missing nodes from LASAGNA multipartite data!"))
      
      N <- G$nodes
      phenos <- NULL
      hh <- grep("PHENO:", rownames(N))
      if (length(hh) > 0) {
        phenos <- N[N$layer == "PHENO", , drop = FALSE]
      }

      layers <- unique(as.character(N$layer))
      layers <- layers[which(layers != "PHENO")]
      N <- do.call(rbind,
        lapply(layers, function(l) { N0=N[N$layer==l, , drop = FALSE];
          N0=N0[order(N0$fc, decreasing = TRUE), , drop = FALSE]
        })
      )
      N <- rbind(N, phenos)

      N$rho <- round(N$rho, 2)
      N$fc <- round(N$fc, 2)
      colnames(N)[which(colnames(N) == "fc")] <- "log2FC"

      if (all(c("id", "label") %in% colnames(N))) {
        kk <- c("label", setdiff(colnames(N), c("value","label")))
        N <- N[, kk, drop = FALSE]
        if (isTRUE(all.equal(N$id, N$label)))         
          N <- N[, colnames(N) != "id", drop = FALSE]
      }

      rm(res, G, layers, LL); gc()

      return(N)

    })
    
    table.RENDER <- function(full = TRUE) {

      dt <- table_data()
      shiny::req(dt)

      ff0 <- rownames(dt)
      phenos <- NULL
      hh <- grep("PHENO:", ff0)
      if (length(hh) > 0) {
        phenos <- ff0[hh]
        ff0 <- ff0[-hh]
      }

      if (input$labeltype == "feature") ff <- c(ff0, phenos)
      jj <- match(ff0, pgx$genes$feature)
      if (input$labeltype == "symbol") {
        ff <- c(pgx$genes$symbol[jj], phenos)
      } else if (input$labeltype == "title") {
        ff <- c(pgx$genes$gene_title[jj], phenos)
      }      
      nas <- which(is.na(ff))
      if (any(nas)) ff[nas] <- rownames(dt)[nas]
      dt$label <- ff
      
      DTable <- DT::datatable(
        dt,
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
