##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


functional_table_go_table_ui <- function(
  id,
  title,
  info.text,
  caption,
  width,
  height
) {
  ns <- shiny::NS(id)

  TableModuleUI(
    ns("datasets"),
    info.text = info.text,
    width = width,
    height = height,
    title = title,
    caption = caption,
    label = "b"
  )
}

functional_table_go_table_server <- function(id,
                                             pgx,
                                             fa_contrast,
                                             fa_filtertable,
                                             fa_filtertable_value,
                                             selected_gsetmethods) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    matchGOid2gset <- function(id, gsets) {
      gsets.id <- sub("\\)$", "", sub(".*\\(GO_", "GO:", gsets))
      match(id, gsets.id)
    }

    table_data <- shiny::reactive({
      comparison <- fa_contrast()
      if (is.null(pgx$meta.go)) {
        return(NULL)
      }
      if (is.null(comparison)) {
        return(NULL)
      }

      go <- pgx$meta.go$graph
      scores <- pgx$meta.go$pathscore[, comparison]
      scores <- scores[which(!is.na(scores) & !is.infinite(scores))]
      scores <- round(scores, digits = 3)
      scores <- scores[order(-abs(scores))]
      go.term <- igraph::V(go)[names(scores)]$Term

      ## get FC and q-value.  match with enrichment table
      gs.meta <- pgx$gset.meta$meta[[comparison]]
      ii <- matchGOid2gset(names(scores), rownames(gs.meta))
      gs.meta <- gs.meta[ii, , drop = FALSE]
      gs.meta$GO.id <- rownames(scores)
      mm <- selected_gsetmethods()
      mm <- intersect(mm, colnames(gs.meta$q))
      qv <- apply(gs.meta$q[, mm, drop = FALSE], 1, max, na.rm = TRUE) ## meta-q
      fx <- gs.meta$meta.fx
      go.term1 <- substring(go.term, 1, 80)
      dt1 <- round(cbind(score = scores, logFC = fx, meta.q = qv), digits = 4)
      dt <- data.frame(id = names(scores), term = go.term1, dt1, stringsAsFactors = FALSE)
      return(dt)
    })

    table_RENDER <- function() {
      dt <- table_data()
      shiny::req(dt)
      filtertable <- fa_filtertable()
      if (filtertable) {
        filter_value <- as.numeric(fa_filtertable_value())
        dt <- dt[which(dt$meta.q < filter_value), ]
      }

      id2 <- paste0("abc(", sub(":", "_", dt$id), ")") ## to match with wrapHyperLink
      id_link <- playbase::wrapHyperLink(
        rep_len(
          "<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>",
          nrow(dt)
        ),
        id2
      ) |> HandleNoLinkFound(
        NoLinkString = "<i class='fa-solid fa-arrow-up-right-from-square weblink'></i>",
        SubstituteString = "<i class='fa-solid fa-arrow-up-right-from-square blank_icon'></i>"
      )
      dt$term <- paste(dt$term, id_link)

      numeric.cols <- colnames(dt)[which(sapply(dt, is.numeric))]

      DT::datatable(dt,
        rownames = FALSE,
        escape = c(-1, -2),
        extensions = c("Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        fillContainer = TRUE,
        plugins = "scrollResize",
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = 200,
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          columnDefs = list(list(
            targets = 1, ## with no rownames column 1 is column 2
            render = DT::JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 50 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
              "}"
            )
          ))
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = color_from_middle(
            dt[, "score"],
            "lightblue",
            "#f5aeae"
          ),
          backgroundSize = "98% 88%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }

    table_RENDER_modal <- shiny::reactive({
      dt <- table_RENDER()
      dt$x$options$scrollY <- SCROLLY_MODAL
      dt
    })

    TableModuleServer(
      "datasets",
      func = table_RENDER,
      func2 = table_RENDER_modal,
      csvFunc = function() {
        table_RENDER()$x$data[, -1]
      },
      selector = "none"
    )
  }) ## end of moduleServer
} ## end of server
