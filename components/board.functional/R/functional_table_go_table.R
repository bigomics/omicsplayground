##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


functional_table_go_table_ui <- function(id, width, height) {
  ns <- shiny::NS(id)

  info_text <- strwrap("<strong>GO score table.</strong> The scoring of a GO
                         term is performed by considering the cumulative score
                         of all terms from that term to the root node. That
                         means that GO terms that are supported by higher level
                         terms levels are preferentially scored.")

  TableModuleUI(
    ns("datasets"),
    info.text = info_text,
    width = width,
    height = height,
    title = "GO score table",
    label = "b"
  )

}


functional_table_go_table_server <- function(id,
                                             inputData,
                                             fa_contrast,
                                             tabH,
                                             selected_gsetmethods)
{
  moduleServer(id, function(input, output, session)
  {
    ns <- session$ns

    matchGOid2gset <- function(id, gsets) {
      gsets.id <- sub("\\)$", "", sub(".*\\(GO_", "GO:", gsets))
      match(id, gsets.id)
    }

    table_data <- shiny::reactive({
      res <- list(
        pgx = inputData(),
        fa_contrast = fa_contrast()
      )
      return(res)
    })

    table_RENDER <- function() {
      res <- table_data()
      pgx <- res$pgx
      comparison <- res$fa_contrast

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
      id2 <- paste0("abc(", sub(":", "_", dt$id), ")") ## to match with wrapHyperLink
      dt$id <- wrapHyperLink(as.character(dt$id), id2) ## add link

      numeric.cols <- colnames(dt)[which(sapply(dt, is.numeric))]

      DT::datatable(dt,
                    rownames = FALSE, escape = c(-1, -2),
                    class = "compact cell-border stripe hover",
                    extensions = c("Scroller"),
                    selection = list(mode = "single", target = "row", selected = 1),
                    fillContainer = TRUE,
                    options = list(
                      dom = "lfrtip",
                      scrollX = TRUE,
                      scrollY = "15vh", scroller = TRUE, deferRender = TRUE
                    ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
                        background = color_from_middle(dt1[, "score"],
                                                       "lightblue",
                                                       "#f5aeae"),
                                                       backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
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
      selector = "none"
    )

  })  ## end of moduleServer
} ## end of server
