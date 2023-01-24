##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FunctionalBoard <- function(id, inputData, selected_gsetmethods) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 750
    rowH <- 660 ## row height of panel
    tabH <- "70vh" ## row height of panel

    fa_infotext <- paste("This module performs specialized pathway analysis.
    <br><br>", a_KEGG, " is a collection of manually curated pathways
    representing the current knowledge of molecular interactions, reactions and
    relation networks as pathway maps. In the <strong>KEGG pathway</strong>
    panel, each pathway is scored for the selected contrast profile and reported
    in the table. A unique feature of the platform is that it provides an
    activation-heatmap comparing the activation levels of pathways across
    multiple contrast profiles. This facilitates to quickly see and detect the
    similarities between profiles in certain pathways.
    <br><br>In the <strong>GO</strong> panel, users can perform ", a_GO, " (GO)
    analysis. GO defines functional concepts/classes and their relationships as
    a hierarchical graph. The GO database provides a computational representation
    of the current knowledge about roles of genes for many organisms in terms of
    molecular functions, cellular components and biological processes. All the
    features described under the KEGG pathway tab, such as scoring the gene sets
    and drawing an activation-heatmap, can be performed for the GO database under
    the GO graph tab. Instead of pathway maps, an annotated graph structure
    provided by the GO database is potted for every selected gene set.
    <br><br><br><br>
    <center><iframe width='500' height='333'
    src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6'
    frameborder='0' allow='accelerometer; autoplay; encrypted-media;
    gyroscope; picture-in-picture' allowfullscreen></iframe></center>")

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$fa_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Functional Analysis Board</strong>"),
        shiny::HTML(fa_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({
      ngs <- inputData()
      shiny::req(ngs)
      ct <- colnames(ngs$model.parameters$contr.matrix)
      ct <- sort(ct)
      shiny::updateSelectInput(session, "fa_contrast", choices = ct)
    })

    ## =========================================================================
    ## KEGG pathways
    ## =========================================================================
    getKeggTable <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs, input$fa_contrast)

      ## ----- get comparison
      comparison <- input$fa_contrast
      if (!(comparison %in% names(ngs$gset.meta$meta))) {
        return(NULL)
      }

      ## ----- get KEGG id
      xml.dir <- file.path(FILES, "kegg-xml")
      kegg.available <- gsub("hsa|.xml", "", dir(xml.dir, pattern = "*.xml"))
      kegg.ids <- getKeggID(rownames(ngs$gsetX))
      ## sometimes no KEGG in genesets...
      if (length(kegg.ids) == 0) {
        shinyWidgets::sendSweetAlert(
          session = session,
          title = "No KEGG terms in enrichment results",
          text = "",
          type = "warning"
        )
        df <- data.frame()
        return(df)
      }

      jj <- which(!is.na(kegg.ids) &
        !duplicated(kegg.ids) &
        kegg.ids %in% kegg.available)
      kegg.gsets <- rownames(ngs$gsetX)[jj]
      kegg.ids <- kegg.ids[jj]

      meta <- ngs$gset.meta$meta[[comparison]]
      meta <- meta[kegg.gsets, ]
      mm <- selected_gsetmethods()
      mm <- intersect(mm, colnames(meta$q))
      meta.q <- apply(meta$q[, mm, drop = FALSE], 1, max, na.rm = TRUE)

      df <- data.frame(
        kegg.id = kegg.ids, pathway = kegg.gsets,
        logFC = meta$meta.fx, meta.q = meta.q,
        check.names = FALSE
      )
      df <- df[!duplicated(df$kegg.id), ] ## take out duplicated gene sets...
      df <- df[order(-abs(df$logFC)), ]
      return(df)
    })

    getFilteredKeggTable <- shiny::reactive({
      df <- getKeggTable()
      do.filter <- FALSE
      do.filter <- input$fa_filtertable
      if (do.filter) df <- df[which(df$meta.q < 0.999), ]
      return(df)
    })

    ## There is a bug in pathview::geneannot.map so we have to override
    ## "Error in pathview::mol.sum(gene.data, gene.idmap) : no ID can be mapped!"
    my.geneannot.map <- function(in.ids, in.type, out.type, org = "Hs", pkg.name = NULL,
                                 unique.map = TRUE, na.rm = TRUE, keep.order = TRUE) {
      if (is.null(pkg.name)) {
        data(bods)
        ridx <- grep(tolower(paste0(org, "[.]")), tolower(bods[, 1]))
        if (length(ridx) == 0) {
          ridx <- grep(tolower(org), tolower(bods[, 2:3])) %% nrow(bods)
          if (length(ridx) == 0) {
            stop("Wrong org value!")
          }
          if (any(ridx == 0)) {
            ridx[ridx == 0] <- nrow(bods)
          }
        }
        pkg.name <- bods[ridx, 1]
      }
      pkg.on <- try(requireNamespace(pkg.name), silent = TRUE)
      if (!pkg.on) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg.name, suppressUpdates = TRUE)
        pkg.on <- try(requireNamespace(pkg.name), silent = TRUE)
        if (!pkg.on) {
          stop(paste("Fail to install/load gene annotation package ",
            pkg.name, "!",
            sep = ""
          ))
        }
      }
      db.obj <- eval(parse(text = paste0(pkg.name, "::", pkg.name)))
      id.types <- AnnotationDbi::columns(db.obj)
      in.type <- toupper(in.type)
      out.type <- toupper(out.type)
      eii <- in.type == toupper("entrez") | in.type == toupper("eg")
      if (any(eii)) {
        in.type[eii] <- "ENTREZID"
      }
      eio <- out.type == toupper("entrez") | out.type == toupper("eg")
      if (any(eio)) {
        out.type[eio] <- "ENTREZID"
      }
      if (in.type == out.type) {
        stop("in.type and out.type are the same, no need to map!")
      }
      nin <- length(in.type)
      if (nin != 1) {
        stop("in.type must be of length 1!")
      }
      out.type <- out.type[!out.type %in% in.type]
      nout <- length(out.type)
      msg <- paste0(
        "must from: ", paste(id.types, collapse = ", "),
        "!"
      )
      if (!in.type %in% id.types) {
        stop("'in.type' ", msg)
      }
      if (!all(out.type %in% id.types)) {
        stop("'out.type' ", msg)
      }
      in.ids0 <- in.ids
      in.ids <- unique(as.character(in.ids))
      out.ids <- character(length(in.ids))

      res <- try(suppressWarnings(
        AnnotationDbi::select(db.obj,
          keys = in.ids,
          keytype = in.type,
          columns = c(in.type, out.type)
        )
      ))
      if (class(res) == "data.frame") {
        res <- res[, c(in.type, out.type)]
        if (nout == 1) {
          na.idx <- is.na(res[, 2])
        } else {
          na.idx <- apply(res[, -1], 1, function(x) all(is.na(x)))
        }
        if (sum(na.idx) > 0) {
          n.na <- length(unique(res[na.idx, 1]))
          if (na.rm) {
            res <- res[!na.idx, ]
          }
        }
        cns <- colnames(res)
        if (unique.map) {
          if (length(out.type) == 1) {
            umaps <- tapply(res[, out.type], res[, in.type],
              paste,
              sep = "", collapse = "; "
            )
          } else {
            umaps <- apply(res[, out.type], 2, function(x) {
              tapply(x, res[, in.type], function(y) {
                paste(unique(y),
                  sep = "", collapse = "; "
                )
              })
            })
          }
          umaps <- cbind(umaps)
          res.uniq <- cbind(rownames(umaps), umaps)
          res <- res.uniq
          colnames(res) <- cns
        }
        res <- as.matrix(res)
        if (!keep.order) {
          rownames(res) <- NULL
          return(res)
        } else {
          res1 <- matrix(NA, ncol = length(cns), nrow = length(in.ids0))
          res1[, 1] <- in.ids0
          rns <- match(in.ids0, res[, 1])
          res1[, -1] <- res[rns, -1]
          colnames(res1) <- cns
          return(res1)
        }
      } else {
        res <- cbind(in.ids, out.ids)
        colnames(res) <- c(in.type, out.type)
        return(res)
      }
    }

    # random global server actions ..
    suppressMessages(require(pathview))
    unlockBinding("geneannot.map", as.environment("package:pathview"))
    assignInNamespace("geneannot.map", my.geneannot.map, ns = "pathview", as.environment("package:pathview"))
    assign("geneannot.map", my.geneannot.map, as.environment("package:pathview"))
    lockBinding("geneannot.map", as.environment("package:pathview"))

    ##############
    ## KEGG TAB ##
    ##############
    functional_plot_kegg_graph_server("kegg_graph",
                                      inputData,
                                      getFilteredKeggTable,
                                      kegg_table,
                                      reactive(input$fa_contrast))

    functional_plot_kegg_actmap_server("kegg_actmap",
                                       inputData,
                                       getKeggTable)

    kegg_table <- functional_table_kegg_table_server("kegg_table",
                                                     inputData,
                                                     getFilteredKeggTable,
                                                     reactive(input$fa_contrast),
                                                     tabH)

    ## ================================================================================
    ## GO graph
    ## ================================================================================



    #GO_network.RENDER <- shiny::reactive({
    #  ngs <- inputData()
    #  shiny::req(ngs)
    #  require(igraph)
#
    #  comparison <- 1
    #  methods <- c("fisher", "gsva", "camera")
    #  comparison <- input$fa_contrast
    #  shiny::req(input$fa_contrast)
    #  if (is.null(comparison)) {
    #    return(NULL)
    #  }
#
    #  sub2 <- go <- ngs$meta.go$graph
    #  if (is.null(go)) {
    #    shinyWidgets::sendSweetAlert(
    #      session = session,
    #      title = "No GO graph in enrichment results",
    #      text = "",
    #      type = "warning"
    #    )
    #    return(NULL)
    #  }
#
    #  score <- ngs$meta.go$pathscore[, comparison]
    #  score[is.na(score) | is.infinite(score)] <- 0
    #  score <- (score / (1e-8 + max(abs(score), na.rm = TRUE)))
    #  igraph::V(sub2)$value <- score
    #  igraph::V(sub2)$color <- BLUERED(32)[16 + round(15 * score)]
    #  igraph::V(sub2)$label <- igraph::V(sub2)$Term
    #  igraph::V(sub2)$label[which(is.na(score) | score == 0)] <- ""
    #  pos <- sub2$layout
#
    #  all.zero <- all(score == 0)
#
    #  if (!all.zero && input$GO_prunetree) {
    #    vv <- igraph::V(sub2)[which(!is.na(score) & abs(score) > 0)]
    #    sp <- igraph::shortest_paths(sub2, from = "all", to = vv, mode = "all", output = "vpath")
    #    sp.vv <- unique(unlist(sp$vpath))
    #    sub2 <- igraph::induced.subgraph(sub2, sp.vv)
    #    pos <- igraph::layout_with_fr(sub2)
    #    score <- score[igraph::V(sub2)$name]
    #  }
#
    #  ## remove root?
    #  removeroot <- TRUE
    #  if (removeroot) {
    #    sub2 <- igraph::induced_subgraph(sub2, which(igraph::V(sub2)$name != "all"))
    #    if (input$GO_prunetree) pos <- igraph::layout_with_fr(sub2)
    #    score <- score[igraph::V(sub2)$name]
    #  }
    #  roots <- c("all", neighbors(go, igraph::V(go)["all"], mode = "all")$name)
    #  roots <- intersect(roots, igraph::V(sub2)$name)
#
    #  astree <- TRUE
    #  if (astree) {
    #    if ("all" %in% igraph::V(sub2)$name) {
    #      pos <- igraph::layout_as_tree(sub2, root = "all", mode = "all")
    #    } else {
    #      pos <- igraph::layout_as_tree(sub2, root = roots, mode = "all")
    #    }
    #    pos[, 2] <- -pos[, 2]
    #  }
#
    #  ## color clusters
    #  if (input$GO_colorclusters) {
    #    clust <- igraph::cluster_louvain(igraph::as.undirected(go))$membership
    #    names(clust) <- igraph::V(go)$name
    #    cc <- c(
    #      RColorBrewer::brewer.pal(12, "Set3"),
    #      RColorBrewer::brewer.pal(8, "Set2"),
    #      RColorBrewer::brewer.pal(8, "Set1")
    #    )
    #    igraph::V(sub2)$color <- rep(cc, 99)[clust[igraph::V(sub2)$name]]
    #    jj <- which(is.na(score) | score == 0)
    #    if (length(jj) > 0) igraph::V(sub2)$color[jj] <- NA
    #  }
#
    #  gr <- visNetwork::toVisNetworkData(sub2)
    #  gr$nodes$color[is.na(gr$nodes$color)] <- "#F9F9F9"
    #  gr$nodes$value <- pmax(abs(gr$nodes$value), 0.001)
    #  gr$nodes$x <- pos[, 1] * 60
    #  gr$nodes$y <- pos[, 2] * 90
    #  gr$nodes$label <- gr$nodes$Term
    #  no.score <- (is.na(score) | score == 0)
    #  gr$nodes$label[which(no.score)] <- NA
#
    #  gr$nodes$shape <- c("box", "circle")[1 + 1 * no.score]
    #  gr$nodes$label <- sapply(gr$nodes$label, breakstring, n = 25, nmax = 95, force = TRUE, brk = "\n")
#
    #  gr.def <- sapply(gr$nodes$Definition, breakstring, n = 50, brk = "<br>")
    #  gr$nodes$title <- paste0(
    #    gr$nodes$Term, "  (", gr$nodes$id, ")<br>",
    #    "<small>", gr.def, "</small>"
    #  )
#
    #  ## rendering
    #  font.size <- 20
    #  cex <- 1
    #  if (input$GO_prunetree) {
    #    font.size <- 20
    #    cex <- 0.6
    #  }
#
    #  visNetwork::visNetwork(gr$nodes, gr$edges) %>%
    #    visNetwork::visEdges(
    #      smooth = FALSE, hidden = FALSE, arrows = list(enabled = TRUE),
    #      scaling = list(min = 10 * cex, max = 30 * cex), width = 5 * cex
    #    ) %>%
    #    visNetwork::visNodes(
    #      font = list(size = font.size * cex, vadjust = 0),
    #      scaling = list(min = 1 * cex, max = 80 * cex)
    #    ) %>%
    #    visNetwork::visPhysics(stabilization = FALSE) %>%
    #    visNetwork::visOptions(highlightNearest = list(enabled = T, degree = 1, hover = TRUE)) %>%
    #    visNetwork::visPhysics(enabled = FALSE)
    #})

    matchGOid2gset <- function(id, gsets) {
      gsets.id <- sub("\\)$", "", sub(".*\\(GO_", "GO:", gsets))
      match(id, gsets.id)
    }

    GO_table.RENDER <- shiny::reactive({
      ngs <- inputData()
      shiny::req(ngs, input$fa_contrast)
      if (is.null(ngs$meta.go)) {
        return(NULL)
      }

      comparison <- input$fa_contrast
      if (is.null(comparison)) {
        return(NULL)
      }

      go <- ngs$meta.go$graph
      scores <- ngs$meta.go$pathscore[, comparison]

      scores <- scores[which(!is.na(scores) & !is.infinite(scores))]
      scores <- round(scores, digits = 3)
      scores <- scores[order(-abs(scores))]
      go.term <- igraph::V(go)[names(scores)]$Term

      ## get FC and q-value.  match with enrichment table
      gs.meta <- ngs$gset.meta$meta[[comparison]]
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
          scrollY = tabH, scroller = TRUE, deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 4) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%") %>%
        DT::formatStyle("score",
          background = color_from_middle(dt1[, "score"], "lightblue", "#f5aeae"),
          backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    })

    #plotGOactmap <- function(score, go, normalize, maxterm, maxfc) {
    #  rownames(score) <- igraph::V(go)[rownames(score)]$Term
#
    #  ## avoid errors!!!
    #  score[is.na(score) | is.infinite(score)] <- 0
    #  score[is.na(score)] <- 0
#
    #  ## reduce score matrix
    #  score <- score[head(order(-rowSums(score**2, na.rm = TRUE)), maxterm), , drop = FALSE] ## max number terms
    #  score <- score[, head(order(-colSums(score**2, na.rm = TRUE)), maxfc), drop = FALSE] ## max comparisons/FC
    #  score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))
#
    #  ## normalize colums
    #  if (normalize) {
    #    ## column scale???
    #    score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2, na.rm = TRUE))))
    #  }
    #  score <- score / max(abs(score), na.rm = TRUE) ## global normalize
    #  score <- sign(score) * abs(score)**0.5 ## fudging for better colors
#
    #  d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
    #  d2 <- as.dist(1 - cor(score, use = "pairwise"))
    #  d1 <- dist(score)
    #  d2 <- dist(t(score))
    #  d1[is.na(d1)] <- 1
    #  d2[is.na(d2)] <- 1
    #  ii <- 1:nrow(score)
    #  jj <- 1:ncol(score)
    #  if (NCOL(score) == 1) {
    #    score <- score[order(-score[, 1]), 1, drop = FALSE]
    #  } else {
    #    ii <- hclust(d1)$order
    #    jj <- hclust(d2)$order
    #    score <- score[ii, jj, drop = FALSE]
    #  }
#
    #  colnames(score) <- substring(colnames(score), 1, 30)
    #  rownames(score) <- substring(rownames(score), 1, 50)
    #  colnames(score) <- paste0(colnames(score), " ")
#
    #  bmar <- 0 + pmax((50 - nrow(score)) * 0.25, 0)
    #  par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), oma = c(0, 1.5, 0, 0.5))
#
    #  corrplot::corrplot(score,
    #    is.corr = FALSE, cl.pos = "n", col = BLUERED(100),
    #    tl.cex = 0.85, tl.col = "grey20", tl.srt = 90,
    #    mar = c(bmar, 0, 0, 0)
    #  )
    #}

    #GO_actmap.RENDER <- shiny::reactive({
    #  ngs <- inputData()
    #  shiny::req(ngs)

    #  if (is.null(ngs$meta.go)) {
    #    return(NULL)
    #  }

    #  score <- ngs$meta.go$pathscore
    #  go <- ngs$meta.go$graph

    #  plotGOactmap(
    #    score = score, go = go,
    #    normalize = input$go_normalize,
    #    maxterm = 50,
    #    maxfc = 25
    #  )
    #})

    #GO_actmap.RENDER2 <- shiny::reactive({
    #  ngs <- inputData()
    #  shiny::req(ngs)
#
    #  if (is.null(ngs$meta.go)) {
    #    return(NULL)
    #  }
#
    #  score <- ngs$meta.go$pathscore
    #  go <- ngs$meta.go$graph
#
    #  plotGOactmap(
    #    score = score, go = go,
    #    normalize = input$go_normalize,
    #    maxterm = 50,
    #    maxfc = 100
    #  )
    #})

    #GO_info1 <- "The <strong>Gene Ontology</strong> (GO) provides a computational representation of the current knowledge about roles of genes for many organisms in terms of molecular functions, cellular components and biological processes. The structure of GO can be described in terms of a graph, where each GO term is a node, and the relationships between the terms are edges between the nodes. GO is loosely hierarchical, with ‘child’ terms being more specialized than their ‘parent’ terms. The graph is interactive. You can move the graph and zoom in using the mouse."

    #GO_network.opts <- shiny::tagList(
    #  withTooltip(
    #    shiny::checkboxInput(ns("GO_prunetree"), "Prune tree", TRUE),
    #    "Prune the tree with only significant branches."
    #  ),
    #  withTooltip(
    #    shiny::checkboxInput(ns("GO_colorclusters"), "Color clusters", FALSE),
    #    "Highlight clusters with different colors."
    #  )
    #)

    #shiny::callModule(
    #  plotModule,
    #  id = "GO_network",
    #  title = "Gene Ontology graph", label = "a",
    #  func = GO_network.RENDER,
    #  plotlib = "visnetwork",
    #  info.text = GO_info1,
    #  download.fmt = c("pdf", "png"), ## no.download=TRUE,
    #  options = GO_network.opts,
    #  pdf.width = 10, pdf.height = 8,
    #  height = c(0.55 * rowH, 750), width = c("100%", 1400),
    #  res = 72,
    #  add.watermark = WATERMARK
    #)

    functional_plot_go_network_server("GO_network",
                                      inputData,
                                      reactive(input$fa_contrast))

    #GO_actmap.opts <- shiny::tagList(
    #  withTooltip(shiny::checkboxInput(ns("go_normalize"), "normalize activation matrix", FALSE), "Click to normalize the columns of the activation matrices.")
    #)
    #go_info <- "The <b>GO activation matrix</b> visualizes the activation of GO terms across conditions. From this figure, you can easily detect GO terms that are consistently up/down across conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile."

    #shiny::callModule(
    #  plotModule,
    #  id = "GO_actmap",
    #  func = GO_actmap.RENDER,
    #  func2 = GO_actmap.RENDER2,
    #  title = "Activation matrix", label = "c",
    #  info.text = go_info,
    #  options = GO_actmap.opts,
    #  pdf.height = 9, pdf.width = 9,
    #  height = c(rowH, 750), width = c("100%", 1400),
    #  res = 72,
    #  add.watermark = WATERMARK
    #)

    functional_plot_go_actmap_server("GO_actmap",
                                     inputData)

    GO_table <- shiny::callModule(
      tableModule,
      id = "GO_table", label = "b",
      func = GO_table.RENDER,
      info.text = "<strong>GO score table.</strong> The scoring of a GO term is performed by considering the cumulative score of all terms from that term to the root node. That means that GO terms that are supported by higher level terms levels are preferentially scored.",
      title = "GO score table",
      height = c(270, 700)
    )

  }) ## end-of-moduleServer
}
