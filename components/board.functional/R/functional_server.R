##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

FunctionalBoard <- function(id, pgx, selected_gsetmethods) {
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
      shiny::req(pgx)
      ct <- colnames(pgx$model.parameters$contr.matrix)
      ct <- sort(ct)
      shiny::updateSelectInput(session, "fa_contrast", choices = ct)
    })

    ## =========================================================================
    ## KEGG pathways
    ## =========================================================================
    getKeggTable <- shiny::reactive({
      shiny::req(pgx, input$fa_contrast)

      ## ----- get comparison
      comparison <- input$fa_contrast
      if (!(comparison %in% names(pgx$gset.meta$meta))) {
        return(NULL)
      }

      ## ----- get KEGG id
      xml.dir <- file.path(FILES, "kegg-xml")
      kegg.available <- gsub("hsa|.xml", "", dir(xml.dir, pattern = "*.xml"))
      kegg.ids <- playbase::getKeggID(rownames(pgx$gsetX))
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
      kegg.gsets <- rownames(pgx$gsetX)[jj]
      kegg.ids <- kegg.ids[jj]

      meta <- pgx$gset.meta$meta[[comparison]]
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
    functional_plot_kegg_graph_server(
      "kegg_graph",
      pgx,
      getFilteredKeggTable,
      kegg_table,
      reactive(input$fa_contrast)
    )

    functional_plot_kegg_actmap_server(
      "kegg_actmap",
      pgx,
      getKeggTable
    )

    kegg_table <- functional_table_kegg_table_server(
      "kegg_table",
      pgx,
      getFilteredKeggTable,
      reactive(input$fa_contrast),
      tabH
    )

    ## ================================================================================
    ## GO Tab
    ## ================================================================================

    functional_plot_go_network_server(
      "GO_network",
      pgx,
      reactive(input$fa_contrast)
    )

    functional_plot_go_actmap_server(
      "GO_actmap",
      pgx
    )

    functional_table_go_table_server(
      "GO_table",
      pgx,
      reactive(input$fa_contrast),
      tabH,
      selected_gsetmethods
    )
  }) ## end-of-moduleServer
}
