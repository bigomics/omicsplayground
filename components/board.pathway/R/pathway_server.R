##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PathwayBoard <- function(id, pgx, selected_gsetmethods = reactive(colnames(pgx$gset.meta$meta[[1]]$fc))) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 750
    rowH <- 660 ## row height of panel
    tabH <- "70vh" ## row height of panel

    fa_infotext <- paste("This module performs specialized pathway analysis.
    <br><br>Reactome and WikiPathways are collections of manually curated pathways
    representing the current knowledge of molecular interactions, reactions and
    relation networks as pathway maps. Each pathway is scored for the selected
    contrast profile and reported
    in the table. A unique feature of the platform is that it provides an
    activation-heatmap comparing the activation levels of pathways across
    multiple contrast profiles. This facilitates to quickly see and detect the
    similarities between profiles in certain pathways.
    <br><br>In the <strong>GO</strong> panel, users can perform ", a_GO, " (GO)
    analysis. GO defines functional concepts/classes and their relationships as
    a hierarchical graph. The GO database provides a computational representation
    of the current knowledge about roles of genes for many organisms in terms of
    molecular functions, cellular components and biological processes. All the
    features described under the Reactome pathway tab, such as scoring the gene sets
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
      shiny::req(pgx$X)
      ct <- colnames(pgx$model.parameters$contr.matrix)
      ct <- sort(ct)
      shiny::updateSelectInput(session, "fa_contrast", choices = ct)
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "WikiPathways" = list(disable = NULL),
      "Reactome" = list(disable = NULL),
      "GO graph" = list(disable = c("fa_filtertable", "fa_filtertable_value"))
    )
    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })

    ## ================================================================================
    ## =========================== FUNCTIONS ==========================================
    ## ================================================================================

    plotActivationMatrix <- function(meta, df, normalize = 1, nterms = 40,
                                     nfc = 10, tl.cex = 1.0, row.nchar = 60) {
      fx <- sapply(meta, function(x) x$meta.fx)
      qv <- sapply(meta, function(x) x$meta.q)
      rownames(fx) <- rownames(qv) <- rownames(meta[[1]])

      kk <- rownames(fx)
      kk <- as.character(df$pathway)

      if (length(kk) < 3) {
        return(NULL)
      }

      if (mean(is.na(qv)) < 0.01) {
        score <- fx[kk, , drop = FALSE] * (1 - qv[kk, , drop = FALSE])**2
      } else {
        score <- fx[kk, , drop = FALSE]
      }

      score <- score[head(order(-rowSums(score**2)), nterms), , drop = FALSE] ## nr gene sets
      score <- score[, head(order(-colSums(score**2)), nfc), drop = FALSE] ## max comparisons/FC
      score <- score + 1e-3 * matrix(rnorm(length(score)), nrow(score), ncol(score))
      d1 <- as.dist(1 - cor(t(score), use = "pairwise"))
      d2 <- as.dist(1 - cor(score, use = "pairwise"))
      d1[is.na(d1)] <- 1
      d2[is.na(d2)] <- 1
      ii <- 1:nrow(score)
      jj <- 1:ncol(score)
      if (NCOL(score) == 1) {
        score <- score[order(-score[, 1]), 1, drop = FALSE]
      } else {
        ii <- hclust(d1)$order
        jj <- hclust(d2)$order
        score <- score[ii, jj, drop = FALSE]
      }

      ## fudged score just for visualization
      score2 <- score
      if (normalize) score2 <- t(t(score2) / apply(abs(score2), 2, max))
      score2 <- sign(score2) * abs(score2 / max(abs(score2)))**1 ## fudging
      rownames(score2) <- tolower(gsub(".*:|wikipathway_|_Homo.*$", "",
        rownames(score2),
        ignore.case = TRUE
      ))
      rownames(score2) <- gsub("(_.*$)", "", rownames(score2))
      #
      rownames(score2) <- playbase::shortstring(rownames(score2), row.nchar)
      colnames(score2) <- playbase::shortstring(colnames(score2), 30)
      colnames(score2) <- paste0(colnames(score2), " ")

      bmar <- 0 + pmax(50 - nrow(score2), 0) * 0.3
      par(mfrow = c(1, 1), mar = c(1, 1, 10, 1), oma = c(0, 1.5, 0, 0.5))

      corrplot::corrplot(
        score2,
        is.corr = FALSE,
        cl.pos = "n",
        col = playdata::BLUERED(100),
        tl.cex = 1.0 * tl.cex,
        tl.col = "grey20",
        tl.srt = 90,
        mar = c(0, 0, 0.5, 0)
      )
    }


    ## =========================================================================
    ## KEGG pathways
    ## =========================================================================

    ## =========================================================================
    ## Get Reactome table
    ## =========================================================================

    getReactomeTable <- shiny::reactive({
      shiny::req(pgx$X, input$fa_contrast)

      ## ----- get comparison
      comparison <- input$fa_contrast
      if (!(comparison %in% names(pgx$gset.meta$meta))) {
        return(NULL)
      }

      ## ----- get REACTOME id
      sbgn.dir <- pgx.system.file("sbgn/", package = "pathway")
      reactome.available <- gsub("^.*reactome_|.sbgn$", "", dir(sbgn.dir, pattern = "*.sbgn"))
      reactome.gsets <- grep("R-HSA", rownames(pgx$gsetX), value = TRUE)
      reactome.ids <- gsub(".*R-HSA", "R-HSA", reactome.gsets)
      ## sometimes no REACTOME in genesets...
      if (length(reactome.ids) == 0) {
        shinyalert::shinyalert(
          title = "No REACTOME terms in enrichment results",
          text = "",
          type = "warning"
        )
        df <- data.frame()
        return(df)
      }

      ## select those of which we have SGBN files
      jj <- which(!is.na(reactome.ids) &
        !duplicated(reactome.ids) &
        reactome.ids %in% reactome.available)
      reactome.gsets <- reactome.gsets[jj]
      reactome.ids <- reactome.ids[jj]

      meta <- pgx$gset.meta$meta[[comparison]]
      meta <- meta[reactome.gsets, ]
      mm <- "fgsea"
      mm <- selected_gsetmethods()
      mm <- intersect(mm, colnames(meta$q))
      meta.q <- apply(meta$q[, mm, drop = FALSE], 1, max, na.rm = TRUE)

      df <- data.frame(
        reactome.id = reactome.ids,
        pathway = reactome.gsets,
        logFC = meta$meta.fx,
        meta.q = meta.q,
        check.names = FALSE
      )
      df <- df[!duplicated(df$reactome.id), ]
      df <- df[order(-abs(df$logFC)), ]
      return(df)
    })

    getFilteredReactomeTable <- shiny::reactive({
      df <- getReactomeTable()
      do.filter <- FALSE
      do.filter <- input$fa_filtertable
      if (do.filter) {
        filter_value <- as.numeric(input$fa_filtertable_value)
        df <- df[which(df$meta.q < filter_value), ]
      }
      return(df)
    })

    functional_plot_reactome_graph_server(
      "reactome_graph",
      pgx,
      getFilteredReactomeTable,
      reactome_table,
      reactive(input$fa_contrast),
      WATERMARK
    )

    functional_plot_reactome_actmap_server(
      "reactome_actmap",
      reactive(pgx$gset.meta$meta),
      getReactomeTable,
      plotActivationMatrix,
      WATERMARK
    )

    reactome_table <- functional_table_reactome_server(
      "reactome_table",
      getFilteredReactomeTable,
      fa_contrast = reactive(input$fa_contrast),
      scrollY = 180
    )

    functional_plot_enrichmap_server(
      "enrichment_map",
      pgx,
      reactive(input$fa_contrast),
      WATERMARK
    )

    ## ================================================================================
    ## GO module servers
    ## ================================================================================

    functional_plot_go_network_server(
      "GO_network",
      pgx,
      reactive(input$fa_contrast),
      WATERMARK
    )

    functional_plot_go_actmap_server(
      "GO_actmap",
      pgx,
      WATERMARK
    )

    functional_table_go_table_server(
      "GO_table",
      pgx = pgx,
      fa_contrast = reactive(input$fa_contrast),
      scrollY = 180,
      selected_gsetmethods = selected_gsetmethods
    )

    ## ================================================================================
    ## WikiPathway module servers
    ## ================================================================================

    getWikiPathwayTable <- shiny::reactive({
      shiny::req(pgx$X, input$fa_contrast)

      ## ----- get comparison
      comparison <- input$fa_contrast
      if (!(comparison %in% names(pgx$gset.meta$meta))) {
        return(NULL)
      }

      ## ----- get WIKIPATHWAY id
      wp.gsets <- grep("_WP", rownames(pgx$gsetX), value = TRUE)
      # extract wp.ids from string
      wp.ids <- gsub(".*_WP", "WP", wp.gsets)
      wp.ids <- gsub("(_.*$)", "", wp.ids)
      ## sometimes no WIKIPATHWAY in genesets...
      if (length(wp.ids) == 0) {
        shinyalert::shinyalert(
          title = "Alas...",
          text = "You have no WikiPathway terms in your enrichment results",
          type = "warning"
        )
        df <- data.frame()
        return(df)
      }

      ## select those with ID
      jj <- which(!is.na(wp.ids) & !duplicated(wp.ids))
      wp.gsets <- wp.gsets[jj]
      wp.ids <- wp.ids[jj]

      meta <- pgx$gset.meta$meta[[comparison]]
      meta <- meta[wp.gsets, ]
      mm <- "fgsea"
      mm <- selected_gsetmethods()
      mm <- intersect(mm, colnames(meta$q))
      meta.q <- apply(meta$q[, mm, drop = FALSE], 1, max, na.rm = TRUE)
      df <- data.frame(
        pathway.id = wp.ids,
        pathway = wp.gsets,
        logFC = meta$meta.fx,
        meta.q = meta.q,
        check.names = FALSE
      )
      df <- df[!duplicated(df$pathway.id), ]
      df <- df[order(-abs(df$logFC)), ]

      return(df)
    })

    getFilteredWikiPathwayTable <- shiny::reactive({
      df <- getWikiPathwayTable()
      do.filter <- FALSE
      do.filter <- input$fa_filtertable
      if (do.filter) {
        filter_value <- as.numeric(input$fa_filtertable_value)
        df <- df[which(df$meta.q < filter_value), ]
      }
      return(df)
    })

    functional_plot_wikipathway_graph_server(
      "wikipathway_graph",
      pgx,
      getFilteredWikiPathwayTable,
      wikipathway_table,
      reactive(input$fa_contrast),
      WATERMARK
    )

    functional_plot_wikipathway_actmap_server(
      "wikipathway_actmap",
      pgx,
      getWikiPathwayTable,
      plotActivationMatrix,
      WATERMARK
    )

    wikipathway_table <- functional_table_wikipathway_server(
      "wikipathway_table",
      pgx,
      getFilteredWikiPathwayTable,
      reactive(input$fa_contrast)
    )
  }) ## end-of-moduleServer
}
