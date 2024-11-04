##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

PathwayBoard <- function(id,
                         pgx,
                         selected_gsetmethods = reactive(colnames(pgx$gset.meta$meta[[1]]$fc))) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 750
    rowH <- 660 ## row height of panel
    tabH <- "70vh" ## row height of panel

    fa_infotext <- tspan(paste("This module performs specialized pathway analysis.
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
    gyroscope; picture-in-picture' allowfullscreen></iframe></center>"), js = FALSE)

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
      "GO graph" = list(disable = NULL),
      "Enrichment Map (beta)" = list(disable = NULL)
    )
    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })


    ## ================================================================================
    ## =========================== FUNCTIONS ==========================================
    ## ================================================================================


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
      ## sometimes no REACTOME in genesets...
      if (nrow(df) == 0) {
        shinyalert::shinyalert(
          title = "No REACTOME terms in enrichment results",
          text = "",
          type = "warning"
        )
        df <- data.frame()
        return(df)
      }
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
      getFilteredReactomeTable,
      pgx = pgx,
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
      fa_filtertable = reactive(input$fa_filtertable),
      fa_filtertable_value = reactive(input$fa_filtertable_value),
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
      getFilteredWikiPathwayTable,
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
