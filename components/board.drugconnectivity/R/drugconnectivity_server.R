##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

DrugConnectivityBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 750
    rowH <- 660 ## row height of panel
    tabH <- 200 ## row height of panel
    tabH <- "60vh" ## row height of panel

    dsea_infotext <- strwrap("<b>This module performs drug enrichment analysis</b> to see if certain drug activity or drug
        sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate
        biological functions through mechanism-of-action (MOA) and known drug molecular targets.<br><br>
        In the <a href='https://portals.broadinstitute.org/cmap/'>Drug Connectivity Map</a> panel,
        you can correlate your signature with known drug profiles from the L1000 database.
        An activation-heatmap compares drug activation profiles across multiple contrasts.
        This facilitates to quickly see and detect the similarities between contrasts for certain drugs.<br><br><br><br>
        <center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6'
        frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture'
        allowfullscreen></iframe></center>")

    shiny::observe({
      shiny::req(pgx$X)
      ct <- names(pgx$drugs)
      shiny::updateSelectInput(session, "dsea_method", choices = ct)
    })

    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$dsea_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Drug Connectivity Analysis Board</strong>"),
        shiny::HTML(dsea_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({
      shiny::req(pgx$X)
      ct <- colnames(pgx$model.parameters$contr.matrix)
      shiny::updateSelectInput(session, "dsea_contrast", choices = sort(ct))
    })

    ## =========================================================================
    ## Shared Reactive functions
    ## =========================================================================

    # common getData-esque function for drug connectivity plots / tables
    getActiveDSEA <- shiny::reactive({
      shiny::req(pgx$drugs, input$dsea_contrast, input$dsea_method)

      contr <- input$dsea_contrast
      if (is.null(contr)) {
        return(NULL)
      }

      dmethod <- input$dsea_method
      if (is.null(dmethod)) {
        return(NULL)
      }

      dr <- pgx$drugs[[dmethod]]

      nes <- round(dr$X[, contr], 4)
      pv <- round(dr$P[, contr], 4)
      qv <- round(dr$Q[, contr], 4)
      drug <- rownames(dr$X)
      if (is.null(ncol(dr$stats))) {
        stats <- dr$stats
      } else {
        stats <- dr$stats[, contr]
      }
      annot <- dr$annot
      nes[is.na(nes)] <- 0
      qv[is.na(qv)] <- 1
      pv[is.na(pv)] <- 1

      ## !!!SHOULD MAYBE BE DONE IN PREPROCESSING???
      if (is.null(annot)) {
        warning("[getActiveDSEA] WARNING:: missing drug annotation in PGX file!")
        annot <- read.csv(file.path(FILESX, "cmap/L1000_repurposing_drugs.txt"),
          sep = "\t", comment.char = "#"
        )
        rownames(annot) <- annot$pert_iname
      }

      ## compile results matrix
      jj <- match(toupper(drug), toupper(rownames(annot)))
      annot <- annot[jj, c("moa", "target")]
      dt <- data.frame(drug = drug, NES = nes, pval = pv, padj = qv, annot)
      dt <- dt[order(-dt$NES), ]

      ## sometimes UI is not ready
      if (length(input$dseatable_filter) == 0) {
        return(NULL)
      }

      if (input$dseatable_filter) {
        sel <- which(dt$moa != "" | dt$target != "")
        dt <- dt[sel, , drop = FALSE]
      }
      dsea <- list(table = dt, clust = dr$clust, stats = stats)

      return(dsea)
    })

    getMOA.target <- shiny::reactive({
      ## meta-GSEA on molecular targets
      dsea <- getActiveDSEA()
      dt <- dsea$table
      shiny::req(dt)
      targets.list <- lapply(
        enc2utf8(as.character(dt$target)),
        function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
      )
      names(targets.list) <- rownames(dt)
      targets <- setdiff(unique(unlist(targets.list)), c(NA, "", " "))
      gmt <- lapply(targets, function(g) {
        names(which(sapply(targets.list, function(t) (g %in% t))))
      })
      names(gmt) <- targets

      rnk <- dt$NES
      names(rnk) <- rownames(dt)
      suppressWarnings(
        moa.target <- fgsea::fgsea(gmt, rnk, nperm = 20000)
      )
      moa.target <- moa.target[order(-abs(moa.target$NES)), ]
      return(moa.target)
    })


    getMOA.class <- shiny::reactive({
      ## meta-GSEA on MOA terms
      dsea <- getActiveDSEA()
      dt <- dsea$table
      shiny::req(dt)
      moa.list <- lapply(
        enc2utf8(as.character(dt$moa)),
        function(s) trimws(strsplit(s, split = "[\\|;,]")[[1]])
      )
      names(moa.list) <- rownames(dt)
      moa <- setdiff(unlist(moa.list), c("", NA, " "))
      gmt <- lapply(moa, function(g) names(which(sapply(moa.list, function(t) (g %in% t)))))
      names(gmt) <- moa
      rnk <- dt$NES
      names(rnk) <- rownames(dt)
      suppressWarnings(
        #
        moa.class <- fgsea::fgsea(gmt, rnk)
      )
      moa.class <- moa.class[order(-abs(moa.class$NES)), ]
      return(moa.class)
    })

    ## =========================================================================
    ## DRUG CONNECTIVITY TAB
    ## =========================================================================

    ## -------- DSEA table
    dsea_table <- drugconnectivity_table_dsea_server(
      "dsea_table",
      getActiveDSEA
    )

    ## --------- DSEA enplot plotting module
    drugconnectivity_plot_enplots_server(
      "dsea_enplots",
      pgx,
      reactive(input$dsea_contrast),
      reactive(input$dsea_method),
      dsea_table,
      getActiveDSEA,
      watermark = WATERMARK
    )

    ## ---------- DSEA Activation map plotting module
    drugconnectivity_plot_moa_server(
      "dsea_moaplot",
      pgx,
      getActiveDSEA,
      getMOA.target,
      getMOA.class,
      watermark = WATERMARK
    )

    ## -------- Activation map plotting module
    drugconnectivity_plot_actmap_server(
      "dsea_actmap",
      pgx,
      reactive(input$dsea_contrast),
      reactive(input$dsea_method),
      dsea_table,
      getActiveDSEA,
      watermark = WATERMARK
    )



    ## =======================================================================================
    ## CONNECTIVITY MAP TAB
    ## =======================================================================================

    drugconnectivity_plot_cmap_enplot_server(
      "cmap_enplot",
      pgx,
      getActiveDSEA,
      cmap_table,
      watermark = WATERMARK
    )

    drugconnectivity_plot_cmap_dsea_server(
      "cmap_dsea",
      pgx = pgx,
      getActiveDSEA = getActiveDSEA,
      cmap_table = cmap_table,
      getMOA.class = getMOA.class,
      getMOA.target = getMOA.target,
      dsea_method = reactive(input$dsea_method),
      dsea_contrast = reactive(input$dsea_contrast),
      watermark = WATERMARK
    )

    cmap_table <- drugconnectivity_table_cmap_server(
      "cmap_table",
      getActiveDSEA
    )
  })
}
