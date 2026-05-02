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

    infotext <- strwrap("<b>This module performs drug enrichment analysis</b> to see if certain drug activity or drug
        sensitivity signatures matches your experimental signatures. Matching drug signatures to your experiments may elicudate
        biological functions through mechanism-of-action (MOA) and known drug molecular targets.<br><br>
        In the <a href='https://portals.broadinstitute.org/cmap/'>Drug Connectivity Map</a> panel,
        you can correlate your signature with known drug profiles from the L1000 database.
        An activation-heatmap compares drug activation profiles across multiple contrasts.
        This facilitates to quickly see and detect the similarities between contrasts for certain drugs.<br><br><br><br>
        <center><iframe width='560' height='315' src='https://www.youtube.com/embed/BtMQ7Y0NoIA?si=3T61_k_onEqsTMcr&amp;start=91' title='YouTube video player' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share' referrerpolicy='strict-origin-when-cross-origin' allowfullscreen></iframe></center>")

    ## ================================================================================
    ## ============================== OBSERVERS  ======================================
    ## ================================================================================

    shiny::observe({
      shiny::req(pgx$X)
      ct <- names(pgx$drugs)
      shiny::updateSelectInput(session, "method", choices = ct)
    })

    shiny::observeEvent(input$dsea_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Drug Connectivity Analysis Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    shiny::observe({
      shiny::req(pgx$X)
      ct <- playbase::pgx.getContrasts(pgx)
      ct <- sort(ct[!grepl("^IA:", ct)])
      shiny::updateSelectInput(session, "contrast", choices = ct)
    })

    # Observe tabPanel change to update Settings visibility
    tab_elements <- list(
      "Drug enrichment" = list(disable = c("aiui")),
      "Connectivity map (beta)" = list(disable = c("aiui")),
      "AI Summary✨" = list(disable = c("filter_table","contrast"))
    )

    shiny::observeEvent(input$tabs, {
      bigdash::update_tab_elements(input$tabs, tab_elements)
    })
    
    ## =========================================================================
    ## Shared Reactive functions
    ## =========================================================================

    get_pgx_drugs <- eventReactive( pgx$drugs, {

      ## check if we have precomputed MOA slots
      if(!is.null(pgx$drugs[[1]]$moa)) {
        return(pgx$drugs)
      }
      
      ## need to recompute MOA (if not done), NOTE: should be done in
      ## pgx computation
      dbg("[DrugConnectivityBoard::get_pgx_drugs] Computing MOA...")
      pgx.showSmallModal("Calculating MOA<br>Please wait...")      
      pgx$drugs$report <- NULL
      db=names(pgx$drugs)[1]
      for(db in names(pgx$drugs)) {
        res <- pgx$drugs[[db]]
        if(is.null(res$moa)) {        
          moa <- metaLINCS::computeMoaEnrichment(res)
          pgx$drugs[[db]][['moa']] <- moa
        }
      }
      shiny::removeModal(session)
      return(pgx$drugs)
    })
    
    # common getData-esque function for drug connectivity plots / tables
    getActiveDSEA <- shiny::reactive({
      
      contr <- input$contrast
      dmethod <- input$method
      shiny::req(contr, dmethod)

      pgxdrugs <- get_pgx_drugs()
      shiny::req(pgxdrugs)

      ## sometimes UI is not ready
      if (length(input$filter_table) == 0) {
        return(NULL)
      }
      do.filter <- input$filter_table
      
      dsea <- playbase::pgx.getDrugConnectivityTable(
        pgx=NULL, contrast=contr, db=dmethod,
        drugs=pgxdrugs, filter=do.filter)
      
      return(dsea)
    })

    getMOA.target <- shiny::reactive({
      contr <- input$contrast
      db <- input$method
      shiny::req(contr,db)
      pgxdrugs <- get_pgx_drugs()
      shiny::req(pgxdrugs)
      moa <- playbase::pgx.getDrugMOATable(
        pgx=NULL, contrast=contr, db=db,
        drugs=pgxdrugs, type="targetGene")       
      return(moa)
    })

    getMOA.class <- shiny::reactive({
      contr <- input$contrast
      db <- input$method
      shiny::req(contr,db)
      pgxdrugs <- get_pgx_drugs()
      shiny::req(pgxdrugs)      
      moa <- playbase::pgx.getDrugMOATable(
        pgx=NULL, contrast=contr, db=db,
        drugs=pgxdrugs, type="drugClass")       
      return(moa)
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
      reactive(input$contrast),
      reactive(input$method),
      dsea_table,
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
      reactive(input$contrast),
      reactive(input$method),
      dsea_table,
      getActiveDSEA,
      watermark = WATERMARK
    )

    ## ==================================================================================
    ## Module servers
    ## ==================================================================================
    
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
      dsea_method = reactive(input$method),
      dsea_contrast = reactive(input$contrast),
      watermark = WATERMARK
    )

    cmap_table <- drugconnectivity_table_cmap_server(
      "cmap_table",
      getActiveDSEA
    )

    drugconnectivity_report_server(
      "cmap_report",
      pgx = pgx,
      drugs = get_pgx_drugs,
      rdb = reactive(input$method)
    )
    


  })
}
