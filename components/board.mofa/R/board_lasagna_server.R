##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

LasagnaBoard <- function(id, pgx) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE
    fullH <- 700 ## full height of page
    rowH1 <- 250 ## row 1 height
    rowH2 <- 440 ## row 2 height

    infotext <- tspan("<b>Weighted gene co-expression network analysis (WGCNA)</b> is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets.

<p>References:<br>
<ol>
<li>Langfelder, P. and Horvath, S., 2008. WGCNA: an R package for weighted correlation network analysis. BMC bioinformatics, 9(1), p.559.
<li>Zhang, B. and Horvath, S., 2005. A general framework for weighted gene co-expression network analysis. Statistical applications in genetics and molecular biology, 4(1).
</ol>
", js = FALSE)


    ## ============================================================================
    ## ======================= OBSERVE FUNCTIONS ==================================
    ## ============================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>LASAGNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    data <- shiny::eventReactive( pgx$mofa, {

      ##shiny::req(pgx$mofa)
      shiny::validate( shiny::need( !is.null(pgx$mofa), "missing MOFA slot"))

      ## if(is.null(pgx$mofa) || !"mofa" %in% names(pgx)) {
      ##   shinyalert::shinyalert(
      ##     title = "Error",
      ##     text = "Please compute MOFA first"
      ##   )
      ##   return(NULL)
      ## }
      
      ##shiny::req(pgx$mofa$lasagna)
      shiny::validate( shiny::need( !is.null(pgx$mofa$lasagna), "missing LASAGNA slot"))
      
      dbg("[LasagnaBoard] names(pgx$mofa) = ", names(pgx$mofa))      
      dbg("[LasagnaBoard] names(pgx$mofa$lasagna) = ", names(pgx$mofa$lasagna))
      
      res <- pgx$mofa$lasagna
      res$posx <- pgx$mofa$posx
      res$posf <- pgx$mofa$posf        
      
      dbg("[LasagnaBoard] len(res$posx) = ", length(res$posx))
      dbg("[LasagnaBoard] len(res$posf) = ", length(res$posf))      
      dbg("[LasagnaBoard] dim(res$posx[1]) = ", dim(res$posx[[1]]))
      dbg("[LasagnaBoard] dim(res$posf[1]) = ", dim(res$posf[[1]]))      
      
      ## update factors in selectInput
      ## contrasts <- colnames(pgx$contrasts)
      contrasts <- colnames(pgx$mofa$contrasts)      
      updateSelectInput(session, "contrast", choices = contrasts,
                        selected = contrasts[1])

      shiny::removeModal()      

      return(res)
    }, ignoreNULL=FALSE)

    
    ## ==========================================================================
    ## ========================== BOARD FUNCTIONS ===============================
    ## ==========================================================================

    
    ## ==========================================================================
    ## =========================== MODULES ======================================
    ## ==========================================================================
    
    mofa_plot_lasagna_server(
      "lasagna",
      data = data,
      ##input_factor = reactive(input$selected_factor),      
      watermark = WATERMARK
    )

    mofa_lasagnaSP_server(
      "lasagnaSP",
      data = data,
      input_contrast = reactive(input$contrast),
      watermark = WATERMARK
    )

    mofa_plot_clustering_server(
      "clusters",
      data = data,
      type = "features",
      input_contrast = reactive(input$contrast),
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
