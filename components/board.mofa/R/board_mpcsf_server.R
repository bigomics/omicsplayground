##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

MPCSF_Board <- function(id, pgx) {
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


    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'

    shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    mofa <- shiny::eventReactive( pgx$mofa, {

      shiny::validate( shiny::need( !is.null(pgx$mofa), "missing MOFA slot"))

      data <- pgx$mofa
      snf <- playbase::snf.cluster(data$xx, pheno=NULL, plot=FALSE) 
      snf$samples <- data$samples
      
      ## update factors in selectInput
      pheno <- colnames(data$samples)
      updateSelectInput(session, "selected_pheno", choices = pheno,
                        selected = pheno[1])

      mofa <- pgx$mofa
      mofa$snf <- snf
      
      return(mofa)
    }, ignoreNULL=FALSE)


    
    ## ================================================================================
    ## ========================== BOARD FUNCTIONS =====================================
    ## ================================================================================

    
    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    mofa_plot_dummy_server(
      "dummy1",
      mofa = mofa,
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
