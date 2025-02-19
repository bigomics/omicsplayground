##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


SNF_Board <- function(id, pgx, board_observers = NULL) {
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


    ## =============================================================================
    ## ======================= OBSERVE FUNCTIONS ===================================
    ## =============================================================================
    
    my_observers <- list()

    infotext <-
      '<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/rRIRMW_RRS4"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write;
        encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe></center>'
    
    my_observers[[1]] <- shiny::observeEvent(input$info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>WGCNA Analysis Board</strong>"),
        shiny::HTML(infotext),
        size = "xl",
        easyClose = TRUE
      ))
    })

    ## add to list global of observers. suspend by default.
    my_observers <- my_observers[!sapply(my_observers,is.null)]
    # lapply( my_observers, function(b) b$suspend() )
    if(!is.null(board_observers)) board_observers[[id]] <- my_observers
    
    ## ========================================================================
    ## ============================= REACTIVES ================================
    ## ========================================================================

    mofa <- shiny::eventReactive( pgx$mofa, {

      shiny::validate( shiny::need(!is.null(pgx$mofa), "missing MOFA slot"))      

      mofa <- pgx$mofa
      if(!"snf" %in% names(pgx$mofa)) {
        snf <- playbase::snf.cluster(mofa$xx, pheno=NULL, plot=FALSE) 
        mofa$snf <- snf
      }
      
      ## update factors in selectInput
      pheno <- colnames(mofa$samples)
      updateSelectInput(session, "selected_pheno", choices = pheno,
                        selected = pheno[1])
      
      return(mofa)
    }, ignoreNULL=FALSE)
  
    
    
    ## ========================================================================
    ## =========================== MODULES ====================================
    ## ========================================================================

    mofa_plot_snf_server(
      "snf_affinity",
      mofa = mofa,
      #type = "affinity",
      #input_pheno = reactive(input$selected_pheno),       
      watermark = WATERMARK
    )
    
    mofa_plot_snf_heatmap_server(
      "snf_heatmap",
      mofa = mofa,
      watermark = WATERMARK
    )

    mofa_plot_snfgraph_server(
      "snf_cluster",
      mofa = mofa,
      watermark = WATERMARK
    )

    return(NULL)
  })
} ## end of Board
