##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


MGseaBoard <- function(id, pgx) {
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

      has.mofa <- "mofa" %in% names(pgx) && !is.null(pgx$mofa)
      shiny::validate( shiny::need( has.mofa, "missing MOFA slot"))      
      
      ## update factors in selectInputs
      enr <- pgx$mofa$fc.gsea      
      contrasts <- names(enr)
      updateSelectInput(session, "contrast", choices = contrasts)
      
      return( pgx$mofa )
    })

    
    ## =====================================================================
    ## ======================= PLOTTING FUNCTIONS ==========================
    ## =====================================================================

    ## =====================================================================
    ## ===================== MODULES =======================================
    ## =====================================================================

    mgsea_table_selected <- reactive({

      shiny::req(mgsea_table$data())
      
      has.selection <- length(mgsea_table$rownames_selected())>0 
      search_key <- mgsea_table$search()
      has.search <- length(search_key)>0 && search_key[1]!=""
      
      if(has.search && !has.selection) {
        sel <- mgsea_table$rownames_all()
      } else if(has.selection) {
        sel <- mgsea_table$rownames_selected()
      } else {
        ##sel <- NULL
        sel <- head( mgsea_table$rownames_all(), 20)
      }
      sel
    })
    
    mofa_plot_enrichment_server(
      "menrichment",
      pgx = pgx,
      gsea = reactive({ mofa()$fc.gsea }),
      input_k = reactive(input$contrast),
      select = mgsea_table_selected,
      req.selection = TRUE,
      watermark = WATERMARK
    )

    mofa_plot_mgsea_server(
      "mgsea_plot",
      gsea = reactive(mofa()$fc.gsea),
      input_k = reactive(input$contrast),
      select = mgsea_table_selected,
      watermark = WATERMARK
    )

    mofa_plot_pathbank_server(
      "pathbank_pathway",
      pgx = pgx,
      sel_pathway = mgsea_table_selected,
      sel_contrast = reactive(input$contrast),
      watermark = WATERMARK
    )

    
    # Table Modules
    mgsea_table <- mofa_table_mgsea_server(
      "mgsea_table",
      gsea = reactive(mofa()$fc.gsea),
      datatypes = NULL,
      input_k = reactive(input$contrast)            
    )

    return(NULL)
  })
} ## end of Board
