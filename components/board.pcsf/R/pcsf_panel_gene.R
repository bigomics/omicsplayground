##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Expression plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
pcsf_genepanel_networkplot_ui <- function(id, caption, info.text, height, width) {
  ns <- shiny::NS(id)

  plot_opts <- tagList(
    withTooltip(
      shiny::radioButtons(ns("layout"), "Layout algorithm:",
        choiceNames = c("Barnes-Hut", "Hierarchical", "Kamada-Kawai"),
        choiceValues = c("BH", "hierarchical", "KK"),
        selected = "KK",
        inline = FALSE
      ),
      "Select graph layout algorithm. Barnes-Hut is a physics-based force-directed layout that is interactive. The Kamada-Kawai layout is based on a physical model of springs but is static. The hierachical layout places nodes as a hierarchical tree."
    ),
    hr(),
    withTooltip(
      radioButtons(
        ns("highlightby"),
        "Highlight labels by:",
        choices = c("centrality", "foldchange" = "prize"),
        selected = "prize",
        inline = FALSE
      ),
      "Highlight labels by scaling label size with selection."
    ),
    hr(),
    withTooltip(
      shiny::selectInput(
        ns("numlabels"), "Number of labels:", choices = c(0,3,5,10,99),
        selected=3 ),
      "Numer of labels to show."
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    title = "PCSF network analysis",
    label = "a",
    plotlib = "visnetwork",
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    options = plot_opts,
    download.fmt = c("png", "pdf", "svg"),
  )

}

#' UI code for table code: expression board
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
pcsf_genepanel_table_ui <- function(
    id,
    title,
    info.text,
    caption,
    width,
    height) {
  ns <- shiny::NS(id)

  table_opts <- shiny::tagList()

  TableModuleUI(
    ns("table"),
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    ##    options = table_opts,
    title = title,
  )
}

pcsf_genepanel_seriesplot_ui <- function(
    id,
    title,
    info.text = "",
    caption = "",
    width = 400,
    height = 400) {

  ns <- shiny::NS(id)

  PlotModuleUI(
    id = ns("plotmodule2"),
    title = "PCSF network analysis",
    plotlib = "base",
    caption = caption,
    info.text = info.text,
    height = height,
    width = width,
    ## options = plot_opts,
    download.fmt = c("png", "pdf", "svg"),
  )
}

pcsf_genepanel_settings_ui <- function(id) {
  ns <- shiny::NS(id)
  ## ??????????
}


#' PCSF network function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
pcsf_genepanel_server <- function(id,
                                  pgx,
                                  r_contrast = shiny::reactive(NULL),
                                  r_ntop = shiny::reactive(500),
                                  r_beta = shiny::reactive(0),
                                  r_cut = shiny::reactive(FALSE),
                                  r_nclust = shiny::reactive(4),
                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$physics, {
      # Update the network
      do.physics <- input$physics
      visNetwork::visNetworkProxy(ns("plotmodule-renderfigure")) %>%
        visNetwork::visPhysics(
          enabled = do.physics,
          barnesHut = list(
            centralGravity = 0  ## pull force to center [def: 0.3]
          )
        )
    })
    
    solve_pcsf <- shiny::eventReactive(
      {
        list( pgx$X, r_contrast(), r_ntop(), r_beta() )
      },
      {
        shiny::req(pgx$X)
        comparisons <- colnames(pgx$model.parameters$contr.matrix)        
        shiny::req(r_contrast() %in% comparisons)

        beta=1;ntop=400;contrast=2
        beta <- 10^as.numeric(r_beta())
        ntop <- as.integer(r_ntop())
        contrast <- r_contrast()
        
        if(pgx$datatype == "multi-omics") {
          ## Multi-omics PCSF
          info("[PcsfBoard:pcsf_compute] computing multi-omics PCSF...")
          datatypes <- unique(playbase::mofa.get_prefix(rownames(pgx$X)))
          if(all(c("gx","px") %in% datatypes)) {
            datatypes <- setdiff(datatypes, c("gx"))
          }
          info("[PcsfBoard:pcsf_compute] datatypes =", datatypes)
          pcsf <- playbase::pgx.computePCSF_multiomics(
            pgx,
            contrast,
            datatypes = datatypes,
            ntop = ntop,
            ncomp = 3,
            beta = beta,
            rm.negedge = TRUE,
            highcor = 0.8, 
            dir = "both",
            ppi = c("STRING", "GRAPHITE"),
            as.name = c("mx")
          ) 
        } else {
          ## Single-omics PCSF
          info("[PcsfBoard:pcsf_compute] computing single-omics PCSF...")
          pcsf <- playbase::pgx.computePCSF(
            pgx,
            contrast,
            level = "gene",
            ntop = ntop,
            ncomp = 2,
            beta = beta,
            dir = "both",
            rm.negedge = TRUE,
            as.name = c("mx")
          )
        }
        
        if (is.null(pcsf)) {
          validate()
          shiny::validate(
            !is.null(pcsf),
            "No PCSF solution found. Beta value is probably too small. Please adjust beta or increase network size."
          )
          return(NULL)
        }
        
        return(pcsf)
      }
    )

    gene_pcsf <- reactive({
      pcsf <- solve_pcsf()
      if(r_cut()) {
        ncomp <- as.integer(r_nclust())
        res <- playbase::pcsf.cut_and_relayout(pcsf, ncomp=ncomp)
      } else {
        pos <- igraph::layout_with_kk(pcsf, weights=NA)
        rownames(pos) <- igraph::V(pcsf)$name
        res <- list(
          graph = pcsf,
          layout = pos
        )
      }
      return(res)
    })

    
    ##-----------------------------------------------------------------
    ##---------------- visnetwork plot --------------------------------
    ##-----------------------------------------------------------------

    visnetwork.RENDER <- function() {
      sel.layout <- input$layout
      req(sel.layout, input$highlightby)
            
      ## compute PCSF
      res <- gene_pcsf()
      pcsf <- res$graph
      layoutMatrix <- res$layout

      physics <- TRUE
      if (sel.layout == "hierarchical") {
        layout <- "hierarchical"
        physics <- FALSE
        layoutMatrix <- NULL
      } else if (sel.layout == "KK") {
        layout <- "layout_with_kk"
        physics <- FALSE
      } else {
        ## barnesHut
        layout <- "layout_with_kk"
        physics <- TRUE
      }
      
      plt <- playbase::plotPCSF(
        pcsf,
        highlightby = input$highlightby,
        layout = layout,
        layoutMatrix = layoutMatrix,
        physics = physics,
        plotlib = "visnet",
        node_cex = 1.4,
        label_cex = 0.85,
        nlabel = as.integer(input$numlabels),
        border_width = 0.2,
        edge_length = 60,
        cut.clusters = FALSE,
        nlargest = -1
      )
      
      return(plt)
    }

    PlotModuleServer(
      "plotmodule",
      func = visnetwork.RENDER,
      plotlib = "visnetwork",
      pdf.width = 10,
      pdf.height = 10,
      add.watermark = watermark,
      vis.delay = 5 ## important! graph physics needs to settle
    )

    ##-----------------------------------------------------------------
    ##------------------- series plot ---------------------------------
    ##-----------------------------------------------------------------

    series.RENDER <- function() {
      sel.layout <- input$layout
      req(sel.layout, input$highlightby)

      ## compute PCSF
      res <- gene_pcsf()
      graph <- res$graph
      pos <- res$layout

      F <- playbase::pgx.getMetaMatrix(pgx, level="gene")$fc
      F <- playbase::rename_by2(F, pgx$genes, "symbol", keep.prefix=TRUE)
      F <- F[igraph::V(graph)$name,]
      colnames(F)
      
      nc <- ceiling(sqrt(ncol(F)))
      nr <- ncol(F)/nc
      par(mfrow=c(nr,nc), mar=c(1,1,4,1)*0.5)
      i=1
      for(i in 1:ncol(F)) {
        fx <- F[,i]
        playbase::plotPCSF(
          graph,
          plotlib = "igraph",
          colorby = fx,
          highlightby = "prize",
          layoutMatrix = pos,
          physics = TRUE, 
          node_cex = 1,
          nlabel = 1,
          label_cex = 0.5,
          border_width = 0.2,
          cut.clusters = FALSE,
          nlargest = -1
        )
        title( colnames(F)[i], cex.main=1.2 )
      }

    }

    PlotModuleServer(
      "plotmodule2",
      func = series.RENDER,
      plotlib = "base",
      res = c(70, 120), ## resolution of plots      
      pdf.width = 12,
      pdf.height = 6,
      add.watermark = watermark
    )
    
    ##-----------------------------------------------------------------
    ##---------------- TABLE ------------------------------------------
    ##-----------------------------------------------------------------
    
    table_data <- function() {

      shiny::req(r_contrast())
      res <- gene_pcsf()
      df <- playbase::pgx.getPCSFcentrality(
        pgx,
        contrast = r_contrast(),
        pcsf = res$graph,
        level = "gene",
        n = 100,
        plot = FALSE
      )

      ## warning. sometimes NaN
      df$centrality[is.nan(df$centrality)] <- 0
      return(df)
    }

    table.RENDER <- function(full=FALSE) {
      df <- table_data()
      
      if(full==FALSE) {
        cols <- c("symbol","gene_title","centrality","logFC")
        cols <- intersect(cols, colnames(df))
        df <- df[,cols, drop=FALSE]
      }      
      num.cols <- c("centrality", "logFC")

      dt <- ui.DataTable(
        df,
        rownames = FALSE,
        num.cols = num.cols,
        color.cols = num.cols,
        substr.cols = c("gene_title"),
        substr.len = 50
      )
            
      return(dt)
    }

    table.RENDER_modal <- function() {
      dt <- table.RENDER(full=TRUE)
      dt$x$options$scrollY <- SCROLLY_MODAL
      return(dt)
    }

    TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "none"
    )

    
  }) ## end of moduleServer
}
