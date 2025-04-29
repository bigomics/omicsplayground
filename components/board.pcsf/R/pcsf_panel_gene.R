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
      radioButtons(
        ns("highlightby"),
        "Highlight labels by:",
        choices = c("centrality", "foldchange" = "prize"),
        selected = "centrality",
        inline = TRUE
      ),
      "Highlight labels by scaling label size with selection."
    ),
    hr(),
    withTooltip(
      shiny::selectInput(
        ns("numlabels"), "Number of labels:", choices = c(0,3,5,10,20,999),
        selected = 10 ),
      "Numer of labels to show."
    )
  )
  
  PlotModuleUI(
    id = ns("plotmodule"),
    title = "Interactive network",
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
  options <- tagList(
    withTooltip(
      shiny::selectInput(
        ns("series_numlabels"), "Number of labels:", choices = c(0,3,5,10,20,999),
        selected = 0 ),
      "Numer of labels to show."
    )
  )

  PlotModuleUI(
    id = ns("plotmodule2"),
    title = "Network by comparison",
    plotlib = "base",
    options = options,
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
  tagList(
    withTooltip(
      shiny::radioButtons(ns("ntop"), "Network size:",
        choices = c("S" = 250, "M" = 750, "L" = 1500),
        selected = 750, inline = TRUE
      ),
      "Select initial network size (number of top genes) for ."
    ),
    hr(),
    withTooltip(
      shiny::radioButtons(ns("mode"), "Solution mode:",
        choices = c("single","common"),
        selected = "single", inline = TRUE
      ),
      "Select solution mode. 'single' solves for selected comparison, 'common' solves across comparisons."
    ),
    hr(),
    withTooltip(
      shiny::selectInput(ns("layout.method"), "Layout algorithm:",
        choices = c("Kamada-Kawai" = "kk","tree","circle"),        
        selected = "kk"
      ),
      "Select graph layout algorithm. Barnes-Hut is a physics-based force-directed layout that is interactive. The Kamada-Kawai layout is based on a physical model of springs but is static. The hierachical layout places nodes as a hierarchical tree."
    ),
    shiny::checkboxInput(ns("physics"),"enable physics",FALSE),    
    hr(),
    withTooltip(
      shiny::checkboxInput(ns("cut"), "Cut clusters", FALSE),
      "Cut network into smaller clusters"
    ),
    shiny::conditionalPanel(
      "input.cut == true",
      ns = ns,
      withTooltip(
        shiny::selectInput(ns("nclust"), "Number of clusters",
          choices = c(1,4,9,16,25,99), selected=16),
        "Maximum number of components"
      ),
      shiny::selectInput(ns("resolution"), "Resolution",
        choices = c(0.01,0.05,0.1,0.2,0.5,1), selected=0.1)
    )
  )
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
                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$physics, {
      # Update the network
      dbg("[pcsf_gene_server:physics] id =", ns("plotmodule-renderfigure"))
      visNetwork::visNetworkProxy(ns("plotmodule-renderfigure")) %>%
        visNetwork::visPhysics(
          enabled = input$physics,
          barnesHut = list(
            centralGravity = 0  ## pull force to center [def: 0.3]
          )
        )
    })

    singlecontrast <- reactive({
      if(input$mode == "common") return(NULL)
      r_contrast()
    })
    
    solve_pcsf <- shiny::eventReactive(
      {
        list( pgx$X, input$ntop, singlecontrast() )        
      },
      {
        shiny::req(pgx$X)
        ntop <- as.integer(input$ntop)

        contrast <- singlecontrast()
        comparisons <- playbase::pgx.getContrasts(pgx)
        if(input$mode == "single") shiny::req(contrast %in% comparisons)
        
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
            contrast = contrast,
            datatypes = datatypes,
            ntop = ntop,
            ncomp = 10,
            beta = 1,
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
            contrast = contrast,
            level = "gene",
            ntop = ntop,
            ncomp = 10,
            beta = 1,
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

      shiny::req(input$layout.method)
      shiny::req(input$nclust)
      
      pcsf <- solve_pcsf()

      input_cut <- input$cut
      input_nclust <- input$nclust
      input_resolution <- input$resolution
      
      if(input_cut) {
        ncomp <- as.integer(input$nclust)
        res <- playbase::pcsf.cut_and_relayout(
          pcsf,
          ncomp = ncomp,
          layout = input$layout.method,
          cluster.method = c("louvain","leiden")[2],
          leiden.resolution = as.numeric(input$resolution)
        )
      } else {
        res <- playbase::pcsf.cut_and_relayout(
          pcsf,
          cut = FALSE,
          layout = input$layout.method,
          component.wise = FALSE,
          ncomp = -1,
          as_grid = FALSE
        )
      }
      
      return(res)
    })

    
    ##-----------------------------------------------------------------
    ##---------------- visnetwork plot --------------------------------
    ##-----------------------------------------------------------------

    visnetwork.RENDER <- function() {
      
      req(input$highlightby)
      
      ## compute PCSF
      res <- gene_pcsf()
      pcsf <- res$graph
      layoutMatrix <- res$layout
      
      physics <- FALSE
      physics <- input$physics
      
      contrast <- r_contrast()
      comparisons <- colnames(pgx$model.parameters$contr.matrix)        
      shiny::req(contrast %in% comparisons)

      F <- playbase::pgx.getMetaMatrix(pgx, level='gene')$fc
      F <- playbase::rename_by2(F, pgx$genes, "symbol", keep.prefix=TRUE)

      shiny::req(contrast)      
      fx <- F[igraph::V(pcsf)$name, contrast]
      igraph::V(pcsf)$foldchange <- fx
      igraph::V(pcsf)$prize <- abs(fx)
      
      plt <- playbase::plotPCSF(
        pcsf,
        colorby = fx,
        highlightby = input$highlightby,
        layout = "layout.norm",
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
      req(input$highlightby)
      
      ## compute PCSF
      res <- gene_pcsf()
      graph <- res$graph
      layout <- res$layout
      
      F <- playbase::pgx.getMetaMatrix(pgx, level="gene")$fc
      F <- playbase::rename_by2(F, pgx$genes, "symbol", keep.prefix=TRUE)
      F <- F[igraph::V(graph)$name,]
      
      nc <- ceiling(1.3*sqrt(ncol(F)))
      nr <- ceiling(ncol(F)/nc)
      par(mfrow=c(nr,nc), mar=c(1,1,4,1)*0.5)
      
      i=1
      for(i in 1:ncol(F)) {
        fx <- F[,i]
        
        playbase::plotPCSF(
          graph,
          plotlib = "igraph",
          colorby = fx,
          highlightby = input$highlightby,
          layoutMatrix = layout,
          physics = TRUE, 
          node_cex = 1,
          nlabel = as.integer(input$series_numlabels),
          label_cex = 0.7,
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
