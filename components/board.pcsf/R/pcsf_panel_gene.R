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
      shiny::radioButtons(
        ns("highlightby"),
        "Highlight labels by:",
        choices = c("foldchange" = "prize", "centrality"),
        selected = "prize",
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
    shiny::checkboxGroupInput(
      ns("solve_options"),
      "Solve options:",
      c("rho as prize" = "rho_prize",
        "meta solution" = "meta_solution",
        "add VHC edges" = "add_vhce" 
      ),
      selected = c("add_vhce"), inline = TRUE
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
      if("meta_solution" %in% input$solve_options) return(NULL)
      r_contrast()
    })
    
    solve_pcsf <- shiny::eventReactive(
      {
        list( pgx$X, input$ntop, singlecontrast(), input$solve_options)        
      },
      {
        shiny::req(pgx$X)
        ntop <- as.integer(input$ntop)

        contrast <- singlecontrast()
        comparisons <- playbase::pgx.getContrasts(pgx)
        if(!"meta_solution" %in% input$solve_options) shiny::req(contrast %in% comparisons)
        
        if(pgx$datatype == "multi-omics") {
          ## Multi-omics PCSF
          info("[PcsfBoard:pcsf_compute] computing multi-omics PCSF...")
          ## If both gx&px are in datasest, we prefer proteomics (px)
          ## for building the PCSF.
          datatypes <- unique(playbase::mofa.get_prefix(rownames(pgx$X)))
          if(all(c("gx","px") %in% datatypes)) {
            datatypes <- setdiff(datatypes, c("gx"))
          }
        } else {
          info("[PcsfBoard:pcsf_compute] computing single-omics PCSF...")          
            datatypes <- NULL
        }
        
        info("[PcsfBoard:pcsf_compute] datatypes =", datatypes)
        pcsf <- playbase::pgx.computePCSF(
          pgx,
          contrast = contrast,
          datatypes = datatypes,
          as_prize = ifelse("rho_prize" %in% input$solve_options,"rho","fc"),
          use_rank = FALSE,
          ntop = ntop,
          ncomp = 5,
          beta = 1,
          rm.negedge = TRUE,
          highcor = ifelse("add_vhce" %in% input$solve_options, 0.9, Inf),
          dir = "both",
          ppi = c("STRING", "GRAPHITE")
        ) 
        
        if (is.null(pcsf)) {
          validate()
          shiny::validate(
            !is.null(pcsf),
            "No PCSF solution found. Perhaps beta value is probably too small. Please adjust beta or increase network size."
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

      dbg("[gene_pcsf] 1: layout.method = ", input$layout.method)
      
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
      shiny::req(contrast)      
      comparisons <- colnames(pgx$model.parameters$contr.matrix)        
      shiny::req(contrast %in% comparisons)
      
      F <- playbase::pgx.getMetaMatrix(pgx, level='gene')$fc
      fx <- F[igraph::V(pcsf)$name, contrast]

      labels <- names(fx)
      labels <- playbase::probe2symbol(names(fx), pgx$genes, "gene_name", fill_na = TRUE)
      
      igraph::V(pcsf)$foldchange <- fx
      igraph::V(pcsf)$prize <- abs(fx)

      plt <- playbase::plotPCSF(
        pcsf,
        sizeby = fx,
        colorby = fx,        
        highlightby = input$highlightby,
        layout = "layout.norm",
        layoutMatrix = layoutMatrix,
        physics = physics,
        plotlib = "visnet",
        node_cex = 6,
        node_gamma = 0.5,
        labels = labels,
        label_cex = 1,
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
      F <- F[igraph::V(graph)$name,,drop=FALSE]

      labels <- rownames(F)
      labels <- playbase::probe2symbol(rownames(F), pgx$genes, "gene_name", fill_na = TRUE)

      hh <- grep("^IA:", colnames(F))
      if (any(hh)) F <- F[, -hh, drop = FALSE]

      nc <- ceiling(1.3*sqrt(ncol(F)))
      nr <- ceiling(ncol(F)/nc)
      par(mfrow=c(nr,nc), mar=c(1,1,4,1)*0.5)      
      i=1
      for(i in 1:ncol(F)) {
        fx <- F[,i]        
        fx[is.na(fx)] <- 0 # NA no size no color
        playbase::plotPCSF(
          graph,
          plotlib = "igraph",
          sizeby = fx,
          colorby = fx,          
          highlightby = input$highlightby,
          layoutMatrix = layout,
          physics = TRUE, 
          node_cex = 1,
          node_gamma = 0.66,
          labels = labels,
          nlabel = as.integer(input$series_numlabels),
          label_cex = 0.7,
          edge_cex = 0.85,
          border_width = 0.2,
          cut.clusters = FALSE,
          nlargest = -1
        )
        title(colnames(F)[i], cex.main = 1.2)
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

      dbg("[pcsf_genepanel_server:table.RENDER] dim.df = ", dim(df))
      dbg("[pcsf_genepanel_server:table.RENDER] colnames.df = ", colnames(df))
      
      if(full==FALSE) {
        cols <- c("symbol","gene_title","logFC","centrality")
        cols <- intersect(cols, colnames(df))
        df <- df[,cols, drop=FALSE]
      }      
      num.cols <- intersect(c("centrality", "logFC"),colnames(df))
      
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
