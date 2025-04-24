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
pcsf_gsetnetwork_plot_ui <- function(id, caption, info.text, height, width) {
  ns <- shiny::NS(id)

  plot_opts <- tagList(
    withTooltip(
      shiny::radioButtons(ns("layout"), "Layout algorithm:",
        choiceNames = c("Barnes-Hut", "Hierarchical", "Kamada-Kawai"),
        choiceValues = c("BH", "hierarchical", "KK"),
        selected = "BH",
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
        selected = "centrality",
        inline = FALSE
      ),
      "Highlight labels by scaling label size with selection."
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
pcsf_gsetnetwork_table_ui <- function(
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

pcsf_gsetnetwork_settings_ui <- function(id) {
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
pcsf_gsetnetwork_server <- function(id,
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
    
    gset_pcsf <- shiny::eventReactive(
      {
        list( pgx$X, r_contrast(), r_ntop(), r_beta(), r_nclust() )
      },
      {
        shiny::req(pgx$X)
        comparisons <- colnames(pgx$model.parameters$contr.matrix)        
        shiny::req(r_contrast() %in% comparisons)
        
        beta <- 10^as.numeric(r_beta())
        ntop <- as.integer(r_ntop())
        contrast <- r_contrast()
        
        pcsf <- playbase::pgx.computePCSF_gset(
          pgx,
          contrast,
          gset.filter = "ext2",              
          gmt.rho = 0.8,
          highcor = 0.9,
          ntop = ntop,
          ncomp = 5,
          beta = beta,
          dir = "both",
          rm.negedge = TRUE
        )
        
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
    
    visnetwork.RENDER <- function() {
      sel.layout <- input$layout
      req(sel.layout, input$highlightby)
      
      physics <- TRUE
      if (sel.layout == "hierarchical") {
        layout <- "hierarchical"
        physics <- FALSE
      } else if (sel.layout == "KK") {
        layout <- "layout_with_kk"
        physics <- FALSE
      } else {
        ## barnesHut
        layout <- "layout_with_kk"
        physics <- TRUE
      }
      
      ## compute PCSF
      pcsf <- gset_pcsf()
      
      nclust <- as.integer(r_nclust())
      dbg("[gset_pcsf] nclust=",nclust)

      plt <- playbase::plotPCSF(
        pcsf,
        highlightby = input$highlightby,
        layout = layout,
        physics = physics,
        plotlib = "visnet",
        node_cex = 30,
        label_cex = 30,
        nlabel = 10,
        border_width = 0.4,
        cut.clusters = r_cut(),
        nlargest = nclust
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
    ##---------------- TABLE ------------------------------------------
    ##-----------------------------------------------------------------
    
    table_data <- function() {

      shiny::req(r_contrast())
      df <- playbase::pgx.getPCSFcentrality(
        pgx,
        contrast = r_contrast(),
        pcsf = gset_pcsf(),
        level = "geneset",
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
        cols <- c("geneset","centrality","logFC")
        cols <- intersect(cols, colnames(df))
        df <- df[,cols, drop=FALSE]
      }      
      num.cols <- c("centrality", "logFC")

      dt <- ui.DataTable(
        df,
        rownames = FALSE,
        num.cols = num.cols,
        color.cols = num.cols,
        substr.cols = c("geneset"),
        substr.len = 60
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
