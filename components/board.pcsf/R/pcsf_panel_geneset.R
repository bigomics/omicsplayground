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
pcsf_gsetpanel_networkplot_ui <- function(id, caption, info.text, height, width) {
  ns <- shiny::NS(id)

  plot_opts <- tagList(
    withTooltip(
      radioButtons(
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
        ns("numlabels"), "Number of labels:",
        choices = c(0, 3, 5, 10, 20, 50, 999),
        selected = 20
      ),
      "Numer of labels to show."
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    title = "Interactive network",
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
pcsf_gsetpanel_table_ui <- function(
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

pcsf_gsetpanel_seriesplot_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    width = 400,
    height = 400) {
  ns <- shiny::NS(id)

  options <- tagList(
    withTooltip(
      shiny::selectInput(
        ns("series_numlabels"), "Number of labels:",
        choices = c(0, 3, 5, 10, 20, 999),
        selected = 0
      ),
      "Numer of labels to show."
    )
  )

  PlotModuleUI(
    id = ns("plotmodule2"),
    title = "Network by comparison",
    plotlib = "base",
    caption = caption,
    info.text = info.text,
    options = options,
    height = height,
    width = width,
    ## options = plot_opts,
    download.fmt = c("png", "pdf", "svg"),
  )
}

pcsf_gsetpanel_settings_ui <- function(id) {
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
      c(
        "rho as prize" = "rho_prize",
        "meta solution" = "meta_solution",
        "add VHC edges" = "add_vhce"
      ),
      selected = c("add_vhce"), inline = TRUE
    ),
    hr(),
    withTooltip(
      shiny::selectInput(ns("layout.method"), "Layout algorithm:",
        choices = c("Kamada-Kawai" = "kk", "tree", "circle"),
        selected = "kk"
      ),
      "Select graph layout algorithm. Barnes-Hut is a physics-based force-directed layout that is interactive. The Kamada-Kawai layout is based on a physical model of springs but is static. The hierachical layout places nodes as a hierarchical tree."
    ),
    shiny::checkboxInput(ns("physics"), "enable physics", FALSE),
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
          choices = c(1, 4, 9, 16, 25, 99), selected = 16
        ),
        "Maximum number of components"
      ),
      shiny::selectInput(ns("resolution"), "Resolution",
        choices = c(0.01, 0.05, 0.1, 0.2, 0.5, 1), selected = 0.1
      )
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
pcsf_gsetpanel_server <- function(id,
                                  pgx,
                                  r_contrast = shiny::reactive(NULL),
                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observeEvent(input$physics, {
      # Update the network
      dbg("[pcsf_gsetpanel_server:physics] id =", ns("plotmodule-renderfigure"))
      visNetwork::visNetworkProxy(ns("plotmodule-renderfigure")) %>%
        visNetwork::visPhysics(
          enabled = input$physics,
          barnesHut = list(
            centralGravity = 0 ## pull force to center [def: 0.3]
          )
        )
    })

    singlecontrast <- reactive({
      if ("meta_solution" %in% input$solve_options) {
        return(NULL)
      }
      r_contrast()
    })

    solve_pcsf <- shiny::eventReactive(
      {
        list(pgx$X, input$ntop, singlecontrast(), input$solve_options)
      },
      {
        shiny::req(pgx$X)
        ntop <- as.integer(input$ntop)

        contrast <- singlecontrast()
        comparisons <- playbase::pgx.getContrasts(pgx)
        if (!"meta_solution" %in% input$solve_options) shiny::req(contrast %in% comparisons)

        ## If this is metabolomics-PPI, let's use ext2 only
        if (any(grepl("ext2", rownames(pgx$gsetX)))) {
          gset.filter <- "ext2"
        } else {
          gset.filter <- NULL
        }

        pcsf <- playbase::pgx.computePCSF_gset(
          pgx,
          contrast = contrast,
          gset.filter = gset.filter,
          as_prize = ifelse("rho_prize" %in% input$solve_options, "rho", "fc"),
          use_rank = FALSE,
          gmt.rho = 0.8,
          highcor = ifelse("add_vhce" %in% input$solve_options, 0.9, Inf),
          ntop = ntop,
          ncomp = 5,
          beta = 1,
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

    gset_pcsf <- reactive({
      pcsf <- solve_pcsf()

      input_cut <- input$cut
      input_nclust <- input$nclust
      input_resolution <- input$resolution

      if (input_cut) {
        ncomp <- as.integer(input_nclust)
        res <- playbase::pcsf.cut_and_relayout(
          pcsf,
          ncomp = ncomp,
          layout = input$layout.method,
          cluster.method = c("louvain", "leiden")[2],
          leiden.resolution = as.numeric(input_resolution)
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


    ## -----------------------------------------------------------------
    ## ---------------- visnetwork plot --------------------------------
    ## -----------------------------------------------------------------

    visnetwork.RENDER <- function() {
      req(input$highlightby)

      ## compute PCSF
      res <- gset_pcsf()
      pcsf <- res$graph
      layoutMatrix <- res$layout

      physics <- FALSE
      physics <- input$physics

      contrast <- r_contrast()
      comparisons <- colnames(pgx$model.parameters$contr.matrix)
      shiny::req(contrast %in% comparisons)

      F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
      fx <- F[igraph::V(pcsf)$name, contrast]
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
        node_cex = 3.5,
        node_gamma = 0.66,
        label_cex = 0.65,
        nlabel = as.integer(input$numlabels),
        border_width = 0,
        edge_length = 40,
        cut.clusters = FALSE,
        nlargest = -1,
        hilight_single = table$rownames_selected()
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

    ## -----------------------------------------------------------------
    ## ------------------- series plot ---------------------------------
    ## -----------------------------------------------------------------

    series.RENDER <- function() {
      req(input$highlightby)

      ## compute PCSF
      res <- gset_pcsf()
      graph <- res$graph
      layout <- res$layout

      F <- playbase::pgx.getMetaMatrix(pgx, level = "geneset")$fc
      F <- F[igraph::V(graph)$name, , drop = FALSE]

      hh <- grep("^IA:", colnames(F))
      if (length(hh) > 0) F <- F[, -hh, drop = FALSE]

      nc <- ceiling(1.3 * sqrt(ncol(F)))
      nr <- ceiling(ncol(F) / nc)
      par(mfrow = c(nr, nc), mar = c(1, 1, 4, 1) * 0.5)

      i <- 1
      for (i in 1:ncol(F)) {
        fx <- F[, i]
        fx[is.na(fx)] <- 0 # NA no size no color
        playbase::plotPCSF(
          graph,
          sizeby = fx,
          colorby = fx,
          plotlib = "igraph",
          highlightby = input$highlightby,
          layoutMatrix = layout,
          node_cex = 1,
          nlabel = as.integer(input$series_numlabels),
          label_cex = 0.6,
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

    ## -----------------------------------------------------------------
    ## ---------------- TABLE ------------------------------------------
    ## -----------------------------------------------------------------

    table_data <- function() {
      shiny::req(r_contrast())
      res <- gset_pcsf()
      df <- playbase::pgx.getPCSFcentrality(
        pgx,
        contrast = r_contrast(),
        pcsf = res$graph,
        level = "geneset",
        n = 100,
        plot = FALSE
      )

      ## warning. sometimes NaN
      df$centrality[is.nan(df$centrality)] <- 0
      return(df)
    }

    table.RENDER <- function(full = FALSE) {
      df <- table_data()

      if (full == FALSE) {
        cols <- c("geneset", "logFC", "centrality")
        cols <- intersect(cols, colnames(df))
        df <- df[, cols, drop = FALSE]
      }
      num.cols <- intersect(c("centrality", "logFC"), colnames(df))

      dt <- ui.DataTable(
        df,
        rownames = FALSE,
        num.cols = num.cols,
        color.cols = num.cols,
        substr.cols = c("geneset"),
        substr.len = 80
      )

      return(dt)
    }

    table.RENDER_modal <- function() {
      dt <- table.RENDER(full = TRUE)
      dt$x$options$scrollY <- SCROLLY_MODAL
      return(dt)
    }

    table <- TableModuleServer(
      "table",
      func = table.RENDER,
      func2 = table.RENDER_modal,
      selector = "single"
    )
  }) ## end of moduleServer
}
