##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Single cell plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#' @param width
#'
#' @export
singlecell_plot_markersplot_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width,
  parent) {
  ns <- shiny::NS(id)

  markersplot.opts <- shiny::tagList(
    withTooltip(shiny::selectInput(parent("mrk_level"), "Level:", choices = c("gene", "geneset")),
      "Specify the level of the marker analysis: gene or gene set level.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::selectInput(parent("mrk_features"), "Feature set:",
        choices = NULL,
        multiple = FALSE
      ),
      "Select a particular functional group for the analysis.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(shiny::textInput(parent("mrk_search"), "Filter:"),
      "Filter markers by a specific keywords.",
      placement = "top", options = list(container = "body")
    ),
    withTooltip(
      shiny::radioButtons(parent("mrk_sortby"), "Sort by:",
        choices = c("intensity", "name"), inline = TRUE
      ),
      "Sort by name or intensity.",
      placement = "top",
      options = list(container = "body")
    )
  )

  PlotModuleUI(
    id = ns("plotmodule"),
    #    plotlib = "plotly",
    plotlib = "ggplot",          
    label = label,
    info.text = info.text,
    title = title,
    caption = caption,
    options = markersplot.opts,
    download.fmt = c("png", "csv"), #FIXME pdf is not working, to avoid crashing other things, we decided to remove it
    height = height,
    width = width
  )
}

#' Single cell plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @export
singlecell_plot_markersplot_server <- function(id,
                                               pgx,
                                               pfGetClusterPositions,
                                               mrk_level,
                                               mrk_features,
                                               mrk_search,
                                               mrk_sortby,
                                               watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- shiny::reactive({
      ## if(!input$tsne.all) return(NULL)

      shiny::req(pgx)
      
      mrk_level <- mrk_level()
      mrk_features <- mrk_features()
      mrk_search <- mrk_search()
      mrk_sortby <- mrk_sortby()

      clust.pos <- pfGetClusterPositions()
      if (is.null(clust.pos)) {
        return(NULL)
      }
      ## pos <- pgx$tsne2d
      pos <- clust.pos

      ## markers <- pgx$families[["CD family"]]
      if (is.null(mrk_features)) {
        return(NULL)
      }
      if (mrk_features == "") {
        return(NULL)
      }

      term <- ""
      if (mrk_level == "gene") {
        markers <- pgx$families[["Transcription factors (ChEA)"]]
        if (mrk_search != "") {
          term <- mrk_search
          jj <- grep(term, pgx$genes$gene_name, ignore.case = TRUE)
          markers <- pgx$genes$gene_name[jj]
          term <- paste("filter:", term)
        } else if (mrk_features %in% names(pgx$families)) {
          markers <- pgx$families[[mrk_features]]
          term <- mrk_features
        } else {
          markers <- pgx$genes$gene_name
        }
        ## markers <- intersect(markers, rownames(pgx$X))
        markers <- intersect(toupper(markers), toupper(pgx$genes$gene_name))
        jj <- match(markers, toupper(pgx$genes$gene_name))
        pmarkers <- intersect(rownames(pgx$genes)[jj], rownames(pgx$X))
        gx <- pgx$X[pmarkers, rownames(pos), drop = FALSE]
      } else if (mrk_level == "geneset") {
        ## markers <- pgx$families[["Immune checkpoint (custom)"]]
        markers <- COLLECTIONS[[1]]
        if (is.null(mrk_features)) {
          return(NULL)
        }
        ft <- mrk_features
        if (mrk_search == "" && ft %in% names(COLLECTIONS)) {
          markers <- COLLECTIONS[[mrk_features]]
          markers <- intersect(markers, rownames(pgx$gsetX))
          term <- mrk_features
        } else if (mrk_search != "") {
          term <- mrk_search
          jj <- grep(term, rownames(pgx$gsetX), ignore.case = TRUE)
          markers <- rownames(pgx$gsetX)[jj]
          term <- paste("filter:", term)
        } else {
          markers <- rownames(pgx$gsetX)
        }
        gx <- pgx$gsetX[markers, rownames(pos), drop = FALSE]
      } else {
        cat("fatal error")
        return(NULL)
      }

      if (!"group" %in% names(pgx$model.parameters)) {
        stop("[markers.plotFUNC] FATAL: no group in model.parameters")
      }

      ## prioritize gene with large variance (groupwise)
      ## grp <- as.character(pgx$samples[rownames(pos),"group"])
      grp <- pgx$model.parameters$group[rownames(pos)]
      zx <- t(apply(gx, 1, function(x) tapply(x, grp, mean)))
      gx <- gx[order(-apply(zx, 1, sd)), , drop = FALSE]
      gx <- gx - min(gx, na.rm = TRUE) + 0.01 ## subtract background??
      rownames(gx) <- sub(".*:", "", rownames(gx))

      NP <- 25
      if (mrk_level == "geneset") NP <- 16
      top.gx <- head(gx, NP) ## match number of plot below!
      if (mrk_sortby == "name") {
        top.gx <- top.gx[order(rownames(top.gx)), , drop = FALSE]
      } else {
        top.gx <- top.gx[order(-rowMeans(top.gx)), , drop = FALSE]
      }
      top.gx <- pmax(top.gx, 0)
      ## top.gx <- tanh(top.gx/mean(top.gx))

      pd <- list(
          top.gx = top.gx,
          pos = pos,
          mrk_level = mrk_level,
          mrk_features = mrk_features
      )

      return(pd)
    })
    
    get_ggplots <- function() {

      pd <- plot_data()
      shiny::req(pd)
      top.gx <- pd$top.gx        
      pos <- pd$pos
      mrk_level <- pd$mrk_level

      ## make smaller dots when more points
      cex1 <- 1.0
      cex1 <- 0.95 * c(2.2, 1.1, 0.6, 0.3)[cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))]

      ## grey to red colorpalette for absolute expression
      klrpal <- colorRampPalette(c("grey90", "grey80", "grey70", "grey60", "red4", "red3"))(16)
      klrpal <- colorRampPalette(c("grey90", "grey60", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66")

      plt <- list()
      i <- 1
      for (i in 1:nrow(top.gx)) {
        colvar <- pmax(top.gx[i, ], 0)
        colvar <- 1 + round(15 * (colvar / (0.7 * max(colvar) + 0.3 * max(top.gx))))
        klr0 <- klrpal[colvar]
        
        if (mrk_level == "gene") {
          label <- sub(".*:", "", rownames(top.gx)[i])
        } else {
          gset <- sub(".*:", "", rownames(top.gx)[i])
          label <- playbase::breakstring(substring(gset, 1, 80), 24, force = TRUE)
          label <- tolower(label)
        }

        ## base::plot(pos[, ],
        ##   pch = 19, cex = cex1, col = klr0,
        ##   xlim = 1.1 * range(pos[, 1]), ylim = 1.1 * range(pos[, 2]),
        ##   fg = gray(0.8), bty = "o",
        ##   xaxt = "n", yaxt = "n", xlab = "tSNE1", ylab = "tSNE2"
        ## )
p
        tt <- rownames(top.gx)[i]
        
        ## ------- start plot ----------       
        p <- playbase::pgx.scatterPlotXY.GGPLOT(
          pos,
          var = colvar,
          col = klrpal,
          cex = 0.5*cex1,
          xlab = "",
          ylab = "",
          xlim = 1.2*range(pos[,1]),
          ylim = 1.2*range(pos[,2]),
          axis = FALSE,
          title = tt,
          cex.title = 0.95,
          ##title.y = 0.85,
          ##cex.clust = cex1*0.8,
          label.clusters = FALSE,
          legend = FALSE,
          gridcolor = "#ffffff",
          bgcolor = "#f8f8f8",          
          box = TRUE
        ) 

        plt[[i]] <- p
      }
      return(plt)
    }

    
    get_plotly_plots <- function() {

      pd <- plot_data()
      shiny::req(pd)
      top.gx <- pd$top.gx        
      pos <- pd$pos
      mrk_level <- pd$mrk_level

      ## make smaller dots when more points
      cex1 <- 1.0
      cex1 <- 0.6 * c(2.2, 1.1, 0.6, 0.3)[cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))]

      ## grey to red colorpalette for absolute expression
      klrpal <- colorRampPalette(c("grey90", "grey80", "grey70", "grey60", "red4", "red3"))(16)
      klrpal <- colorRampPalette(c("grey90", "grey60", "red3"))(16)
      klrpal <- paste0(gplots::col2hex(klrpal), "66")

      plt <- list()
      i <- 1
      for (i in 1:nrow(top.gx)) {
        colvar <- pmax(top.gx[i, ], 0)
        colvar <- 1 + round(15 * (colvar / (0.7 * max(colvar) + 0.3 * max(top.gx))))
        klr0 <- klrpal[colvar]
        
        if (mrk_level == "gene") {
          label <- sub(".*:", "", rownames(top.gx)[i])
        } else {
          gset <- sub(".*:", "", rownames(top.gx)[i])
          label <- breakstring(substring(gset, 1, 80), 24, force = TRUE)
          label <- tolower(label)
        }

        ## base::plot(pos[, ],
        ##   pch = 19, cex = cex1, col = klr0,
        ##   xlim = 1.1 * range(pos[, 1]), ylim = 1.1 * range(pos[, 2]),
        ##   fg = gray(0.8), bty = "o",
        ##   xaxt = "n", yaxt = "n", xlab = "tSNE1", ylab = "tSNE2"
        ## )
p
        tt <- rownames(top.gx)[i]
        
        ## ------- start plot ----------       
        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = colvar,
          col = klrpal,
          cex = 0.6*cex1,
          xlab = "",
          ylab = "",
          xlim = 1.2*range(pos[,1]),
          ylim = 1.2*range(pos[,2]),
          axis = FALSE,
          title = tt,
          cex.title = 0.85,
          title.y = 0.85,
#         cex.clust = cex1*0.8,
          label.clusters = FALSE,
          legend = FALSE,
          gridcolor = "fff",
          bgcolor = "#f8f8f8",
          tooltip = FALSE
        ) %>%
          plotly::style(
            hoverinfo = 'none'
          )

        plt[[i]] <- p
      }
      return(plt)
    }
    
    
    plotly.RENDER <- function() {
      pd <- plot_data()
      plt <- get_plotly_plots()
      shiny::req(plt)        
      nr  <- ceiling(sqrt(length(plt)))
      title <- pd$mrk_features
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = 0.01
      ) %>%
        plotly_default() %>%
        plotly::layout(
          title = list(text=title, size=12),
          margin = list(l=0,r=0,b=0,t=30) # lfbt            
        )
      
      return(fig)
    }

    plotly_modal.RENDER <- function() {
      fig <- plotly.RENDER() %>%
        plotly_modal_default() %>%
        plotly::layout(
          margin = list(l=0,r=0,b=0,t=50), # lfbt  
          title = list(size=18)
        ) 
      return(fig)
    }

    ggplot.RENDER <- function() {
      pd <- plot_data()  
      plt <- get_ggplots()
      shiny::req(plt)              
      nr  <- ceiling(sqrt(length(plt)))
      title <- pd$mrk_features
      fig <- gridExtra::grid.arrange(
        grobs = plt,
        nrow = nr,
        ncol = nr,
        padding = unit(0.01,"line"),
        top = textGrob(title,gp=gpar(fontsize=15))
      )
      return(fig)
    }
    
    PlotModuleServer(
      "plotmodule",
      func = ggplot.RENDER,
      #      func = plotly.RENDER,
      #      func2 = plotly_modal.RENDER,
      plotlib = "ggplot",
      res = c(85, 90),
      pdf.width = 10, pdf.height = 10,
      add.watermark = watermark
    )
  }) ## end of moduleServer
}
