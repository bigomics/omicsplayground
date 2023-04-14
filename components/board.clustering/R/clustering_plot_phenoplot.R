##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


clustering_plot_phenoplot_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width) {
  ns <- shiny::NS(id)
  
  phenoplot.opts <- shiny::tagList(
    shiny::checkboxInput(ns("showlabels"), "Show group labels", TRUE)
  )

  PlotModuleUI(
    ns("pltmod"),
    title = title,
    label = label,
    #    plotlib = "base",
    plotlib = "plotly",
    info.text = info,
    caption = caption,
    options = phenoplot.opts,
    download.fmt = c("png", "pdf", "csv"),
    width = width,
    height = height
  )
}

clustering_plot_phenoplot_server <- function(id,
                                             pgx,
                                             selected_phenotypes,
                                             hm_getClusterPositions,
                                             watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    plot_data <- reactive({
      pgx <- pgx
      shiny::req(pgx$Y)

      ## get t-SNE positions
      clust <- hm_getClusterPositions()
      ## pos = pgx$tsne2d
      pos <- clust$pos
      colnames(pos) <- c("x","y")
      
      Y <- pgx$Y[rownames(pos), , drop = FALSE]
      pheno <- selected_phenotypes()

      # ## don't show these...
      # removed the code below because it was removing the batch and sample, overwritting user wishes
      # on selected_phenotypes
      # pheno <- grep("batch|sample|donor|repl|surv", pheno,
      #   invert = TRUE, ignore.case = TRUE, value = TRUE
      # )
      Y <- Y[,pheno] 
      
      ## complete dataframe for downloading
      df <- data.frame( pos, Y)
      
      return(
        list(
          df = df,
          pheno = pheno,
          showlabels = input$showlabels
        )
      )
    })
    
    render_plotly <- function(pd, pheno, cex=1) {

      pheno <- pd[["pheno"]]
      Y <- pd[["df"]][,pheno]
      showlabels <- pd[["showlabels"]]
      pos <- pd[["df"]][,c("x","y")]

      ## points size depending on how many points we have
      ncex <- cut(nrow(pos), breaks = c(-1, 40, 200, 1000, 1e10))
      cex1 <- 0.8 * cex * c(1.8, 1.3, 0.8, 0.5)[ncex]
      cex1 <- cex1 * ifelse(length(pheno) > 6, 0.8, 1)
      cex1 <- cex1 * ifelse(length(pheno) > 12, 0.8, 1)

      plt <- list()
      i=1
      for (i in 1:min(20, length(pheno))) {
        ## ------- set colors
        colvar <- factor(Y[, 1])
        colvar <- factor(Y[, pheno[i]])
        colvar[which(colvar %in% c(NA, "", " ", "NA", "na"))] <- NA
        colvar <- factor(as.character(colvar))
        klrpal <- COLORS
        klr1 <- klrpal[colvar]
        klr1 <- paste0(gplots::col2hex(klr1), "99")
        jj <- which(is.na(klr1))
        if (length(jj)) klr1[jj] <- "#AAAAAA22"
        tt <- tolower(pheno[i])

        ## ------- start plot        
        p <- playbase::pgx.scatterPlotXY.PLOTLY(
          pos,
          var = colvar,
          col = klrpal,
          cex = cex1,
          xlab = "",
          ylab = "",
          title = tt,
          cex.title = cex*1.1,
          cex.clust = cex*1.1,
          label.clusters = showlabels
        ) %>% plotly::layout(
          ## showlegend = TRUE,
          plot_bgcolor = "#f8f8f8"
        )
        
        plt[[i]] <- p
      } 
      return(plt)
    }

    plotly.RENDER <- function() {

      pd <- plot_data()
      pheno <- pd[["pheno"]]      
      plt <- render_plotly(pd, pheno, cex=0.85) 
      
      nr = min(3,length(plt))
      if (length(plt) >= 6)  nr = 4
      if (length(plt) >= 12)  nr = 5
      
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = 0.04
      ) %>% plotly_default()

      return(fig)
    }

    plotly_modal.RENDER <- function() {

      pd <- plot_data()
      pheno <- pd[["pheno"]]      
      plt <- render_plotly(pd, pheno, cex=1.3) 
      
      nc = min(3,length(plt))
      if (length(plt) >= 6)  nc = 4
      if (length(plt) >= 12)  nc = 5
      nr <- ceiling(length(plt)/nc)
      
      fig <- plotly::subplot(
        plt,
        nrows = nr,
        margin = 0.06
      ) %>% plotly_modal_default()
      
      return(fig)
    }
    
    PlotModuleServer(
      "pltmod",
      ##plotlib = "base",
      plotlib = "plotly",      
      func = plotly.RENDER,
      func2 = plotly_modal.RENDER,
      csvFunc = plot_data, ##  *** downloadable data as CSV
      res = c(85), ## resolution of plots
      pdf.width = 6, pdf.height = 9,
      add.watermark = watermark
    )
  })
}
