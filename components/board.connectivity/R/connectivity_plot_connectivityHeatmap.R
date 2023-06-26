##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Importance plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
connectivity_plot_connectivityHeatmap_ui <- function(
  id,
  title,
  info.text,
  caption,
  label = "",
  height,
  width) {
  ns <- shiny::NS(id)
  
  plot_opts <- shiny::tagList(
    withTooltip(shiny::radioButtons(ns("ngenes"), "Number of genes",
      choices=c(50,200,500), inline=TRUE),
      "Number of genes to show."
    ),
    hr(),
    withTooltip(shiny::radioButtons(ns("nsig"), "Number of signatures",
      choices=c(10,20,40,100), selected=20, inline=TRUE),
      "Number of nearest signatures to show."
    ),
    hr(),
    withTooltip(shiny::checkboxInput(ns("clusterx"), "Cluster genes"),
      "Cluster genes or sort by expression.."
      ),
    hr(),
    withTooltip(shiny::checkboxInput(ns("cumFCplot_absfc"), "Use absolute foldchange", FALSE),
      "Take the absolute foldchange for calculating the cumulative sum.",
      placement = "right", options = list(container = "body")
    ),
    hr(),    
    withTooltip(shiny::checkboxInput(ns("cumFCplot_absfc"), "Use absolute foldchange", FALSE),
      "Take the absolute foldchange for calculating the cumulative sum.",
      placement = "right", options = list(container = "body")
    )
  )
  PlotModuleUI(ns("plotmodule"),
    title = title,
    label = label,
    ## plotlib = "base",
    plotlib = "plotly",
    ##plotlib = "iheatmapr",
    info.text = info.text,
    options = plot_opts,
    download.fmt = c("pdf", "png", "csv"),    
    height = height,
    width = width,
    caption = caption
  )
}

#' Importance plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
connectivity_plot_connectivityHeatmap_server <- function(id,
                                                         getTopProfiles,
                                                         getConnectivityScores,
                                                         getCurrentContrast,
                                                         watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      
      plot_data <- shiny::reactive({
        F <- getTopProfiles()
        F[is.na(F)] <- 0

        ## get correlation
        df <- getConnectivityScores()
        rho1 <- df$rho[match(colnames(F), df$pathway)]
        
        ## add current contrast
        cc <- getCurrentContrast()
        shiny::req(cc)
        fc <- cc$fc[match(rownames(F),names(cc$fc))]
        names(fc) <- rownames(F)
        ## fc[is.na(fc)] <- 0
        F <- cbind(fc[rownames(F)], F)
        colnames(F)[1] <- "thisFC"
        colnames(F)[1] <- cc$name
        colnames(F)[1] <- paste("********",cc$name,"********")
        rho2 <- c(1, rho1)
        names(rho2) <- colnames(F)        
        if (input$cumFCplot_absfc) {
          F <- abs(F)
        }
        F <- F[order(-rowMeans(F**2,na.rm=TRUE)),, drop = FALSE]
        F <- scale(F, center=FALSE)
        
        list(
          F = F,
          score = rho2
        )
        
      })

      plot_heatmap <- function(F, maxfc, maxgenes=60) {
        F <- F[, 1:min(NCOL(F), maxfc), drop = FALSE]
        F1 <- head(F, maxgenes)
        par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
        playbase::gx.splitmap(t(F1),
          split = 1,
          ## cluster_columns = FALSE,
          cluster_columns = TRUE,
          cluster_rows = TRUE,
          rowlab.maxlen = 80,
          symm.scale = TRUE,
          mar = c(15, 0, 0, 110),
          key.offset = c(0.85, 0.15),
          cexRow = 0.9,
          cexCol = 0.75
        )
      }

      create_iheatmap <- function(F, score, maxfc=20, maxgenes=60) {

        sel <- 1:min(NCOL(F), maxfc)
        F <- F[,sel, drop = FALSE]
        score <- score[colnames(F)]          
        F <- head(F, maxgenes)        
        ii <- order(rowMeans(F,na.rm=TRUE))
        F <- F[ii,]

        plt <- iheatmapr::main_heatmap(
          data = t(F),
          layout = list(margin = list(r=0))) %>% 
          iheatmapr::add_row_clustering() %>% 
          iheatmapr::add_row_labels(size=0.5) %>%
          iheatmapr::add_col_labels() 

        if(input$clusterx) {
          plt <- plt %>%
            iheatmapr::add_col_clustering() 
        }

        ## add average logFC barplot on top
        avgF <- rowMeans(F,na.rm=TRUE)
        plt <- plt %>%
          iheatmapr::add_col_barplot(
                       y = avgF,
                       layout = list(title="avg logFC"),
                       buffer = 0.10,
                       size = 0.30
                     ) %>%
          iheatmapr::add_row_barplot(
                       x = score,
                       layout = list(title="similarity"),
                       size = 0.10
                     )

        plt <- plt %>% iheatmapr::to_plotly_list()
        plt <- plotly::as_widget(plt) %>%
          plotly::layout(
                    margin = list(l=0,r=0,t=0,b=40)
                  )
        plt
      }

      plot_RENDER <- function(){
        pd <- plot_data()
        F <- pd$F
        score <- pd$score
        shiny::req(F)
        ngenes <- as.numeric(input$ngenes)
        nsig <- as.numeric(input$nsig)
        create_iheatmap(F, score, maxfc=nsig, maxgenes=ngenes)
      }

      plot_RENDER2 <- function(){
        pd <- plot_data()
        F <- pd$F
        score <- pd$score
        shiny::req(F)
        ngenes <- as.numeric(input$ngenes) * 1.5
        nsig <- as.numeric(input$nsig) * 2
        create_iheatmap(F, score, maxfc=nsig, maxgenes=ngenes)
      }
      
      PlotModuleServer(
        "plotmodule",
        ##plotlib = "base",
        plotlib = "plotly",
        ##plotlib = "iheatmapr",
        func = plot_RENDER,
        func2 = plot_RENDER2,
        csvFunc = plot_data,
        pdf.width = 14, pdf.height = 5.5,
        res = c(90, 90),
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
