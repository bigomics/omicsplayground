##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' AUC plot UI input function
#'
#' @description A shiny Module for plotting (UI code).
#'
#' @param id
#' @param label
#' @param height
#'
#' @export
biomarker_plot_auc_ui <- function(
    id,
    title,
    info.text,
    info.methods,
    info.extra_link,
    caption,
    label = "",
    height,
    width) {
  ns <- shiny::NS(id)

  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    plotlib = "base",
    info.text = info.text,
    info.methods = info.methods,
    info.extra_link = info.extra_link,
    options = NULL,
    caption = caption,
    download.fmt = c("png", "pdf"),
    width = width,
    height = height
  )
}

#' AUC plot Server function
#'
#' @description A shiny Module for plotting (server code).
#'
#' @param id
#'
#' @return
#' @export
biomarker_plot_auc_server <- function(id,
                                      calcVariableImportance,
                                      is_computed,
                                      watermark = FALSE) {
  moduleServer(
    id, function(input, output, session) {
      plot_data <- shiny::reactive({
        res <- calcVariableImportance()
        shiny::req(res)
        shiny::req(is_computed())
        return(res)
      })

      plot.RENDER <- function() {
        res <- plot_data()

        shiny::req(res)
        par(mfrow = c(1, 1), mar = c(1, 0, 2, 0))
        rf <- partykit::as.party(res$rf)
        require(tidyr)
        require(pROC)
        bestpreds <- unique(res$rf$frame$var)
        bestpreds <- bestpreds[which(bestpreds != "<leaf>")]

        D0 <- rf$data
        y <- D0[, 1]
        D <- as.data.frame(t(D0[,-1]))
        genes <- rep(rownames(D), each=ncol(D))
        D <- tidyr::pivot_longer(
          D,
          cols = everything(), 
          names_to = "samples",
          values_to = "value"
        )
        D <- cbind(D, gene=genes)
        D <- D[which(D$gene %in% bestpreds), ] 
        rownames(D) <- 1:nrow(D)
        jj <- match(D$samples, rownames(D0))
        D$pheno <- D0$y          

        is.multinomial <- length(table(res$y)) > 2
        if (is.multinomial) {
          res.roc <- pROC::multiclass.roc(
            response = D$pheno,
            predictor = D$value,
            levels = levels(D$pheno)
          )
          ## legends <- cols <- c()
          ## i=1
          ## for(i in 1:length(res.roc$rocs)) {
          ##   pl <- res.roc$rocs[[i]]
          ##   legends <- c(legends, paste0(pl$levels, collapse=";"))
          ## if(i==1) { col="black"; auc.adj=c(0,3) }
          ## if(i==2) { col="red"; auc.adj=c(0,5) }
          ## if(i==3) { col="blue"; auc.adj=c(0,7) }
          ## cols <- c(cols, col)
          ## add <- ifelse(i>1, TRUE, FALSE)          
          pROC::plot.roc(res.roc$rocs[[1]], add = FALSE,
            main = paste0(bestpreds, collapse="; "),
            print.auc = TRUE, legacy.axes = TRUE, las = 1)
          pROC::plot.roc(res.roc$rocs[[2]], add = TRUE, col = "red",
            print.auc = TRUE, legacy.axes = TRUE, print.auc.adj = c(0,3))
          pROC::plot.roc(res.roc$rocs[[3]], add = TRUE, col = "blue",
            print.auc = TRUE, legacy.axes = TRUE, print.auc.adj = c(0,5))
          cols <- c("black","red","blue")
          legend("bottomright", legend = c("A","B","C"), col = cols, lwd = 2)
        } else {
          res.roc <- pROC::roc(D, response = "pheno", predictor = "value")
          pROC::plot.roc(res.roc, legacy.axes = TRUE,
            print.auc = TRUE, las = 1,
            main = paste0(bestpreds, collapse="; "))
        }
      }
      
      PlotModuleServer(
        "plot",
        plotlib = "base", # does not use plotly
        func = plot.RENDER,
        func2 = plot.RENDER, # no separate modal plot render
        res = c(60, 100),
        pdf.width = 10, pdf.height = 6,
        add.watermark = watermark
      )
    } ## end of moduleServer
  )
}
