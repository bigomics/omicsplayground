##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_lasagnaSP_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("plot"),
    ## plotlib = "plotly",
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_lasagnaSP_frequency_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("freq_plot"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_lasagnaSP_scores_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  PlotModuleUI(
    ns("scores"),
    title = title,
    label = label,
    info.text = info.text,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_table_lasagnaSP_ui <- function(
    id,
    label = "",
    title = "title",
    info.text = "info",
    caption = "caption",
    width = 400,
    height = 400
    ) {
  ns <- shiny::NS(id)

  options = tagList(
    checkboxInput(ns("showpath"),"Show path details",FALSE)
  )
  
  TableModuleUI(
    ns("table"),
    title = title,
    options = options,
    info.text = info.text,
    caption = caption,
    width = width,
    height = height,
    label = label
  )
}


mofa_lasagnaSP_server <- function(id,
                                  data,
                                  input_contrast = reactive(NULL),      
                                  watermark = FALSE) {
  moduleServer(id, function(input, output, session) {
    
    sp_data <- reactive({
      data <- data()
      pheno <- input_contrast()
      shiny::req(pheno)

      dbg("[mofa_lasagnaSP_server] pheno = ", pheno)      
      dbg("[mofa_lasagnaSP_server] colnames.Y = ", colnames(data$Y))
      dbg("[mofa_lasagnaSP_server] pheno in colnames.Y = ", pheno %in% colnames(data$Y))
      
      shiny::validate( shiny::need( pheno %in% colnames(data$Y), "Invalid PHENO"))

      px <- paste0("PHENO:",pheno)
      sp <- playbase::lasagna.solve_SP(data, pheno, vtop=200) 
      sp.score <- sp$score
      sp.score[which(rownames(sp)==px)] <- 9999  ## match at 1
      sp <- sp[order(-sp.score),]
      head(sp, 1000)
    })
    
    plot.RENDER <- function() {
      sp <- sp_data()
      shiny::req(sp)
      
      sel.row <- rownames(sp)[1]
      sel.row <- sp_table$rownames_selected()
      if(length(sel.row)==0 || sel.row[1]=='') sel.row <- NULL

      playbase::lasagna.plot_SP(
        sp,
        ntop = 1000,
        hilight = sel.row,
        labcex = 0.8,
        colorby = "pheno",
        plotlib = "ggplot") 

    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      ##plotlib = "plotly",
      pdf.width = 8, pdf.height = 8,
      res = c(85, 125),
      add.watermark = watermark
    )

    
    ##----------------------------------------------
    ##----------------- scores plot ----------------
    ##----------------------------------------------

    scoreplot.RENDER <- function(ntop=18) {
      sp <- sp_data()
      shiny::req(sp)
      par(mar=c(4,4,1,0.5))
      plot( sort(sp$score, decreasing=FALSE),
        ylab="path score", xlab="ordered paths")
      abline(h=0, lty=2)
    }
    
    PlotModuleServer(
      "scores",
      func = scoreplot.RENDER,
      ##plotlib = "plotly",
      pdf.width = 12, pdf.height = 8,
      res = c(70, 110),
      add.watermark = watermark
    )

    ##----------------------------------------------
    ##----------------- freq plot ------------------
    ##----------------------------------------------

    freqplot.RENDER <- function(ntop=20) {
      sp <- sp_data()
      shiny::req(sp)
      
      pp <- table(unlist(strsplit(sp$path, split="->")))
      dt <- sub(":.*","",names(pp))
      table(dt)
      ## take top 100 of each type
      nn <- tapply(pp, dt, function(a) head(sort(a,decreasing=TRUE),100),
        simplify=FALSE)
      names(nn) <- NULL
      
      nn$PHENO <- NULL
      nn$SOURCE <- NULL
      nn$SINK <- NULL
      
      names(nn) <- NULL
      nn <- sort(unlist(nn))
      nn <- nn * sign(sp[names(nn),"score"])
      nn <- nn[!is.na(nn)]
      nn.top <- head(nn[order(-abs(nn))], ntop)
      
      par(mfrow=c(1,2),mar=c(4,4,0,0.5))
      plot.new()
      barplot( sort(nn.top), horiz=TRUE, las=1, xlab="frequency", main="")
    }

    freqplot.RENDER2 <- function() {
      freqplot.RENDER(ntop=35) 
    }
    
    PlotModuleServer(
      "freq_plot",
      func = freqplot.RENDER,
      func2 = freqplot.RENDER2,      
      ##plotlib = "plotly",
      pdf.width = 12, pdf.height = 8,
      res = c(70, 100),
      add.watermark = watermark
    )



    ##----------------------------------------------
    ##----------------- TABLE ----------------------
    ##----------------------------------------------
    
    table.RENDER <- function() {
      df <- sp_data()
      df <- cbind( feature=rownames(df), df)
      df$feature <- stringr::str_trunc(df$feature, 100)
      if(!input$showpath) {
        df$path <- NULL
      }
      
      numeric.cols <- grep("score|rho", colnames(df))

      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        selection = list(mode = "single", target = "row", selected = 1),
        class = "compact cell-border stripe hover",
        plugins = "scrollResize",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "70vh",
          scroller = TRUE,
          scrollResize = TRUE,
          deferRender = TRUE
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }
    
    sp_table <- TableModuleServer(
      "table",
      func = table.RENDER,
      selector = "single"
    )
    
    
  })
}
