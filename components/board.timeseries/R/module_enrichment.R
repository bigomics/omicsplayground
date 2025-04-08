TimeSeriesBoard.enrichment_table_ui <- function(
    id,
    label = "label",
    title = "title",
    info.text = "Table of enriched genesets. Table reporting the genesets along with the correlation value (rho) and the p.value. ",
    caption = "Table of enriched genesets. Table reporting the genesets along with the correlation value (rho) and the p.value. ",
    height = c("40%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)
  
  ## options <- tagList(
  ##   withTooltip(
  ##     shiny::checkboxInput(
  ##      ns("show_statdetails"),
  ##      "Show detailed statistical methods"),
  ##      title = "Show detailed statistical methods."
  ##   )
  ## )
  
  TableModuleUI(
    ns("table"),
    info.text = info.text,
    #options = options,
    height = height,
    caption = caption,
    width = width,
    title = title,
    label = "b"
  )
}


TimeSeriesBoard.enrichment_lolliplot_ui <- function(
    id,
    label = "label",
    title = "title",
    caption = "Top-correlated gene sets with the selected gene module. X-axis shows the Pearson's correlation coefficient. Correlation coefficient values and p-values can be retrieved from the 'Enriched gene sets table' below'.",
    info.text = "Top-correlated gene sets with the selected gene module. X-axis shows the Pearson's correlation coefficient. Correlation coefficient values and p-values can be retrieved from the 'Enriched gene sets table' below'.",
    info.methods = "First, the log2-transformed and normalized data matrix is feature-scaled. Feature scaling is performed by subtracting the feature average value from each feature value and then dividing by the feature standard deviation. For each sample, the module-wise average expression is calculated. Pearson's correlation between gene set average expression and gene module average expression is calculated and mostly significantly correlated gene sets can be identified.",
    info.references = list(),
    info.extra_link = "extra.link",
    height = c("calc(100vh - 310px)", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)
    
  PlotModuleUI(ns("plot"),
    title = title,
    label = label,
    ## plotlib = "plotly",
    info.text = info.text,
    info.methods = info.methods,
    info.references = info.references,
    info.extra_link = info.extra_link,
    caption = caption,
    download.fmt = c("png", "pdf", "csv", "svg"),
    width = width,
    height = height
  )
}


TimeSeriesBoard.enrichment_server <- function(id,
                                              pgx,
                                              data,
                                              select_module,
                                              watermark = FALSE) {
  
  moduleServer(id, function(input, output, session) {

    ##-----------------------------------------------------
    ##----------------------- Data -----------------------
    ##-----------------------------------------------------

    gset_data <- function() {
      k <- select_module()
      shiny::req(k)
      res <- data()
      gset.rho <- res$gset.rho
      if(is.null(gset.rho)) {
        df <- data.frame( geneset = "", rho = 1, p.value=1 )[0,]        
        return(df)
      }
      
      rho <- res$gset.rho[,k]
      names(rho) <- rownames(gset.rho)

      rho <- sort(rho, decreasing=TRUE)      
      pv <- playbase::cor.pvalue(rho, ncol(pgx$X))      
      df <- data.frame( geneset = names(rho), rho = rho, p.value=pv )
      return(df)
    }
    
    ##-----------------------------------------------------
    ##----------------------- Plot -----------------------
    ##-----------------------------------------------------
    
    render_plot <- function(nshort=60, ntop=12, cex=0.85) {
      library(ggplot2)
      library(plotly)
      df <- gset_data()
      shiny::validate(shiny::need(nrow(df)>0, "no genesets"))
      
      sel <- table_module$rows_all()      
      shiny::req(sel)
      
      order1 <- order(df$rho[sel])
      sel <- sel[c(head(order1,ntop), tail(order1,ntop))]
      top <- df[sel,]
      values <- top$rho
      names(values) <- playbase::shortstring(top$geneset,nshort)
      ## sizes <- rowSums(pgx$GMT[,names(values)]!=0)
      playbase::ggLollipopPlot(values, sizes = NULL, xlab = "value",
        cex.text = cex) 
    }

    render_plot.modal <- function() {
      render_plot(nshort=120, ntop=18, cex=1)
    }
    
    PlotModuleServer(
      "plot",
      func = render_plot,
      func2 = render_plot,      
      plotlib = "base",
      ##csvFunc = gset_data, ##  *** downloadable data as CSV
      res = c(72, 110), ## resolution of plots
      pdf.width = 14,
      pdf.height = 3.5,
      add.watermark = watermark
    )

    ##-----------------------------------------------------
    ##----------------------- Table -----------------------
    ##-----------------------------------------------------

    render_table <- function() {

      df <- gset_data()
      shiny::req(df)

      dbg("[TimeSeriesBoard.enrichment_server:render_table] dim(df) =", dim(df))
      shiny::validate(shiny::need(nrow(df)>0, "no genesets"))      
      
      numeric.cols <- c("rho","p.value")
      DT::datatable(
        df,
        rownames = FALSE,
        extensions = c("Buttons", "Scroller"),
        plugins = "scrollResize",
        selection = list(mode = "single", target = "row", selected = NULL),
        class = "compact hover",
        fillContainer = TRUE,
        options = list(
          dom = "lfrtip",
          scrollX = TRUE,
          scrollY = "23vh",
          scrollResize = TRUE,
          scroller = TRUE,
          deferRender = TRUE,
          columnDefs = list(
            list(
              targets = "geneset", ## with no rownames column 1 is column 2
              render = DT::JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 60 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 60) + '...</span>' : data;",
                "}"
              )
            )
          )
        ) ## end of options.list
      ) %>%
        DT::formatSignif(numeric.cols, 3) %>%
        DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
    }

    table_module <- TableModuleServer(
      "table",
      func = render_table,
      selector = "none"
    )
    
  }) ## end of moduleServer
}
