TimeSeriesBoard.enrichment_table_ui <- function(
    id,
    label = "label",
    title = "title",
    info.text = "info.text",
    caption = "caption",
    height = c("40%", TABLE_HEIGHT_MODAL),
    width = c("auto", "100%")
    ) {
  ns <- shiny::NS(id)


  options <- tagList(
    withTooltip(
      shiny::checkboxInput(ns("show_statdetails"), "Show detailed statistical methods"),
      title = "Show detailed statistical methods."
    )
  )
  
  TableModuleUI(
    ns("table"),
    info.text = info.text,
    options = options,
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
    caption = "caption",
    info.text = "info.text",
    info.methods = "info.methods",
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
      
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] names(res) =", names(res))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] dim(gset.rho) =", dim(gset.rho))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] rownames.gset.rho =", head(rownames(gset.rho)))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] colnames.gset.rho =", head(colnames(gset.rho)))
      
      rho <- res$gset.rho[,k]
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] 0: names.rho =", head(names(rho)))

      
      names(rho) <- rownames(gset.rho)
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] 1: names.rho =", head(names(rho)))

      rho <- sort(rho, decreasing=TRUE)
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] 2: names.rho =", head(names(rho)))
      
#      gset.filter = "HALLMARK"
#      gset.filter = "GOBP|GO_BP"
#      rho <- rho[grepl(gset.filter,names(rho))]
      
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] length(rho) =", length(rho))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] names.rho =", head(names(rho)))
      
      pv <- playbase::cor.pvalue(rho, ncol(pgx$X))

      dbg("[TimeSeriesBoard.enrichment_server:gset_data] length(pv) =", length(pv))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] head(rho) =", head(rho))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] length.names.rho =", length(names(rho)))
      dbg("[TimeSeriesBoard.enrichment_server:gset_data] 3: names.rho =", head(names(rho)))
      
      df <- data.frame( geneset = names(rho), rho = rho, p.value=pv )

      dbg("[TimeSeriesBoard.enrichment_server:gset_data] dim(df) =", dim(df))
      
      return(df)
    }
    
    ##-----------------------------------------------------
    ##----------------------- Plot -----------------------
    ##-----------------------------------------------------
    
    render_plot <- function() {
      library(ggplot2)
      library(plotly)
      df <- gset_data()
      shiny::req(df)
      dbg("[TimeSeriesBoard.enrichment_server:render_plot] dim(df) =", dim(df))
      
      sel <- table_module$rows_all()      
      dbg("[TimeSeriesBoard.enrichment_server:render_plot] sel =", head(sel))

      sel <- head(sel,24)
      top <- df[sel,]
      values <- top$rho
      names(values) <- top$geneset
      ## sizes <- rowSums(pgx$GMT[,names(values)]!=0)
      playbase::ggLollipopPlot(values, sizes = NULL, xlab = "value",
        cex.text = 0.9) 
    }

    PlotModuleServer(
      "plot",
      func = render_plot,
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
      sum.na <- sum(rownames(df)=="" | is.na(rownames(df)))
      dbg("[TimeSeriesBoard.enrichment_server:render_table] sum.na =", sum.na)
      
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
                "return type === 'display' && data.length > 80 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 80) + '...</span>' : data;",
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
