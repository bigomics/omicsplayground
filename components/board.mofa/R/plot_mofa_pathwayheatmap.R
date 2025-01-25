##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

mofa_plot_pathwayheatmap_ui <- function(
    id,
    title = "",
    info.text = "",
    caption = "",
    label = "",
    height = 400,
    width = 400) {
  ns <- shiny::NS(id)

  options <- shiny::tagList(
    shiny::checkboxInput(
      inputId = ns("split"),
      label = "Split heatmap by data types",
      value = TRUE
    )
  )

  PlotModuleUI(
    ns("plot"),
    title = title,
    label = label,
    info.text = info.text,
    options = options,
    caption = caption,
    height = height,
    width = width,
    download.fmt = c("png", "pdf")
  )
}

mofa_plot_pathwayheatmap_server <- function(id,
                                           mofa,
                                           input_factor = reactive(1),
                                           selected = reactive(1),
                                           watermark = FALSE) {
  moduleServer(id, function(input, output, session) {

    plot.RENDER <- function() {
      mofa <- mofa()      
      k <- input_factor()
      factors <- colnames(mofa$F)
      if(!is.null(k)) shiny::req(k %in% factors)

      pw <- selected()
      shiny::validate(need(length(pw)==1, "Please select a pathway"))

      dbg("[plot_pathwayheatmap_server] selected_pathway =", pw)

      features <- NULL
      ## if(!is.null(mofa$GMT)) {
      ##   if(pw %in% colnames(mofa$GMT)) {
      ##     features <- names(which(mofa$GMT[,pw] != 0))
      ##   }
      ## }
      
      playbase::mofa.plot_heatmap(
        mofa, k=k, ## main=k,
        ##features = features,
        pathway = pw,        
        ntop = 40,
        split = input$split,
        type = "splitmap",
        annot = "pheno",
        maxchar = 40,
        show_types = NULL,
        mar = c(3,0,0,0),
        annot.ht = 0.9,
        cexRow = 0.9)
    }

    PlotModuleServer(
      "plot",
      func = plot.RENDER,
      pdf.width = 8, pdf.height = 12,
      res = c(80, 100),
      add.watermark = watermark
    )

    
  })
}



