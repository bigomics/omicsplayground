##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##



## Annotate clusters ############

clustering_plot_hm_parcoord_ui <- function(id,
                                           label='',
                                            height,
                                            width
                                           )
{
  ns <- shiny::NS(id)

  info_text = "The <strong>Parallel Coordinates</strong> panel
displays the expression levels of selected genes across all conditions in the analysis. On the x-axis the experimental conditions are plotted. The y-axis shows the expression level of the genes grouped by condition. The colors correspond to the gene groups as defined by the hierarchical clustered heatmap."


  hm_parcoord_opts = shiny::tagList(
            withTooltip( shiny::checkboxInput(ns('hm_pcscale'),'Scale values',TRUE),
                   "Scale expression values to mean=0 and SD=1.",
                   placement="right",options = list(container = "body"))
        )

  PlotModuleUI(
    ns("pltmod"),
    title = "Parallel coordinates",
    label = label,
    plotlib = "plotly",
    info.text = info_text,
    options = hm_parcoord_opts,
    download.fmt=c("png","pdf","csv"),
    width = width,
    height = height
  )
}

clustering_plot_hm_parcoord_server <- function(id,
                                               pgx,
                                               hm_parcoord.matrix,
                                               watermark=FALSE
                                               )
{
  moduleServer( id, function(input, output, session) {
    ns <- session$ns

    hm_parcoord.RENDER <- shiny::reactive({

        pc <- hm_parcoord.matrix()
        shiny::req(pc)
        zx <- pc$mat
        ## build dimensions
        dimensions <- list()
        for(i in 1:ncol(zx)) {
            dimensions[[i]] <-  list(
                range = c(min(zx[,i]),max(zx[,i])),
                ## constraintrange = c(100000,150000),
                ## tickvals = c(0,0.5,1,2,3),
                ## ticktext = c('A','AB','B','Y','Z'),
                visible = TRUE,
                label = colnames(zx)[i],
                values = zx[,i]
            )
        }

        clust.id <- as.integer(factor(pc$clust))
        table(clust.id)

        df <- data.frame(clust.id=clust.id, zx)
        klrpal = rep(RColorBrewer::brewer.pal(8,"Set2"),99)
        ##klrpal = rep(c("red","blue","green","yellow","magenta","cyan","black","grey"),99)
        klrpal = klrpal[1:max(clust.id)]
        ##klrpal <- setNames(klrpal, sort(unique(clust.id)))
        klrpal2 <- lapply(1:length(klrpal),function(i) c((i-1)/(length(klrpal)-1),klrpal[i]))

        plt <-  plotly::plot_ly(df, source = "pcoords") %>%
            plotly::add_trace(type = 'parcoords',
                      line = list(color = ~clust.id,
                                  ## colorscale = list(c(0,'red'),c(0.5,'green'),c(1,'blue'))
                                  ##colorscale = 'Jet',
                                  colorscale = klrpal2,
                                  cmin = min(clust.id), cmax = max(clust.id),
                                  showscale = FALSE
                                  ##reversescale = TRUE
                                  ),
                      dimensions = dimensions)
        plt <- plt %>%
            plotly::layout(margin = list(l=60, r=60, t=0, b=30)) %>%
            ##config(displayModeBar = FALSE) %>%
            ##config(modeBarButtonsToRemove = setdiff(all.plotly.buttons,"toImage") ) %>%
            plotly::config(toImageButtonOptions = list(format='svg', width=900, height=350, scale=1.2)) %>%
            plotly::config(displaylogo = FALSE) %>%
            plotly::event_register("plotly_restyle")

        plt

    })


    PlotModuleServer(
      "pltmod",
      plotlib = "plotly",
      ##plotlib2 = "plotly",
      func = hm_parcoord.RENDER,
      ##renderFunc = plotly::renderPlotly,
      ##renderFunc2 = plotly::renderPlotly,
      res = c(90,170),                ## resolution of plots
      pdf.width = 8, pdf.height = 8,
      add.watermark = watermark
    )

    # shiny::callModule(
    #     plotModule,
    #     ## hm_parcoord_module <- plotModule(
    #     "hm_parcoord",
    #     func = hm_parcoord.RENDER, ## ns = ns,
    #     plotlib = "plotly", ## renderFunc="renderPlotly",
    #     ## download.fmt = c("png","pdf","html"),  ## PNG & PDF do not work!!!
    #     ## download.fmt = c("html"),
    #     options = hm_parcoord_opts,
    #     height = c(0.45*fullH,600), width = c("100%",1000),
    #     pdf.width=10, pdf.height=6, info.width="350px",
    #     title = "Parallel coordinates", label = "a",
    #     info.text = hm_parcoord_text,
    #     add.watermark = WATERMARK
    #     ## caption = hm_parcoord_text,
    # )
    })


}
