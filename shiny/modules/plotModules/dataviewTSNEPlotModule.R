dataviewtSNEModuleUI <- function(id) {
  ns <- NS(id)
  
    tagList(
      shiny::fillRow(
        flex = c(NA,NA,NA,NA,NA,1), height = "10%",
        shiny::HTML(paste0("<span class='module-label'>d</span>")),
        shinyWidgets::dropdownButton(
            downloadButton(ns('downloadPlotPNG'),'PNG'),
            downloadButton(ns('downloadPlotPDF'),'PDF'),
            circle = TRUE, size = "xs", ## status = "danger",
            icon = shiny::icon("download"), width = "40px", right=FALSE,
            tooltip = shinyWidgets::tooltipOptions(title = "Download", placement = "right")
        ),
        shinyWidgets::dropdownButton(
          shiny::tags$p(shiny::HTML("Figure")),
          shiny::br(),
          circle = TRUE, size = "xs", inline = TRUE,
          icon = shiny::icon("info"), width = "300px",
          inputId = ns("info"), right=FALSE,
          tooltip = shinyWidgets::tooltipOptions(title = "Info", placement = "right")
        ), 
        shiny::actionButton(inputId=ns("zoombutton"),label=NULL,
                            icon=icon("window-maximize"),
                            class="btn-circle-xs"),
        shiny::HTML(paste("<center>t-SNE clustering</center>")),
      ),
      plotOutput(ns("plot")),
      shiny::div(class="popup-plot",
                shinyBS::bsModal(ns("plotPopup"), title, size="l",
                        ns("zoombutton"),
                        plotOutput(ns("modal_plot"))
                        )
      )
    )
}

dataviewtSNEModuleServer <- function(id, filterStates, data) {
  moduleServer(
    id,
    function(input, output, session) {
      
      plot_dl <- reactiveValues()

      plot_data <- shiny::reactive({
        
        ngs <- data
        
        gene <- ngs$genes$gene_name[1]
        
        if(!is.null(filterStates$search_gene) && filterStates$search_gene!="") gene <- filterStates$search_gene
        samples <- colnames(ngs$X)
        if(!is.null(filterStates$data_samplefilter)) {
          samples <- selectSamplesFromSelectedLevels(ngs$Y, filterStates$data_samplefilter)
        }
        nsamples = length(samples)
        
        ## precompute
        pp <- rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]
        
        gx <- NULL
        ylab <- NULL
        
        if(filterStates$data_type == "counts") {
          gx <- ngs$counts[pp,samples]
          ylab <- "expression (counts)"
        } else if(filterStates$data_type == "CPM") {
          gx <- 2**ngs$X[pp,samples]
          ylab <- "expression (CPM)"
        } else if(filterStates$data_type == "logCPM") {
          gx <- ngs$X[pp,samples]
          ylab <- "expression (log2CPM)"
        }
        
        pos <- ngs$tsne2d[samples,]
        
        fc1 <- tanh(0.99 * scale(gx)[,1])
        fc1 <- tanh(0.99 * scale(gx, center = FALSE)[,1])
        ##fc1 <- tanh(0.99 * gx/sd(gx))
        fc2 <- (fc1 - min(fc1))
        
        jj2 <- order(abs(fc1))
        ## determine how to do grouping for group labels
        
        data <- data.frame(pos[jj2,])
        data$fc2 <- fc2
        
        grp <- NULL
        
        if(!is.null(filterStates$data_groupby) && filterStates$data_groupby != "<ungrouped>") {
          grp <- factor(ngs$samples[samples, filterStates$data_groupby])
          data$grp <- grp
        }
        
        return(data)
      })

      output$plot <-  renderCachedPlot({
        
        fig <- ggplot(plot_data(), aes(tSNE.x, tSNE.y)) +
               labs(x = "tSNE1", y = "tSNE2") +
               scale_color_continuous(name = "Expression") +
               guides(color = guide_colorbar(barwidth = unit(.7, "lines"))) +
               theme_omics(base_size = 10, axis_num = "xy", legendnum = TRUE)
        
        if (!is.null(plot_data()$grp)) {
          fig <- fig +
            ggforce::geom_mark_hull(
              aes(fill = stage(grp, after_scale = colorspace::desaturate(fill, 1)), label = grp),
              color = "grey33", 
              size = .8,
              alpha = .33 / length(unique(data$grp)),
              expand = unit(3.4, "mm"), 
              con.cap = unit(.01, "mm"), 
              con.colour = "grey33", 
              label.buffer = unit(3, "mm"),
              label.fontsize = 22, 
              label.fontface = "plain"
            ) +
            geom_point(aes(color = fc2), size = 3.5) +
            scale_x_continuous(expand = c(.15, .15)) +
            scale_y_continuous(expand = c(.15, .15)) +
            scale_fill_discrete(guide = "none")
          
        } else {
          fig <- fig +
            geom_point(aes(color = fc2), size = 4.5)
        }
        
        gridExtra::grid.arrange(fig)
        
        plot_dl$plot <- fig

      }, res = 96, cacheKeyExpr = { list(plot_data()) },)


      output$modal_plot <- renderCachedPlot({
        
        fig <-
        ggplot(plot_data(), aes(tSNE.x, tSNE.y)) +
        labs(x = "tSNE1", y = "tSNE2") +
        scale_color_continuous(name = "Expression") +
        guides(color = guide_colorbar(barwidth = unit(.7, "lines"))) +
        theme_omics(base_size = 20, axis_num = "xy", legendnum = TRUE)

      if (!is.null(plot_data()$grp)) {
        fig <- fig +
          ggforce::geom_mark_hull(
            aes(fill = stage(grp, after_scale = colorspace::desaturate(fill, 1)), label = grp),
            color = "grey33", 
            size = .8,
            alpha = .33 / length(unique(data$grp)),
            expand = unit(3.4, "mm"), 
            con.cap = unit(.01, "mm"), 
            con.colour = "grey33", 
            label.buffer = unit(3, "mm"),
            label.fontsize = 22, 
            label.fontface = "plain"
          ) +
          geom_point(aes(color = fc2), size = 3.5) +
          scale_x_continuous(expand = c(.15, .15)) +
          scale_y_continuous(expand = c(.15, .15)) +
          scale_fill_discrete(guide = "none")
        
      } else {
        fig <- fig +
          geom_point(aes(color = fc2), size = 4.5)
      }

      gridExtra::grid.arrange(fig)

      }, res = 96, cacheKeyExpr = { list(plot_data()) },)
     
     output$downloadPlotPNG <- downloadHandler(
         filename = function(){paste("bigomics-tSNE",'.png',sep='')},
         content = function(file){
          ggsave(file, plot = plot_dl$plot)
        }
      )
     output$downloadPlotPDF <- downloadHandler(
       filename = function(){paste("bigomics-tSNE",'.pdf',sep='')},
       content = function(file){
         pdf(file, width = 6, height = 6)
         print(plot_dl$plot)
         dev.off()
       }
     ) 
    }
  )
}