
dataviewtSNEplotModuleUI3 <- function(id, height=c(600,800)) {

    ns <- shiny::NS(id)    
    options <- tagList(
        actionButton(ns("redraw"),"redraw")
    )
    
    PlotModuleUI(
        NS(id,"pltsrv"),  ## nested NS
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,        
        info.text = "Information: This is a T-SNE figure.",
        caption = "T-SNE figure caption.",        
        options = options,
        download.fmt=c("png","pdf","csv"),         
        width = c("auto","1200"),
        height = height,
        label = "A",
        title = "t-SNE clustering"
    )
    
}

dataviewtSNEplotModuleServer3 <- function(id, filterStates, data, label="", watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {
        
        message("[dataviewtSNEModuleServer3] ***called***")

        plot_dl <- reactiveValues()
        
        plot_data <- shiny::reactive({
            
            if(class(filterStates)[1]=="reactiveExpr") {
                filterStates <- filterStates()
            }
            
            ##input$redraw
            
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
            colnames(data) <- c("pos_x", "pos_y")
            data$fc2 <- fc2
            
            grp <- NULL
            
            if(!is.null(filterStates$data_groupby) && filterStates$data_groupby != "<ungrouped>") {
                grp <- factor(ngs$samples[samples, filterStates$data_groupby])
            }
            data$grp <- grp
            
            return(data)
        })

        ##      output$plot <-  renderCachedPlot({
        ##plot.RENDER <- shiny::reactive({
        plot.RENDER <- function() {        
            
            message("[dataviewtSNEModuleServer3] plot.RENDER called")

            data <- plot_data()
            fig_base <- 
                ggplot(data, aes(pos_x, pos_y)) +
                labs(x = "tSNE1", y = "tSNE2") +
                scale_color_viridis_c(
                    option = "rocket", 
                    direction = -1, 
                    begin = .05, end = .97,
                    limits = c(0, 1),
                    labels = function(x) sprintf("%1.2f", x),
                    name = "Expression"
                ) +
                guide_continuous(aes = "color", type = "steps", width = .4) +
                theme_omics(base_size = 12, axis_num = "xy", legendnum = TRUE)
            
            plot_dl$base <- fig_base
            
            if (!is.null(plot_data()$grp)) {
                fig <- fig_base +
                    ggforce::geom_mark_hull(
                                 aes(fill = stage(grp, after_scale = colorspace::desaturate(fill, 1)), label = grp),
                                 color = "grey33", 
                                 size = .4,
                                 alpha = .33 / length(unique(plot_data()$grp)),
                                 expand = unit(2.7, "mm"), 
                                 con.cap = unit(.01, "mm"), 
                                 con.colour = "grey33", 
                                 label.buffer = unit(2, "mm"),
                                 label.fontsize = 12.5, 
                                 label.fontface = "plain"
                             ) +
                    geom_point(
                        aes(color = stage(fc2, after_scale = colorspace::darken(color, .35)), 
                            fill = after_scale(color)), 
                        size = 1.8, 
                        shape = 21, 
                        stroke = .5
                    ) +
                    scale_x_continuous(expand = c(.4, .4)) +
                    scale_y_continuous(expand = c(.4, .4)) +
                    scale_fill_discrete(guide = "none")
                
            } else {
                fig <- fig_base +
                    geom_point(
                        aes(color = stage(fc2, after_scale = colorspace::darken(color, .35)), 
                            fill = after_scale(color)), 
                        size = 2.3, 
                        shape = 21, 
                        stroke = .5
                    )
            }
            
            ##Sys.sleep(5)
            
            plot_dl$plot <- fig
            fig
            ##gridExtra::grid.arrange(fig)
            ##cowplot::as_grob(fig)
            ##}, res = 96, cacheKeyExpr = { list(plot_data()) },)
        }
        
        ##modal_plot.RENDER <- function() {}
        ##output$modal_plot <- renderCachedPlot({
        ##modal_plot.RENDER <- shiny::reactive({
        modal_plot.RENDER <- function() {
            
            message("[dataviewtSNEplotModuleServer3] modal_plot.RENDER called")
            message("[tsne3::modal_plot.RENDER] names(plot_dl) = ",names(plot_dl))
            message("[tsne3::modal_plot.RENDER] length(plot_dl) = ",length(plot_dl))
            message("[tsne3::modal_plot.RENDER] class(plot_dl$base) = ",class(plot_dl$base))
            if(is.null(class(plot_dl$base))) {
                message("[tsne3::modal_plot.RENDER] call plot.RENDER ")
                tmp <- plot.RENDER()
            }
            
            fig_base <- plot_dl$base +
                ##fig_base <- plot_dl$plot +         
                guide_continuous(aes = "color", type = "steps", width = .7) +
                theme_omics(base_size = 20, axis_num = "xy", legendnum = TRUE)
            
            if (!is.null(plot_data()$grp)) {
                fig <- fig_base #+
                                        # ggforce::geom_mark_hull(
                                        #   aes(fill = stage(grp, after_scale = colorspace::desaturate(fill, 1)), label = grp),
                                        #   color = "grey33", 
                                        #   size = .8,
                                        #   alpha = .33 / length(unique(plot_data()$grp)),
                                        #   expand = unit(3.4, "mm"), 
                                        #   con.cap = unit(.01, "mm"), 
                                        #   con.colour = "grey33", 
                                        #   label.buffer = unit(2, "mm"),
                                        #   label.fontsize = 22, 
                                        #   label.fontface = "plain"
                                        # ) +
                                        # geom_point(
                                        #   aes(color = stage(fc2, after_scale = colorspace::darken(color, .35)), 
                                        #       fill = after_scale(color)), 
                                        #   size = 1.8, 
                                        #   shape = 21, 
                                        #   stroke = .5
                                        # ) +
                                        # scale_x_continuous(expand = c(.15, .15)) +
                                        # scale_y_continuous(expand = c(.15, .15)) +
                                        # scale_fill_discrete(guide = "none")
                
            } else {
                fig <- fig_base #+
                                        # geom_point(
                                        #   aes(color = stage(fc2, after_scale = colorspace::darken(color, .35)), 
                                        #       fill = after_scale(color)), 
                                        #   size = 4.7, 
                                        #   shape = 21, 
                                        #   stroke = .5
                                        # )
            }
            fig
            ##gridExtra::grid.arrange(fig)
        }
        ##})        
        ##}, res = 96, cacheKeyExpr = { list(plot_data()) },)

        message("[dataviewtSNEModuleServer3] creating PlotModuleServer... ")          
        
        PlotModuleServer(
            "pltsrv",
            plotlib = "ggplot",
            plotlib2 = "ggplot",
            ##plotlib = "base",            
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            ##func2 = plot.RENDER,            
            csvFunc = plot_data,             ##  *** downloadable data as CSV
            ##renderFunc = shiny::renderPlot,
            ##renderFunc = shiny::renderCachedPlot,
            ##renderFunc2 = shiny::renderPlot,        
            ##renderFunc2 = shiny::renderCachedPlot,        
            res = c(96,120),                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            ##label = label, title = "t-SNE clustering",
            add.watermark = watermark
        )
        
        ## output$downloadPlotPNG <- downloadHandler(
        ##     filename = function(){paste("bigomics-tSNE",'.png',sep='')},
        ##     content = function(file){
        ##      ggsave(file, plot = plot_dl$plot)
        ##    }
        ##  )
        ## output$downloadPlotPDF <- downloadHandler(
        ##   filename = function(){paste("bigomics-tSNE",'.pdf',sep='')},
        ##   content = function(file){
        ##     pdf(file, width = 6, height = 6)
        ##     print(plot_dl$plot)
        ##     dev.off()
        ##   }
        ## )
        ##  })

    })  ## end of moduleServer
}
