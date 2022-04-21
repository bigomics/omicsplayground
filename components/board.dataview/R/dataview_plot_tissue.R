##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_plot_tissue_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)
    info_text = paste("Tissue expression for the selected gene in the tissue expression ",a_GTEx," dataset. Colors corresponds to 'tissue clusters' as computed by unsupervised clustering.")
    
    PlotModuleUI(
        ns("pltmod"),
        title = "Tissue expression (GTEX)",
        label = label,
        outputFunc = plotOutput,
        outputFunc2 = plotOutput,        
        info.text = info_text,
        caption = NULL,
        caption2 = NULL,        
        options = NULL,
        download.fmt=c("png","pdf","csv"),         
        width = c("auto","1200"),
        height = height
    )
    
}

dataview_plot_tissue_server <- function(id, pgxdata, parent.input, watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        plot_data  <- shiny::reactive({

            dbg("[dataview_tissueplot_server:plot_data] reacted!")
            
            ngs <- pgxdata()
            shiny::req(ngs)
            if(is.null(parent.input$data_type)) return(NULL)
            
            dbg("[dataview_tissueplot_server:plot_data] called..")
            
            gene <- parent.input$search_gene
            pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]
            hgnc.gene = toupper(as.character(ngs$genes[pp,"gene_name"]))

            tx = tissue.klr = NULL
            if( hgnc.gene %in% rownames(TISSUE)) {
                tx = TISSUE[hgnc.gene,]
                grp = TISSUE.grp[names(tx)]
                tissue.klr = COLORS[grp]
                ylab="expression (TPM)"
                if(parent.input$data_type=="logCPM") {
                    ylab = "expression (log2TPM)"
                    tx = log(1 + tx)
                }
                jj <- 1:length(tx)
                sorting="no"
                if(sorting=="decr") jj <- order(-tx)
                if(sorting=="inc") jj <- order(tx)

                tx <- tx[jj]
                tissue.klr <- tissue.klr[jj]
            }

            list(
                df = data.frame(
                    x = tx,
                    color = tissue.klr
                ),
                gene = hgnc.gene,
                ylab = ylab
            )
                

        })
        
        plot.RENDER <- function() {
            pdat <- plot_data()
            df   <- pdat$df
            ylab <- pdat$ylab
            gene <- pdat$gene
            if(is.null(df$x)) return(NULL)

            par(mar=c(6,4,1,1), mgp=c(1.5,0.5,0))
            barplot(df$x, las=3, main=gene, cex.main=1, col.main="#7f7f7f",
                    col = df$color, border=NA,ylab=ylab, cex.names=0.9,
                    names.arg=rep(NA,length(df$x)))

            text((1:length(df$x)-0.5)*1.2, -0.04*max(df$x), names(df$x), las=3,
                 cex=0.85, pos=2, adj=0, offset=0, srt=55, xpd=TRUE)

        }
        
        modal_plot.RENDER <- function() {
            plot.RENDER()
        }
        
        PlotModuleServer(
            "pltmod",
            plotlib = "base",
            plotlib2 = "base",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            ##renderFunc = shiny::renderPlot,
            ##renderFunc2 = shiny::renderPlot,        
            renderFunc = shiny::renderCachedPlot,
            renderFunc2 = shiny::renderCachedPlot,        
            res = c(96,120)*1,                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


