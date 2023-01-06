##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_plot_tissue_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)
    info_text = paste("Tissue expression for the selected gene in the tissue expression ",
                      a_GTEx," dataset. Colors corresponds to 'tissue clusters' as computed by unsupervised clustering.")

    PlotModuleUI(
        ns("pltmod"),
        title = "Tissue expression (GTEX)",
        label = label,
        plotlib = "plotly",
        #outputFunc = plotOutput,
        #outputFunc2 = plotOutput,
        info.text = info_text,
        options = NULL,
        download.fmt=c("png","pdf","csv"),
        width = c("auto","100%"),
        height = height
    )

}

dataview_plot_tissue_server <- function(id, pgx, r.gene, r.data_type, watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {

        plot_data  <- shiny::reactive({

            shiny::req(pgx$X)
            shiny::req(r.gene(), r.data_type())

            ## dereference reactive
            gene <- r.gene()
            data_type <- r.data_type()

            pp <- rownames(pgx$genes)[match(gene,pgx$genes$gene_name)]
            hgnc.gene = toupper(as.character(pgx$genes[pp,"gene_name"]))

            tx = tissue.klr = grp = NULL
            if( hgnc.gene %in% rownames(TISSUE)) {
                tx = TISSUE[hgnc.gene,]
                grp = TISSUE.grp[names(tx)]
                tissue.klr = COLORS[grp]
                ylab="expression (TPM)"
                if(data_type=="logCPM") {
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
                    tissue = names(tx),
                    x = tx,
                    group = grp,
                    color = tissue.klr
                ),
                gene = hgnc.gene,
                ylab = ylab
            )


        })

        plot.RENDER.base <- function() {
            pdat <- plot_data()
            shiny::req(pdat)

            df   <- pdat$df
            ylab <- pdat$ylab
            gene <- pdat$gene

            par(mar=c(6,4,1,0), mgp=c(1.5,0.5,0))
            barplot(df$x, las=3, main=gene, cex.main=1, col.main="#7f7f7f",
                    col = df$color, border=NA,ylab=ylab, cex.names=0.9,
                    names.arg=rep(NA,length(df$x)))

            text((1:length(df$x)-0.5)*1.2, -0.04*max(df$x), df$tissue, las=3,
                 cex=0.85, pos=2, adj=0, offset=0, srt=55, xpd=TRUE)

        }

        plot.RENDER <- function() {
            pdat <- plot_data()
            shiny::req(pdat)

            df   <- pdat$df
            ylab <- stringr::str_to_sentence(pdat$ylab)
            gene <- pdat$gene

            ## plot as regular bar plot
            df <- dplyr::mutate(df, tissue = forcats::fct_reorder(stringr::str_to_title(paste(tissue, " ")), x))

            #df$tissue <- factor(df$tissue, levels = df$tissue)

            plotly::plot_ly(
              data = df,
              ## name = pd$gene
              y = ~tissue,
              x = ~x,
              type = 'bar',
              orientation = 'h',
              color = ~color, ## TODO: use variable that encodes grouping
              colors = omics_pal_d()(length(unique(df$color)))
            ) %>%
            plotly::layout(
              yaxis = list(title = FALSE),
              xaxis = list(title = ylab),
              font = list(family = "Lato"),
              showlegend = FALSE,
              bargap = .4,
              margin = list(l = 10, r = 10, b = 10, t = 10)
            ) %>%
            plotly_default1()
        }

        modal_plot.RENDER <- function() {
            plot.RENDER() %>%
                plotly::layout(
                    showlegend = TRUE, ## TODO: I guess a legend makes sense here?
                    font = list(
                        size = 18
                    )
                )
        }

        PlotModuleServer(
            "pltmod",
            plotlib = "plotly",
            func = plot.RENDER,
            func2 = modal_plot.RENDER,
##          renderFunc = shiny::renderCachedPlot,
##          renderFunc2 = shiny::renderCachedPlot,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            res = c(90,170),                ## resolution of plots
            pdf.width = 8, pdf.height = 4,
            add.watermark = watermark
        )

    })  ## end of moduleServer
}


