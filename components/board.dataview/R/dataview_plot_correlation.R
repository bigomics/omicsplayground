##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

dataview_plot_correlation_ui <- function(id, label='', height=c(600,800)) {

    ns <- shiny::NS(id)
    info_text = "Barplot of the top positively and negatively correlated genes with the selected gene. Absolute expression levels of genes are colored in the barplot, where the low and high expressions range between the light and dark colors, respectively."
    
    PlotModuleUI(
        ns("pltsrv"),
        title = "Top correlated genes",
        label = label,
        plotlib = "plotly",
        ##outputFunc = plotly::plotlyOutput,
        ##outputFunc2 = plotly::plotlyOutput,        
        info.text = info_text,
        options = NULL,
        download.fmt=c("png","pdf","csv"),         
        width = c("auto","100%"),
        height = height
    )
  
}

dataview_plot_correlation_server <- function(id,
                                             pgx,
                                             r.gene    = reactive(""),
                                             r.samples = reactive(NULL),
                                             watermark=FALSE)
{
    moduleServer( id, function(input, output, session) {
    
        ##source("../../app/R/global.R",chdir=TRUE)
        ##load("../../../data/example-data.pgx")

        getTopCorrelatedGenes <- function(pgx, gene, n=30, samples=NULL) {
            
            ## precompute
            if(is.null(samples)) samples = colnames(pgx$X)
            samples <- intersect(samples, colnames(pgx$X))
            pp <- rownames(pgx$genes)[match(gene,pgx$genes$gene_name)]
            
            ## corr always in log.scale and restricted to selected samples subset
            ## should match exactly the rawtable!!
            if(pp %in% rownames(pgx$X)) {
                rho <- cor(t(pgx$X[,samples]), pgx$X[pp,samples], use="pairwise")[,1]
            } else if(pp %in% rownames(pgx$counts)) {
                x0 <- logCPM(pgx$counts[,samples])
                x1 <- x0[pp,]
                rho <- cor(t(x0), x1, use="pairwise")[,1]
            } else {
                rho <- rep(0, nrow(pgx$genes))
                names(rho) <- rownames(pgx$genes)
            }
            
            rho[is.na(rho)] <- 0
            jj = head(order(-abs(rho)),n)
            jj <- c(head(order(-rho),n/2), tail(order(-rho),n/2))
            top.rho = rho[jj]
            
            gx1 <- sqrt(rowSums(pgx$X[names(top.rho),samples]**2,na.rm=TRUE))
            gx1 <- (gx1 / max(gx1))
            klr1 <- rev(colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)),
                alpha = TRUE)(16))[1+round(15*gx1) ]
            klr1[which(is.na(klr1))] <- rgb(0.2,0.5,0.8,0.1)
            
            names(top.rho) = sub(".*:","",names(top.rho))

            getGeneAnnot <- function(ngs, genes) {
                ann <- ngs$genes[,grep("name|title|chr|map",colnames(ngs$genes),ignore.case=TRUE)]
                apply(ann[genes,],1, function(x)
                    paste(mapply(paste0,colnames(ann),": ",x),collapse='\n'))
            }
            annot <- getGeneAnnot(pgx, names(top.rho))
            
            df <- data.frame(
                genes = names(top.rho),
                rho = top.rho,
                color = klr1,
                value = gx1,
                annot = annot
            )

            res <- list(df, gene)
            
            return(res)
        }
        
        plot_data <- shiny::reactive({
            shiny::req(pgx$X,pgx$Y)
            shiny::req(r.gene())             
            samples <- r.samples()
            gene <- r.gene()
            pd <- getTopCorrelatedGenes(pgx, gene=gene, n=40, samples=samples)
            pd
        })
        
        plot.RENDER <- function() {
            pd <- plot_data()
            shiny::req(pd)
            
            par(mar=c(4.3,3.0,1,0), mgp=c(2.0,0.6,0))
            barplot(pd$top.rho, col=pd$klr1, ## horiz=TRUE,
                las = 3,  
                main = paste(pd$gene),
                ylab="correlation (r)",
                cex.names = 0.85, cex.main = 1, cex.axis = 0.9,
                col.main="#7f7f7f", border=NA)
            legend("topright", legend=c("expr high","expr low"),
                fill=c("#3380CCCC","#3380CC40"),
                cex=0.8, y.intersp=0.85)
        }

        plotly.RENDER <- function() {
            pd <- plot_data()
            shiny::req(pd)
            
            df <- pd[[1]]
            df$genes <- factor(df$genes, levels=df$genes)

            ## plot as regular bar plot            
            plotly::plot_ly(data = df,
                type = 'bar', ## name = pd$gene
                x = ~genes,
                y = ~rho,
                ## hoverinfo = "text",
                hovertext = ~annot,
                marker = list(color = ~color)
            )
        }
        
        modal_plotly.RENDER <- function() {
            plotly.RENDER() %>%
                plotly::layout(
                    ## showlegend = TRUE,
                    font = list(
                        size = 16
                    )
                )
        }
        
        PlotModuleServer(
            "pltsrv",
            plotlib = "plotly",
            func = plotly.RENDER,
            func2 = modal_plotly.RENDER,
            csvFunc = plot_data,   ##  *** downloadable data as CSV
            res = c(80,170),                ## resolution of plots
            pdf.width = 6, pdf.height = 6,
            add.watermark = watermark
        )
        
    })  ## end of moduleServer
}


