##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##
DataViewBoard <- function(input, output, session, inputData)
{
    ns <- session$ns ## NAMESPACE
    rowH = 355  ## row height of panels
    imgH = 315  ## height of images
    fullH = 750 ## full height of panel
    tabH = 600  ## height of tables

    ##----------------------------------------------------------------------
    ## More Info (pop up window)
    ##----------------------------------------------------------------------
    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>Group by</code>'

    data_infotext =paste0(
        'The <strong>DataView module</strong> provides information and visualisations of the dataset to quickly lookup a gene,
        check the counts, or view the data tables.<br><br>
        The <strong>Plots</strong> panel displays figures related to the expression level of the selected gene,
        correlation, and average expression ranking within the dataset.
        More information about the gene and hyperlinks to external databases are provided. Furthermore,
        it displays the correlation and tissue expression for a selected gene in external reference datasets.
        In the <strong>Counts</strong> panel, the total number of counts (abundance) per sample and their distribution among the samples are displayed.
        This is most useful to check the technical quality of the dataset, such as total read counts or abundance of ribosomal genes.
        In <strong>Gene Table</strong> panel, the exact expression values across the samples can be looked up,
        where genes are ordered by the correlation with respect to the first gene. Gene-wise average expression of a phenotype sample grouping
        is also presented in this table. In the <strong>Samples</strong> panel, more complete information about samples can be found.
        Finally, the <strong>Contrasts</strong> panel, shows information about the phenotype comparisons.
        <br><br><br>
        <center><iframe width="560" height="315" src="https://www.youtube.com/embed/S32SPINqO8E"
        title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
        allowfullscreen></iframe></center>
    ')

    ## ------- observe functions -----------
    shiny::observeEvent( input$data_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Data View Board</strong>"),
            shiny::HTML(data_infotext),
            easyClose = TRUE, size="l"))
    })

    ## update filter choices upon change of data set
    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)

        ## levels for sample filter
        levels = getLevels(ngs$Y)
        shiny::updateSelectInput(session, "data_samplefilter", choices=levels)
        genes <- sort(ngs$genes[rownames(ngs$X),]$gene_name)
        listGenes <- data.table::as.data.table(genes)
        names(listGenes) = '(type gene if missing!)'
        fc2 = rowMeans(pgx.getMetaFoldChangeMatrix(ngs)$fc**2)
        selgene = names(sort(-fc2))[1] ## most var gene??
        
        sel1 = which(genes==selgene)
        shiny::updateSelectInput(session,'search_gene', choices=listGenes, selected=selgene)
        grps <- pgx.getCategoricalPhenotypes(ngs$samples, min.ncat=2, max.ncat=999)
        grps <- sort(grps)
        selgrp <- grps[1]
        grps <- c("<ungrouped>",grps)
        if("group" %in% grps) selgrp = "group"
        if(nrow(ngs$samples)<=20) selgrp = "<ungrouped>"
        shiny::updateSelectInput(session,'data_groupby', choices=grps, selected=selgrp)
    })


    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================


    ##----------------------------------------------------------------------
    ##                     Info messages
    ##----------------------------------------------------------------------

    genePlots_tsne_text=paste0('<b>T-SNE clustering</b> of samples (or cells) colored by an expression of the gene selected in the ',dropdown_search_gene, ' dropdown menu. The red color represents an over-expression of the selected gene across samples (or cells).')
    genePlots_barplot_text=paste0('Expression barplot of grouped samples (or cells) for the gene selected in the ',dropdown_search_gene,'. Samples (or cells) in the barplot can be ungrouped by setting the ',menu_grouped, ' under the main <i>Options</i>.')
    genePlots_correlationplot_text=paste0('Barplot of the top positively and negatively correlated genes with the selected gene. Absolute expression levels of genes are colored in the barplot, where the low and high expressions range between the light and dark colors, respectively.')
    genePlots_averageRankPlot_text=paste0('Ranking of the average expression of the selected gene.')
    data_geneInfo_text = paste0('For more information about the the selected gene, follow the hyperlinks to public databases, including ', a_OMIM,', ', a_KEGG, ' and ',a_GO,'.')
    data_tissueplot_text = paste0('Tissue expression for the selected gene in the tissue expression ',a_GTEx,' dataset. Colors corresponds to "tissue clusters" as computed by unsupervised clustering.')


    ##----------------------------------------------------------------------
    ##                     Average Rank plot
    ##----------------------------------------------------------------------
    MARGINS1 = c(7,3.5,2,1)

    genePlots_averageRankPlot.RENDER <- shiny::reactive({


        ngs <- inputData()
        alertDataLoaded(session,ngs)
        shiny::req(ngs)

        dbg("[genePlots_averageRankPlot.RENDER] reacted")
        
        gene <- ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples <- colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples <- length(samples)
        
        if(input$data_type=="counts") {
            mean.fc <- sort(rowMeans(ngs$counts[,samples,drop=FALSE]),decreasing=TRUE)
            ylab="expression (counts)"
        }
        if(input$data_type=="logCPM") {
            mean.fc <- sort(rowMeans(ngs$X[,samples,drop=FALSE]),decreasing=TRUE)
            ylab="expression (log2CPM)"
        }

        j <- which(sub(".*:","",names(mean.fc))==gene)

        mar <- MARGINS1
        mar[4] <- 0
        par(mar=mar, mgp=c(2.1,0.8,0))
        par(mar=c(2.3,3.0,2,2), mgp=c(2.0,0.6,0))
        ##MARGINS1
        base::plot( mean.fc, type="h", lwd=0.4,
                   col="#bbd4ee", cex.axis=0.9,
                   ylab=ylab, xlab="ordered genes", xaxt="n")
        points( j, mean.fc[j], type="h", lwd=2, col="black")
        text( j, mean.fc[j], gene, pos=3, cex=0.9)

        dbg("[genePlots_averageRankPlot.RENDER] done")

    })
   
    shiny::callModule(
        plotModule, id="genePlots_averageRankPlot",
        func = genePlots_averageRankPlot.RENDER,
        func2 = genePlots_averageRankPlot.RENDER,
        info.text = genePlots_averageRankPlot_text,
        height = imgH, ## width = '100%',
        pdf.width=6, pdf.height=6,
        label="c", title="Average rank",
        add.watermark = WATERMARK
    )
    
    ##----------------------------------------------------------------------
    ##                     Correlation plot
    ##----------------------------------------------------------------------

    getTopCorrelatedGenes <- function(ngs, gene, n=30, samples=NULL) {

        ## precompute
        if(is.null(samples)) samples = colnames(ngs$X)
        samples <- intersect(samples, colnames(ngs$X))
        pp=rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]

        ## corr always in log.scale and restricted to selected samples subset
        ## should match exactly the rawtable!!
        if(pp %in% rownames(ngs$X)) {
            rho <- cor(t(ngs$X[,samples]), ngs$X[pp,samples], use="pairwise")[,1]
        } else if(pp %in% rownames(ngs$counts)) {
            x0 <- logCPM(ngs$counts[,samples])
            x1 <- x0[pp,]
            rho <- cor(t(x0), x1, use="pairwise")[,1]
        } else {
            rho <- rep(0, nrow(ngs$genes))
            names(rho) <- rownames(ngs$genes)
        }

        rho[is.na(rho)] <- 0
        jj = head(order(-abs(rho)),30)
        jj <- c(head(order(-rho),15), tail(order(-rho),15))
        top.rho = rho[jj]

        gx1 <- sqrt(rowSums(ngs$X[names(top.rho),samples]**2,na.rm=TRUE))
        gx1 <- (gx1 / max(gx1))
        klr1 <- rev(colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)),
                                     alpha = TRUE)(16))[1+round(15*gx1) ]
        klr1[which(is.na(klr1))] <- rgb(0.2,0.5,0.8,0.1)

        names(top.rho) = sub(".*:","",names(top.rho))
        offset = min(top.rho)*0.95
        offset = 0

        dbg("[genePlots_correlationplot_data()] done")

        res = list(top.rho=top.rho, offset=offset, klr1=klr1)
        return(res)
    }

    genePlots_correlationplot.RENDER <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)

        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene

        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }

        res <- getTopCorrelatedGenes(ngs, gene=gene, n=30, samples=samples)
        ## res = genePlots_correlationplot_data()

        mar = MARGINS1
        par(mar=mar, mgp=c(2.1,0.8,0))
        barplot(res$top.rho - res$offset, col=res$klr1, ## horiz=TRUE,
                las = 3,#main=paste("top correlated genes\nwith",gene),
                main = paste(gene),
                offset = res$offset, ylab="correlation (r)",
                ## names.arg=rep(NA,length(top.rho)),
                cex.names = 0.85, cex.main = 1, cex.axis = 0.9,
                col.main="#7f7f7f", border=NA)
        ##text( (1:length(top.rho) - 0.5)*1.2, offset, names(top.rho),
        ##col="black", cex=0.75, srt=90, pos=3, offset=0.4, font=1)
        legend("topright", legend=c("expr high","expr low"),
               ##fill=c("grey30","grey80"),
               fill=c("#3380CCCC","#3380CC40"),
               cex=0.8, y.intersp=0.85)

        dbg("[genePlots_correlationplot.RENDER] done")

    })

    ##genePlots_correlationplot_module <- plotModule(
    ##    id="genePlots_correlationplot", ns=ns,
    shiny::callModule(
        plotModule, "genePlots_correlationplot",
        func = genePlots_correlationplot.RENDER,
        func2 = genePlots_correlationplot.RENDER,
        info.text=genePlots_correlationplot_text,
        height = imgH, pdf.width=6, pdf.height=6,
        label="e", title="Top correlated genes",
        add.watermark = WATERMARK
    )
    ##output <- attachModule(output, genePlots_correlationplot_module)

    ##----------------------------------------------------------------------
    ##                     Bar/box plot
    ##----------------------------------------------------------------------
    
    genePlots_barplot.RENDER <- shiny::reactive({


        cat("[dataview] genePlots_barplot.RENDER reacted\n")

        ngs <- inputData()
        shiny::req(ngs)
        shiny::req(input$data_groupby,input$search_gene,input$data_type)

        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples = length(samples)

        grpvar=1
        grpvar <- input$data_groupby
        ##grp  = factor(ngs$Y[samples,grpvar])
        grp  = factor(as.character(ngs$Y[samples,grpvar]))
        klr0 = COLORS
        klr  = klr0[as.integer(grp)]

        ## precompute
        pp=rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]

        gx = NULL
        ylab = NULL
        if(input$data_type=="counts") {
            gx = ngs$counts[pp,samples]
            ylab="expression (counts)"
        } else if(input$data_type=="CPM") {
            gx = 2**ngs$X[pp,samples]
            ylab="expression (CPM)"
        } else if(input$data_type=="logCPM") {
            gx = ngs$X[pp,samples]
            ylab="expression (log2CPM)"
        }

        mar=MARGINS1
        par(mar=mar, mgp=c(2.1,0.8,0))

        BLUE = col=rgb(0.2,0.5,0.8,0.8)
        bee.cex = ifelse(length(gx)>500,0.1,0.2)
        bee.cex = c(0.3,0.1,0.05)[cut(length(gx),c(0,100,500,99999))]

        ##if(input$data_grouped) {
        if(input$data_groupby != "<ungrouped>") {
            nnchar = nchar(paste(unique(grp),collapse=''))
            srt = ifelse(nnchar < 20, 0, 35)
            srt
            ngrp <- length(unique(grp))
            cx1 = ifelse( ngrp < 10, 1, 0.8)
            cx1 = ifelse( ngrp > 20, 0.6, cx1)
            ##cx1 = ifelse( ngrp > 10, 0.6, 0.9)
            if(input$geneplot_type == 'bar') {
                gx.b3plot(
                    gx, grp, las=3, main=gene, ylab=ylab,
                    cex.main=1, col.main="#7f7f7f",
                    bar=TRUE, border=NA, ## bee = ifelse(length(gx) < 500,TRUE,FALSE),
                    bee.cex = bee.cex, ## sig.stars=TRUE, max.stars=5,
                    xlab="", names.cex=cx1, srt=srt,
                    ## col=klr0[ii],
                    col = rgb(0.4,0.6,0.85,0.85)
                )
            } else if(input$geneplot_type == 'violin') {
                ##vioplot::vioplot( gx ~ grp, main = gene, cex.main=1.0,
                ##                 ylab = ylab, xlab='',
                ##                 col = rgb(0.2,0.5,0.8,0.8))
                pgx.violinPlot(gx, grp, main = gene, cex.main=1,
                               xlab = '', ylab = ylab,
                               ##vcol = rgb(0.2,0.5,0.8,0.8),
                               vcol = rgb(0.4,0.6,0.85,0.85),
                               srt = srt)

            } else {
                boxplot(
                    gx ~ grp, main = gene, cex.main=1.0,
                    ylab = ylab, xlab='', xaxt='n',
                    col =  rgb(0.4,0.6,0.85,0.85)
                )
                yy <- sort(unique(grp))
                text(x = 1:length(yy),
                     y = par("usr")[3] - 0.03*diff(range(gx)),
                     labels = yy,
                     xpd = NA,
                     srt = srt,
                     adj = ifelse(srt==0, 0.5, 0.965),
                     cex = cx1)
            }



        }  else {
            jj <- 1:length(gx)
            sorting="no"
            if(sorting == "decr")  jj <- order(-gx)
            if(sorting == "inc")  jj <- order(gx)
            tt=""
            barplot(gx[jj], col=BLUE, ##col=klr[jj],
                    las = 3, cex.names = 0.8,
                    ylab = ylab, xlab = tt,
                    main = gene, cex.main=1, col.main="#7f7f7f", border=NA,
                    names.arg = rep(NA,length(gx)) )
            if(length(gx)<100) {
                cx1 = ifelse(length(gx) > 20, 0.8, 0.9)
                cx1 = ifelse(length(gx) > 40, 0.6, cx1)
                cx1 = ifelse(length(gx) < 10, 1, cx1)
                text((1:length(gx)-0.5)*1.2, -0.04*max(gx), names(gx)[jj],
                     las=3, cex=cx1, pos=2, adj=0, offset=0, srt=45, xpd=TRUE)
            }
        }
    })

    genePlots_barplot.opts <- shiny::tagList(
        shiny::radioButtons(ns('geneplot_type'),'plot type (grouped)', c('bar','violin','box'),
                     inline=TRUE)
    )

    shiny::callModule(
        plotModule, "genePlots_barplot",
        func = genePlots_barplot.RENDER,
        func2 = genePlots_barplot.RENDER,
        options = genePlots_barplot.opts,
        info.text = genePlots_barplot_text,
        height = imgH,
        pdf.width = 7, pdf.height = 5,
        label="b", title="Abundance/expression",
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ## t-SNE
    ##----------------------------------------------------------------------

    genePlots_tsne.RENDER <- shiny::reactive({


        ngs <- inputData()
        shiny::req(ngs)

        dbg("[genePlots_tsne.RENDER] reacted")

        gene = "KCNN4"
        gene = ngs$genes$gene_name[1]
        if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
        samples = colnames(ngs$X)
        if(!is.null(input$data_samplefilter)) {
            samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        }
        nsamples = length(samples)

        ## precompute
        pp=rownames(ngs$genes)[1]
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]

        gx = NULL
        ylab = NULL
        if(input$data_type=="counts") {
            gx = ngs$counts[pp,samples]
            ylab="expression (counts)"
        } else if(input$data_type=="CPM") {
            gx = 2**ngs$X[pp,samples]
            ylab="expression (CPM)"
        } else if(input$data_type=="logCPM") {
            gx = ngs$X[pp,samples]
            ylab="expression (log2CPM)"
        }

        ##par(mar=c(12,2,2,1), mgp=c(2.1,0.8,0), oma=c(3,0.5,1.5,0.3))

        pos <- ngs$tsne2d[samples,]

        cex1 <- 1.8*c(1.6,1.0,0.6,0.3)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]
        klrpal = colorRampPalette(c("blue3", "aliceblue", "grey85", "lavenderblush", "red3"))(16)
        klrpal = colorRampPalette(c("grey80", "grey50", "red3"))(16)

        fc1 <- tanh(0.99 * scale(gx)[,1])
        fc1 <- tanh(0.99 * scale(gx,center=FALSE)[,1])
        ##fc1 <- tanh(0.99 * gx/sd(gx))
        fc2 <- (fc1 - min(fc1))
        klr1 = klrpal[1 + round(15*fc2/max(abs(fc2)))]
        klr1 = paste0(gplots::col2hex(klr1),"88")

        jj2 <- order(abs(fc1))

        ## determine how to do grouping for group labels
        groupby <- input$data_groupby
        grp <- NULL
        if(groupby != "<ungrouped>") {
            grp <- factor(ngs$samples[samples, groupby])
        }

        data <- data.frame(pos[jj2,])
        data$grp <- grp

        ## TODO: does currently not render in app, throws `need finite 'xlim' values` error
        ## NOTE: for now I have removed the individual color to use colors for the ellipses
       fig <-
         ggplot(data, aes(tSNE.x, tSNE.y)) +
          labs(x = "tSNE1", y = "tSNE2") +
          theme_bw(base_size = 13)

       if (!is.null(grp)) {
         if(input$show_cluster == 'yes') {
         fig <- fig +
           ggforce::geom_mark_hull(
             aes(fill = grp, label = grp, color = grp,
                 color = after_scale(colorspace::desaturate(color, .3)),
                 fill = after_scale(colorspace::desaturate(color, .5))),
             expand = unit(2.7, "mm"), con.cap = unit(.01, "mm"),
             label.buffer = unit(2, "mm"), alpha = .15,
             label.fontsize = 12.5, label.fontface = "plain"
           ) +
           geom_point(aes(color = grp), size = 1.5) +
           scale_x_continuous(expand = c(.4, .4)) +
           scale_y_continuous(expand = c(.4, .4)) +
           scale_color_discrete(guide = "none") +
           scale_fill_discrete(guide = "none")
         } else {
           fig <- fig +
             geom_point(aes(color = grp), size = 2) +
             scale_color_discrete(name = NULL)
         }
       } else {
         fig <- fig +
           #geom_point(aes(color = expression), size = 2)
           geom_point(size = 2)
       }

       gridExtra::grid.arrange(fig)

        dbg("[genePlots_tsne.RENDER] done")

    })

    genePlots_tsne_max.RENDER <- shiny::reactive({


      ngs <- inputData()
      shiny::req(ngs)

      dbg("[genePlots_tsne.RENDER] reacted")

      gene = "KCNN4"
      gene = ngs$genes$gene_name[1]
      if(!is.null(input$search_gene) && input$search_gene!="") gene <- input$search_gene
      samples = colnames(ngs$X)
      if(!is.null(input$data_samplefilter)) {
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
      }
      nsamples = length(samples)

      ## precompute
      pp=rownames(ngs$genes)[1]
      pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]

      gx = NULL
      ylab = NULL
      if(input$data_type=="counts") {
        gx = ngs$counts[pp,samples]
        ylab="expression (counts)"
      } else if(input$data_type=="CPM") {
        gx = 2**ngs$X[pp,samples]
        ylab="expression (CPM)"
      } else if(input$data_type=="logCPM") {
        gx = ngs$X[pp,samples]
        ylab="expression (log2CPM)"
      }

      ##par(mar=c(12,2,2,1), mgp=c(2.1,0.8,0), oma=c(3,0.5,1.5,0.3))

      pos <- ngs$tsne2d[samples,]

      cex1 <- 1.8*c(1.6,1.0,0.6,0.3)[cut(nrow(pos),breaks=c(-1,40,200,1000,1e10))]
      klrpal = colorRampPalette(c("blue3", "aliceblue", "grey85", "lavenderblush", "red3"))(16)
      klrpal = colorRampPalette(c("grey80", "grey50", "red3"))(16)

      fc1 <- tanh(0.99 * scale(gx)[,1])
      fc1 <- tanh(0.99 * scale(gx,center=FALSE)[,1])
      ##fc1 <- tanh(0.99 * gx/sd(gx))
      fc2 <- (fc1 - min(fc1))
      klr1 = klrpal[1 + round(15*fc2/max(abs(fc2)))]
      klr1 = paste0(gplots::col2hex(klr1),"88")

      jj2 <- order(abs(fc1))

      ## determine how to do grouping for group labels
      groupby <- input$data_groupby
      grp <- NULL
      if(groupby != "<ungrouped>") {
        grp <- factor(ngs$samples[samples, groupby])
      }

      data <- data.frame(pos[jj2,])
      data$grp <- grp

      ## TODO: does currently not render in app, throws `need finite 'xlim' values` error
      ## NOTE: for now I have removed the individual color to use colors for the ellipses
      fig <-
        ggplot(data, aes(tSNE.x, tSNE.y)) +
        labs(x = "tSNE1", y = "tSNE2") +
        theme_bw(base_size = 13)

      if (!is.null(grp)) {
        if(input$show_cluster == 'yes') {
          fig <- fig +
            ggforce::geom_mark_hull(
              aes(fill = grp, label = grp, color = grp,
                  color = after_scale(colorspace::desaturate(color, .3)),
                  fill = after_scale(colorspace::desaturate(color, .5))),
              expand = unit(3.4, "mm"), con.cap = unit(.01, "mm"),
              label.buffer = unit(3, "mm"), alpha = .15,
              label.fontsize = 22, label.fontface = "plain"
            ) +
            geom_point(aes(color = grp), size = 3.5) +
            scale_x_continuous(expand = c(.15, .15)) +
            scale_y_continuous(expand = c(.15, .15)) +
            scale_color_discrete(guide = "none") +
            scale_fill_discrete(guide = "none")
        } else {
          fig <- fig +
            geom_point(aes(color = grp), size = 4.5) +
            scale_color_discrete(name = NULL)
        }
      } else {
        fig <- fig +
          #geom_point(aes(color = expression), size = 4.5)
          geom_point(size = 4.5)
      }

      gridExtra::grid.arrange(fig)

      dbg("[genePlots_tsne.RENDER] done")

    })

    genePlots_tsne.opts <- shiny::tagList(
      shiny::radioButtons(ns('show_cluster'), 'show cluster?', c('yes', 'no'),
                          inline = TRUE)
    )

    shiny::callModule(
        plotModule, "genePlots_tsne",
        #plotlib = "ggplot",
        func = genePlots_tsne.RENDER,
        func2 = genePlots_tsne_max.RENDER,
        options = genePlots_tsne.opts,
        info.text = genePlots_tsne_text,
        height = imgH, pdf.width = 6, pdf.height = 6,
        label = "d", title = "t-SNE clustering",
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ##  Tissue expression plot
    ##----------------------------------------------------------------------

    data_tissueplot.RENDER  <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)
        if(is.null(input$data_type)) return(NULL)

        dbg("[data_tissueplot.RENDER] reacted")

        gene <- input$search_gene
        pp <- rownames(ngs$genes)[match(gene,ngs$genes$gene_name)]
        hgnc.gene = toupper(as.character(ngs$genes[pp,"gene_name"]))


        par(mar=c(6,4,1,1), mgp=c(2.2,0.8,0))
        mar=MARGINS1
        par(mar=mar, mgp=c(1.5,0.5,0))

        if( hgnc.gene %in% rownames(TISSUE)) {
            tx = TISSUE[hgnc.gene,]
            grp = TISSUE.grp[names(tx)]
            tissue.klr = COLORS[grp]
            ylab="expression (TPM)"
            if(input$data_type=="logCPM") {
                ylab = "expression (log2TPM)"
                tx = log(1 + tx)
            }
            jj <- 1:length(tx)
            sorting="no"
            if(sorting=="decr") jj <- order(-tx)
            if(sorting=="inc") jj <- order(tx)

            barplot(tx[jj], las=3, main=hgnc.gene, cex.main=1, col.main="#7f7f7f",
                    col = tissue.klr[jj], border=NA,ylab=ylab, cex.names=0.9,
                    names.arg=rep(NA,length(tx)))
                                        #title(main="tissue expression", line=-0.3, cex.main=1.1)
            text((1:length(tx)-0.5)*1.2, -0.04*max(tx), names(tx)[jj], las=3,
                 cex=0.85, pos=2, adj=0, offset=0, srt=55, xpd=TRUE)

        } else {
            frame()
        }

        dbg("[data_tissueplot.RENDER] done")

    })

    shiny::callModule(
        plotModule, "data_tissueplot",
        func = data_tissueplot.RENDER,
        func2 = data_tissueplot.RENDER,
        info.text = data_tissueplot_text,
        height = imgH, pdf.width=9, pdf.height=6,
        label="f", title="Tissue expression",
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ## Gene information
    ##----------------------------------------------------------------------

    data_geneInfo.RENDER  <- shiny::reactive({
        require(org.Hs.eg.db)
        gene = "A1BG-AS1"
        gene = "CD4"
        gene <- input$search_gene
        gene = toupper(sub(".*:","",gene))

        eg = "1017"
        eg = names(which(as.list(org.Hs.egSYMBOL)==gene))
        eg <- mget(gene, envir=org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
        if(is.na(eg)) eg <- mget(gene, envir=org.Hs.egALIAS2EG, ifnotfound=NA)[[1]]
        eg
        eg = eg[1]
        if(is.null(eg) || length(eg)==0) return(NULL)

        output = "(gene info not available)"
        if(length(eg)>0 && !is.na(eg)) {
            ##as.list(org.Hs.eg.db::org.Hs.egSYMBOL)[[eg]]
            info <- getHSGeneInfo(eg)  ## defined in pgx-functions.R
            info$summary <- '(no info available)'
            if(FALSE && input$data_geneinfo) {
                suppressWarnings(suppressMessages(
                    info2 <- try(getMyGeneInfo(eg, fields="summary"))
                ))
                if(!any(class(info2) == 'try-error')) {
                    info <- c(info, info2)  ## defined in pgx-functions.R
                }
            }
            if(gene %in% names(GENE.SUMMARY)) {
                info$summary <- GENE.SUMMARY[gene]
                info$summary <- gsub('Publication Note.*|##.*','',info$summary)
            }

            ## reorder
            nn <- intersect(c("symbol","name","map_location","summary",names(info)),names(info))
            info <- info[nn]
            info$symbol <- paste0(info$symbol,'<br>')

            output <- c()
            for(i in 1:length(info)) {
                xx <- paste(info[[i]], collapse=", ")
                output[[i]] <- paste0("<b>",names(info)[i],"</b>: ",xx)
            }
            output <- paste(output, collapse="<p>")
        }    
        shiny::wellPanel(shiny::HTML(output))
    })

    shiny::callModule(
        plotModule, "data_geneInfo",
        title = "Gene info", label="a",
        plotlib = "generic",
        func = data_geneInfo.RENDER,
        func2 = data_geneInfo.RENDER,
        renderFunc = "renderUI", outputFunc = "htmlOutput",
        just.info = FALSE, no.download = TRUE,
        info.text = data_geneInfo_text,
        ## options = shiny::tagList(),
        height = c(fullH,600), width=c('auto',800),
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ##                     Interface
    ##----------------------------------------------------------------------
    dataview_caption1 = "<b>Gene plots.</b> <b>(a)</b> Further information about the selected gene from public databases. <b>(b)</b> Abundance/expression of selected gene across groups. <b>(c)</b> Average rank of the selected gene compared to other genes. <b>(d)</b> t-SNE of samples colored by expression of selected gene. <b>(e)</b> Top correlated genes. Darker color corresponds to higher expression of the gene. <b>(f)</b> Tissue expression of selected gene."
    output$plotsUI <- shiny::renderUI({
        shiny::fillCol(
            height = fullH,
            flex = c(NA,0.03,1),
            shiny::div(shiny::HTML(dataview_caption1), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,0.06,5),
                plotWidget(ns("data_geneInfo")),
                shiny::br(),
                shiny::fillCol(
                    flex = c(1,0.2,1),
                    shiny::fillRow(
                        flex = c(1.5,1,1), id = "genePlots_row1",
                        height = rowH, ## width=1600,
                        plotWidget(ns("genePlots_barplot")),
                        plotWidget(ns("genePlots_averageRankPlot")),
                        plotWidget(ns("genePlots_tsne"))
                    ),
                    shiny::br(),
                    shiny::fillRow(
                        flex = c(1.5,2), id = "genePlots_row2",
                        height = rowH, ## width=1600,
                        plotWidget(ns("genePlots_correlationplot")),
                        ##plotWidget(ns("data_corplot")),
                        plotWidget(ns("data_tissueplot"))
                    )
                )
            )
        )
    })

    ##dragula("genePlots_row1")
    ##dragula("genePlots_row2")
    ##dragula(c("genePlots_row1","genePlots_row2"))

    ##----------------------------------------------------------------------
    ##                     Info messages for Counts
    ##----------------------------------------------------------------------

    counts_tab_barplot_text=paste0('Barplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_boxplot_text=paste0('Boxplot of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_histplot_text=paste0('Histogram of the total number of counts (abundance) for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_abundanceplot_text=paste0('Barplot showing the percentage of counts in terms of major gene types such as CD molecules, kinanses or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')
    counts_tab_average_countplot_text=paste0('Barplot showing the average count levels of major gene types such as CD molecules, kinanses or RNA binding motifs for each group. The samples (or cells) can be grouped/ungrouped in the ',menu_grouped, ' setting uder the main <i>Options</i>.')


    ##----------------------------------------------------------------------
    ##                     Count information barplot
    ##----------------------------------------------------------------------
    MARGINS2 = c(9,3.5,2,0.5)
    MARGINS2 = c(8,3.5,2,0.5)

    counts_tab_barplot.RENDER <- shiny::reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
        shiny::req(input$data_groupby)
                                        #par(mar=c(20,6,10,6), mgp=c(2.2,0.8,0))
        par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))

        ylab = "counts (million)"
        if(input$data_groupby!="<ungrouped>") {
            ylab = "average group counts (million)"
        }
        subtt0=""
        if(!is.null(res$subtt)) subtt0 <- paste0("\n(",paste(res$subtt,collapse=","),")")
        ## ---- xlab ------ ###
        names.arg = names(res$total.counts)
        if( length(names.arg) > 20){ names.arg = "" }
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        ## ---- xlab ------ ###
        barplot(res$total.counts/1e6, las=3, border = NA,
                ##main=paste("total counts",subtt0), cex.main=1.7,
                col=rgb(0.2,0.5,0.8,0.8), #col="grey40",
                                        #cex.names=res$cx1+0.01, cex.lab=1, ylab=ylab,
                cex.names=cex.names, cex.lab=1, ylab=ylab,
                ylim=c(0,max(res$total.counts)/1e6)*1.1,
                names.arg=names.arg)
    })

    shiny::callModule(
        plotModule, "counts_tab_barplot",
        func = counts_tab_barplot.RENDER,
        func2 = counts_tab_barplot.RENDER,
        info.text = counts_tab_barplot_text,
        height=imgH, pdf.width=7, pdf.height=6, ## res=45,
        label="a",title='Total counts',
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ##                     Count information boxplot
    ##----------------------------------------------------------------------

    counts_tab_boxplot.RENDER <- shiny::reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
        ##par(mar=c(3,3,3,3), mgp=c(2.4,0.7,0), oma=c(1,1,1,1)*0.2 )
        par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))
        ## ---- xlab ------ ###
        xaxt="l"
        names.arg = colnames(res$log2counts[res$jj,])
        if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)); xaxt="n"}
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        boxplot(res$log2counts[res$jj,], col=rgb(0.2,0.5,0.8,0.4),
                ##col=rgb(0.2,0.5,0.8,0.3), #col="grey70",
                ##main="counts distribution", cex.main=1.6,
                ## cex.names=res$cx1+0.1,
                names = names.arg, cex.axis=cex.names,#border=rgb(0.2,0.5,0.8,0.8),
                border = 	rgb(0.824,0.824,0.824,0.9),xaxt=xaxt,
                las=3, cex.lab=1, ylab="counts (log2)", outline=FALSE, varwidth = FALSE)
    })

    shiny::callModule(
        plotModule, "counts_tab_boxplot",
        func = counts_tab_boxplot.RENDER,
        func2 = counts_tab_barplot.RENDER,
        info.text = counts_tab_boxplot_text,
        height=imgH, pdf.width=7, pdf.height=6, ## res=50,
        label="b", title='Counts distribution',
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ##                     Count information histogram
    ##----------------------------------------------------------------------

    counts_tab_histplot.RENDER <- shiny::reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
                                        #par(mfrow=c(2,3), mar=c(9,4,3,1.5), mgp=c(2.4,0.7,0), oma=c(1,1,1,1)*0.2 )
        par(mar=c(8,4,1,2), mgp=c(2.2,0.8,0))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))

        n=1000
        gx.hist <- function(gx, n=1000, main="",ylim=NULL) {
            jj <- 1:nrow(gx)
            if(length(jj)>n) jj <- sample(jj,n,replace=TRUE)
            h0 <- hist(as.vector(c(gx[jj],min(gx),max(gx))),
                       breaks=120, main=main, border=FALSE,
                                        #col=rgb(0.2,0.5,0.8,0.5),
                       col="grey",
                       freq=FALSE, ## ylim=ylim,
                       xlim=c(min(gx),max(gx)),
                       xlab="expression (log2)",
                       ##cex.names=res$cx1+0.1,
                       cex.lab=1)
            i = 1
            for(i in 1:ncol(gx)) {
                h1 <- hist(gx[jj,i], breaks=h0$breaks,plot=FALSE)
                ##lines( h0$mids, h1$density, col=i+1 )
                lines( h0$mids, h1$density, col="black", lwd=0.5 )
            }
        }
        gx.hist(gx=res$log2counts[res$jj,], n=2000) #, main="histogram")
    })

    shiny::callModule(
        plotModule, "counts_tab_histplot",
        func = counts_tab_histplot.RENDER,
        func2 = counts_tab_histplot.RENDER,
        info.text = counts_tab_histplot_text,
        height=imgH, pdf.width=7, pdf.height=6, ## res=50,
        label="c", title='Counts histogram',
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ##  Count information abundance of major gene types
    ##----------------------------------------------------------------------

    counts_tab_abundanceplot.RENDER <- shiny::reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
                                        #par(mar=c(6,4,0,4), mgp=c(2.2,0.8,0))
        par(mar=c(6,4,0,4))
        par(mar=c(7,4,0,2))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))

        klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)), alpha = TRUE)(nrow(res$prop.counts))

        ## ----------- gene type counts
        ymax = max(colSums(res$prop.counts, na.rm=TRUE))
        ## ---- xlab ------ ###
        names.arg = colnames(res$prop.counts)
        if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)) }
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        ## ---- xlab ------ ###
        barplot(res$prop.counts, las=3, #main="abundance of major gene types", cex.main=1.6,
                                        #cex.names=res$cx1+0.04,
                cex.lab=1.0, border = NA,
                ylim = c(0,ymax)*1.6, ylab = "abundance (%)",
                names.arg = names.arg, cex.names = cex.names,
                col = klr)
        leg <- legend("topleft", legend=rev(rownames(res$prop.counts)),
                                        #fill=rev(grey.colors(nrow(res$prop.counts))),
                      fill=rev(klr),cex=1, y.intersp=0.75, bty="n", plot = FALSE)
        leftx <- leg$rect$left*0.9
        rightx <- leg$rect$right*0.9
        topy <- leg$rect$top
        bottomy <- leg$rect$bottom
        legend(x = c(leftx, rightx), y = c(topy, bottomy),
               legend=rev(rownames(res$prop.counts)),
                                        #fill=rev(grey.colors(nrow(res$prop.counts))),
               fill=rev(klr), bty="n", cex=0.9, y.intersp=0.75)

    })

    shiny::callModule(
        plotModule, "counts_tab_abundanceplot",
        func = counts_tab_abundanceplot.RENDER,
        func2 = counts_tab_abundanceplot.RENDER,
        info.text = counts_tab_abundanceplot_text,
        height=imgH, pdf.width=10, pdf.height=6, ## res=50,
        label="d", title='Abundance of major gene types',
        add.watermark = WATERMARK
    )

    ##----------------------------------------------------------------------
    ## Count information average count by gene type
    ##----------------------------------------------------------------------
    counts_tab_average_countplot.RENDER <- shiny::reactive({
        res = getCountsTable()
        if(is.null(res)) return(NULL)
                                        #par(mar=c(6,4,0,4), mgp=c(2.2,0.8,0))
        par(mar=c(6,4,0,4))
        par(mar=c(7,4,0,2))
        par(mar=MARGINS2, mgp=c(2.2,0.8,0))

        klr <- colorRampPalette(c(rgb(0.2,0.5,0.8,0.8), rgb(0.2,0.5,0.8,0.1)),
                                alpha = TRUE)(nrow(res$avg.counts))
        ## klr <- grey.colors(nrow(res$avg.counts))

        ## ---- xlab ------ ###
        names.arg = colnames(res$avg.counts)
        if( length(names.arg) > 20){ names.arg = rep("",length(names.arg)) }
        cex.names <- ifelse(length(names.arg)>10,0.8,0.9)
        ## ---- xlab ------ ###
        ymax = max(colSums(res$avg.counts, na.rm=TRUE))
        barplot(res$avg.counts, las=3, #main="average counts by gene type", cex.main=1.6,
                                        #cex.names=res$cx1+0.04,
                border=NA, cex.lab=1.0,
                names.arg=names.arg, cex.names=cex.names,
                ylim=c(0,ymax)*1.6, ylab="average count", col=klr)
        leg <- legend("topleft", legend=rev(rownames(res$avg.counts)),
                                        #fill=rev(grey.colors(nrow(res$avg.counts))),
                      fill=rev(klr), cex=1, y.intersp=0.75, bty="n", plot = FALSE)
        leftx <- leg$rect$left*0.9
        rightx <- leg$rect$right*0.9
        topy <- leg$rect$top
        bottomy <- leg$rect$bottom
        legend(x = c(leftx, rightx), y = c(topy, bottomy),
               legend=rev(rownames(res$avg.counts)),
                                        #fill=rev(grey.colors(nrow(res$avg.counts))),
               fill=rev(klr), cex=0.9, y.intersp=0.75, bty="n")
    })

    shiny::callModule(
        plotModule, "counts_tab_average_countplot",
        func = counts_tab_average_countplot.RENDER,
        func2 = counts_tab_average_countplot.RENDER,
        info.text = counts_tab_average_countplot_text,
        height=imgH, pdf.width=10, pdf.height=6, ##res=50,
        label="e", title='Average count by gene type',
        add.watermark = WATERMARK
    )

    getCountsTable <- shiny::reactive({
        ngs = inputData()
        shiny::req(ngs)

        shiny::validate(shiny::need("counts" %in% names(ngs), "no 'counts' in object."))
        subtt=NULL

        samples = colnames(ngs$X)
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        nsamples = length(samples)
        if("counts" %in% names(ngs)) {
            counts = ngs$counts[,samples,drop=FALSE]
        } else {
            cat("WARNING:: no counts table. estimating from X\n")
            counts = pmax(2**ngs$X - 1,0)
            k = grep("lib.size",colnames(ngs$samples))[1]
            if(length(k)>0) {
                libsize = ngs$samples[colnames(counts),k]
                libsize
                counts = t(t(counts) * libsize)
            }
            ##counts <- round(counts)
        }
        if(sum(is.na(counts))>0) {
            cat("WARNING:: plot counts: counts has missing values!\n")
        }

        ##if(input$data_sampling=="grouped") {
        grpvar <- input$data_groupby
        gr = ngs$Y[samples,grpvar]
        grps = sort(unique(gr))
        ##if(input$data_grouped && length(grps)>1 ) {
        if(input$data_groupby != "<ungrouped>" && length(grps)>1) {
            newx = c()
            for(g in grps) {
                mx = rowMeans(counts[,which(gr==g),drop=FALSE], na.rm=TRUE)
                ## mx = rowSums(counts[,which(gr==g),drop=FALSE], na.rm=TRUE)  ## SUM or MEAN???
                newx = cbind(newx, mx)
            }
            if(NCOL(newx)==1) newx <- matrix(newx,ncol=1)
            rownames(newx) = rownames(counts)
            colnames(newx) = grps
            counts = newx
        }

        ## if too many samples (like scRNA-seq do subsampling...)
        if(ncol(counts) > 500) {
            kk <- sample(ncol(counts),400,replace=TRUE)
            counts <- counts[,kk,drop=FALSE]
            subtt=c(subtt,"random subset")
        }
        colnames(counts) <- substring(colnames(counts),1,24)

        gset <- list()
        gg = ngs$genes[rownames(counts),]$gene_name
        tt = ngs$genes[rownames(counts),]$gene_title
        g1 <- gg[grep("^rpl|^rps",gg,ignore.case=TRUE)]
        g2 <- gg[grep("^mrpl|^mrps",gg,ignore.case=TRUE)]
        g3 <- gg[grep("^MT-",gg,ignore.case=TRUE)]
        g4 <- gg[grep("mitochondr",tt,ignore.case=TRUE)]
        gset[["Ribosomal (RPL/RPS)"]] = g1
        gset[["Mitochondrial ribosomal (MRPL/MRPS)"]] = g2
        gset[["Mitochondrial (MT)"]] = g3
        gset[["Other mitochondrial"]] = setdiff(g4,g3)
        jj <- grep("mitochondr|ribosom",names(FAMILIES),invert=TRUE,ignore.case=TRUE)
        gset.other <- lapply( FAMILIES[jj], function(x) setdiff(x, c(g1,g2,g3,g4)))
        gset <- c(gset, gset.other)
        gset <- gset[grep("<all>",names(gset),invert=TRUE)]
        gset <- gset[sapply(gset,length) > 10]

        ## Counts per samples, by category
        total.counts = Matrix::colSums(counts,na.rm=TRUE)
        summed.counts = t(sapply(gset, function(f)
            Matrix::colSums(counts[which(gg %in% f),,drop=FALSE], na.rm=TRUE)))
        avg.counts   = t(sapply(gset, function(f)
            Matrix::colMeans(counts[which(gg %in% f),,drop=FALSE], na.rm=TRUE)))
        prop.counts = 100 * t(t(summed.counts) / total.counts)

        head(sort(rowSums(prop.counts,na.rm=TRUE),decreasing=TRUE),10)
        head(sort(rowSums(avg.counts,na.rm=TRUE),decreasing=TRUE),10)
        ##jj <- order(-apply(avg.counts,1,sd,na.rm=TRUE))

        jj <- head(order(-rowSums(prop.counts,na.rm=TRUE)),6)
        prop.counts <- prop.counts[jj,,drop=FALSE]
        jj <- head(order(-rowSums(avg.counts,na.rm=TRUE)),6)
        avg.counts  <- avg.counts[jj,,drop=FALSE]
        sorting="no"
        if(sorting=="decr") {
            total.counts <- sort(total.counts, decreasing=TRUE)
            prop.counts  <- prop.counts[,order(-colMeans(prop.counts))]
            avg.counts   <- avg.counts[,order(-colMeans(avg.counts))]
            counts <- counts[,order(-colMeans(counts))]
        }
        if(sorting=="inc") {
            total.counts <- sort(total.counts, decreasing=FALSE)
            prop.counts  <- prop.counts[,order(colMeans(prop.counts))]
            avg.counts   <- avg.counts[,order(colMeans(avg.counts))]
            counts <- counts[,order(colMeans(counts))]
        }

        ss <- names(total.counts)
        prop.counts <- prop.counts[,ss,drop=FALSE]
        avg.counts <- avg.counts[,ss,drop=FALSE]
        counts <- counts[,ss,drop=FALSE]

        log2counts <- log2(1 + counts)
        ##log2counts[which(log2counts==0)] <- NA
        jj <- sample(nrow(counts),100,replace=TRUE)
        jj <- sample(nrow(counts),1000,replace=TRUE)

        ## create the plots
        par(mar=c(3,3,3,3), mgp=c(2.4,0.7,0), oma=c(1,1,1,1)*0.2 )
        cx1=1
        if(length(total.counts)>12) cx1=0.9
        if(length(total.counts)>30) cx1=0.8
        if(length(total.counts)>50) cx1=0.7
        if(length(total.counts)>80) cx1=0.6

        if(1) {
            names(total.counts) <- substring(names(total.counts),1,30)
            colnames(log2counts) <- substring(colnames(log2counts),1,30)
            colnames(prop.counts) <- substring(colnames(prop.counts),1,30)
            colnames(avg.counts) <- substring(colnames(avg.counts),1,30)
        }

        res = list(total.counts=total.counts,subtt=subtt, cx1=cx1,
                   log2counts=log2counts, jj=jj, prop.counts=prop.counts,
                   avg.counts=avg.counts)

    })

    dataview_counts_caption = "<b>Counts distribution</b>. Plots associated with the counts, abundance or expression levels across the samples/groups.  <b>(a)</b> Total counts per sample or average per group.  <b>(b)</b> Distribution of total counts per sample/group. The center horizontal bar correspond to the median.  <b>(c)</b> Histograms of total counts distribution per sample/group. <b>(d)</b> Abundance of major gene types per sample/group. <b>(e)</b> Average count by gene type per sample/group."

    output$countsUI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.04,1,1),
            height = fullH,
            shiny::div(shiny::HTML(dataview_counts_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(1,1,1), id = "counts_tab_row1", height=rowH,
                plotWidget(ns("counts_tab_barplot")),
                plotWidget(ns("counts_tab_boxplot")),
                plotWidget(ns("counts_tab_histplot"))
            ),
            shiny::fillRow(
                flex = c(1,1), id = "counts_tab_row2", height=rowH,
                plotWidget(ns("counts_tab_abundanceplot")),
                plotWidget(ns("counts_tab_average_countplot"))
            )
        )
    })
    ##dragula(c("counts_tab_row1","counts_tab_row2"))

    ##================================================================================
    ##======================  Raw counts/abundance table =============================
    ##================================================================================

    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>grouped</code>'
    menu_options='<code>Options</code>'

    data_rawdataTable_text = paste0('Under the <strong>gene table </strong>, the average expression values of genes across the groups can be read. The samples (or cells) can be ungrouped by unclicking the ',menu_grouped, ' in the main <i>Options</i> to see the exact expression values per sample (or cell).', 'The genes in the table are ordered by the correlation (<b>rho</b> column) with respect to the gene selected by users from the ',dropdown_search_gene, ' setting. <b>SD</b> column reports the standard deviation of expression across samples (or cells).')

    shiny::observeEvent( input$data_type, {
        ngs = inputData()
        if(input$data_type %in% c("counts","CPM")) {
            pp <- rownames(ngs$counts)
        } else {
            ## log2CPM
            pp <- rownames(ngs$X)
        }
        ## levels for sample filter
        genes <- sort(ngs$genes[pp,]$gene_name)
        sel = genes[1]  ## most var gene
        shiny::updateSelectizeInput(session,'search_gene', choices=genes, selected=sel, server=TRUE)

    })

    data_rawdataTable.RENDER <- shiny::reactive({
        ## get current view of raw_counts
        ngs = inputData()
        shiny::req(ngs)
        shiny::req(input$data_groupby)

        dbg("[data_rawdataTable.RENDER] reacted")

        pp <- rownames(ngs$X)
        if(input$data_type=="counts") {
            ##x <- ngs$counts[pp,]
            x <- ngs$counts
        } else if(input$data_type=="CPM") {
            ##x <- 2**ngs$X
            ##x <- edgeR::cpm(ngs$counts[pp,], log=FALSE)
            x <- edgeR::cpm(ngs$counts, log=FALSE)
        } else {
            ## log2CPM
            x <- ngs$X
        }
        x0=x

        ##------------------ select samples
        dbg("[data_rawdataTable.RENDER] select samples")
        samples <- colnames(ngs$X)
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        samples <- intersect(colnames(x),samples)
        x <- x[,samples,drop=FALSE]

        gene = "CD3E"
        gene = "CCR6"
        gene = input$search_gene
        if(is.null(gene) | gene=="" | is.na(gene)) return(NULL)
        ##xgene = sub(".*:","",rownames(x))

        ## Quickly (?) calculated correlation to selected gene
        dbg("[data_rawdataTable.RENDER] calculate rho")
        rho = sdx = avg = NULL
        if(input$data_type == "logCPM") {
            k=1
            logx <- ngs$X[rownames(x),]
            xgenes <- ngs$genes[rownames(x),"gene_name"]
            k <- which(xgenes==gene)
            rho = cor( t(logx[,samples]), logx[k,samples], use="pairwise")[,1]
            rho = round(rho[rownames(x)], digits=3)
            sdx = round(apply(logx[,samples],1,sd),digits=3)
        }
        avg <- round(rowMeans(x),digits=3)

        ##if(input$data_sampling=="grouped") {
        ##do.grouped <- input$data_grouped
        grpvar = "group"
        group <- NULL
        grpvar <- input$data_groupby

        if(grpvar %in% colnames(ngs$Y)) {
            group = ngs$Y[colnames(x),grpvar]
        }
        if(length(samples)>500 && grpvar=="<ungrouped>") {
            ##grpvar="group"
            group <- ngs$model.parameters$group
        }
        do.grouped <- (grpvar!="<ungrouped>")
        if(do.grouped && !is.null(group) ) {
            ##group = ngs$Y[colnames(x),grpvar]
            allgroups = sort(unique(group))
            newx = c()
            for(gr in allgroups) {
                mx = rowMeans(x[,which(group==gr),drop=FALSE],na.rm=TRUE)
                newx = cbind(newx, mx)
            }
            rownames(newx) = rownames(x)
            ##colnames(newx) = paste0("[",allgroups,"]")
            colnames(newx) = paste0("avg.",allgroups,"")
            x = newx
        }

        x = round( as.matrix(x), digits=3)
        x95 = quantile(as.vector(x0[which(x0>0)]),probs=0.95)
        x99 = quantile(as.vector(x0[which(x0>0)]),probs=0.99)

        if(NCOL(x)==0 || nrow(x)==0) return(NULL)

        dbg("[data_rawdataTable.RENDER] create dataframe")
        ##rownames(x) = sub(".*:","",rownames(x))
        xgenes <- ngs$genes[rownames(x),"gene_name"]
        gene.title <- GENE.TITLE[toupper(xgenes)]
        gene.title <- substring(gene.title,1,50)
        if(is.null(rho)) {
            x = data.frame( gene=xgenes, title=gene.title,
                           AVG=avg,
                           as.matrix(x), check.names=FALSE)
        } else {
            x = data.frame( gene=xgenes, title=gene.title,
                           rho=rho, SD=sdx, AVG=avg,
                           as.matrix(x), check.names=FALSE)
        }
        x = x[order(x$gene),,drop=FALSE]

        ## put selected gene always on top
        j1 <- which(x$gene==gene)
        jj <- c(j1, setdiff(1:nrow(x),j1))
        x <- x[jj,,drop=FALSE]

        if(ncol(x) > 100) {
            max.row <- 1e5 / ncol(x)
            max.row <- 100*ceiling(max.row/100)
            max.row
            x <- head(x, max.row)
        }
        numcols <- grep('gene|title',colnames(x),value=TRUE,invert=TRUE)

        DT::datatable( x, rownames=FALSE,
                      class = 'compact cell-border stripe hover',
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip',
                          ##pageLength = 60,##  lengthMenu = c(20, 30, 40, 60, 100, 250),
                          scroller=TRUE, scrollX = TRUE, scrollY = tabH,
                          deferRender=TRUE
                      )  ## end of options.list
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
                DT::formatStyle(numcols,
                                background = DT::styleColorBar(c(0,x99), 'lightblue'),
                                ##background = color_from_middle(x99, 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')
    }) %>%
    bindCache(input$search_gene, input$data_type, input$data_groupby)

    data_rawdataTable_caption = "<b>Gene table.</b> The table shows the gene expression values per sample, or average expression values across the groups. The column 'rho' reports the correlation with the gene selected in 'Search gene' in the left side bar."

    data_rawdataTable <- shiny::callModule(
        tableModule, "data_rawdataTable",
        func = data_rawdataTable.RENDER,
        title = "Gene expression table",
        filename = "counts.csv",
        info.text = data_rawdataTable_text
    )

    output$genetableUI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.025,1),
            height = fullH,
            shiny::div(shiny::HTML(data_rawdataTable_caption), class="caption"),
            shiny::br(),
            tableWidget(ns("data_rawdataTable"))
        )
    })

    ##================================================================================
    ##================================= Samples ======================================
    ##================================================================================

    data_phenoHeatmap.RENDER <- shiny::reactive({
        ngs = inputData()
        shiny::req(ngs)
        dbg("[data_phenoHeatmap.RENDER] reacted")

        annot <- ngs$samples
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        annot <- annot[samples,,drop=FALSE]
        annot.ht <- ifelse( ncol(annot) > 10, 5, 6)
        annot.ht <- ifelse( ncol(annot) > 20, 4, annot.ht)
        annot.ht <- ifelse( ncol(annot) > 30, 3, annot.ht)

        do.clust <- input$data_phenoclustsamples
        plt <- pgx.plotPhenotypeMatrix0(
            annot, annot.ht=annot.ht, cluster.samples=do.clust)
        ## plt <- plt %>% plotly::config(displayModeBar = FALSE)
        dbg("[data_phenoHeatmap.RENDER] reacted] done!")
        plt
    })

    data_phenoHeatmap_opts <- shiny::tagList(
        shinyBS::tipify( shiny::checkboxInput(ns('data_phenoclustsamples'),'cluster samples',TRUE),
               "Cluster samples.", placement="top")
    )

    data_phenoHeatmap_caption = "<b>Phenotype clustering.</b> Clustered heatmap of sample information (i.e. phenotype data)."
    data_phenoHeatmap_info = "<b>Phenotype clustering.</b> Clustered heatmap of sample information (i.e. phenotype data). Column ordering has been performed using hierarchical clustering on a one-hot encoded matrix."

    shiny::callModule(
        plotModule,
        id = "data_phenoHeatmap", label="a",
        func = data_phenoHeatmap.RENDER,
        func2 = data_phenoHeatmap.RENDER,
        ## plotlib = "iheatmapr",
        title = "Phenotype clustering",
        info.text = data_phenoHeatmap_info,
        options = data_phenoHeatmap_opts,
        height = c(360,600), width = c('auto',1200),
        res=c(68,75), pdf.width=10, pdf.height=6,
        add.watermark = WATERMARK
    )

    data_phenotypeAssociation.RENDER <- shiny::reactive({

        ngs = inputData()
        shiny::req(ngs)
        dbg("[data_phenotypeAssociation.RENDER] reacted")
        annot <- ngs$samples
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        annot <- annot[samples,,drop=FALSE]
        pq <- pgx.testPhenoCorrelation(annot, plot=TRUE)
        dbg("[data_phenotypeAssociation.RENDER] done")
    })

    data_phenotypeAssociation_opts <- shiny::tagList(
        shinyBS::tipify( shiny::checkboxInput(ns('data_phenoclustsamples'),'cluster samples',TRUE),
               "Cluster samples.", placement="top")
    )

    data_phenotypeAssociation_caption = "<b>Phenotype association matrix.</b> Clustered heatmap of phenotype association. The values corresponds to the -log10(p) value of the corresponding statistical test between two phenotype variables. A higher value corresponds to stronger 'correlation'."
    data_phenotypeAssociation_info = "<b>Phenotype clustering.</b> Clustered heatmap of sample information (i.e. phenotype data). The values corresponds to the -log10(p) value of the corresponding statistical test between two phenotype variables. A higher value corresponds to stronger 'correlated' variables. For discrete-discrete pairs the Fisher's exact test is used. For continuous-discrete pairs, the Kruskal-Wallis test is used. For continuous-continous pairs, Pearson's correlation test is used."

    shiny::callModule(
        plotModule,
        id = "data_phenotypeAssociation", label="b",
        func = data_phenotypeAssociation.RENDER,
        func2 = data_phenotypeAssociation.RENDER,
        ## plotlib = "iheatmapr",
        info.text = data_phenotypeAssociation_info,
        title = "Phenotype association",
        ##info.text = "Sample information table with information about phenotype of samples.",
        options = data_phenotypeAssociation_opts,
        height = c(360,700), width = c('auto',900), res=c(72,75),
        pdf.width=8, pdf.height=6,
        add.watermark = WATERMARK
    )

    data_sampleTable.RENDER <- shiny::reactive({
        ## get current view of raw_counts
        ngs = inputData()
        shiny::req(ngs)

        dbg("[data_sampleTable.RENDER] reacted")
        
        dt <- NULL
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        dt <- ngs$samples[samples,,drop=FALSE]
        DT::datatable( dt,
                      class = 'compact cell-border stripe hover',
                      rownames = TRUE,
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip',
                          scroller=TRUE, scrollX = TRUE, scrollY = 190,
                          deferRender=TRUE
                      )) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')         
    }) %>%
    bindCache(input$data_samplefilter)
    
    data_sampleTable_caption="<b>Sample information table.</b> Phenotype information about the samples. Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."
    data_sampleTable_info = "<b>Sample information table.</b> Phenotype information about the samples. Phenotype variables starting with a 'dot' (e.g. '.cell cycle' and '.gender' ) have been estimated from the data."

    data_sampleTable <- shiny::callModule(
        tableModule, "data_sampleTable", label="c",
        func = data_sampleTable.RENDER,
        func2 = data_sampleTable.RENDER,
        title = "Sample information",
        filename = "samples.csv",
        info.text = data_sampleTable_info,
        height = c(280,750), width=c('auto',1280)
    )

    sampletableUI_caption <- paste(
        ## "<b>Phenotype clustering and sample information table.</b>",
        "<b>(a)</b>",data_phenoHeatmap_caption,
        "<b>(b)</b>",data_phenotypeAssociation_caption,
        "<b>(c)</b>",data_sampleTable_caption
    )

    output$sampletableUI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.04,1.2,1),
            height = fullH,
            shiny::div(shiny::HTML(sampletableUI_caption), class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(2,0.07,1),
                shiny::div(plotWidget(ns("data_phenoHeatmap")), style="overflow-y: auto;"),
                shiny::br(),
                plotWidget(ns("data_phenotypeAssociation"))
            ),
            tableWidget(ns("data_sampleTable"))
        )
    })

    ##================================================================================
    ##================================= CONTRASTS ====================================
    ##================================================================================

    data_contrastTable.RENDER <- shiny::reactive({
        ## get current view of raw_counts
        ngs = inputData()
        shiny::req(ngs)

        dbg("[data_contrastTable.RENDER] reacted")

        ##if(is.null(input$data_samplefilter)) return(NULL)
        dt <- NULL
        samples <- selectSamplesFromSelectedLevels(ngs$Y, input$data_samplefilter)
        names(ngs$model.parameters)
        if(input$data_ctbygroup=="group") {
            ct <- ngs$model.parameters$contr.matrix
            ##kk <- which(rownames(ct) %in% ngs$samples[samples,"group"])
            kk <- which(rownames(ct) %in% ngs$model.parameters$group[samples])
            dt <- ct[kk,,drop=FALSE]
        } else {
            dt <- ngs$model.parameters$exp.matrix[samples,,drop=FALSE]
        }
        dt <- sign(dt)
        colnames(dt) <- sub("[_. ]vs[_. ]","\nvs ",colnames(dt))
        dt[dt==0] <- NA

        DT::datatable( dt,
                      class = 'compact cell-border stripe hover',
                      rownames = TRUE,
                      extensions = c('Buttons','Scroller'),
                      selection = list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip',
                          scroller=TRUE, scrollX = TRUE, scrollY = tabH,
                          deferRender=TRUE,
                          autoWidth = TRUE
                      )) %>%
            DT::formatStyle(0, target='row', fontSize='12px', lineHeight='70%') %>%
                DT::formatStyle(colnames(dt),
                                background = color_from_middle(c(-1,1), 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%',
                                backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center')

    })

    data_contrastTable_info = "<b>Contrast table.</b> Table summarizing the contrasts of all comparisons. Here, you can check which samples belong to which groups for the different comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison."

    data_contrastTable_caption = "<b>Contrast table.</b> summarizing the contrasts of all comparisons. Non-zero entries '+1' and '-1' correspond to the group of interest and control group, respectively. Zero or empty entries denote samples not use for that comparison."

    data_contrastTable_opts = shiny::tagList(
        shinyBS::tipify( shiny::radioButtons(ns('data_ctbygroup'),
                             "Show by:", choices=c("sample","group")),
               "Show contrasts by group or by samples.",
               placement="right", options = list(container = "body"))
    )

    data_contrastTable <- shiny::callModule(
        tableModule, "data_contrastTable",
        func = data_contrastTable.RENDER,
        options = data_contrastTable_opts,
        title = "Contrast table",
        filename = "contrasts.csv",
        info.text = data_contrastTable_info
        ##caption = data_contrastTable_caption
    )

    output$contrasttableUI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.03,1), height = fullH,
            shiny::div(shiny::HTML(data_contrastTable_caption),class="caption"),
            shiny::br(),
            tableWidget(ns("data_contrastTable"))
        )
    })

    ##================================================================================
    ## Resource info (dev)
    ##================================================================================

    datatable_timings.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)

        dbg("[datatable_timings.RENDER] reacted")

        ##if(is.null(ngs$timings)) return(NULL)
        D <- data.frame()
        if(!is.null(ngs$timings)) {
            D <- round(ngs$timings[,1:3],digits=3)
            D <- apply(D, 2, function(x) tapply(x, rownames(D), sum))
            catg <- gsub("^\\[|\\].*","",rownames(D))
            metd <- gsub("^.*\\]","",rownames(D))
            D <- data.frame(category=catg, method=metd, D)
        }
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='tp', pageLength = 100),
                      class = 'compact cell-border stripe hover' ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    datatable_timings_text = 'The <b>timings</b> table reports more detailed information about the object dimensions, object sizes and execution times of the methods.'

    datatable_timings <- shiny::callModule(
        tableModule, "datatable_timings",
        func = datatable_timings.RENDER,
        info.text = datatable_timings_text,
        options = NULL, title='Timings'
    )

    datatable_objectdims.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)

        dims1 <- lapply( ngs, dim)
        lens <- sapply( ngs, length)
        dims2 <- t(sapply( ngs[which(!sapply(dims1,is.null)) ], dim))
        kk <- which(sapply(dims1,is.null))
        dims2 <- rbind(dims2, cbind(lens[kk],0))
        colnames(dims2) = c("nrows","ncols")
        D = data.frame( object=rownames(dims2), dims2, check.names=FALSE)
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='t', pageLength = 50),
                      class = 'compact cell-border stripe hover') %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    datatable_objectdims_text = 'This table provides details about the data dimensions of objects.'

    datatable_objectdims <- shiny::callModule(
        tableModule, "datatable_objectdims",
        func = datatable_objectdims.RENDER,
        info.text = datatable_objectdims_text,
        options = NULL, title='Object dimensions'
    )

    datatable_objectsize.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs)
        objsize <- sapply(ngs,object.size)
        objsize <- round( objsize/1e6, digits=2)
        D = data.frame( object=names(ngs), "size.Mb"=objsize, check.names=FALSE)
        DT::datatable( D, rownames=FALSE,
                      options = list(dom='t', pageLength = 50),
                      class = 'compact cell-border stripe hover') %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%')
    })

    datatable_objectsize_text = "This table provides information about  about the memory sizes of objects"

    datatable_objectsize <- shiny::callModule(
        tableModule, "datatable_objectsize",
        func = datatable_objectsize.RENDER,
        options = NULL, title='Object sizes',
        info.text = datatable_objectsize_text
        ## caption = datatable_objectsize_caption
    )

    resourceinfo_caption="<b>Resource information.</b> Details about the execution times of the methods, dimensions and memory sizes of objects."

    output$resourceinfoUI <- shiny::renderUI({
        shiny::fillCol(
            flex = c(NA,0.02,1),
            height = fullH,
            shiny::div(shiny::HTML(resourceinfo_caption),class="caption"),
            shiny::br(),
            shiny::fillRow(
                flex = c(5,1, 2,1, 1.5, 2), ## width = 600,
                tableWidget(ns("datatable_timings")),
                shiny::br(),
                tableWidget(ns("datatable_objectdims")),
                shiny::br(),
                tableWidget(ns("datatable_objectsize")),
                shiny::br()
            )
        )
    })
}