##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

WordCloudBoard <- function(id, pgx)
{
  moduleServer(id, function(input, output, session)
  {

    ns <- session$ns ## NAMESPACE
    # fullH = 750
    rowH = 660  ## row height of panel
    # tabH = 200  ## row height of panel
    # tabH = '70vh'  ## row height of panel

    wc_infotext = paste("This module performs WordCloud analysis or 'keyword enrichment', i.e. it computes the enrichment of keywords for the contrasts. Frequently appearing words in the top ranked gene sets form an unbiased description of the contrast.
<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
")

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================

    shiny::observeEvent( input$wc_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>WordCloud Analysis Board</strong>"),
            shiny::HTML(wc_infotext),
            easyClose = TRUE, size="l" ))
    })

    shiny::observe({
        shiny::req(pgx$gset.meta)
        ct <- names(pgx$gset.meta$meta)
        ct <- sort(ct)
        shiny::updateSelectInput(session, "wc_contrast", choices=ct )
    })

    ##---------------------------------------------------------------
    ##------------- Functions for WordCloud -------------------------
    ##---------------------------------------------------------------

    getWordFreqResults <- shiny::reactive({
        shiny::req(pgx$gset.meta)
        if("wordcloud" %in% names(pgx)) {
            res <- pgx$wordcloud
        } else {
            dbg("**** CALCULATING WORDCLOUD ****\n")
            progress <- shiny::Progress$new()
            res <- pgx.calculateWordCloud(pgx, progress=progress, pg.unit=1)
            on.exit(progress$close())

            ## save in object??
            pgx[["wordcloud"]] <- res
        }
        return(res)
    })

    getCurrentWordEnrichment <- shiny::reactive({

        res <- getWordFreqResults()
        shiny::req(res, input$wc_contrast)

        contr=1
        contr <- input$wc_contrast
        gsea1 <- res$gsea[[ contr ]]
        topFreq <- data.frame( gsea1, tsne=res$tsne, umap=res$umap)
        topFreq <- topFreq[order(-topFreq$NES),]

        return(topFreq)
    })

    wordtsne.RENDER <- shiny::reactive({

        topFreq <- getCurrentWordEnrichment()

        df <- topFreq
        klr = ifelse( df$padj<=0.05, "red", "grey")
        ps1 = 0.5 + 3*(1-df$padj)*(df$NES/max(df$NES))**3

        ## label top 20 words
        df$label <- rep(NA, nrow(df))
        jj <- head(order(-abs(df$NES)),20)
        df$label[jj] <- as.character(df$word[jj])
        cex=1

        if(input$wordtsne_algo=="tsne") {
            p <- ggplot2::ggplot(df, ggplot2::aes(tsne.x, tsne.y, label=label))
        } else {
            p <- ggplot2::ggplot(df, ggplot2::aes(umap.x, umap.y, label=label))
        }
        p <- p +
            ggplot2::geom_point( size=cex*ps1, color=klr) +
            ggrepel::geom_text_repel(size=4*cex) +
            ggplot2::theme_bw() +
            ggplot2::theme( axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank())

        return(p)
    })


    plotWcActmap <- function(score, normalize, nterm, nfc) {

        ## reduce score matrix
        score = score[head(order(-rowMeans(score[,]**2)),nterm),,drop=FALSE] ## max terms
        score = score[,head(order(-colSums(score**2)),nfc),drop=FALSE] ## max comparisons/FC

        cat("<wordcloud_actmap> dim(score)=",dim(score),"\n")
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        jj=1;ii=1:nrow(score)
        ii <- hclust(d1)$order
        jj <- hclust(d2)$order
        score <- score[ii,jj,drop=FALSE]

        colnames(score) = substring(colnames(score),1,30)
        rownames(score) = substring(rownames(score),1,50)
        colnames(score) <- paste0(colnames(score)," ")
        cex2=0.85

        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,2,0,1))

        score2 <- score
        if(normalize) score2 <- t(t(score2) / apply(abs(score2),2,max)) ## normalize cols???
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**3   ## fudging
        bmar <- 0 + pmax((50 - nrow(score2))*0.25,0)
        corrplot::corrplot( score2, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex = 0.9*cex2, tl.col = "grey20", tl.srt = 90,
                 mar=c(bmar,0,0,0) )
    }

    wordcloud_actmap.RENDER <- shiny::reactive({

        res <- getWordFreqResults()
        score <- sapply(res$gsea, function(x) x$NES)
        rownames(score) <- res$gsea[[1]]$word

        plotWcActmap(
            score = score,
            normalize = input$wc_normalize,
            nterm = 50,
            nfc = 20
        )

    })

    wordcloud_actmap.RENDER2 <- shiny::reactive({

        res <- getWordFreqResults()
        score <- sapply(res$gsea, function(x) x$NES)
        rownames(score) <- res$gsea[[1]]$word

        plotWcActmap(
            score = score,
            normalize = input$wc_normalize,
            nterm = 50,
            nfc = 100
        )

    })

    ##-------- Activation map plotting module
    wordcloud_actmap.opts = shiny::tagList()
    shiny::callModule(
        plotModule,
        id="wordcloud_actmap",
        func = wordcloud_actmap.RENDER,
        func2 = wordcloud_actmap.RENDER2,
        title = "Activation matrix", label="d",
        info.text = "The <strong>Activation Matrix</strong> visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = wordcloud_actmap.opts,
        pdf.width=6, pdf.height=10,
        height = c(rowH,750), width=c("100%",1400),
        res=72,
        add.watermark = WATERMARK
    )

    ## ================================================================================
    ## =========================== MODULES ============================================
    ## ================================================================================

    WATERMARK <- FALSE

    # Enrichment plots

    wordcloud_plot_enrichment_server(
      "gseaplots",
      getWordFreqResults = getWordFreqResults,
      getCurrentWordEnrichment = getCurrentWordEnrichment,
      wordcloud_enrichmentTable = wordcloud_enrichmentTable,
      watermark = WATERMARK
    )

    # Word cloud

    wordcloud_plot_wordcloud_server(
      "wordcloud",
      getCurrentWordEnrichment = getCurrentWordEnrichment,
      watermark = WATERMARK
    )

    # Word t-SNE

    wordcloud_plot_wordtsne_server(
      "wordtsne",
      getCurrentWordEnrichment = getCurrentWordEnrichment,
      watermark = WATERMARK
    )

    # Enrichment table

    wordcloud_enrichmentTable <- wordcloud_table_enrichment_server(
      "wordcloud_enrichmentTable",
      getCurrentWordEnrichment = getCurrentWordEnrichment
    )

    # Leading-edge table

    wordcloud_leadingEdgeTable <- wordcloud_table_leading_edge_server(
      "wordcloud_leadingEdgeTable",
      pgx = pgx,
      wc_contrast = shiny::reactive(input$wc_contrast),
      wordcloud_enrichmentTable = wordcloud_enrichmentTable,
      getCurrentWordEnrichment = getCurrentWordEnrichment
    )

  })
}
