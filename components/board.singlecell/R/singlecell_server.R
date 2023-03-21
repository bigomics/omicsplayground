##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


SingleCellBoard <- function(id, inputData) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    fullH <- 750 ## full height of panel
    imgH <- 680 ## row height of panel
    tabH <- 200 ## row height of panel

    infotext <-
      "The <strong>Cell Profiling Board</strong> infers the type of cells using computational deconvolution methods and reference datasets from the literature. Currently, we have implemented a total of 8 methods and 9 reference datasets to predict immune cell types (4 datasets), tissue types (2 datasets), cell lines (2 datasets) and cancer types (1 dataset). However, we plan to expand the collection of methods and databases and to infer other cell types.

<br><br>The <strong>Proportions tab</strong> visualizes the interrelationships between two categorical variables (so-called cross tabulation). Although this feature is very suitable for a single-cell sequencing data, it provides useful information about the proportion of different cell types in samples obtained by the bulk sequencing method.

<br><br>For each combination of gene pairs, the platform can generate a cytometry-like plot of samples under the <strong>Cytoplot</strong> tab. The aim of this feature is to observe the distribution of samples in relation to the selected gene pairs. For instance, when applied to single-cell sequencing data from immunological cells, it can mimic flow cytometry analysis and distinguish T helper cells from the other T cells by selecting the CD4 and CD8 gene combination.

<br><br>The <strong>Markers</strong> section provides potential marker genes, which are the top N=36 genes with the highest standard deviation within the expression data across the samples. For every gene, it produces a t-SNE plot of samples, with samples colored in red when the gene is overexpressed in corresponding samples. Users can also restrict the marker analysis by selecting a particular functional group in which genes are divided into 89 groups, such as chemokines, transcription factors, genes involved in immune checkpoint inhibition, and so on.

<br><br>It is also possible to perform a copy number variation analysis under the <strong>CNV tab</strong>. The copy number is estimated from gene expression data by computing a moving average of the relative expression along the chromosomes. CNV generates a heatmap of samples versus chromosomes, where samples can be annotated further with a phenotype class provided in the data."


    ## ================================================================================
    ## ======================= OBSERVE FUNCTIONS ======================================
    ## ================================================================================

    shiny::observeEvent(input$infotext, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Single Cell Board</strong>"),
        shiny::HTML(infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    ## update filter choices upon change of data set
    shiny::observe({
      ngs <- inputData()
      shiny::req(ngs)
      ## levels for sample filter
      levels <- getLevels(ngs$Y)
      shiny::updateSelectInput(session, "samplefilter", choices = levels)

      ## update cluster methods if available in object
      if ("cluster" %in% names(ngs)) {
        clustmethods <- names(ngs$cluster$pos)
        clustmethods <- c("default", clustmethods)
        shiny::updateSelectInput(session, "clustmethod",
          choices = clustmethods
        )
      }
    })

    shiny::observe({
      ngs <- inputData()
      shiny::req(ngs)
      refsets <- "LM22"
      refsets <- sort(names(ngs$deconv))
      refsel <- unique(c(grep("LM22", refsets, value = TRUE), refsets))[1]
      shiny::updateSelectInput(session, "refset", choices = refsets, selected = refsel)
      shiny::updateSelectInput(session, "refset2", choices = refsets, selected = refsel)

      ## dcmethods <- names(ngs$deconv[[1]])
      ## dcsel <- intersect(c("meta.prod","meta"),dcmethods)[1]
      ## shiny::updateSelectInput(session, "dcmethod", choices=dcmethods, selected=dcsel)
      ## shiny::updateSelectInput(session, "dcmethod2", choices=dcmethods, selected=dcsel)

      grpvars <- c("<ungrouped>", colnames(ngs$samples))
      sel <- grpvars[1]
      if (ncol(ngs$X) > 30) sel <- grpvars[2]
      shiny::updateSelectInput(session, "group2", choices = grpvars, selected = sel)
    })

    shiny::observeEvent(input$refset, {
      shiny::req(input$refset)
      ngs <- inputData()
      dcmethods <- names(ngs$deconv[[input$refset]])
      dcsel <- intersect(c("meta.prod", "meta"), dcmethods)[1]
      shiny::updateSelectInput(session, "dcmethod", choices = dcmethods, selected = dcsel)
    })

    shiny::observeEvent(input$refset2, {
      shiny::req(input$refset2)
      ngs <- inputData()
      dcmethods <- names(ngs$deconv[[input$refset2]])
      dcsel <- intersect(c("meta.prod", "meta"), dcmethods)[1]
      shiny::updateSelectInput(session, "dcmethod2", choices = dcmethods, selected = dcsel)
    })

    shiny::observe({
      ngs <- inputData()
      ## if(is.null(ngs)) return(NULL)
      shiny::req(ngs)

      ## if(is.null(input$crosstaboptions)) return(NULL)
      pheno0 <- grep("group|sample|donor|id|batch", colnames(ngs$samples), invert = TRUE, value = TRUE)
      pheno0 <- grep("sample|donor|id|batch", colnames(ngs$samples), invert = TRUE, value = TRUE)
      kk <- selectSamplesFromSelectedLevels(ngs$Y, input$samplefilter)
      nphenolevel <- apply(ngs$samples[kk, pheno0, drop = FALSE], 2, function(v) length(unique(v)))
      pheno0 <- pheno0[which(nphenolevel > 1)]
      genes <- sort(as.character(ngs$genes$gene_name))
      pheno1 <- c("<cell type>", pheno0) # pheno1 <- c("<cell type>", pheno0)
      genes1 <- c("<none>", genes)
      shiny::updateSelectInput(session, "crosstabvar", choices = pheno1)
      shiny::updateSelectInput(session, "crosstabpheno", choices = pheno1,,selected = pheno1[1])
      shiny::updateSelectizeInput(session, "crosstabgene", choices = genes1, server = TRUE, selected = genes1[2])
    })

    shiny::observe({
      ngs <- inputData()
      shiny::req(ngs, input$mrk_level)

      choices <- names(ngs$families)
      selected <- grep("^CD", choices, ignore.case = TRUE, value = TRUE)[1]
      if (input$mrk_level == "geneset") {
        nn <- sapply(COLLECTIONS, function(k) sum(k %in% rownames(ngs$gsetX)))
        choices <- names(COLLECTIONS)[nn >= 5]
        selected <- grep("HALLMARK", names(COLLECTIONS), ignore.case = TRUE, value = TRUE)
      }
      shiny::updateSelectInput(session, "features", choices = choices, selected = selected)
      shiny::updateSelectInput(session, "mrk_features", choices = choices, selected = selected)
    })

    shiny::observe({
      ngs <- inputData()
      ## if(is.null(ngs)) return(NULL)
      shiny::req(ngs)
      ## just at new data load
      genes <- NULL
      g1 <- g2 <- NULL

      F <- pgx.getMetaFoldChangeMatrix(ngs)$fc
      F <- F[order(-apply(F, 1, sd)), ]
      genes <- rownames(F)
      g1 <- rownames(F)[1]
      g2 <- rownames(F)[2]

      if (length(g1) == 0) g1 <- genes[1]
      if (length(g2) == 0) g2 <- genes[2]

      shiny::updateSelectizeInput(session, "cytovar1", choices = genes, selected = g1, server = TRUE)
      shiny::updateSelectizeInput(session, "cytovar2", choices = genes, selected = g2, server = TRUE)
    })



    # REACTIVE FUNCTIONS #########

    pfGetClusterPositions <- shiny::reactive({ # used by many plots
      ngs <- inputData()
      shiny::req(ngs)

      ## zx <- filtered_matrix1()
      zx <- ngs$X
      kk <- colnames(zx)
      kk <- selectSamplesFromSelectedLevels(ngs$Y, input$samplefilter)
      if (length(kk) == 0) {
        return(NULL)
      }
      zx <- zx[, kk, drop = FALSE]
      zx <- head(zx[order(-apply(zx, 1, sd)), ], 1000)
      zx <- t(scale(t(zx))) ## scale??

      pos <- NULL
      m <- "tsne"
      m <- input$clustmethod
      has.clust <- ("cluster" %in% names(ngs) && m %in% names(ngs$cluster$pos))
      has.clust
      if (!has.clust && m == "pca") {
        pos <- irlba::irlba(zx, nv = 3)$v
        rownames(pos) <- colnames(zx)
      } else if (has.clust) {
        pos <- ngs$cluster$pos[[m]][, 1:2]
      } else {
        pos <- ngs$tsne2d
      }
      dim(pos)
      pos <- pos[colnames(zx), ]
      pos <- scale(pos) ## scale
      colnames(pos) <- paste0("dim", 1:ncol(pos))
      rownames(pos) <- colnames(zx)


      # code snipped from pfGetClusterPositions2, pfGetClusterPositions2 is currently never called

      # dbg("[pfGetClusterPositions2] computing distances and clusters...")
      # dbg("[pfGetClusterPositions2] dim(pos) = ",dim(pos))
      #
      # ##dist = as.dist(dist(pos))
      # dist = 0.001+dist(pos)**2
      #
      # dbg("[pfGetClusterPositions2] creating graph")
      #
      # gr = igraph::graph_from_adjacency_matrix(
      #   1.0/dist, diag=FALSE, mode="undirected")
      #
      # dbg("[pfGetClusterPositions2] cluster louvain")
      #
      # clust <- igraph::cluster_louvain(gr)$membership
      #
      # dbg("pfGetClusterPositions2:: done!")
      # return( list(pos=pos, clust=clust) )

      return(pos)
    })

    # Type mapping (heatmap) reactivity ##########

    getDeconvResults2 <- shiny::reactive({ # used by many functions
      ngs <- inputData()
      shiny::req(ngs)

      method <- "meta"
      method <- input$dcmethod2
      if (is.null(method)) {
        return(NULL)
      }
      shiny::req(input$refset2)

      refset <- "LM22"
      refset <- input$refset2
      if (!("deconv" %in% names(ngs))) {
        return(NULL)
      }
      results <- ngs$deconv[[refset]][[method]]
      ## threshold everything (because DCQ can be negative!!!)
      results <- pmax(results, 0)


      return(results)
    })


    # plots -------------------------------------------------------------------

    singlecell_plot_icpplot_server(
      id = "icpplot",
      inputData = inputData,
      pfGetClusterPositions = pfGetClusterPositions,
      method = shiny::reactive(input$dcmethod),
      refset = shiny::reactive(input$refset),
      lyo = shiny::reactive(input$layout),
      sortby = shiny::reactive(input$sortby)
    )

    singlecell_plot_phenoplot_server(
      id = "phenoplot",
      inputData = inputData,
      pfGetClusterPositions = pfGetClusterPositions
    )

    singlecell_plot_mappingplot_server(
      id = "mappingplot",
      inputData = inputData,
      getDeconvResults2 = getDeconvResults2,
      pfGetClusterPositions = pfGetClusterPositions,
      grpvar = shiny::reactive(input$group2),
      refset = shiny::reactive(input$refset2),
      group = shiny::reactive(input$group2),
      view = shiny::reactive(input$view2)
    )

    singlecell_plot_crosstabPlot_server(
      id = "crosstabPlot",
      inputData = inputData,
      samplefilter = shiny::reactive(input$samplefilter),
      crosstabvar = shiny::reactive(input$crosstabvar),
      pheno = shiny::reactive(input$crosstabpheno),
      gene = shiny::reactive(input$crosstabgene),
      getDeconvResults2 = getDeconvResults2,
      watermark = FALSE
    )

    singlecell_plot_markersplot_server(
      id = "markersplot",
      inputData = inputData,
      pfGetClusterPositions = pfGetClusterPositions,
      mrk_level = shiny::reactive(input$mrk_level),
      mrk_features = shiny::reactive(input$mrk_features),
      mrk_search = shiny::reactive(input$mrk_search),
      mrk_sortby = shiny::reactive(input$mrk_sortby),
      watermark = FALSE
    )

    singlecell_plot_cytoplot_server(
      id = "cytoplot",
      inputData = inputData,
      pfGetClusterPositions = pfGetClusterPositions,
      samplefilter = shiny::reactive(input$samplefilter),
      cytovar1 = shiny::reactive(input$cytovar1),
      cytovar2 = shiny::reactive(input$cytovar2),
      selectSamplesFromSelectedLevels = selectSamplesFromSelectedLevels,
      watermark = FALSE
    )





    # CNV #######

    # getCNAfromExpression <- shiny::reactive({ # Currently not used
    #     ngs <- inputData()
    #     shiny::req(ngs)
    #
    #     dbg("[SingleCellBoard:getCNAfromExpression] calculating CNV with SMA40 ...")
    #
    #     ##source("../R/pgx-cna.R");source("../R/gx-heatmap.r")
    #     shiny::withProgress( message='calculating CNV (sma40)...', value=0.33, {
    #         res <- pgx.CNAfromExpression(ngs, nsmooth=40)
    #     })
    #     return(res)
    # })
    #
    # getCNAfromExpression.inferCNV <- shiny::reactive({ # Currently not used
    #     ngs <- inputData()
    #     shiny::req(ngs)
    #
    #     dbg("[SingleCellBoard:getCNAfromExpression] calculating CNV using inferCNV...")
    #
    #
    #     shiny::withProgress( message='calculating CNV (inferCNV)...', value=0.33, {
    #         res <- pgx.inferCNV(ngs, refgroup=NULL)
    #     })
    #     return(res)
    # })


    # Currently not used Stefan

    # cna.plotFUNC <- shiny::reactive({

    #     ##return(NULL)
    #     ngs <- inputData()
    #     shiny::req(ngs,input$cna_method,input$cna_annotvar,input$cna_orderby)

    #     if(input$cna_method=="inferCNV") {
    #         res <- getCNAfromExpression.inferCNV()
    #         par(mfrow=c(1,1))
    #         grid::grid.raster(res$png)
    #     } else {
    #         res <- getCNAfromExpression()
    #         annotvar=NA
    #         annotvar <- input$cna_annotvar
    #         if(annotvar=="<none>") annotvar <- NULL
    #         order.by <- input$cna_orderby
    #         pgx.plotCNAHeatmap(
    #             ngs, res, annot=annotvar, order.by=order.by,
    #             downsample=10 )
    #     }

    # })

    # cna.opts = shiny::tagList(
    #     withTooltip(shiny::radioButtons(ns("cna_method"),label="Method:", choices=c("sma40","inferCNV"),inline=TRUE),
    #            "Select the computational method for CNV inference. The <tt>sma40</tt> method uses a fast moving average of the relative expression values with a window of 40 genes. <tt>inferCNV</tt> uses the method inferCNV of the Trinity CTAT Project (warning: this method is very slow!)."),
    #     withTooltip(shiny::selectInput(ns("cna_annotvar"),label="Annotate with:", choices=NULL, multiple=FALSE),
    #            "Select what annotation variable to show together with the heatmap", placement = "top"),
    #     ##checkboxGroupInput('cnaoptions','',c("bin20"), inline=TRUE),
    #     ##radioButtons('cnaplottype',NULL,c("image","heatmap","splitmap"), inline=TRUE),
    #     withTooltip(shiny::radioButtons(ns("cna_orderby"),"Order samples by:",c("clust","pc1"), inline=TRUE),
    #            "Select how to order the vertical (sample) axis: clustering or according the loading of the first principal component.")
    # )

    # cnaModule_info = "<strong>Copy number variation (CNV)</strong> analysis. The copy number is estimated from gene expression data by computing a moving average of the relative expression along the chromosomes. CNV generates a heatmap of samples versus chromosomes, where samples can be annotated further with a phenotype class provided in the data."

    # shiny::callModule(
    #     plotModule,
    #     id = "cnaplot",
    #     func = cna.plotFUNC,
    #     func2 = cna.plotFUNC,
    #     options = cna.opts,
    #     title = "Inferred copy number variation (CNV)",
    #     info.text = cnaModule_info,
    #     pdf.width=10, pdf.height=8,
    #     height = c(fullH - 60,700), width = c('100%',1000),
    #     res = 110,
    #     add.watermark = WATERMARK
    # )

    # shiny::observe({
    #     ngs <- inputData()
    #     ##if(is.null(ngs)) return(NULL)
    #     shiny::req(ngs)
    #     ## levels for sample filter

    #     annotvar <- c(colnames(ngs$Y),"<none>")
    #     shiny::updateSelectInput(session, "cna_annotvar", choices=annotvar)

    # })

    # Currently not used Stefan 22.03.22

    # iTALK ######

    #     italk_getResults <- shiny::reactive({
    #         ngs <- inputData()
    #         ## if(is.null(ngs)) return(NULL)
    #         shiny::req(ngs)
    #         shiny::req(input$italk_groups)

    #         db <- iTALK::database
    #         db.genes <- unique(c(db$Ligand.ApprovedSymbol,db$Receptor.ApprovedSymbol))
    #         length(db.genes)
    #         ##genes <- intersect(genes, rownames(ngs$X))
    #         xgenes <- toupper(ngs$genes[rownames(ngs$X),"gene_name"])
    #         db.genes <- intersect(db.genes, xgenes)
    #         length(db.genes)

    #         ## make groups
    #         ph <- "group"
    #         ph <- "cell.type"
    #         ph <- input$italk_groups
    #         ct <- as.character(ngs$samples[,ph])
    #         table(ct)

    #         ##data <- data.frame(cell_type=ct, t(log2(1 + ngs$counts[genes,])))
    #         pp1 <- rownames(ngs$X)[match(db.genes, toupper(xgenes))]
    #         gx <- t(ngs$X[pp1,,drop=FALSE])

    #         colnames(gx) <- db.genes ## UPPERCASE
    #         gx0 <- apply(gx,2,function(x) tapply(x,ct,mean))
    #         dim(gx0)

    #         top_genes <- 50
    #         top_genes <- input$italk_netview_topgenes

    #         colnames(gx) <- toupper(colnames(gx))
    #         data1 <- data.frame(cell_type=ct, gx)
    #         dim(data1)

    #         ## find the ligand-receptor pairs from highly expressed genes
    #         cell_col <- rep(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b','#9ef49a'),99)
    #         ct.names = unique(as.character(data1$cell_type))
    #         cell_col <- cell_col[1:length(ct.names)]
    #         names(cell_col) <- ct.names

    #         comm_type='cytokine'
    #         comm_type <- input$italk_category

    #         mode = "absolute"
    #         ##mode <- input$italk_mode
    #         if(mode=="absolute") {
    #             highly_exprs_genes <- iTALK::rawParse(data1, top_genes=50, stats='mean')
    #             res_cat <- iTALK::FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
    #             dim(res_cat)
    #             xx <- res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs
    #             res_cat <- res_cat[order(xx,decreasing=TRUE),]
    #         } else {
    #             ## contrast <- input$fa_contrast
    #             ## group <- ngs$model.parameters$exp.matrix[,contrast]
    #             ## data1$compare_groups <- group
    #             ## data1 <- data1[which(data1$compare_groups!=0),]
    #             ## ## find DEGenes of regulatory T cells and NK cells between these 2 groups
    #             ## deg_t  <- DEG(data %>% plotly::filter(cell_type=='CD4Tcells'),method='Wilcox',contrast=c(2,1))
    #             ## deg_nk <- DEG(data %>% plotly::filter(cell_type=='NKcells'),method='Wilcox',contrast=c(2,1))
    #             ## ## find significant ligand-receptor pairs and do the plotting
    #             ## res_cat <- FindLR(deg_t,deg_nk,datatype='DEG',comm_type=comm_type)
    #             ## res_cat <- res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
    #         }

    #         gx0.genes <- toupper(colnames(gx0))
    #         is.ligand = (gx0.genes %in% db[,"Ligand.ApprovedSymbol"])
    #         is.receptor = (gx0.genes %in% db[,"Receptor.ApprovedSymbol"])
    #         lr.type = c("L","R","LR")[ 1*is.ligand + 2*is.receptor]
    #         names(lr.type) <- colnames(gx0)
    #         table(lr.type)
    #         dim(gx0)
    #         res <- list( table=res_cat, exprs=gx0, cell_col=cell_col, lr.type=lr.type)
    #         return(res)
    #     })

    #     italk_netview.RENDER <- shiny::reactive({
    #         res <- italk_getResults()
    #         shiny::req(res)
    #         ##if(is.null(res)) return(NULL)
    #         res_cat <- res$table
    #         ## Communication graph
    #         iTALK::NetView(res_cat, col=res$cell_col, vertex.label.cex=1, arrow.width=1, edge.max.width=5)
    #     })

    #     italk_LRPlot.RENDER <- shiny::reactive({
    #         ## Circos plot
    #         res <- italk_getResults()
    #         shiny::req(res)
    #         ##if(is.null(res)) return(NULL)
    #         res_cat <- res$table
    #         if(nrow(res_cat)<1) return(NULL)

    #         ntop=25
    #         ntop = as.integer(input$italk_LRPlot_ntop)
    #         res_top <- head(res_cat,ntop)
    #         iTALK::LRPlot(res_top, datatype='mean count', cell_col=res$cell_col,
    #                       link.arr.lwd = head(res_cat$cell_from_mean_exprs,ntop),
    #                       link.arr.width = head(res_cat$cell_to_mean_exprs,ntop))
    #         comm_type <- shiny::isolate(input$italk_category)
    #         title((paste(comm_type,"genes     ")), line=0.5)
    #     })

    #     italk_heatmap.RENDER <- shiny::reactive({
    #         ## Expression heatmap
    #         ngs <- inputData()
    #         res <- italk_getResults()
    #         shiny::req(ngs,res)
    #         ##if(is.null(res)) return(NULL)

    #         res_cat <- res$table
    #         ntop=50
    #         ntop = as.integer(input$italk_LRPlot_ntop)
    #         res_top <- head(res_cat,ntop)
    #         genes_top <- sort(unique(c(res_top$ligand,res_top$receptor)))
    #         if(length(genes_top)==0) return(NULL)
    #         gx0  <- t(res$exprs[,genes_top,drop=FALSE])
    #         rownames(gx0) <- paste0(rownames(gx0)," (",res$lr.type[rownames(gx0)],")")

    #         par(oma=c(3,2,3,0))
    #         gx.heatmap(gx0, scale="none", mar=c(15,8),
    #                    cexRow=1, cexCol=1.3, col=BLUERED(64),
    #                    key=FALSE, keysize=0.6)
    #     })



    #     italk_LRPlot_info = "The Ligand-Receptor plot visualizes the communication structure of ligand-receptor genes as
    # a circle plot. The width of the arrow represents the expression level/log fold change of the ligand; while the width of arrow head represents the expression level/log fold change of the receptor. Different color and the type of the arrow stands for whether the ligand and/or receptor are upregulated or downregulated. For further information, see iTALK R package (Wang et al., BioRxiv 2019)."

    #     shiny::callModule(
    #         plotModule,
    #         id = "italk_LRPlot",
    #         func = italk_LRPlot.RENDER,
    #         func2 = italk_LRPlot.RENDER,
    #         title = "Ligand-Receptor plot", label="a",
    #         info.text = italk_LRPlot_info,
    #         options = shiny::tagList(
    #             withTooltip( shiny::selectInput(ns("italk_LRPlot_ntop"),"ntop pairs",
    #                                 choices=c(10,15,25,50,75,100),selected=25),
    #                    "Select the maximum number of LR pairs to include in the LR plot.",
    #                    placement="top")
    #         ),
    #         pdf.width=6, pdf.height=8,
    #         height = imgH, res=80,
    #         add.watermark = WATERMARK
    #     )

    #     italk_heatmap_info = "The heatmap visualizes the expression level/log fold change of the ligand/receptor genes. For further information, see iTALK R package (Wang et al., BioRxiv 2019)."

    #     shiny::callModule(
    #         plotModule,
    #         id = "italk_heatmap",
    #         func = italk_heatmap.RENDER,
    #         func2 = italk_heatmap.RENDER,
    #         title = "Expression heatmap", label="b",
    #         info.text = italk_heatmap_info,
    #         options = shiny::tagList(),
    #         pdf.width=6, pdf.height=8,
    #         height = imgH, res=80,
    #         add.watermark = WATERMARK
    #     )

    #     italk_netview_info = "The NetView plot visualizes the communication structure of ligand-receptor genes as a graph. The colors represent different types of cells as a structure and the width of edges represent the strength of the communication. Labels on the edges show exactly how many interactions exist between two types of cells. For further information, see iTALK R package (Wang et al., BioRxiv 2019)."

    #     shiny::callModule(
    #         plotModule,
    #         id = "italk_netview",
    #         func = italk_netview.RENDER,
    #         func2 = italk_netview.RENDER,
    #         title = "NetView", label="c",
    #         info.text = italk_netview_info,
    #         options = shiny::tagList(
    #             withTooltip( shiny::selectInput(ns("italk_netview_topgenes"),"top genes",
    #                                 choices=c(10,25,50,75,100),selected=50),
    #                    "Select the number of topgenes to search for ligand-receptor pairs.",
    #                    placement="top" )
    #         ),
    #         pdf.width=6, pdf.height=8,
    #         height = imgH, res=80,
    #         add.watermark = WATERMARK
    #     )

    #     shiny::observe({
    #         ngs <- inputData()
    #         shiny::req(ngs)
    #         ph <- sort(colnames(ngs$samples))
    #         sel = ph[1]
    #         ct <- grep("cell.fam|cell.type|type|cluster",ph,value=TRUE)
    #         if(length(ct)>0) sel <- ct[1]
    #         shiny::updateSelectInput(session, "italk_groups", choices=ph, selected=sel)
    #     })

    # Currently not used Stefan 22.03.22

    # # Trajectory (dev) ####

    # monocle_getResults <- shiny::reactive({

    #     ngs <- inputData()
    #     shiny::req(ngs)
    #     ##if(is.null(ngs)) return(NULL)

    #     ## Create a Progress object
    #     progress <- shiny::Progress$new()
    #     on.exit(progress$close())
    #     progress$set(message = "Calculating trajectories", value = 0)

    #      ## Step 0: Make Monocle object from ngs #######

    #     NGENES=2000
    #     jj <- head(order(-apply(log(1+ngs$counts),1,sd,na.rm=TRUE)),NGENES)
    #     ngs$counts <- ngs$counts[jj,]
    #     ngs$genes  <- ngs$genes[rownames(ngs$counts),]
    #     ngs$genes$gene_short_name <- ngs$genes$gene_name
    #     ngs$samples$.cluster <- ngs$samples$cluster  ## save to avoid name clash
    #     ngs$samples$cluster <- NULL

    #     X <- ngs$counts
    #     G <- as.data.frame(ngs$genes)
    #     rownames(X) <- toupper(rownames(X))
    #     sum(duplicated(rownames(X)))

    #     ii <- which(!duplicated(rownames(X)))
    #     X <- X[ii,]
    #     X <- X[ii,]
    #     G <- G[ii,]
    #     rownames(G) <- rownames(X)
    #     G$gene_name <- rownames(G)
    #     G$gene_short_name <- rownames(G)

    #     cds <- monocle3::new_cell_data_set(
    #                          expression_data = X,
    #                          cell_metadata = ngs$samples,
    #                          gene_metadata = G)

    #     expression_data = X
    #     cell_metadata = ngs$samples
    #     gene_metadata = G
    #     sce <- SingleCellExperiment::SingleCellExperiment(
    #         list(counts = as(expression_data, "dgCMatrix")),
    #         rowData = gene_metadata,
    #         colData = cell_metadata )

    #     cds <- methods::new("cell_data_set",
    #                         assays = SummarizedExperiment::Assays(
    #                                                            list(counts = as(expression_data,
    #                                                                             "dgCMatrix"))),
    #                         colData = SingleCellExperiment::colData(sce),
    #                         int_elementMetadata = sce@int_elementMetadata,
    #                         int_colData = sce@int_colData,
    #                         int_metadata = sce@int_metadata,
    #                         metadata = sce@metadata,
    #                         NAMES = sce@NAMES,
    #                         elementMetadata = sce@elementMetadata,
    #                         rowRanges = sce@rowRanges)

    #     metadata(cds)$cds_version <- Biobase::package.version("monocle3")
    #     clusters <- stats::setNames( S4Vectors::SimpleList(), character(0))
    #     cds <- monocle3::estimate_size_factors(cds)
    #     cds

    #      ## Step 1: Normalize and pre-process the data ####

    #     progress$inc(0.1, detail = "preprocessing")
    #     ##cds <- monocle3::preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ batch")
    #     cds <- monocle3::preprocess_cds(cds, num_dim = 40)
    #     ##plot_pc_variance_explained(cds)


    #      ## Step 2: Reduce the dimensions using UMAP ######

    #     progress$inc(0.1, detail = "reduce dimensions")
    #     cds <- monocle3::reduce_dimension(cds)  ## default is UMAP
    #     ##cds <- monocle3::reduce_dimension(cds, reduction_method="tSNE")


    #      ## Step 3: Cluster the cells #####

    #     progress$inc(0.2, detail = "clustering cells")
    #     cds <- cluster_cells(cds, reduction_method="UMAP")


    #     ## Step 4: Learn trajectory graph #######

    #     progress$inc(0.2, detail = "learning graph")
    #     cds <- learn_graph(cds)


    #     ## Step 5: Order cells ######

    #     ##cds <- order_cells(cds)
    #     ##plot_cells(cds)


    #     ## After: update selectors ######

    #     GENES = rownames(cds)
    #     shiny::updateSelectizeInput(session, "monocle_plotgene", choices=GENES, server=TRUE)
    #     grps = setdiff(colnames(cds@colData),c("Size_Factor"))
    #     shiny::updateSelectInput(session, "monocle_groupby", choices=grps, selected=".cluster" )

    #     progress$inc(0.2, detail = "done")

    #     dbg("monocle_getResults: DONE!")

    #     return(cds)
    # })

    # ## monocle_topMarkers <- shiny::reactive({ })

    # monocle_plotTopMarkers.RENDER <- shiny::reactive({

    #     cds <- monocle_getResults()

    #     shiny::req(cds,input$monocle_groupby)

    #     ## Find marker genes expressed by each cluster
    #     pheno1 = "cluster"
    #     pheno1 = input$monocle_groupby
    #     if(pheno1=="cluster") pheno1 <- ".cluster"
    #     pheno1
    #     marker_test_res = top_markers(cds, group_cells_by=pheno1, cores=4)

    #     NTOP = 1
    #     NTOP = 3
    #     NTOP = as.integer(input$monocle_ntop)
    #     top_specific_markers = marker_test_res %>%
    #         plotly::filter(fraction_expressing >= 0.10) %>%
    #         plotly::group_by(cell_group) %>%
    #         top_n(NTOP, specificity)
    #     ##top_n(1, marker_test_q_value)
    #     ##top_n(1, pseudo_R2)

    #     dbg("[SingleCellBoard:monocle_plotTopMarkers.RENDER] 3")

    #     top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
    #     scale_max1 = 0.8 * log(max(assay(cds))+0.1)
    #     ##par(mar=c(0,0,0,0))
    #     g <- monocle3::plot_genes_by_group(cds,
    #                                        top_specific_marker_ids,
    #                                        group_cells_by = pheno1,
    #                                        ##ordering_type="maximal_on_diag",
    #                                        scale_max = scale_max1, max.size = 5)
    #     ## g <- ggplot2::qplot(seq(0,4*pi,0.1), sin(seq(0,4*pi,0.1)))
    #     g <- g + ggplot2::theme(plot.margin=unit(c(1,1,1,1)*0.5,"cm"))

    #     dbg("[SingleCellBoard:monocle_plotTopMarkers.RENDER] done!")

    #     return(g)
    # })

    # monocle_plotTrajectory.RENDER <- shiny::reactive({


    #     cds <- monocle_getResults()
    #     shiny::req(cds,input$monocle_groupby)

    #     pheno1 = ".cluster"
    #     pheno1 = input$monocle_groupby
    #     if(pheno1=="cluster") pheno1 <- ".cluster"

    #     if(pheno1=="pseudotime") {
    #         root_cells <- input$monocle_rootcells
    #         cds = order_cells(cds, root_cells=root_cells)
    #     }

    #     size1 <- ifelse( ncol(cds) < 100, 2, 1)

    #     g <- monocle3::plot_cells(cds,
    #                               cell_size = size1,
    #                               group_label_size = 5,
    #                               color_cells_by = pheno1,
    #                               label_groups_by_cluster=FALSE,
    #                               label_leaves=FALSE,
    #                               label_branch_points=FALSE)
    #     g <- g + ggplot2::theme(plot.margin=unit(c(1,1,1,1)*0.5,"cm"))
    #     return(g)
    # })

    # monocle_plotGene.RENDER <- shiny::reactive({

    #     cds <- monocle_getResults()
    #     shiny::req(cds, input$monocle_plotgene)

    #     genes = "CD14"
    #     genes = input$monocle_plotgene

    #     size1 <- ifelse( ncol(cds) < 100, 2, 1)
    #     g <- monocle3::plot_cells(cds, genes = genes[1],
    #                               cell_size = size1,
    #                               show_trajectory_graph=FALSE,
    #                               label_cell_groups=FALSE)
    #     g <- g + ggplot2::theme(plot.margin=unit(c(1,1,1,1)*0.5,"cm"))
    #     return(g)
    # })


    # ##------- plotTopMarkers module -------
    # shiny::callModule(
    #     plotModule,
    #     id = "monocle_plotTopMarkers",
    #     func = monocle_plotTopMarkers.RENDER,
    #     func2 = monocle_plotTopMarkers.RENDER,
    #     plotlib = "ggplot",
    #     title = "Top markers by group", label="a",
    #     info.text = "The heatmap visualizes the expression of group-specific markers.",
    #     options = shiny::tagList(
    #         shiny::selectInput(ns("monocle_groupby"),"Group by:",choices=NULL),
    #         shiny::selectInput(ns("monocle_ntop"),"ntop:",choices=c(1,2,3,4,5,10,25),selected=5)
    #     ),
    #     pdf.width=5, pdf.height=8,
    #     height=imgH, res=c(72,95),
    #     add.watermark = WATERMARK
    # )

    # ##------- plotTrajectory module -------
    # monocle_plotTrajectory.opts = shiny::tagList()
    # shiny::callModule(
    #     plotModule,
    #     id = "monocle_plotTrajectory",
    #     plotlib = "ggplot",
    #     func = monocle_plotTrajectory.RENDER,
    #     func2 = monocle_plotTrajectory.RENDER,
    #     title = "Single-cell trajectory", label="b",
    #     info.text = "Single-cell trajectory analysis how cells choose between one of several possible end states. Reconstruction algorithms can robustly reveal branching trajectories, along with the genes that cells use to navigate these decisions.",
    #     options = monocle_plotTrajectory.opts,
    #     pdf.width=8, pdf.height=8,
    #     height = 0.45*imgH, res=c(80,95),
    #     add.watermark = WATERMARK
    # )

    # ##------- plotGene module -------
    # monocle_plotGene.opts = shiny::tagList(
    #     shiny::selectInput(ns("monocle_plotgene"),"Gene:",choices=NULL)
    # )
    # shiny::callModule(
    #     plotModule,
    #     id = "monocle_plotGene", plotlib="ggplot",
    #     func = monocle_plotGene.RENDER,
    #     func2 = monocle_plotGene.RENDER,
    #     title = "Gene expression", label="c",
    #     info.text = ".",
    #     options = monocle_plotGene.opts,
    #     pdf.width=8, pdf.height=8,
    #     height = 0.45*imgH, res = c(80,95),
    #     add.watermark = WATERMARK
    # )

    return(NULL)
  })
}
