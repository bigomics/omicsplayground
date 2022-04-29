##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

#' DataView module server function
#'
#' @description A shiny Module (server code).
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param pgx Reactive expression that provides the input pgx data object 
#'
#' @export 
DataViewBoard <- function(id, pgx)
{
  moduleServer(id, function(input, output, session)
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

      shiny::req(pgx$Y, pgx$samples)

      ## levels for sample filter
      levels = getLevels(pgx$Y)
      shiny::updateSelectInput(session, "data_samplefilter", choices=levels)

      grps <- pgx.getCategoricalPhenotypes(pgx$samples, min.ncat=2, max.ncat=999)
      grps <- sort(grps)
      grps <- c(grep("^[.]",grps,value=TRUE,invert=TRUE),grep("^[.]",grps,value=TRUE))
      selgrp <- grps[1]
      grps <- c("<ungrouped>",grps)
      if("group" %in% grps) selgrp = "group"
      if("condition" %in% grps) selgrp = "condition"      
      if(nrow(pgx$samples)<=20) selgrp = "<ungrouped>"
      shiny::updateSelectInput(session,'data_groupby', choices=grps, selected=selgrp)
    })


    shiny::observeEvent({
      input$data_type
      pgx$X
      pgx$counts
    }, {

      if(input$data_type %in% c("counts","CPM")) {
        pp <- rownames(pgx$counts)
      } else {
        ## log2CPM
        pp <- rownames(pgx$X)
      }

      ## gene filter. 
      genes <- sort(pgx$genes[pp,]$gene_name)
      fc2 = rowMeans(pgx.getMetaFoldChangeMatrix(pgx)$fc**2)
      genes = intersect(names(sort(-fc2)),genes) ## most var gene??
      selgene <- genes[1]
      genes1 <- unique(c(selgene,sort(genes)))
      if(length(genes1)>1000) {
        genes1 <- c(sort(genes1[1:1000]),"(type SYMBOL for more genes...)",
                    genes1[1001:length(genes1)])
      }
      shiny::updateSelectizeInput(session,'search_gene', choices=genes1, selected=selgene,
                                  ##options = list(maxOptions = 9999999),
                                  options = list(maxOptions = 1001),                                    
                                  server = TRUE)
    })

    last_search_gene <- reactiveVal()
    
    input_search_gene <- reactive({
      if( input$search_gene %in% c("(type SYMBOL for more genes...)","")) {
        gene1 <- last_search_gene() 
        return(gene1)
      }
      last_search_gene(input$search_gene)
      return(input$search_gene)
    })

    
    ##================================================================================
    ##=========================== MODULES ============================================
    ##================================================================================

    WATERMARK = FALSE
    
    ## get selected samples after sample filtering
    selected_samples <- reactive({
      samples <- colnames(pgx$X)
      if(!is.null(input$data_samplefilter)) {
        samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
      }
      samples
    })
    
    ## dbg("[***dataview_server] names.input = ",names(input))    
    dataview_module_geneinfo_server(
      "geneinfo", 
      r.gene  = reactive(input$search_gene),
      watermark = WATERMARK      
    )

    ## first tab ---------------------------------------
    dataview_plot_tsne_server(
      "tsneplot",
      pgx,
      r.gene         = reactive(input$search_gene),
      r.samples      = selected_samples,
      r.data_type    = reactive(input$data_type),            
      r.groupby = reactive(input$data_groupby),
      watermark = WATERMARK      
    )        

    dataview_plot_averagerank_server(
      "averagerankplot",
      pgx,
      r.gene         = reactive(input$search_gene),
      r.samples      = selected_samples,
      r.data_type    = reactive(input$data_type),
      watermark = WATERMARK      
    )            

    dataview_plot_correlation_server(
      "correlationplot",
      pgx,
      r.gene         = reactive(input$search_gene),
      r.samples      = selected_samples,
      watermark = WATERMARK      
    )

    dataview_plot_tissue_server(
      "tissueplot",
      pgx,
      r.gene         = reactive(input$search_gene),
      r.data_type    = reactive(input$data_type),
      watermark = WATERMARK      
    )            

    dataview_plot_expression_server(
      "expressionplot",
      pgx,
      r.gene         = reactive(input$search_gene),
      r.samples      = selected_samples,
      r.data_type    = reactive(input$data_type),            
      r.data_groupby = reactive(input$data_groupby),
      watermark = WATERMARK
    )

    ## second tab -----------------------------------
    dataview_plot_totalcounts_server(
      "counts_total",
      getCountStatistics,
      watermark = WATERMARK
    )
    
    dataview_plot_boxplot_server(
      "counts_boxplot",
      input,
      getCountStatistics,
      watermark = WATERMARK
    )

    dataview_plot_histogram_server(
      "counts_histplot",
      getCountStatistics,
      watermark = WATERMARK      
    )

    dataview_plot_abundance_server(
      "counts_abundance",
      getCountStatistics,
      watermark = WATERMARK
    )

    dataview_plot_averagecounts_server(
      "counts_averagecounts",
      getCountStatistics,
      watermark = WATERMARK
    )

    ## fourth tab
    dataview_plot_phenoheatmap_server(
      "phenoheatmap",
      pgx,
      r.samples = selected_samples,
      watermark = WATERMARK      
    )

    dataview_plot_phenoassociation_server(
      "phenoassociation",
      pgx,
      r.samples = selected_samples,
      watermark = WATERMARK      
    )    
    
    ##================================================================================
    ##===============================  TABLES ========================================
    ##================================================================================

    dataview_table_rawdata_server(
      "rawdatatable", pgx, 
      r.gene         = reactive(input$search_gene),
      r.data_type    = reactive(input$data_type),            
      r.samples      = selected_samples,
      r.groupby = reactive(input$data_groupby)      
    )    

    dataview_table_samples_server(
      "sampletable", pgx,
      r.samples      = selected_samples
    )

    dataview_table_resources_server(
      "resources", pgx
    )

    dataview_table_contrasts_server(
      "contrastTable", pgx,
      r.samples      = selected_samples      
    )
    
    ##================================================================================
    ##========================= FUNCTIONS ============================================
    ##================================================================================
    
    getCountStatistics <- shiny::reactive({
      shiny::req(pgx$X,pgx$Y,pgx$samples)

      shiny::validate(shiny::need("counts" %in% names(pgx), "no 'counts' in object."))
      subtt=NULL

      samples = colnames(pgx$X)
      samples <- selectSamplesFromSelectedLevels(pgx$Y, input$data_samplefilter)
      nsamples = length(samples)
      if("counts" %in% names(pgx)) {
        counts = pgx$counts[,samples,drop=FALSE]
      } else {
        cat("WARNING:: no counts table. estimating from X\n")
        counts = pmax(2**pgx$X - 1,0)
        k = grep("lib.size",colnames(pgx$samples))[1]
        if(length(k)>0) {
          libsize = pgx$samples[colnames(counts),k]
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
      gr = pgx$Y[samples,grpvar]
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
      gg = pgx$genes[rownames(counts),]$gene_name
      tt = pgx$genes[rownames(counts),]$gene_title
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

      res = list(
        total.counts = total.counts,
        subtt = subtt,
        cx1 = cx1,
        log2counts = log2counts,
        jj = jj,
        prop.counts = prop.counts,
        avg.counts = avg.counts
      )
    })
    
    
    ##================================================================================
    ##================================= END ====================================
    ##================================================================================



  })
}
