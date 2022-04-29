##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##


dataview_table_rawdata_ui <- function(id) {
  ns <- shiny::NS(id)
  tableWidget(ns("tbl"))  
}

dataview_table_rawdata_server <- function(id,
                                          pgx,
                                          r.gene = reactive(""),
                                          r.data_type = reactive("counts"),
                                          r.samples = reactive(""),
                                          r.groupby = reactive("")
                                          )
{
  moduleServer(id, function(input, output, session) {
    
    table_data <- shiny::reactive({
      ## get current view of raw_counts
      
      shiny::req(pgx$X,pgx$Y,pgx$genes,pgx$model.parameters)
      shiny::req(r.gene(),r.data_type())
      
      dbg("[dataview_rawdata:table_data] reacted!")

      ## dereference reactives
      gene <- r.gene()
      data_type <- r.data_type()
      samples <- r.samples()
      groupby <- r.groupby()
      
      if(is.null(gene) || gene=="" || is.na(gene)) {
        gene <- rownames(pgx$X)[1]
      }

      if(data_type=="counts") {
        x <- pgx$counts
      } else if(data_type=="CPM") {
        x <- edgeR::cpm(pgx$counts, log=FALSE)
      } else {
        ## log2CPM
        x <- pgx$X
      }
      x0=x

      ##------------------ select samples
      if(samples[1]=="") samples <- colnames(pgx$X)
      samples <- intersect(colnames(x),samples)
      x <- x[,samples,drop=FALSE]

      ## Quickly (?) calculated correlation to selected gene
      dbg("[dataview_rawdata:table_data] calculate rho")
      dbg("[dataview_rawdata:table_data] data_type = ",data_type)

      ## compute correlation (always in logCPM)
      rho = sdx = avg = NULL
      logx <- pgx$X[rownames(x),]
      xgenes <- pgx$genes[rownames(x),"gene_name"]
      k <- which(xgenes==gene)
      rho = cor( t(logx[,samples]), logx[k,samples], use="pairwise")[,1]
      rho = round(rho[rownames(x)], digits=3)
      sdx = round(apply(logx[,samples],1,sd),digits=3)
      avg <- round(rowMeans(x),digits=3)

      dbg("[dataview_rawdata:table_data] compute groupings")
      
      group <- NULL
      if(groupby %in% colnames(pgx$Y)) {
        group = pgx$Y[colnames(x),groupby]
      }
      if(length(samples)>500 && groupby=="<ungrouped>") {
        group <- pgx$model.parameters$group
      }
      do.grouped <- (groupby!="<ungrouped>")
      if(do.grouped && !is.null(group) ) {
        allgroups = sort(unique(group))
        newx = c()
        for(gr in allgroups) {
          mx = rowMeans(x[,which(group==gr),drop=FALSE],na.rm=TRUE)
          newx = cbind(newx, mx)
        }
        rownames(newx) = rownames(x)
        colnames(newx) = paste0("avg.",allgroups,"")
        x = newx
      }

      x = round(as.matrix(x), digits=3)
      x95 = quantile(as.vector(x0[which(x0>0)]),probs=0.95)
      x99 = quantile(as.vector(x0[which(x0>0)]),probs=0.99)

      if(NCOL(x)==0 || nrow(x)==0) return(NULL)

      dbg("[dataview_rawdata:table_data] create dataframe")
      ##rownames(x) = sub(".*:","",rownames(x))
      xgenes <- pgx$genes[rownames(x),"gene_name"]
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
      ##x = x[order(x$gene),,drop=FALSE]
      x = x[order(-x$rho,-x$SD),,drop=FALSE]      

      dbg("[dataview_rawdata:table_data] table_data() done!")
      
      list(
        x = x,
        x95 = x95,
        x99 = x99
      )

    }) ## %>% bindCache(pgx$Y, r.gene(), r.data_type(), r.groupby())

    rawdataTable.RENDER <- function() {      

      dt <- table_data()
      req(dt, dt$x)
      
      numcols <- grep('gene|title',colnames(dt$x),value=TRUE,invert=TRUE)
      tabH = 700  ## height of table
            
      DT::datatable(
        dt$x, rownames=FALSE,
        class = 'compact cell-border stripe hover',
        extensions = c('Buttons','Scroller'),
        selection = list(mode='single', target='row', selected=1),
        options=list(
          dom = 'lfrtip',
          pageLength = 25,
          lengthMenu = c(25, 40, 100, 250),
          scroller=FALSE, scrollY = FALSE,                          
          deferRender=TRUE
        )  ## end of options.list
      ) %>%
        DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>%
        DT::formatStyle(numcols,
                        background = DT::styleColorBar(c(0,dt$x99), 'lightblue'),
                        ##background = color_from_middle(x99, 'lightblue', '#f5aeae'),
                        backgroundSize = '98% 88%',
                        backgroundRepeat = 'no-repeat',
                        backgroundPosition = 'center')
    }

    rawdataTable_modal.RENDER <- function() {
      rawdataTable.RENDER() %>%
        DT::formatStyle(0, target='row', fontSize='20px', lineHeight='70%') 
    }
    
    dropdown_search_gene='<code>Search gene</code>'
    menu_grouped='<code>grouped</code>'
    menu_options='<code>Options</code>'
    info_text = paste0('Under the <strong>gene table </strong>, the average expression values of genes across the groups can be read. The samples (or cells) can be ungrouped by unclicking the ',menu_grouped, ' in the main <i>Options</i> to see the exact expression values per sample (or cell).', 'The genes in the table are ordered by the correlation (<b>rho</b> column) with respect to the gene selected by users from the ',dropdown_search_gene, ' setting. <b>SD</b> column reports the standard deviation of expression across samples (or cells).')
    
    shiny::callModule(
      tableModule, "tbl",
      func = rawdataTable.RENDER,
      csvFunc = table_data,
      title = "Gene expression table",
      filename = "counts.csv",
      info.text = info_text,
      caption2 = info_text
    )

  })  ## end of moduleServer
} ## end of server

