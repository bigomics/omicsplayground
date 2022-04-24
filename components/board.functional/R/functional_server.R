##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

FunctionalBoard <- function(id, inputData, selected_gsetmethods)
{
  moduleServer(id, function(input, output, session)
  {

    ns <- session$ns ## NAMESPACE
    fullH = 750
    rowH = 660  ## row height of panel
    tabH = 200  ## row height of panel
    tabH = '70vh'  ## row height of panel    

    fa_infotext = paste("This module performs specialized pathway analysis. <br><br>",a_KEGG," is a collection of manually curated pathways representing the current knowledge of molecular interactions, reactions and relation networks as pathway maps. In the <strong>KEGG pathway</strong> panel, each pathway is scored for the selected contrast profile and reported in the table. A unique feature of the platform is that it provides an activation-heatmap comparing the activation levels of pathways across multiple contrast profiles. This facilitates to quickly see and detect the similarities between profiles in certain pathways.

<br><br>In the <strong>GO</strong> panel, users can perform ",a_GO," (GO) analysis. GO defines functional concepts/classes and their relationships as a hierarchical graph. The GO database provides a computational representation of the current knowledge about roles of genes for many organisms in terms of molecular functions, cellular components and biological processes. All the features described under the KEGG pathway tab, such as scoring the gene sets and drawing an activation-heatmap, can be performed for the GO database under the GO graph tab. Instead of pathway maps, an annotated graph structure provided by the GO database is potted for every selected gene set.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
")

    

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    shiny::observeEvent( input$fa_info, {
        shiny::showModal(shiny::modalDialog(
            title = shiny::HTML("<strong>Functional Analysis Board</strong>"),
            shiny::HTML(fa_infotext),
            easyClose = TRUE, size="l" ))
    })
    
    shiny::observe({
        ngs <- inputData()
        shiny::req(ngs)
        ct <- colnames(ngs$model.parameters$contr.matrix)
        ##ct <- c(ct,"<sd>")
        ct <- sort(ct)
        shiny::updateSelectInput(session, "fa_contrast", choices=ct )
    })
    

    ##================================================================================
    ## KEGG pathways
    ##================================================================================

    getKeggTable <- shiny::reactive({

        ngs = inputData()
        shiny::req(ngs)
        shiny::req(input$fa_contrast)
        
        ## ----- get comparison
        comparison=2
        comparison <- input$fa_contrast
        if(!(comparison %in% names(ngs$gset.meta$meta))) return(NULL)

        ## ----- get KEGG id
        xml.dir <- file.path(FILES,"kegg-xml")
        kegg.available <- gsub("hsa|.xml","",dir(xml.dir, pattern="*.xml"))
        kegg.available
        kegg.ids <- getKeggID(rownames(ngs$gsetX))
        kegg.ids    
        ## sometimes no KEGG in genesets...
        if(length(kegg.ids)==0) {
            shinyWidgets::sendSweetAlert(
                session=session,
                title = "No KEGG terms in enrichment results",
                text="",
                type = "warning")
            dbg("[FunctionalBoard::getKeggTable] no KEGG terms in gsetX")
            df <- data.frame()
            return(df)
        }
        
        jj <- which(!is.na(kegg.ids) &
                    !duplicated(kegg.ids) &
                    kegg.ids %in% kegg.available )
        kegg.gsets <- rownames(ngs$gsetX)[jj]
        kegg.ids <- kegg.ids[jj]

        meta <- ngs$gset.meta$meta[[comparison]]
        meta <- meta[kegg.gsets,]
        mm <- selected_gsetmethods()
        dbg("[FunctionalBoard::getKeggTable] 1: gset methods = ",mm)
        mm <- intersect(mm, colnames(meta$q))
        dbg("[FunctionalBoard::getKeggTable] 2: gset methods = ",mm)
        meta.q <- apply(meta$q[,mm,drop=FALSE],1,max,na.rm=TRUE)
        
        df <- data.frame( kegg.id=kegg.ids, pathway=kegg.gsets,
                         logFC = meta$meta.fx, meta.q = meta.q,
                         check.names=FALSE)
        ##df <- df[order(-fx),]
        df <- df[!duplicated(df$kegg.id), ]  ## take out duplicated gene sets...
        df <- df[order(-abs(df$logFC)),]    
        return(df)
    })
    
    getFilteredKeggTable <- shiny::reactive({
        df <- getKeggTable()
        do.filter = FALSE
        do.filter <- input$fa_filtertable
        if(do.filter) df <- df[which(df$meta.q < 0.999),]
        return(df)
    })

    ## There is a bug in pathview::geneannot.map so we have to override
    ## "Error in pathview::mol.sum(gene.data, gene.idmap) : no ID can be mapped!"
    my.geneannot.map <- function(in.ids, in.type, out.type, org = "Hs", pkg.name = NULL, 
                                 unique.map = TRUE, na.rm = TRUE, keep.order = TRUE) 
    {
        if (is.null(pkg.name)) {
            data(bods)
            ridx = grep(tolower(paste0(org, "[.]")), tolower(bods[,1]))
            if (length(ridx) == 0) {
                ridx = grep(tolower(org), tolower(bods[, 2:3]))%%nrow(bods)
                if (length(ridx) == 0) 
                    stop("Wrong org value!")
                if (any(ridx == 0)) 
                    ridx[ridx == 0] = nrow(bods)
            }
            pkg.name = bods[ridx, 1]
        }
        pkg.on = try(requireNamespace(pkg.name), silent = TRUE)
        if (!pkg.on) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) 
                install.packages("BiocManager")
            BiocManager::install(pkg.name, suppressUpdates = TRUE)
            pkg.on = try(requireNamespace(pkg.name), silent = TRUE)
            if (!pkg.on) 
                stop(paste("Fail to install/load gene annotation package ", 
                           pkg.name, "!", sep = ""))
        }
        db.obj <- eval(parse(text = paste0(pkg.name, "::", pkg.name)))
        id.types <- AnnotationDbi::columns(db.obj)
        in.type = toupper(in.type)
        out.type = toupper(out.type)
        eii = in.type == toupper("entrez") | in.type == toupper("eg")
        if (any(eii)) 
            in.type[eii] = "ENTREZID"
        eio = out.type == toupper("entrez") | out.type == toupper("eg")
        if (any(eio)) 
            out.type[eio] = "ENTREZID"
        if (in.type == out.type) 
            stop("in.type and out.type are the same, no need to map!")
        nin = length(in.type)
        if (nin != 1) 
            stop("in.type must be of length 1!")
        out.type = out.type[!out.type %in% in.type]
        nout = length(out.type)
        msg <- paste0("must from: ", paste(id.types, collapse = ", "), 
                      "!")
        if (!in.type %in% id.types) 
            stop("'in.type' ", msg)
        if (!all(out.type %in% id.types)) 
            stop("'out.type' ", msg)
        in.ids0 = in.ids
        in.ids <- unique(as.character(in.ids))
        out.ids = character(length(in.ids))
###res <- try(suppressWarnings( plotly::select(db.obj, keys = in.ids, 
        res <- try(suppressWarnings(
            AnnotationDbi::select(db.obj, keys = in.ids, 
                                  keytype = in.type,
                                  columns = c(in.type, out.type))))
        if (class(res) == "data.frame") {
            res <- res[, c(in.type, out.type)]
            if (nout == 1) 
                na.idx <- is.na(res[, 2])
            else na.idx <- apply(res[, -1], 1, function(x) all(is.na(x)))
            if (sum(na.idx) > 0) {
                n.na <- length(unique(res[na.idx, 1]))
                if (n.na > 0) {
                    print(paste("Note:", n.na, "of", length(in.ids), 
                                "unique input IDs unmapped."))
                }
                if (na.rm) 
                    res <- res[!na.idx, ]
            }
            cns = colnames(res)
            if (unique.map) {
                if (length(out.type) == 1) 
                    umaps = tapply(res[, out.type], res[, in.type], 
                                   paste, sep = "", collapse = "; ")
                else umaps = apply(res[, out.type], 2, function(x)
                    tapply(x, res[, in.type], function(y) paste(unique(y),
                                                                sep = "", collapse = "; ")))
                umaps = cbind(umaps)
                res.uniq = cbind(rownames(umaps), umaps)
                res = res.uniq
                colnames(res) = cns
            }
            res = as.matrix(res)
            if (!keep.order) {
                rownames(res) = NULL
                return(res)
            }
            else {
                res1 = matrix(NA, ncol = length(cns), nrow = length(in.ids0))
                res1[, 1] = in.ids0
                rns = match(in.ids0, res[, 1])
                res1[, -1] = res[rns, -1]
                colnames(res1) = cns
                return(res1)
            }
        }
        else {
            res = cbind(in.ids, out.ids)
            colnames(res) = c(in.type, out.type)
            return(res)
        }
    }

    if(1) {
        suppressMessages(require(pathview))
        unlockBinding("geneannot.map", as.environment("package:pathview"))
        assignInNamespace("geneannot.map", my.geneannot.map, ns="pathview", as.environment("package:pathview"))
        assign("geneannot.map", my.geneannot.map, as.environment("package:pathview"))
        lockBinding("geneannot.map", as.environment("package:pathview"))
    }

    kegg_graph.RENDER <- shiny::reactive({

        ngs <- inputData()
        alertDataLoaded(session,ngs)
        ##NULL.IMG <- list(src=NULL, contentType = 'image/png')
        ##NULL.IMG <- list(src=NA, contentType = 'image/png')
        NULL.IMG <- list(src="", contentType = 'image/png')
        if(is.null(ngs)) return(NULL.IMG)


        comparison=1
        comparison = input$fa_contrast
        if(is.null(comparison) || length(comparison)==0) return(NULL.IMG)
        if(comparison=="") return(NULL.IMG)
        
        ## get fold-change vector
        fc <- ngs$gx.meta$meta[[comparison]]$meta.fx
        pp <- rownames(ngs$gx.meta$meta[[comparison]])

        if("hgnc_symbol" %in% colnames(ngs$genes)) {
            names(fc) <- ngs$genes[pp,"hgnc_symbol"]
        } else {
            names(fc) <- toupper(ngs$genes[pp,"gene_name"])
        }
        fc <- fc[order(-abs(fc))]
        fc <- fc[which(!duplicated(names(fc)) & names(fc)!="")]

        ## get selected KEGG id
        df <- getFilteredKeggTable()
        if(is.null(df)) return(NULL.IMG)
        
        sel.row=1
        sel.row <- kegg_table$rows_selected()
        ##if(is.null(sel.row) || length(sel.row)==0) return(NULL)
        if(is.null(sel.row) || length(sel.row)==0) return(NULL.IMG)
        sel.row <- as.integer(sel.row)
        
        pathway.id = "05213" 
        pathway.id = "04110" ## CELL CYCLE
        pathway.name = pw.genes = "x"
        if(is.null(sel.row) || length(sel.row)==0) return(NULL.IMG)

        if(!is.null(sel.row) && length(sel.row)>0) {
            pathway.id <- df[sel.row,"kegg.id"]
            pathway.name <- df[sel.row,"pathway"]
            ##pw.genes <- GSETS[[as.character(pathway.name)]]
            pw.genes <- unlist(getGSETS(as.character(pathway.name)))
        }

        ## folder with predownloaded XML files
        xml.dir <- file.path(FILES,"kegg-xml")
        xml.dir <- normalizePath(xml.dir)  ## absolute path

        ## We temporarily switch the working directory to always readable
        ## TMP folder
        curwd <- getwd()
        curwd
        tmpdir <- tempdir()
        tmpdir
        dbg("[FunctionalBoard::kegg_graph] switch to /tmp folder=",tmpdir,"\n")        
        setwd(tmpdir)
        pv.out <- pathview::pathview(
            gene.data = fc, pathway.id=pathway.id, gene.idtype="SYMBOL",
            gene.annotpkg = "org.Hs.eg.db", species = "hsa",
            out.suffix="pathview", limit = list(gene=2, cpd=1),
            low = list(gene = "dodgerblue2", cpd = "purple"),
            high = list(gene = "firebrick2", cpd = "yellow"), 
            kegg.dir=xml.dir, kegg.native=TRUE, same.layer=FALSE )
        Sys.sleep(0.2) ## wait for graph

        ## back to previous working folder
        setwd(curwd)   
        dbg("[FunctionalBoard::kegg_graph] back to working folder=",getwd(),"\n")
        
        ##width  <- session$clientData$output_kegg_graph_width
        ##height <- session$clientData$output_kegg_graph_height    
        outfile = file.path(tmpdir,paste0("hsa",pathway.id,".pathview.png"))
        dbg("[FunctionalBoard::kegg_graph] outfile=",outfile,"\n")
        dbg("[FunctionalBoard::kegg_graph] file.exists(outfile)=",file.exists(outfile),"\n")
        file.exists(outfile)
        if(!file.exists(outfile)) return(NULL.IMG)
                
        list(src = outfile, 
             contentType = 'image/png',
             ##width = 1040*0.8, height = 800*0.9, ## actual size: 1040x800
             ##width = 900, height = 600, ## actual size: 1040x800
             width = "100%", height = "100%", ## actual size: 1040x800
             ##width = img[2], height = img[1], ## actual size: 1040x800         
             alt = "pathview image")        
        ##    }, deleteFile = TRUE)
    })        

    ##output$kegg_table <- shiny::renderDataTable({
    kegg_table.RENDER <- shiny::reactive({
        
        ngs <- inputData()
        shiny::req(ngs)
        if(is.null(ngs$meta.go)) return(NULL)
        
        comparison=1
        comparison = input$fa_contrast

        if(is.null(comparison)) return(NULL)
        df <- getFilteredKeggTable()
        if(is.null(df)) return(NULL)
        if(nrow(df)==0) return(NULL)

        ## add hyperlink
        url = paste0("https://www.genome.jp/kegg-bin/show_pathway?map=hsa",df$kegg.id,"&show_description=show")
        df$kegg.id <- paste0("<a href='",url,"' target='_blank'>",df$kegg.id,"</a>")    
        ##df$pathway <- wrapHyperLink(df$pathway, df$pathway)
        
        numeric.cols <- colnames(df)[which(sapply(df, is.numeric))]
        
        DT::datatable( df, rownames=FALSE, escape = c(-1,-2),
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "logFC",
                                ##background = DT::styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle( df[,"logFC"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })
    

    normalize=1;nterms=40;nfc=10
    plotKEGGactmap <- function(meta, df, normalize, nterms, nfc)
    {
        
        fx <- sapply(meta, function(x) x$meta.fx)
        qv <- sapply(meta, function(x) x$meta.q)    
        rownames(fx) <- rownames(qv) <- rownames(meta[[1]])

        kk <- rownames(fx)
        kk <- as.character(df$pathway)
        ##if(is.na(kk) || kk=="" || !(kk %in% rownames(fx))) return(NULL)
        if(length(kk) < 3) return(NULL)

        if(mean(is.na(qv))<0.01) {
            score <- fx[kk,,drop=FALSE] * (1 - qv[kk,,drop=FALSE])**2
        } else {
            score <- fx[kk,,drop=FALSE] 
        }
        dim(score)
        ## if(NCOL(score)==1) score <- cbind(score,score)  ## UGLY...
        
        score <- score[head(order(-rowSums(score**2)),nterms),,drop=FALSE]  ## nr gene sets
        score <- score[,head(order(-colSums(score**2)),nfc),drop=FALSE]  ## max comparisons/FC    
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        ii=1:nrow(score);jj=1:ncol(score)
        if(NCOL(score)==1) {
            score <- score[order(-score[,1]),1,drop=FALSE]
        } else {
            ii <- hclust(d1)$order
            ## ii <- order(-rowMeans(score))
            jj <- hclust(d2)$order
            score <- score[ii,jj,drop=FALSE]
        }
        dim(score)
        rownames(score) = substring(rownames(score),1,50)
        
        score2 <- score
        if(normalize) score2 <- t(t(score2) / apply(abs(score2),2,max))
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**3   ## fudging
        rownames(score2) <- tolower(gsub(".*:|kegg_|_Homo.*$","",
                                         rownames(score2),ignore.case=TRUE))
        rownames(score2) <- substring(rownames(score2), 1, 40)
        colnames(score2) <- shortstring(colnames(score2), 30)
        colnames(score2) <- paste0(colnames(score2)," ")
        
        ## heatmap(score2, scale="none", mar=c(8,20))
        bmar <- 0 + pmax(50 - nrow(score2),0)*0.3
        par(mfrow=c(1,1), mar=c(1,1,10,1), oma=c(0,1.5,0,0.5) )

        corrplot::corrplot( score2, is.corr=FALSE, cl.pos="n", col=BLUERED(100),
                 tl.cex = 0.85, tl.col="grey20", tl.srt = 90,
                 mar=c(bmar,0,0.5,0) )


    }
    
    kegg_actmap.RENDER <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)
        df <- getKeggTable()
        if(is.null(df) || nrow(df)==0) {
            dbg("[FunctionalBoard::kegg_actmap.RENDER] emtpy KEGG table")
            ##par(mfrow=c(1,1), mar=c(1,1,1,1)*0, oma=c(0,2,0,1)*0 )
            ##frame()
            return(NULL)
        }
        meta <- ngs$gset.meta$meta
        plotKEGGactmap(meta, df, normalize = input$kegg_normalize, nterms=50, nfc=25)
        
    })    

    kegg_actmap.RENDER2 <- shiny::reactive({


        ngs <- inputData()
        shiny::req(ngs)
        df <- getKeggTable()
        if(is.null(df) || nrow(df)==0) {
            dbg("[FunctionalBoard::kegg_actmap.RENDER] emtpy KEGG table")
            ##par(mfrow=c(1,1), mar=c(1,1,1,1)*0, oma=c(0,2,0,1)*0 )
            ##frame()
            return(NULL)
        }
        meta <- ngs$gset.meta$meta
        plotKEGGactmap(meta, df, normalize = input$kegg_normalize, nterms=50, nfc=100)        
    })    

    kegg_info1 = "<strong>KEGG pathways</strong> are a collection of manually curated pathways representing the current knowledge of molecular interactions, reactions and relation networks as pathway maps. In the pathway map, genes are colored according to their upregulation (red) or downregulation (blue) in the contrast profile. Each pathway is scored for the selected contrast profile and reported in the table below. "

    kegg_graph.opts <- shiny::tagList(
    )
    
    shiny::callModule(
        plotModule,
        id = "kegg_graph", label="a",
        title = "Kegg pathway map",
        func = kegg_graph.RENDER,
        plotlib = "image",
        ##renderFunc = "renderImage", outputFunc = "imageOutput",
        options = kegg_graph.opts,
        download.fmt = "png", just.info=TRUE,
        info.text = kegg_info1, info.width="350px",
        height = c(0.53*rowH,700), width = c("100%",1280),
        add.watermark = WATERMARK
    )

    kegg_table_info = "<strong>Enrichment table.</strong> The table is interactive; enabling user to sort on different variables and select a pathway by clicking on the row in the table. The scoring is performed by considering the total number of genes in the pathway (n), the number of genes in the pathway supported by the contrast profile (k), the ratio of k/n, and the ratio of |upregulated or downregulated genes|/k. Additionally, the table contains the list of the upregulated and downregulated genes for each pathway and a q value from the Fisher’s test for the overlap."

    kegg_table_opts <- shiny::tagList(
        ##selectInput(ns("kegg_table_logfc"),"logFC threshold for Fisher-test",c(0.2,0.5,1))
    )
    
    kegg_table <- shiny::callModule(
        tableModule,
        id = "kegg_table", label="b",
        func = kegg_table.RENDER,
        options = kegg_table_opts,
        info.text = kegg_table_info, 
        info.width = "350px",
        title = "Enrichment table",
        height = c(270,700)
    )

    kegg_actmap.opts = shiny::tagList(
        withTooltip(shiny::checkboxInput(ns('kegg_normalize'),'normalize activation matrix',FALSE),"Click to normalize the columns of the activation matrices.")
    )
    shiny::callModule(
        plotModule,
        id = "kegg_actmap",
        func  = kegg_actmap.RENDER,
        func2 = kegg_actmap.RENDER2, 
        title = "Activation matrix", label="c",
        info.text = "The <strong>KEGG activation matrix</strong> visualizes the activation levels of pathways (or pathway keywords) across multiple contrast profiles. This facilitates to quickly see and detect the similarities of certain pathways between contrasts. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = kegg_actmap.opts,
        pdf.height = 9, pdf.width = 9 , 
        height = c(rowH,750), width = c("100%",1400),
        res=72,
        add.watermark = WATERMARK
    )
    
    

    ##================================================================================
    ## GO graph
    ##================================================================================
    


    GO_network.RENDER <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)
        require(igraph)
        
        comparison=1;methods=c("fisher","gsva","camera")
        comparison = input$fa_contrast
        shiny::req(input$fa_contrast)
        if(is.null(comparison)) return(NULL)
        
        ##sub2 <- getSigGO(comparison, methods, nterms=250, ntop=25, ngs=ngs)
        sub2 <- go <- ngs$meta.go$graph
        if(is.null(go)) {
            shinyWidgets::sendSweetAlert(
                session=session,
                title = "No GO graph in enrichment results",
                text="",
                type = "warning")
            dbg("[GO_network.RENDER] ***ERROR*** no META.GO in pgx object!")
            return(NULL)
        }
        
        score = ngs$meta.go$pathscore[,comparison]        
        score[is.na(score) | is.infinite(score)] = 0
        score = (score / (1e-8+max(abs(score),na.rm=TRUE)))
        igraph::V(sub2)$value <- score
        ##igraph::V(sub2)$color <- gplots::bluered(32)[16 + round(15*score)]
        igraph::V(sub2)$color <- BLUERED(32)[16 + round(15*score)]        
        igraph::V(sub2)$label <- igraph::V(sub2)$Term
        igraph::V(sub2)$label[which(is.na(score)|score==0)] = ""
        pos = sub2$layout

        dbg("[FunctionalBoard::GO_network.RENDER] 1: sum(is.na(score)))=",sum(is.na(score)))
        dbg("[FunctionalBoard::GO_network.RENDER] 1: all(score=0)=",all(score==0))        
        all.zero <- all(score==0)
        
        if(!all.zero && input$GO_prunetree) {
            ##cat("pruning GO graph\n")
            vv = igraph::V(sub2)[which(!is.na(score) & abs(score)>0)]
            sp = igraph::shortest_paths(sub2, from="all", to=vv, mode="all", output="vpath")
            sp.vv = unique(unlist(sp$vpath))
            sub2 = igraph::induced.subgraph(sub2, sp.vv)
            pos = igraph::layout_with_fr(sub2)
            score = score[igraph::V(sub2)$name]
        }

        dbg("[FunctionalBoard::GO_network.RENDER] 2: len.V(sub2)=",length(igraph::V(sub2)))
        if(length(igraph::V(sub2))==0) {
            ## return(NULL)
        }
        
        ## remove root?
        removeroot=TRUE
        if(removeroot) {
            sub2 <- igraph::induced_subgraph(sub2, which(igraph::V(sub2)$name!="all"))
            if(input$GO_prunetree) pos = igraph::layout_with_fr(sub2)
            score <- score[igraph::V(sub2)$name]
            ##pos <- pos[igraph::V(sub2)$name,]        
        }
        roots <- c("all",neighbors(go, igraph::V(go)["all"], mode="all" )$name)
        roots <- intersect(roots, igraph::V(sub2)$name)
        
        astree = TRUE
        if(astree) {
            if("all" %in% igraph::V(sub2)$name) {
                pos = igraph::layout_as_tree(sub2, root="all", mode="all")
            } else {
                pos = igraph::layout_as_tree(sub2, root=roots, mode="all")
            }
            pos[,2] = -pos[,2]
        }

        ## color clusters
        if(input$GO_colorclusters) {
            ##clust = igraph::cluster_louvain(igraph::as.undirected(sub2))$membership
            clust = igraph::cluster_louvain(igraph::as.undirected(go))$membership
            names(clust) = igraph::V(go)$name
            cc = c(RColorBrewer::brewer.pal(12,"Set3"),
                   RColorBrewer::brewer.pal(8,"Set2"),
                   RColorBrewer::brewer.pal(8,"Set1"))
            igraph::V(sub2)$color = rep(cc,99)[clust[igraph::V(sub2)$name]]
            jj = which(is.na(score) | score==0)
            if(length(jj)>0) igraph::V(sub2)$color[jj] = NA
        }
        

        ##pos <- pos[igraph::V(sub2)$name,]
        gr = visNetwork::toVisNetworkData(sub2)
        gr$nodes$color[is.na(gr$nodes$color)] = "#F9F9F9"
        gr$nodes$value = pmax(abs(gr$nodes$value),0.001)
        gr$nodes$x = pos[,1]*60
        gr$nodes$y = pos[,2]*90
        gr$nodes$label = gr$nodes$Term
        no.score <- (is.na(score)|score==0)
        gr$nodes$label[which(no.score)] = NA

        ##if(input$fa_boxnode) {
        gr$nodes$shape <- c("box","circle")[1 + 1*no.score]
        gr$nodes$label <- sapply(gr$nodes$label,breakstring,n=25,nmax=95,force=TRUE,brk="\n")
        
        ##gr.def <- breakstring(gr$nodes$Definition,80)
        gr.def <- sapply(gr$nodes$Definition,breakstring,n=50,brk="<br>")
        gr$nodes$title = paste0(gr$nodes$Term,"  (",gr$nodes$id,")<br>",
                                "<small>",gr.def,"</small>")
        ##gr$edges$title = edge.info   

        ## rendering
        font.size=20; cex=1
        if(input$GO_prunetree) {
            font.size=20; cex=0.6
        }

        visNetwork::visNetwork(gr$nodes, gr$edges) %>%
            visNetwork::visEdges(smooth = FALSE, hidden=FALSE, arrows=list(enabled=TRUE),
                     scaling=list(min=10*cex, max=30*cex), width=5*cex )  %>%
            visNetwork::visNodes(font = list( size=font.size*cex, vadjust=0),
                     scaling=list(min=1*cex, max=80*cex))  %>%
            visNetwork::visPhysics(stabilization = FALSE)  %>%
            ##visInteraction(hideEdgesOnDrag = TRUE) %>%
            ##visInteraction(navigationButtons = TRUE) %>%
            ##visOptions(nodesIdSelection = TRUE) %>%
            ##visOptions(selectedBy="component") %>%
            visNetwork::visOptions(highlightNearest = list(enabled=T, degree=1, hover=TRUE)) %>%
            ##visEvents(select="function(nodes){Shiny.onInputChange('current_node_id',nodes.nodes);;}") %>%
            visNetwork::visPhysics(enabled=FALSE) 
    })    

    matchGOid2gset <- function(id,gsets) {
        gsets.id <- sub("\\)$","",sub(".*\\(GO_","GO:",gsets))
        match(id,gsets.id)
    }

    GO_table.RENDER <- shiny::reactive({
        ngs <- inputData()
        shiny::req(ngs, input$fa_contrast)
        if(is.null(ngs$meta.go)) {
            dbg("[FunctionalBoard::GO_table.RENDER] no META.GO in pgx object!")
            return(NULL)
        }
        
        comparison=1
        comparison = input$fa_contrast
        ##req(input$fa_contrast)
        if(is.null(comparison)) return(NULL)

        go = ngs$meta.go$graph
        scores <- ngs$meta.go$pathscore[,comparison]

        dbg("[GO_table.RENDER] 1: len.scores = ",length(scores))
        
        scores <- scores[which(!is.na(scores) & !is.infinite(scores))]
        scores <- round(scores, digits=3)
        ##scores <- sort(scores, decreasing=TRUE)
        scores <- scores[order(-abs(scores))]
        go.term = igraph::V(go)[names(scores)]$Term

        dbg("[GO_table.RENDER] 2: len.scores = ",length(scores))
        
        ## get FC and q-value.  match with enrichment table
        gs.meta <- ngs$gset.meta$meta[[comparison]]
        ii <- matchGOid2gset(names(scores),rownames(gs.meta))
        gs.meta <- gs.meta[ii,,drop=FALSE]
        gs.meta$GO.id <- rownames(scores)        
        mm <- selected_gsetmethods()
        mm <- intersect(mm, colnames(gs.meta$q))
        qv <- apply(gs.meta$q[,mm,drop=FALSE],1,max,na.rm=TRUE) ## meta-q
        fx <- gs.meta$meta.fx
        
        dbg("[GO_table.RENDER] 2: dim(gs.meta) = ",paste(dim(gs.meta),collapse="x"))
        
        go.term1 = substring(go.term, 1, 80)
        dt1 = round( cbind(score=scores, logFC=fx, meta.q=qv), digits=4)    
        dt = data.frame( id=names(scores), term=go.term1, dt1, stringsAsFactors=FALSE)    
        id2 = paste0("abc(",sub(":","_",dt$id),")")  ## to match with wrapHyperLink
        dt$id <- wrapHyperLink(as.character(dt$id), id2)  ## add link
        
        numeric.cols <- colnames(dt)[which(sapply(dt, is.numeric))]
        numeric.cols

        DT::datatable( dt, rownames=FALSE, escape = c(-1,-2),
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      fillContainer = TRUE,
                      options=list(
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "score",
                                ##background = DT::styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle( dt1[,"score"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })

    ## normalize=1;maxterm=50;maxfc=10
    plotGOactmap <- function(score, go, normalize, maxterm, maxfc)
    {

        rownames(score) = igraph::V(go)[rownames(score)]$Term
        
        dbg("[plotGOactmap] 0: sum(isna(score)) = ", sum(is.na(score)))
        dbg("[plotGOactmap] 0: min(score,na.rm=0) = ", min(score,na.rm=FALSE))
        dbg("[plotGOactmap] 0: max(score,na.rm=0) = ", max(score,na.rm=FALSE))
        dbg("[plotGOactmap] 0: min(score,na.rm=1) = ", min(score,na.rm=TRUE))
        dbg("[plotGOactmap] 0: max(score,na.rm=1) = ", max(score,na.rm=TRUE))

        ## avoid errors!!!
        score[is.na(score) | is.infinite(score)] <- 0
        score[is.na(score)] = 0
        
        ## reduce score matrix
        ##score = head(score[order(-rowSums(abs(score))),],40)
        score = score[head(order(-rowSums(score**2,na.rm=TRUE)),maxterm),,drop=FALSE] ## max number terms    
        score = score[,head(order(-colSums(score**2,na.rm=TRUE)),maxfc),drop=FALSE] ## max comparisons/FC
        ##if(NCOL(score)==1) score <- cbind(score,score)   
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))

        ## normalize colums
        if(normalize) {
            ## column scale???
            score <- t(t(score) / (1e-8 + sqrt(colMeans(score**2,na.rm=TRUE))))
        }
        score <- score / max(abs(score),na.rm=TRUE) ## global normalize
        score <- sign(score) * abs(score)**0.5   ## fudging for better colors       

        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1 <- dist(score)
        d2 <- dist(t(score))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        ii=1:nrow(score); jj=1:ncol(score)
        if(NCOL(score)==1) {
            score <- score[order(-score[,1]),1,drop=FALSE]
        } else {
            ii <- hclust(d1)$order
            ## ii <- order(-rowMeans(score))
            jj <- hclust(d2)$order
            score <- score[ii,jj,drop=FALSE]
        }
        
        colnames(score) = substring(colnames(score),1,30)
        rownames(score) = substring(rownames(score),1,50)
        colnames(score) <- paste0(colnames(score)," ")
        
        bmar <- 0 + pmax((50 - nrow(score))*0.25,0)
        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,1.5,0,0.5))

        corrplot::corrplot( score, is.corr=FALSE, cl.pos="n", col=BLUERED(100),
                           tl.cex = 0.85, tl.col="grey20", tl.srt = 90,
                           mar=c(bmar,0,0,0) )


    }

    GO_actmap.RENDER <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)

        if(is.null(ngs$meta.go)) {
            dbg("[FunctionalBoard:GO_actmap.RENDER] no META.GO in pgx object!")
            return(NULL)
        }
        
        score = ngs$meta.go$pathscore
        go = ngs$meta.go$graph

        plotGOactmap(
            score = score, go = go,
            normalize = input$go_normalize,
            maxterm = 50,
            maxfc = 25
        )
                    
    })    

    GO_actmap.RENDER2 <- shiny::reactive({

        ngs <- inputData()
        shiny::req(ngs)

        if(is.null(ngs$meta.go)) {
            dbg("[FunctionalBoard:GO_actmap.RENDER] no META.GO in pgx object!")
            return(NULL)
        }

        score = ngs$meta.go$pathscore
        ##if(is.null(score)) return(NULL)
        go = ngs$meta.go$graph

        plotGOactmap(
            score = score, go = go,
            normalize = input$go_normalize,
            maxterm = 50,
            maxfc = 100
        )
                    
    })    

    GO_info1 = "The <strong>Gene Ontology</strong> (GO) provides a computational representation of the current knowledge about roles of genes for many organisms in terms of molecular functions, cellular components and biological processes. The structure of GO can be described in terms of a graph, where each GO term is a node, and the relationships between the terms are edges between the nodes. GO is loosely hierarchical, with ‘child’ terms being more specialized than their ‘parent’ terms. The graph is interactive. You can move the graph and zoom in using the mouse."

    GO_network.opts = shiny::tagList(
        withTooltip( shiny::checkboxInput(ns("GO_prunetree"), "Prune tree", TRUE),
               "Prune the tree with only significant branches."),
        withTooltip( shiny::checkboxInput(ns("GO_colorclusters"), "Color clusters", FALSE),
               "Highlight clusters with different colors.")
    )

    shiny::callModule(
        plotModule,
        id = "GO_network", 
        title = "Gene Ontology graph", label="a",
        func = GO_network.RENDER,
        plotlib = "visnetwork",
        info.text = GO_info1,
        download.fmt = c("pdf","png"), ## no.download=TRUE,
        options = GO_network.opts,
        pdf.width = 10, pdf.height = 8,
        height = c(0.55*rowH,750), width = c("100%",1400),        
        res=72,
        add.watermark = WATERMARK
    )
    ##output <- attachModule(output, GO_network_module)

    GO_actmap.opts = shiny::tagList(
        withTooltip(shiny::checkboxInput(ns('go_normalize'),'normalize activation matrix',FALSE),"Click to normalize the columns of the activation matrices.")
    )
    go_info = "The <b>GO activation matrix</b> visualizes the activation of GO terms across conditions. From this figure, you can easily detect GO terms that are consistently up/down across conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile."
    
    shiny::callModule(
        plotModule,
        id = "GO_actmap",
        func = GO_actmap.RENDER,
        func2 = GO_actmap.RENDER2, 
        title = "Activation matrix", label="c",
        info.text = go_info,
        options = GO_actmap.opts,
        pdf.height = 9, pdf.width = 9, 
        height = c(rowH,750), width = c("100%",1400),
        res=72,
        add.watermark = WATERMARK
    )
    
    GO_table <- shiny::callModule(
        tableModule,
        id = "GO_table", label="b",
        func = GO_table.RENDER,
        info.text="<strong>GO score table.</strong> The scoring of a GO term is performed by considering the cumulative score of all terms from that term to the root node. That means that GO terms that are supported by higher level terms levels are preferentially scored.", 
        title = "GO score table",
        height = c(270,700)        
    )

    ##================================================================================
    ## Fire plot (dev)
    ##================================================================================

    pgx.firePlot <- function(ngs, cmp, gsets, shownames=TRUE) {
        k0 <- which(ngs$model.parameters$exp.matrix[,cmp] <0)
        gx <- rowMeans(ngs$X[,k0])
        fx <- ngs$gx.meta$meta[[cmp]]$meta.fx
        names(fx) <- rownames(ngs$gx.meta$meta[[1]])
        names(gx) <- toupper(sub(".*:","",names(gx)))
        names(fx) <- toupper(sub(".*:","",names(fx)))
        
        ## gene set expression as mean of members
        gsets0 <- gsets
        gsets <- lapply(gsets, intersect, names(fx))    
        zx <- sapply( gsets, function(gg) mean(fx[gg],na.rm=TRUE))
        zx <- zx * (1 - exp( -(sapply(gsets,length)/10)**1 ))  ## moderate by set size
        
        gsets <- head(gsets[order(-abs(zx[names(gsets)]))],60)
        gsets <- gsets[order(zx[names(gsets)])] ## order by gset expression
        
        yy  <- lapply(gsets, function(g) gx[g])
        dx <- unlist(lapply(gsets, function(g) fx[g]))
        dx <- sign(dx) * abs(dx/max(abs(dx)))**0.66
        xx <- mapply(rep, 1:length(yy), sapply(yy,length))
        
        BLUERED1 <- colorRampPalette(
            rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#DDDDDD",
                  "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))    
        klr <- paste0(BLUERED1(15)[8 + 7*dx],"88")
        
        xx0 <- unlist(xx)
        yy0 <- unlist(yy)
        gg0 <- as.vector(unlist(gsets))
        
        ## plot the genesets with gene members
        ii <- order(abs(dx))
        par(mfrow=c(2,1), mar=c(0,3.5,2,2), mgp=c(2.1,0.8,0))
        dy <- 0.02*diff(range(yy0))
        base::plot( xx0[ii], yy0[ii] + dy, col=klr[ii],
             xaxs="i", xlim=c(-2, max(xx0)+3), 
             xaxt="n", xlab="", ylab="expression (log2.CPM)", 
             pch=19, cex=1.66*abs(dx)[ii], ylim=c(-4,max(yy0)*1.05) )
        title(cmp, cex.main=1.0, line=0.9)
        abline(h=0, lty=2)
        
        gset.name <- names(gsets)
        if(0 && shownames) {    
            n0 <- sapply(gsets0[names(gsets)], length)
            n1 <- sapply(gsets, length)
            gset.name <- paste0( gset.name," (",n1,"/",n0,")")
        }
        gset.name <- sapply(gset.name, function(s) shortstring(s,48,dots=0.5))
        mtext( gset.name, side=1, at=1:length(gsets), 
              las=3, cex=0.75, line=0.8, col="grey35", xpd=TRUE)
        
        ## label some top FC genes
        if(shownames) {
            top.genes <- head(names(sort(-abs(fx[unique(gg0)]))),25)
            idx <- 1:length(gsets)
            ii <- order(-abs(fx[gg0]))
            g <- top.genes[1]
            jj <- c()
            for(g in top.genes) {
                i0 <- ii[match(g, gg0[ii])]
                jj <- c(jj, i0)
                idx <- setdiff(idx, xx0[i0])
                ii <- ii[ which(xx0[ii] %in% idx)]
            }

            jj <- setdiff(jj,NA)
            boxes = sapply(nchar(gg0[jj]),function(n) paste(rep("█",n),collapse=""))
            text( xx0[jj], yy0[jj], labels=boxes, cex=0.8, pos=3, offset=0.4, col="#FFFFFF99")        
            text( xx0[jj], yy0[jj], gg0[jj], cex=0.75, pos=3, offset=0.4 )
            for(i in 1:2) {
                points( xx0[jj], yy0[jj], col=klr[jj], pch=19, cex=1.66*abs(dx)[jj] )
            }
        }
        
        ## plot the "missing" genes
        gsets1 <- mapply( setdiff, gsets0[names(gsets)], gsets)
        tt <- sort(table(unlist(gsets1)),decreasing=TRUE)
        tt <- tt / max(tt)
        head(tt,10)
        top.missing <- head(names(tt),10)
        m = top.missing[1]
        for(i in 1:length(top.missing)) {
            m <- top.missing[i]
            jj <- which(sapply(gsets1, function(gs) (m %in% gs)))
            y0 <- -0.2 - 0.40*i
            points( jj, rep(y0,length(jj)), col="grey60", pch=21, cex=0.7 )
            if(shownames) {
                mtext(m, side=4, at=y0, las=1, cex=0.55,
                      col="grey35", line=0.2)
            }
            ##text(sample(jj,1), y0, m, pos=4, cex=0.5)    }
        }
    } ## end of firePlot()

    output$fa_fireplot <- shiny::renderPlot({

        ngs <- inputData()
        shiny::req(ngs)
        ##df <- getFilteredKeggTable()
        
        cmp=1
        cmp <- input$fa_contrast

        sel="B-cell related"
        sel="KEGG metabolic pathways"
        sel="Hallmark collection"
        sel=5
        sel <- input$fire_xpcollection    
        if(is.null(sel)) return(NULL)

        ##gsets <- GSETS[ COLLECTIONS[[sel]] ]
        gsets <- getGSETS( COLLECTIONS[[sel]] )     
        ##gsets0 <- gsets0[ intersect( names(gsets0), names(zx)) ]
        shownames <- input$fire_shownames
        pgx.firePlot(ngs, cmp, gsets, shownames=shownames)     
        
    }, res=85)

    output$fireplot_UI <- shiny::renderUI({
        shiny::fillRow(
            height = fullH,
            flex = c(2,1),
            shiny::fillCol(
                flex = c(NA,1), 
                shiny::inputPanel(
                    shiny::selectInput(ns("fire_xpcollection"),NULL,
                                choices=setdiff(names(COLLECTIONS),"<all>") ),
                    shiny::checkboxInput(ns('fire_shownames'),'shownames',TRUE),            
                    cellArgs=list(width='100%')
                ),
                shiny::plotOutput('fa_fireplot')
            ),
            shiny::br()
        )
    })

  }) ## end-of-moduleServer

}
