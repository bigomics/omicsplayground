FunctionalInputs <- function(id) {
    ns <- NS(id)  ## namespace
    tagList(
        uiOutput(ns("description")),
        uiOutput(ns("inputsUI"))
    )
}

FunctionalUI <- function(id) {
    ns <- NS(id)  ## namespace
    fillCol(
        flex = c(1),
        height = 780,
        tabsetPanel(
            id = ns("tabs"),
            tabPanel("KEGG",uiOutput(ns("kegg_analysis_UI"))),
            tabPanel("GO graph",uiOutput(ns("GO_analysis_UI"))),
            tabPanel("Drug CMap",uiOutput(ns("DSEA_analysis_UI"))),
            tabPanel("WordCloud",uiOutput(ns("wordcloud_UI")))
            ## tabPanel("Fire plot (dev)",uiOutput(ns("fireplot_UI")))            
        )
    )
}

FunctionalModule <- function(input, output, session, env)
{
    ns <- session$ns ## NAMESPACE
    inputData <- env[["load"]][["inputData"]]
    rowH = 780  ## row height of panel
    tabH = 200  ## row height of panel    
    description = "<b>Functional analysis</b>. <br> Perform specialized functional analysis
to understand biological functions including GO, KEGG, and drug connectivity mapping."
    output$description <- renderUI(HTML(description))
    
    description = "<b>Functional analysis</b>. <br> Perform specialized functional analysis
to understand biological functions including GO, KEGG, and drug connectivity mapping."
    output$description <- renderUI(HTML(description))

    fa_infotext = paste("This module performs specialized pathway and drug enrichment analysis. <br><br>",a_KEGG," is a collection of manually curated pathways representing the current knowledge of molecular interactions, reactions and relation networks as pathway maps. In the <strong>KEGG pathway</strong> panel, each pathway is scored for the selected contrast profile and reported in the table. A unique feature of the platform is that it provides an activation-heatmap comparing the activation levels of pathways across multiple contrast profiles. This facilitates to quickly see and detect the similarities between profiles in certain pathways.

<br><br>In the <strong>GO</strong> panel, users can perform ",a_GO," (GO) analysis. GO defines functional concepts/classes and their relationships as a hierarchical graph. The GO database provides a computational representation of the current knowledge about roles of genes for many organisms in terms of molecular functions, cellular components and biological processes. All the features described under the KEGG pathway tab, such as scoring the gene sets and drawing an activation-heatmap, can be performed for the GO database under the GO graph tab. Instead of pathway maps, an annotated graph structure provided by the GO database is potted for every selected gene set.

<br><br> In the <a href='https://portals.broadinstitute.org/cmap/'>Drug Connectivity Map</a> panel, you can correlate your signature with more than 5000 known drug profiles from the L1000 database. An activation-heatmap compares drug activation profiles across multiple contrasts. This facilitates to quickly see and detect the similarities between contrasts for certain drugs.

<br><br><br><br>
<center><iframe width='500' height='333' src='https://www.youtube.com/embed/watch?v=qCNcWRKj03w&list=PLxQDY_RmvM2JYPjdJnyLUpOStnXkWTSQ-&index=6' frameborder='0' allow='accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe></center>
")

    
    ##================================================================================
    ##========================= INPUTS UI ============================================
    ##================================================================================

    output$inputsUI <- renderUI({
        ui <- tagList(
            tipify( actionLink(ns("fa_info"), "Info", icon = icon("info-circle") ),
                   "Show more information about this module."),
            hr(), br(),             
            tipify( selectInput(ns("fa_contrast"),"Contrast:", choices=NULL),
                   "Select the contrast corresponding to the comparison of interest.",
                   placement="top"),
            tipify( actionLink(ns("fa_options"), "Options", icon=icon("cog", lib = "glyphicon")),
                   "Show/hide advanced options", placement="top"),
            br(),
            conditionalPanel(
                "input.fa_options % 2 == 1", ns=ns,
                tagList(
                    tipify(checkboxInput(ns('fa_normalize'),'normalize activation matrix',TRUE),
                           "Click to fine-tune the coloring of an activation matrices."),
                    tipify(checkboxInput(ns('fa_filtertable'),'filter signficant (tables)',FALSE),
                           "Click to filter the significant entries in the tables.")
                )
            )
        )
        ui
    })
    outputOptions(output, "inputsUI", suspendWhenHidden=FALSE) ## important!!!

    ##================================================================================
    ##======================= OBSERVE FUNCTIONS ======================================
    ##================================================================================
    
    observeEvent( input$fa_info, {
        showModal(modalDialog(
            title = HTML("<strong>Functional Analysis Module</strong>"),
            HTML(fa_infotext),
            easyClose = TRUE, size="l" ))
    })

    observe({
        ngs <- inputData()
        req(ngs)
        ct <- colnames(ngs$model.parameters$contr.matrix)
        ##ct <- c(ct,"<sd>")
        updateSelectInput(session, "fa_contrast", choices=ct )
    })
    

    ##================================================================================
    ## KEGG pathways
    ##================================================================================

    getKeggTable <- reactive({
        ngs = inputData()
        req(ngs)
        req(input$fa_contrast)
        
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
            sendSweetAlert(
                session=session,
                title = "No KEGG terms in enrichment results",
                text="",
                type = "warning")
            dbg("[functional:getKeggTable] no KEGG terms in gsetX")
            df <- data.frame()
            return(df)
        }
        
        jj <- which(!is.na(kegg.ids) &
                    !duplicated(kegg.ids) &
                    kegg.ids %in% kegg.available )
        kegg.gsets <- rownames(ngs$gsetX)[jj]
        kegg.ids <- kegg.ids[jj]

        ## get gene set FC and q-value
        names(ngs$gset.meta$meta)
        meta <- ngs$gset.meta$meta[[comparison]]
        fx <- meta$meta.fx
        qv <- meta$meta.q
        names(fx) <- names(qv) <- rownames(meta)
        fx <- round(fx[kegg.gsets], digits=3)
        qv <- round(qv[kegg.gsets], digits=5)

        ## get gene set FC and q-value
        gene.meta <- ngs$gx.meta$meta[[comparison]]
        gene.fc <- gene.meta$meta.fx
        ##gene.qv <- gene.meta$meta.q
        gene.qv <- unclass(gene.meta$q)
        gene.qv <- apply( gene.qv[,setdiff(colnames(gene.qv),c("ttest","t.test"))],1,max)  ## no t-test...
        genesUPPERCASE <- toupper(sub(".*:","",rownames(gene.meta)))
        names(gene.fc) <- names(gene.qv) <- genesUPPERCASE
        
        ## calculate set size
        sig.genes <- names(gene.fc)
        sig.genes <- names(which(abs(gene.fc) >= 1))
        sig.genes <- names(which(gene.qv < 0.25 & abs(gene.fc) >= 1))
        pw.genes <- lapply(GSETS[kegg.gsets], intersect, names(gene.fc))
        ngene0 <- sapply(pw.genes, length)
        ngene1 <- sapply(lapply(pw.genes, intersect, sig.genes), length)

        sig.genes <- sig.genes[order(-abs(gene.fc[sig.genes]))]
        top.genes <- sapply(GSETS[kegg.gsets], function(x)
            paste0(head(intersect(sig.genes,x),6),collapse=","))
        top.genes <- paste0(top.genes, c("",",...")[1 + 1*(ngene1>6)])
        top.genes    
        
        sig.up <- names(which(gene.qv < 0.25 & gene.fc >= 1))
        sig.dn <- names(which(gene.qv < 0.25 & gene.fc <= -1))
        delta <- sapply( pw.genes, function(gg) mean(gg %in% sig.up) - mean(gg %in% sig.dn) )
        delta <- round(100*delta, digits=2)
        
        ## fast Fisher test
        require(corpora)
        ft.pv <- rep(1, length(ngene1))
        jj <- which(ngene0>0)  ## error for ngene0==0
        if(length(jj)>0 && length(sig.genes)>0) {
            ft.jj <- corpora::fisher.pval(
                                  k1=ngene1[jj], n1=ngene0[jj],
                                  k2=length(sig.genes), n2=length(gene.fc),
                                  alternative="greater")
            ft.pv[jj] <- ft.jj
        }
        ft.r <- round(ngene1/ngene0,digits=2)
        ft.qv <- round(p.adjust(ft.pv, method="fdr"),digits=4)
        ft.res <- data.frame( ratio=ft.r, "k/n"=paste0(ngene1,"/",ngene0), q=ft.qv, check.names=FALSE)

        ##
        kegg.ids2 <- paste0("KEGG:",kegg.ids)
        df <- data.frame( kegg.id=kegg.ids, pathway=kegg.gsets,
                         ##"delta.pct" = delta,
                         ##fx=fx, ## meta.q=qv,
                         ##n0=ngene0, n1=ngene1,
                         ft.res, 
                         top.genes=top.genes, check.names=FALSE)
        ##df <- df[order(-fx),]
        df <- df[order(-df$ratio),]    
        df <- df[!duplicated(df$kegg.id), ]  ## take out duplicated gene sets...
        return(df)
    })

    getFilteredKeggTable <- reactive({
        df <- getKeggTable()
        do.filter = FALSE
        do.filter <- input$fa_filtertable
        if(do.filter) df <- df[which(df$q < 0.999),]
        return(df)
    })

    ## There is a bug in pathview::geneannot.map so we have to override
    ## "Error in mol.sum(gene.data, gene.idmap) : no ID can be mapped!"
    my.geneannot.map <- function(in.ids, in.type, out.type, org = "Hs", pkg.name = NULL, 
                                 unique.map = TRUE, na.rm = TRUE, keep.order = TRUE) 
    {
        if (is.null(pkg.name)) {
            data(bods)
            ridx = grep(tolower(paste0(org, "[.]")), tolower(bods[, 
                                                                  1]))
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
        id.types <- columns(db.obj)
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
###res <- try(suppressWarnings( select(db.obj, keys = in.ids, 
        res <- try(suppressWarnings( AnnotationDbi::select(db.obj, keys = in.ids, 
                                                           keytype = in.type, columns = c(in.type, out.type))))
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
                else umaps = apply(res[, out.type], 2, function(x) tapply(x, 
                                                                          res[, in.type], function(y) paste(unique(y), 
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
        library(pathview)
        unlockBinding("geneannot.map", as.environment("package:pathview"))
        assignInNamespace("geneannot.map", my.geneannot.map, ns="pathview", as.environment("package:pathview"))
        assign("geneannot.map", my.geneannot.map, as.environment("package:pathview"))
        lockBinding("geneannot.map", as.environment("package:pathview"))
    }

    KEGGIMG = paste0(tempfile(),".png")
    output$kegg_graph <- renderImage({

        ngs <- inputData()
        ##NULL.IMG <- list(src=NULL, contentType = 'image/png')
        ##NULL.IMG <- list(src=NA, contentType = 'image/png')
        NULL.IMG <- list(src="", contentType = 'image/png')
        if(is.null(ngs)) return(NULL.IMG)

        require(KEGGgraph)
        require(KEGG.db)
        require(pathview)

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
        sel.row <- input$kegg_table_rows_selected
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
            pw.genes <- GSETS[[as.character(pathway.name)]]
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
        setwd(tmpdir)
        pv.out <- pathview(
            gene.data = fc, pathway.id=pathway.id, gene.idtype="SYMBOL",
            gene.annotpkg = "org.Hs.eg.db", species = "hsa",
            out.suffix="pathview", limit = list(gene=2, cpd=1),
            low = list(gene = "dodgerblue2", cpd = "purple"),
            high = list(gene = "firebrick2", cpd = "yellow"), 
            kegg.dir=xml.dir, kegg.native=TRUE, same.layer=FALSE )
        Sys.sleep(0.2) ## wait for graph

        ## back to previous working folder
        setwd(curwd)   
        dbg("kegg_graph: back to working folder=",getwd(),"\n")
        
        ##width  <- session$clientData$output_kegg_graph_width
        ##height <- session$clientData$output_kegg_graph_height    
        outfile = file.path(tmpdir,paste0("hsa",pathway.id,".pathview.png"))
        dbg("kegg_graph: outfile=",outfile,"\n")
        file.exists(outfile)
        if(!file.exists(outfile)) return(NULL.IMG)

        dbg("kegg_graph: copy",outfile,"to",KEGGIMG,"\n")
        unlink(KEGGIMG,force=TRUE)
        file.copy(outfile,KEGGIMG)
        
        list(src = outfile,
             contentType = 'image/png',
             ##width = 1040*0.8, height = 800*0.9, ## actual size: 1040x800
             ##width = 900, height = 600, ## actual size: 1040x800
             width = "100%", height = "100%", ## actual size: 1040x800         
             alt = "pathview image")
        
    }, deleteFile = TRUE)        


    ##output$kegg_table <- renderDataTable({
    kegg_table.RENDER <- reactive({
        
        ngs <- inputData()
        req(ngs)
        if(is.null(ngs$meta.go)) return(NULL)
        
        comparison=1
        comparison = input$fa_contrast

        if(is.null(comparison)) return(NULL)
        df <- getFilteredKeggTable()
        if(is.null(df)) return(NULL)

        ## add hyperlink
        url = paste0("https://www.genome.jp/kegg-bin/show_pathway?map=hsa",df$kegg.id,"&show_description=show")
        df$kegg.id <- paste0("<a href='",url,"' target='_blank'>",df$kegg.id,"</a>")    
        df$pathway <- wrapHyperLink(df$pathway, df$pathway)
        
        numeric.cols <- which(sapply(df, is.numeric))
        numeric.cols
        
        DT::datatable( df, rownames=FALSE, escape = c(-1,-2),
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "ratio",
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle( df[,"ratio"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })

    kegg_actmap.RENDER <- reactive({
        require(igraph)
        ngs <- inputData()
        req(ngs)
        df <- getKeggTable()
        if(is.null(df) || nrow(df)==0) {
            dbg("[functional:kegg_actmap.RENDER] emtpy KEGG table")
            ##par(mfrow=c(1,1), mar=c(1,1,1,1)*0, oma=c(0,2,0,1)*0 )
            ##frame()
            return(NULL)
        }
        meta <- ngs$gset.meta$meta
        fx <- sapply(meta, function(x) x$meta.fx)
        qv <- sapply(meta, function(x) x$meta.q)    
        rownames(fx) <- rownames(qv) <- rownames(meta[[1]])
        
        kk <- as.character(df$pathway)
        ##if(is.na(kk) || kk=="" || !(kk %in% rownames(fx))) return(NULL)
        if(length(kk) < 3) return(NULL)

        if(mean(is.na(qv))<0.01) {
            score <- fx[kk,,drop=FALSE] * (1 - qv[kk,,drop=FALSE])**2
        } else {
            score <- fx[kk,,drop=FALSE] 
        }
        dim(score)
        if(NCOL(score)==1) score <- cbind(score,score)
        
        score <- score[head(order(-rowSums(score**2)),50),,drop=FALSE]  ## maximum nr of gene sets
        score <- score[,head(order(-colSums(score**2)),25),drop=FALSE]  ## max comparisons/FC    
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        ii=1:nrow(score);jj=1
        jj <- hclust(d2)$order
        ii <- hclust(d1)$order
        score <- score[ii,jj,drop=FALSE]
        
        cex2=0.9
        rownames(score) = substring(rownames(score),1,50)
        if(ncol(score)>15) {
            rownames(score) = substring(rownames(score),1,40)
            cex2=0.8
        }
        if(ncol(score)>25) {
            rownames(score) = substring(rownames(score),1,30)
            colnames(score) <- rep("",ncol(score))
            cex2=0.65
        }
        
        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,2,0,1) )
        require(corrplot)
        score2 <- score
        if(input$fa_normalize) score2 <- t(t(score2) / apply(abs(score2),2,max))
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**3   ## fudging
        rownames(score2) <- tolower(gsub(".*:|kegg_|_Homo.*$","",
                                         rownames(score2),ignore.case=TRUE))
        rownames(score2) <- substring(rownames(score2), 1, 40)
        ## heatmap(score2, scale="none", mar=c(8,20))
        bmar <- 0 + pmax(50 - nrow(score2),0)*0.3
        corrplot( score2, is.corr=FALSE, cl.pos="n", col=BLUERED(100),
                 tl.cex=cex2, tl.col="grey20", mar=c(bmar,0,0,0) )
        
    })    

    kegg_info1 = "<strong>KEGG pathways</strong> are a collection of manually curated pathways representing the current knowledge of molecular interactions, reactions and relation networks as pathway maps. In the pathway map, genes are colored according to their upregulation (red) or downregulation (blue) in the contrast profile. Each pathway is scored for the selected contrast profile and reported in the table below. "
    kegg_buttons1 = plotModuleButtons(
        "kegg_graph", ns=ns, text=kegg_info1, just.info=TRUE, 
        download.fmt="png", ##no.download=TRUE,
        title = "KEGG pathway map", label="a"
    )
    output$kegg_graph_png  <- downloadHandler(
        filename = "plot.png",
        content = function(file) {
            file.copy(KEGGIMG,file)
        }
    )

    kegg_actmap.opts = tagList()
    kegg_actmap_module <- plotModule(
        id="kegg_actmap", ns=ns, func=kegg_actmap.RENDER,
        title = "Activation matrix", label="c",
        info.text = "The <strong>KEGG activation matrix</strong> visualizes the activation levels of pathways (or pathway keywords) across multiple contrast profiles. This facilitates to quickly see and detect the similarities of certain pathways between contrasts. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        ##options = kegg_actmap.opts,
        pdf.width=10, pdf.height=10, res=72
    )
    output <- attachModule(output, kegg_actmap_module)

    kegg_table_info = "<strong>Enrichment table.</strong> The table is interactive; enabling user to sort on different variables and select a pathway by clicking on the row in the table. The scoring is performed by considering the total number of genes in the pathway (n), the number of genes in the pathway supported by the contrast profile (k), the ratio of k/n, and the ratio of |upregulated or downregulated genes|/k. Additionally, the table contains the list of the upregulated and downregulated genes for each pathway and a q value from the Fisher’s test for the overlap."
    
    kegg_table_module <- tableModule(
        id = "kegg_table", ns=ns, label="b",
        func = kegg_table.RENDER,
        info.text = kegg_table_info, 
        info.width="350px",
        title = "Enrichment table"
    )
    output <- attachModule(output, kegg_table_module)

    ##----------------------------------------------------------------------

    kegg_analysis_caption = "<b>KEGG pathway analysis.</b> KEGG pathways are a collection of manually curated pathways representing the current knowledge of molecular interactions, reactions and relation networks as pathway maps. <b>(a)</b> Colored KEGG pathway map. Genes are colored according to their upregulation (red) or downregulation (blue) in the contrast profile. <b>(b)</b> Table reporting (Fisher's exact) enrichment score for each pathway for the selected contrast profile. <b>(c)</b> Activation matrix visualizing activation levels of pathways across multiple contrast profiles." 

    output$kegg_analysis_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1,NA),
            fillRow(
                flex = c(1.3,0.04,1),
                height = 650,
                fillCol(
                    flex = c(NA,1.8,0.06,1),
                    height = 600,
                    kegg_buttons1,
                    imageOutput(ns("kegg_graph"), height="100%"),
                    br(),
                    moduleWidget(kegg_table_module, outputFunc="dataTableOutput", ns=ns)
                ),
                br(),  ## horizontal space
                moduleWidget(kegg_actmap_module, ns=ns)
            ),
            div(HTML(kegg_analysis_caption), class="caption")
        )
    })

    ##================================================================================
    ## GO graph
    ##================================================================================
    
    require(visNetwork)

    GO_network.RENDER <- reactive({

        require(igraph)
        require(RColorBrewer)
        require(visNetwork)
        ngs <- inputData()
        req(ngs)

        comparison=1;methods=c("fisher","gsva","camera")
        comparison = input$fa_contrast
        req(input$fa_contrast)
        if(is.null(comparison)) return(NULL)
        
        ##sub2 <- getSigGO(comparison, methods, nterms=250, ntop=25, ngs=ngs)
        sub2 <- go <- ngs$meta.go$graph
        if(is.null(go)) {
            sendSweetAlert(
                session=session,
                title = "No GO graph in enrichment results",
                text="",
                type = "warning")
            dbg("[functional:GO_network.RENDER] no META.GO in pgx object!")
            return(NULL)
        }

        score = ngs$meta.go$pathscore[,comparison]        
        score = (score/max(abs(score),na.rm=TRUE))
        score[is.na(score)] = 0
        V(sub2)$value <- score
        V(sub2)$color <- bluered(32)[16 + round(15*score)]
        V(sub2)$label <- V(sub2)$Term
        V(sub2)$label[which(is.na(score)|score==0)] = ""
        pos = sub2$layout
        
        ##if("prune" %in% input$GO_options) {
        if(input$GO_prunetree) {
            ##cat("pruning GO graph\n")
            vv = V(sub2)[which(!is.na(score) & score!=0)]
            sp = shortest_paths(sub2, from="all", to=vv, mode="all", output="vpath")
            sp.vv = unique(unlist(sp$vpath))
            sub2 = induced.subgraph(sub2, sp.vv)
            pos = layout_with_fr(sub2)
            score = score[V(sub2)$name]
        }

        ## remove root?
        removeroot=TRUE
        if(removeroot) {
            sub2 <- induced_subgraph(sub2, which(V(sub2)$name!="all"))
            ##if("prune" %in% input$GO_options) pos = layout_with_fr(sub2)
            if(input$GO_prunetree) pos = layout_with_fr(sub2)
            score <- score[V(sub2)$name]
            ##pos <- pos[V(sub2)$name,]        
        }
        roots <- c("all",neighbors(go, V(go)["all"], mode="all" )$name)
        roots <- intersect(roots, V(sub2)$name)
        
        astree = TRUE
        if(astree) {
            if("all" %in% V(sub2)$name) {
                pos = layout_as_tree(sub2, root="all", mode="all")
            } else {
                pos = layout_as_tree(sub2, root=roots, mode="all")
            }
            pos[,2] = -pos[,2]
        }

        ## color clusters
        if(input$GO_colorclusters) {
            ##clust = cluster_louvain(as.undirected(sub2))$membership
            clust = cluster_louvain(as.undirected(go))$membership
            names(clust) = V(go)$name
            cc = c(brewer.pal(12,"Set3"),brewer.pal(8,"Set2"),brewer.pal(8,"Set1"))
            V(sub2)$color = rep(cc,99)[clust[V(sub2)$name]]
            jj = which(is.na(score) | score==0)
            if(length(jj)>0) V(sub2)$color[jj] = NA
        }
        
        require(visNetwork)
        ##pos <- pos[V(sub2)$name,]
        gr = toVisNetworkData(sub2)
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
        ##if("prune" %in% input$GO_options) {
        if(input$GO_prunetree) {
            font.size=20; cex=0.6
        }

        visNetwork(gr$nodes, gr$edges) %>%
            visEdges(smooth = FALSE, hidden=FALSE, arrows=list(enabled=TRUE),
                     scaling=list(min=10*cex, max=30*cex), width=5*cex )  %>%
            visNodes(font = list( size=font.size*cex, vadjust=0),
                     scaling=list(min=1*cex, max=80*cex))  %>%
            visPhysics(stabilization = FALSE)  %>%
            ##visInteraction(hideEdgesOnDrag = TRUE) %>%
            ##visInteraction(navigationButtons = TRUE) %>%
            ##visOptions(nodesIdSelection = TRUE) %>%
            ##visOptions(selectedBy="component") %>%
            visOptions(highlightNearest = list(enabled=T, degree=1, hover=TRUE)) %>%
            ##visEvents(select="function(nodes){Shiny.onInputChange('current_node_id',nodes.nodes);;}") %>%
            visPhysics(enabled=FALSE) 
    })    


    GO_table.RENDER <- reactive({
        ngs <- inputData()
        req(ngs, input$fa_contrast)
        if(is.null(ngs$meta.go)) {
            dbg("[functional:GO_table.RENDER] no META.GO in pgx object!")
            return(NULL)
        }
        
        comparison=1
        comparison = input$fa_contrast
        ##req(input$fa_contrast)
        if(is.null(comparison)) return(NULL)

        go = ngs$meta.go$graph
        scores <- ngs$meta.go$pathscore[,comparison]        
        scores <- scores[which(!is.na(scores))]
        scores <- round(scores, digits=3)
        scores <- sort(scores, decreasing=TRUE)
        go.term = V(go)[names(scores)]$Term
        qv=fx=NULL
        if("qvalue" %in% names(ngs$meta.go)) qv = ngs$meta.go$qvalue[names(scores),comparison]
        if("foldchange" %in% names(ngs$meta.go)) fx = ngs$meta.go$foldchange[names(scores),comparison]
        
        go.term = substring(go.term, 1, 80)
        dt1 = round( cbind(score=scores, meta.fx=fx, meta.q=qv), digits=4)    
        dt = data.frame( id=names(scores), term=go.term, dt1, stringsAsFactors=FALSE)    
        id2 = paste0("abc(",sub(":","_",dt$id),")")  ## to match with wrapHyperLink
        dt$id <- wrapHyperLink(as.character(dt$id), id2)  ## add link
        
        numeric.cols <- which(sapply(dt, is.numeric))
        numeric.cols

        DT::datatable( dt, rownames=FALSE, escape = c(-1,-2),
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      options=list(
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "score",
                                ##background = styleColorBar(c(0,3), 'lightblue'),
                                background = color_from_middle( dt1[,"score"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })

    GO_actmap.RENDER <- reactive({
        require(igraph)
        ngs <- inputData()
        req(ngs)

        if(is.null(ngs$meta.go)) {
            dbg("[functional:GO_actmap.RENDER] no META.GO in pgx object!")
            return(NULL)
        }

        score = ngs$meta.go$pathscore
        ##if(is.null(score)) return(NULL)
        go = ngs$meta.go$graph
        score[is.na(score)] = 0
        rownames(score) = V(go)[rownames(score)]$Term
        cat("GO_actmap:: dim(score)=",dim(score),"\n")
        
        ## reduce score matrix
        ##score = head(score[order(-rowSums(abs(score))),],40)
        score = score[head(order(-rowSums(score**2)),50),,drop=FALSE] ## max number of terms    
        score = score[,head(order(-colSums(score**2)),25),drop=FALSE] ## max comparisons/FC
        if(NCOL(score)==1) score <- cbind(score,score)   
        cat("GO_actmap:: dim(score)=",dim(score),"\n")

        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        jj=1;ii=1:nrow(score)
        jj <- hclust(d2)$order
        ii <- hclust(d1)$order
        score <- score[ii,jj,drop=FALSE]
        
        cex2=1
        colnames(score) = substring(colnames(score),1,30)
        rownames(score) = substring(rownames(score),1,50)
        if(ncol(score)>15) {
            rownames(score) = substring(rownames(score),1,40)
            cex2=0.85
        }
        if(ncol(score)>25) {
            rownames(score) = substring(rownames(score),1,30)
            colnames(score) <- rep("",ncol(score))
            cex2=0.7
        }

        ##pdf("module-functional.pdf",w=8,h=12)
        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,2,0,1))
        require(corrplot)
        score2 <- score
        if(input$fa_normalize) score2 <- t( t(score2) / apply(abs(score2),2,max)) ## column scale???
        ##score2 <- abs(score2)**0.8 * sign(score2)  ## fudging
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**0.8   ## fudging
        ## heatmap(score2, scale="none", mar=c(8,20))
        bmar <- 0 + pmax((50 - nrow(score2))*0.25,0)
        corrplot( score2, is.corr=FALSE, cl.pos="n", col=BLUERED(100),
                 tl.cex=cex2, tl.col="grey20", mar=c(bmar,0,0,0) )
        ##dev.off()
        
    })    


    GO_info1 = "The <strong>Gene Ontology</strong> (GO) provides a computational representation of the current knowledge about roles of genes for many organisms in terms of molecular functions, cellular components and biological processes. The structure of GO can be described in terms of a graph, where each GO term is a node, and the relationships between the terms are edges between the nodes. GO is loosely hierarchical, with ‘child’ terms being more specialized than their ‘parent’ terms. The graph is interactive. You can move the graph and zoom in using the mouse."

    GO_network.opts = tagList(
        tipify( checkboxInput(ns("GO_prunetree"), "Prune tree", TRUE),
               "Prune the tree with only significant branches."),
        tipify( checkboxInput(ns("GO_colorclusters"), "Color clusters", TRUE),
               "Highlight clusters with different colors.")
    )
    GO_network_module <- plotModule(
        "GO_network", GO_network.RENDER, plotlib="visnetwork",
        title = "Gene Ontology graph", label="a",
        info.text = GO_info1,
        download.fmt = c("pdf","png"), ## no.download=TRUE,
        options = GO_network.opts,
        pdf.width=10, pdf.height=10, res=80
    )
    output <- attachModule(output, GO_network_module)

    GO_actmap.opts = tagList()
    GO_actmap_module <- plotModule(
        "GO_actmap", GO_actmap.RENDER,
        title = "Activation matrix", label="c",
        info.text = "The Gene Ontology (GO) activation matrix visualizes the activation of GO terms across conditions. From this figure, you can easily detect GO terms that are consistently up/down across conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = GO_actmap.opts,
        pdf.width=10, pdf.height=10, res=72
    )
    output <- attachModule(output, GO_actmap_module)

    GO_table_module <- tableModule(
        id = "GO_table", label="b",
        func = GO_table.RENDER,
        info.text="<strong>GO score table.</strong> The scoring of a GO term is performed by considering the cumulative score of all terms from that term to the root node. That means that GO terms that are supported by higher level terms levels are preferentially scored.", 
        title = "GO score table"
    )
    output <- attachModule(output, GO_table_module)

    GO_analysis_caption = "<b>Gene Ontology (GO) analysis.</b> GO describes our knowledge of the biological domain of genes with respect to three aspects: molecular function, cellular component and biological process. <b>(a)</b> Hierarchical graph representing the enrichment of the GO terms as a tree structure. <b>(b)</b> GO scoring table. The score of a GO term is the cumulative score of all higher order terms. <b>(c)</b> Activation matrix visualizing the enrichment of GO terms across multiple contrast profiles."
    
    output$GO_analysis_UI <- renderUI({
        fillCol(
            flex=c(1,NA),
            height = rowH,
            fillRow(
                height = 650,
                flex = c(1.2,1),
                fillCol(
                    flex = c(2,1),
                    height = 600,
                    moduleWidget(GO_network_module, outputFunc="visNetworkOutput", ns=ns),
                    moduleWidget(GO_table_module, outputFunc="dataTableOutput", ns=ns)
                ),
                moduleWidget(GO_actmap_module, ns=ns)
            ),
            div(HTML(GO_analysis_caption),class="caption")
        )
    })

    
    ##================================================================================
    ## Drug signature enrichment analysis L1000
    ##================================================================================

    getDseaTable <- reactive({
        ngs <- inputData()
        req(ngs)        
        req(input$fa_contrast, input$dsea_monocombo)

        comparison=3
        names(ngs$gx.meta$meta)
        comparison = input$fa_contrast
        if(is.null(comparison)) return(NULL)
        
        dmethod="combo"
        dmethod="mono"
        dmethod <- input$dsea_monocombo
        if(is.null(dmethod)) return(NULL)
        
        fc <- ngs$gx.meta$meta[[comparison]]$meta.fx
        names(fc) <- rownames(ngs$gx.meta$meta[[1]])
        nes <- round(ngs$drugs[[dmethod]]$X[,comparison],4)
        pv  <- round(ngs$drugs[[dmethod]]$P[,comparison],4)
        qv  <- round(ngs$drugs[[dmethod]]$Q[,comparison],4)
        drug <- rownames(ngs$drugs[[dmethod]]$X)
        stats <- ngs$drugs[[dmethod]]$stats
        nes[is.na(nes)] <- 0
        qv[is.na(qv)] <- 1
        pv[is.na(pv)] <- 1
        
        ## SHOULD MAYBE BE DONE IN PREPROCESSING....
        descr0 <- read.csv(file.path(FILES,"L1000_repurposing_drugs.txt"),
                           sep="\t", comment.char="#")

        if(dmethod=="combo") {
            drugs <- strsplit(names(nes),split="[+]")
            drug1 <- sapply(drugs,"[",1)
            drug2 <- sapply(drugs,"[",2)
            j1 <- match( drug1, descr0$pert_iname)
            j2 <- match( drug2, descr0$pert_iname)
            cmoa <- paste( descr0[j1,"moa"],"+",descr0[j2,"moa"])
            ctarget <- paste( descr0[j1,"target"],"+",descr0[j2,"target"])
            descr <- data.frame( MOA=cmoa, target=ctarget)
        } else {
            jj <- match( drug, descr0$pert_iname)
            descr <- descr0[jj,c("moa","target")]
        }
        
        res <- data.frame( drug=drug, NES=nes, pval=pv, padj=qv, descr)
        res <- res[order(-abs(res$NES)),]
        
        return(res)
    })


    dsea_enplots.RENDER <- reactive({

        ngs <- inputData()
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))        
        req(input$fa_contrast, input$dsea_monocombo)

        comparison=1
        comparison = input$fa_contrast
        if(is.null(comparison)) return(NULL)

        res <- getDseaTable()

        dmethod="mono"
        dmethod="combo"
        dmethod <- input$dsea_monocombo

        ## rank vector for enrichment plots
        rnk <- ngs$drugs[[dmethod]]$stats[,comparison]
        dctype <- sub("_.*$","",names(rnk))
        ##table(rownames(res) %in% dctype)
        ##table(sapply(rownames(res), function(g) sum(grepl(g,names(rnk),fixed=TRUE))))
        
        ## ENPLOT TYPE
        itop <- c( head(order(-res$NES),10), tail(order(-res$NES),10))
        par(oma=c(0,1,0,0))
        par(mfrow=c(4,5), mar=c(1,1.5,1.8,1))
        i=1
        for(i in itop) {
            dx <- rownames(res)[i]
            dx
            gmtdx <- grep(dx,names(rnk),fixed=TRUE,value=TRUE)  ## L1000 naming allows this...
            length(gmtdx)
            ##if(length(gmtdx) < 3) { frame(); next }
            gsea.enplot( rnk, gmtdx, main=dx, cex.main=1.25)
            nes <- round(res$NES[i],2)
            qv  <- round(res$padj[i],3)
            tt <- c( paste("NES=",nes), paste("q=",qv) )
            legend("topright", legend=tt, cex=0.9)
        }
        
    })    

    dsea_moaplot.RENDER <- reactive({
        ngs <- inputData()

        req(ngs, input$fa_contrast, input$dsea_monocombo)
        
        if(is.null(ngs$drugs)) return(NULL)
        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        
        comparison=1
        comparison = input$fa_contrast
        if(is.null(comparison)) return(NULL)

        res <- getDseaTable()

        dmethod="mono"
        dmethod="combo"
        dmethod <- input$dsea_monocombo
        
        j1 <- which( res$padj < 0.2 & res$NES > 0)
        j2 <- which( res$padj < 0.2 & res$NES < 0)
        moa.pos <- strsplit(as.character(res$moa[j1]), split="\\|")
        moa.neg <- strsplit(as.character(res$moa[j2]), split="\\|")
        moa <- strsplit(as.character(res$moa), split="\\|")
        fx <- mapply(function(x,n) rep(x,n), res$NES, sapply(moa,length))
        moa.avg <- sort(tapply( unlist(fx), unlist(moa), mean))
        moa.sum <- sort(tapply( unlist(fx), unlist(moa), sum))
        head(moa.pos)
        head(moa.neg)
        moa.pos <- sort(table(unlist(moa.pos)),decreasing=TRUE)
        moa.neg <- sort(table(unlist(moa.neg)),decreasing=TRUE)
        
        dtg.pos <- strsplit(as.character(res$target[j1]), split="\\|")
        dtg.neg <- strsplit(as.character(res$target[j2]), split="\\|")
        dtg <- strsplit(as.character(res$target), split="\\|")
        dx <- mapply(function(x,n) rep(x,n), res$NES, sapply(dtg,length))
        dtg.avg <- sort(tapply( unlist(dx), unlist(dtg), mean))
        dtg.sum <- sort(tapply( unlist(dx), unlist(dtg), sum))    
        dtg.pos <- sort(table(unlist(dtg.pos)),decreasing=TRUE)
        dtg.neg <- sort(table(unlist(dtg.neg)),decreasing=TRUE)
        head(dtg.pos)
        head(dtg.neg)

        NTOP=10
        if(1) {
            moa.top <- sort(c( head(moa.pos,NTOP), -head(moa.neg,NTOP)),decreasing=TRUE)
            dtg.top <- sort(c( head(dtg.pos,NTOP), -head(dtg.neg,NTOP)),decreasing=TRUE)

            ##layout(matrix(1:2,nrow=1),widths=c(1.4,1))
            ##par(mfrow=c(2,1))
            par(mar=c(4,15,5,0.5), mgp=c(2,0.7,0))

            par(mfrow=c(1,1))

            if(input$dsea_moatype=="drug class") {
                par(mar=c(15,4,3,0.5), mgp=c(2,0.7,0))
                barplot(moa.top, horiz=FALSE, las=3, ylab="drugs (n)")
                ##title(main="MOA", line=1 )
            } else {
                par(mar=c(15,4,3,0.5), mgp=c(2,0.7,0))
                barplot(dtg.top, horiz=FALSE, las=3, ylab="drugs (n)")
                ##title(main="target gene", line=1 )
            }
        }
        
    })    

    dsea_table.RENDER <- reactive({
        ngs <- inputData()
        req(ngs)
        if(is.null(ngs$drugs)) return(NULL)
        
        res <- getDseaTable()
        req(res)
        res$moa <- shortstring(res$moa,50)
        res$target <- shortstring(res$target,50)
        
        ## limit number of results??
        ##jj <- unique(c( head(order(-res$NES),250), tail(order(-res$NES),250)))
        jj <- unique(c( head(order(-res$NES),1000), tail(order(-res$NES),1000)))
        jj <- jj[order(-abs(res$NES[jj]))]
        res <- res[jj,]
        
        DT::datatable( res, rownames=FALSE,
                      class = 'compact cell-border stripe hover',                  
                      extensions = c('Scroller'),
                      selection=list(mode='single', target='row', selected=1),
                      options=list(
                          ##dom = 'Blfrtip', buttons = c('copy','csv','pdf'),
                          dom = 'lfrtip', 
                          scrollX = TRUE, ##scrollY = TRUE,
                          scrollY = tabH, scroller=TRUE, deferRender=TRUE
                      )  ## end of options.list 
                      ) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
                DT::formatStyle( "NES",
                                background = color_from_middle( res[,"NES"], 'lightblue', '#f5aeae'),
                                backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                                backgroundPosition = 'center') 
    })

    dsea_actmap.RENDER <- reactive({
        require(igraph)
        ngs <- inputData()
        req(ngs, input$fa_contrast, input$dsea_monocombo)

        shiny::validate(need("drugs" %in% names(ngs), "no 'drugs' in object."))    
        if(is.null(ngs$drugs)) return(NULL)
        
        dmethod="mono"
        dmethod <- input$dsea_monocombo
        comparison=1
        comparison = input$fa_contrast
        if(is.null(comparison)) return(NULL)
        
        nes <- ngs$drugs[[dmethod]]$X
        qv <- ngs$drugs[[dmethod]]$Q
        score <- nes * (1 - qv)**2
        score[is.na(score)] <- 0
        if(NCOL(score)==1) score <- cbind(score,score)
        
        ## reduce score matrix
        ##score = head(score[order(-rowSums(abs(score))),],40)
        ##score = score[head(order(-rowSums(score**2)),50),] ## max number of terms
        score = score[head(order(-score[,comparison]**2),50),,drop=FALSE] ## max number of terms    
        score = score[,head(order(-colSums(score**2)),25),drop=FALSE] ## max comparisons/FC

        cat("dsea_actmap:: dim(score)=",dim(score),"\n")
        score <- score + 1e-3*matrix(rnorm(length(score)),nrow(score),ncol(score))
        d1 <- as.dist(1-cor(t(score),use="pairwise"))
        d2 <- as.dist(1-cor(score,use="pairwise"))
        d1[is.na(d1)] <- 1
        d2[is.na(d2)] <- 1
        jj=1;ii=1:nrow(score)
        ii <- hclust(d1)$order
        jj <- hclust(d2)$order
        score <- score[ii,jj,drop=FALSE]
        
        cex2=1
        colnames(score) = substring(colnames(score),1,30)
        rownames(score) = substring(rownames(score),1,50)
        if(ncol(score)>15) {
            rownames(score) = substring(rownames(score),1,40)
            cex2=0.85
        }
        if(ncol(score)>25) {
            rownames(score) = substring(rownames(score),1,30)
            colnames(score) <- rep("",ncol(score))
            cex2=0.7
        }

        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,2,0,1))
        require(corrplot)
        score2 <- score
        if(input$fa_normalize) score2 <- t( t(score2) / apply(abs(score2),2,max)) 
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**3   ## fudging
        bmar <- 0 + pmax((50 - nrow(score2))*0.25,0)
        corrplot( score2, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex=cex2, tl.col="grey20", mar=c(bmar,0,0,0) )
        
    })    

    
    ##--------- DSEA enplot plotting module
    dsea_enplots.opts = tagList(
        tipify( radioButtons(ns('dsea_monocombo'),"Analysis type:",c("mono","combo"),inline=TRUE),
               "Select type of drug enrichment analysis: mono or combo (if available). ")
    )
    dsea_enplots_module <- plotModule(
        "dsea_enplots", dsea_enplots.RENDER,
        title = "Drug profile enrichment", label="a",
        info.text = "The <strong>Drug Connectivity Map</strong> correlates your signature with more than 5000 known drug profiles from the L1000 database, and shows the top N=10 similar and opposite profiles by running the GSEA algorithm on the contrast-drug profile correlation space.",
        options = dsea_enplots.opts,
        pdf.width=11, pdf.height=7, res=72
    )
    output <- attachModule(output, dsea_enplots_module)


    ##---------- DSEA Activation map plotting module
    dsea_moaplot.opts = tagList(
        tipify( radioButtons(ns('dsea_moatype'),'Plot type:',c("drug class","target gene"),inline=TRUE),
               "Select plot type of MOA analysis: by class description or by target gene.")
    )
    dsea_moaplot_module <- plotModule(
        "dsea_moaplot", dsea_moaplot.RENDER,
        title = "Mechanism of action", label="c",
        info.text = "This plot visualizes the <strong>mechanism of action</strong> (MOA) across the enriched drug profiles. On the vertical axis, the number of drugs with the same MOA are plotted. You can switch to visualize between MOA or target gene.",
        options = dsea_moaplot.opts,
        pdf.width=4, pdf.height=6, res=72
    )
    ##output$dsea_moaplot     <- dsea_moaplot_module$render
    ##output$dsea_moaplot_pdf <- dsea_moaplot_module$pdf
    output <- attachModule(output, dsea_moaplot_module)

    ##-------- Activation map plotting module
    dsea_actmap.opts = tagList()
    dsea_actmap_module <- plotModule(
        "dsea_actmap", dsea_actmap.RENDER,
        title = "Activation matrix", label="d",
        info.text = "The <strong>Activation Matrix</strong> visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = dsea_actmap.opts,
        pdf.width=6, pdf.height=10, res=72
    )
    output <- attachModule(output, dsea_actmap_module)

    ##--------buttons for table
    dsea_table_module <- tableModule(
        id = "dsea_table", label="b",
        func = dsea_table.RENDER,
        info.text="Drug profile enrichment table. Enrichment is calculated by correlating your signature with more than 5000 known drug profiles from the L1000 database. Because the L1000 has multiple perturbation experiment for a single drug, drugs are scored by running the GSEA algorithm on the contrast-drug profile correlation space. In this way, we obtain a single score for multiple profiles of a single drug.", 
        title = "Profile enrichment table"
    )
    output <- attachModule(output, dsea_table_module)


    observe({
        ngs <- inputData()
        req(ngs)
        ct <- c("mono","combo")
        if(is.null(ngs$drugs$combo)) {
            ct <- c("mono")
        }
        updateRadioButtons(session, "dsea_monocombo", choices=ct, inline=TRUE)
    })
       
    ##-----------------------------------------
    ## Page layout
    ##-----------------------------------------

    dsea_analysis_caption = "<b>Drug Connectivity Map.</b> Drug CMap correlates your signature with more than 5000 known drug perturbation profiles from the L1000 database. <b>(a)</b> Figure showing the top N=10 similar and opposite profiles by running the GSEA algorithm on the contrast-drug profile correlation space. <b>(b)</b> Table summarizing the statistical results of the drug enrichment analysis. <b>(c)</b> Mechanism-of-action plot showing the top most frequent drug class (or target genes) having similar or opposite enrichment compared to the query signature. <b>(d)</b> Activation matrix visualizing enrichment levels of drug signatures across multiple contrast profiles." 

    output$DSEA_analysis_UI <- renderUI({
        fillCol(
            flex = c(1,NA),
            height = rowH,
            fillRow(
                height = 650,
                flex = c(2.6,1), 
                fillCol(
                    flex = c(1.3,0.15,1),
                    height = 650,
                    fillRow(
                        flex=c(2.2,1),
                        moduleWidget(dsea_enplots_module, ns=ns),
                        moduleWidget(dsea_moaplot_module, ns=ns)
                    ),
                    br(),  ## vertical space
                    moduleWidget(dsea_table_module, outputFunc="dataTableOutput", ns=ns)        
                ),
                moduleWidget(dsea_actmap_module, ns=ns)        
            ),
            div(HTML(dsea_analysis_caption),class="caption")
        )
    })


    ##---------------------------------------------------------------
    ##------------- Functions for WordCloud ------------------------
    ##---------------------------------------------------------------

    pgx.calculateWordFreq <- function(ngs, progress=NULL, pg.unit=1) {

        if(!is.null(progress)) progress$set(message = "WordCloud", value = 0)
        
        ## get gset meta foldchange-matrix
        S <- sapply( ngs$gset.meta$meta, function(x) x$meta.fx)
        rownames(S) <- rownames(ngs$gset.meta$meta[[1]])
        S <- S[order(-apply(S,1,sd)),]
        S <- S[order(-rowMeans(S**2)),]

        ## exclude down, GSE gene sets??????
        S <- S[grep("dn|down|^gse",rownames(S),ignore.case=TRUE,invert=TRUE),]
        
        if(!is.null(progress)) progress$inc(0.2*pg.unit, detail="calculating word frequencies")
        
        ## Determine top most frequent terms
        sname <- gsub(".*:","",tolower(rownames(S)))
        sname <- gsub("b.cell","bcell",sname)
        sname <- gsub("t.cell","tcell",sname)
        words <- strsplit(sname,split="[-_ ]")    
        names(words) <- rownames(S)
        terms <- names(sort(-table(unlist(words))))
        stopwords = strsplit(
            "with int like strand tconv pid lee and or mouse dn small big human homo sapiens mus musculus drug early late hsa gse culture in the a g for line r up down events anti large targets tissue vitro process cells ctrl regulation processing common pathway months days white pre post process all mice from", split=" ")[[1]]
        terms <- terms[which(!terms %in% stopwords)]
        terms <- terms[sapply(terms,nchar)>2]
        terms <- grep("[0-9]|^\\(",terms,invert=TRUE,value=TRUE)        
        length(terms)
        terms <- head(terms,1000)
        
        ## Calculate incidence matrix
        words2 <- lapply(words, function(w) intersect(w, terms))
        words2 <- words2[sapply(words2,length)>0]
        idx <- lapply(1:length(words2), function(i) cbind(i,match(words2[[i]],terms)))
        idx <- do.call(rbind, idx)
        library(qlcMatrix)
        W <- sparseMatrix(idx[,1], idx[,2], x=1)
        dim(W)
        rownames(W) = names(words2)
        colnames(W) = terms
        
        ## filter on minimal size and maximum ratio
        nn <- Matrix::colSums(W,na.rm=TRUE)
        nr <- nn / nrow(W)
        W <- W[,which(nn >= 10 & nr <= 0.10)]
        dim(W)
        
        if(!is.null(progress)) progress$inc(0.3*pg.unit, detail="computing GSEA")        

        ## compute for average contrast
        require(fgsea)
        rms.FC <- Matrix::rowMeans(S[rownames(W),]**2)**0.5
        rms.FC <- rms.FC + 0.01*rnorm(length(rms.FC))
        gmt <- apply(W,2,function(x) names(which(x!=0)))
        res <- fgsea( gmt, rms.FC, nperm=1000 )
        res$leadingEdge <- sapply(res$leadingEdge,paste,collapse="//")
        ## res$leadingEdge <- NULL
        colnames(res)[1] <- "word"
        
        ## --------- only significant and positive
        ##res <- res[(res$padj < 0.20 & res$NES>0),]
        res <- res[(res$padj < 1 & res$NES>0),]
        res <- res[order(-abs(res$NES)),]
        dim(res)

        ## now compute significant terms for all contrasts
        all.gsea <- list()
        for(i in 1:ncol(S)) {
            fc <- as.vector(S[rownames(W),i])
            names(fc) <- rownames(W)
            fc <- fc + 0.01*rnorm(length(fc))
            gmt1 <- gmt[as.character(res$word)]
            res1 <- fgsea( gmt1, fc, nperm=1000 )
            res1$leadingEdge <- sapply(res1$leadingEdge,paste,collapse="//")
            ## res$leadingEdge <- NULL
            colnames(res1)[1] <- "word"
            all.gsea[[colnames(S)[i]]] <- res1
        }
        all.gsea[["rms.FC"]] <- res
        
        if(!is.null(progress)) progress$inc(0.25*pg.unit, detail="clustering")

        require(Rtsne)
        library(umap)
        pos1 = Rtsne(as.matrix(t(W)),perplexity=10,check_duplicates=FALSE)$Y
        pos2 = umap(as.matrix(t(W)))$layout
        rownames(pos1) = rownames(pos2) = colnames(W)
        colnames(pos1) = colnames(pos2) = c("x","y")
        pos1 = pos1[res$word,]
        pos2 = pos2[res$word,]
        
        res = list(gsea=all.gsea, tsne=pos1, umap=pos2)
        res
    }

    enrich_getWordFreqResults <- reactive({
        ngs <- inputData()
        req(ngs)
        
        ## ------------ Create a Progress object
        progress <- shiny::Progress$new()
        on.exit(progress$close())    

        res <- pgx.calculateWordFreq(ngs, progress=progress, pg.unit=1)    
        return(res)
    })


    enrich_getWordEnrichment <- reactive({

        res <- enrich_getWordFreqResults()
        req(res, input$fa_contrast)

        contr=1
        contr <- input$fa_contrast
        gsea1 <- res$gsea[[ contr ]]
        topFreq <- data.frame( gsea1, tsne=res$tsne, umap=res$umap)
        
        ## update selectors
        words <- sort(res$gsea[[1]]$word)
        updateSelectInput(session, "enrich_wordcloud_exclude", choices=words)
        
        return(topFreq)
    })

    enrich_wordtsne.RENDER <- reactive({

        topFreq <- enrich_getWordEnrichment()

        df <- topFreq
        klr = ifelse( df$padj<=0.05, "red", "grey")    
        ps1 = 0.5 + 3*(1-df$padj)*(df$NES/max(df$NES))**3

        ## label top 20 words
        df$label <- rep(NA, nrow(df))
        jj <- head(order(-abs(df$NES)),20)
        df$label[jj] <- as.character(df$word[jj])
        cex=1
        ##cex=2.5
        
        require(ggrepel)
        if(input$enrich_wordtsne_algo=="tsne") {
            p <- ggplot(df, aes(tsne.x, tsne.y, label=label))
        } else {
            p <- ggplot(df, aes(umap.x, umap.y, label=label))
        }
        p <- p +
            geom_point( size=cex*ps1, color=klr) +
            geom_text_repel(size=4*cex) +
            ##geom_text_repel(point.padding=NA, size=cex) +
            ##scale_x_continuous( expand=c(0,0) ) +
            ##scale_y_continuous( expand=c(0,0) ) +
            ##coord_cartesian( xlim=c(0,1), ylim=c(0,1)) +
            theme_bw() +
            theme( axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks=element_blank()) 

        return(p)
    })

    enrich_wordtsne.PLOTLY <- reactive({

        topFreq <- enrich_getWordEnrichment()
        
        df <- topFreq
        klr = ifelse( df$padj<=0.05, "red", "grey")    
        ps1 = 0.5 + 3*(1-df$padj)*(df$NES/max(df$NES))**3

        ## label top 20 words
        df$label <- rep("", nrow(df))
        jj <- head(order(-abs(df$NES)),20)
        df$label[jj] <- as.character(df$word[jj])
        cex=1
        ##cex=2.5
        df$abs.NES <- abs(df$NES)**2
        
        require(ggrepel)
        if(input$enrich_wordtsne_algo=="tsne") {
            pos = cbind( x=df$tsne.x, y=df$tsne.y)
        } else {
            pos = cbind( x=df$umap.x, y=df$umap.y)
        }
        
        plt <- plot_ly(
            df,
            text = df$word, hoverinfo = 'text'
            ## hovertemplate = paste0("%{text}<br>NES: %{NES}<extra> </extra>")
        ) %>%
            add_markers(
                type="scatter",
                x = pos[,1], y = pos[,2], 
                color = klr,
                size = ~abs.NES,
                ## sizes = c(5,100),
                marker = list(
                    ##size = 16,
                    ##sizes=c(20,400),
                    line = list(color="grey20", width=0.6)
                )) %>%
            add_annotations(
                x = pos[,1], y = pos[,2],
                text = df$label,
                font = list(size=12),
                ##xref = "x", yref = "y",
                showarrow = FALSE)

        ax <- list(
            title = "",
            showticklabels = FALSE,
            showgrid = FALSE
        )

        m <- list(
            l = 0,
            r = 0,
            b = 0,
            t = 0,
            pad = 4
        )

        plt <- plt %>%
            layout(
                xaxis = ax,
                yaxis = ax,
                showlegend = FALSE,
                margin = m
            )
        
        return(plt)
    })

    enrich_wordcloud.RENDER <- reactive({

        topFreq <- enrich_getWordEnrichment()
        df <- topFreq
        
        excl.words <- ""
        excl.words <- input$enrich_wordcloud_exclude
        ##cat("<enrich_wordcloud> 0: excl.words=",excl.words,"\n")
        ##cat("<enrich_wordcloud> 0: len.excl.words=",length(excl.words),"\n")
        if(length(excl.words)>0) {
            df <- df[ which(!df$word %in% excl.words), ]
        }
        
        cex1 <- 1+round((5*rank(abs(df$NES))/nrow(df))**2)    
        cex2 <- (-log10(df$padj))**1.0
        size <- 10*abs(cex1 * cex2)**1
        minsize <- tail(sort(size),250)[1]

        color.pal = input$enrich_wordcloud_colors
        
        ##require(wordcloud2)
        ##require(rWordCloud)
        require(wordcloud)
        ##wordcloud2( data.frame(word=df$word, size=size), size=1)
        ##d3Cloud(text = df$word, size = size)
        
        par(mar=c(1,1,1,1)*0)
        wordcloud(words = df$word, freq = size,
                  ##colors=brewer.pal(8, "Dark2"),
                  colors = brewer.pal(8, color.pal),
                  scale=c(2,0.1)*0.9, min.freq=minsize)
        

    })

    enrich_keyword.RENDER <- reactive({

        ngs <- inputData()
        req(ngs)

        topFreq <- enrich_getWordEnrichment()
        
        ## get gset meta foldchange-matrix
        S <- sapply( ngs$gset.meta$meta, function(x) x$meta.fx)
        S <- S + 0.0001*matrix(rnorm(length(S)),nrow(S),ncol(S))
        rownames(S) <- rownames(ngs$gset.meta$meta[[1]])
        
        keyword = "lipid"
        keyword = "apoptosis"
        keyword = "cell.cycle"

        ##keyword <- input$enrich_wordcloud_clicked_word ## wordcloud2
        keyword <- input$d3word ## rWordcloud
        sel.row <- input$wordcloud_enrichmentTable_rows_selected
        keyword <- topFreq$word[sel.row]
        
        if( length(keyword)==0 || keyword[1] %in% c(NA,"") ) keyword <- "cell.cycle"
        cat("<enrich_keyword> 3: selected keyword=",keyword,"\n")
        
        targets <- grep(keyword, rownames(S), ignore.case=TRUE, value=TRUE)
        length(targets)
        gmt <- list("set1"=targets)
        names(gmt)[1] = keyword
        
        require(fgsea)
        i=1
        res <- c()
        for(i in 1:ncol(S)) {
            res1 <- fgsea( gmt, S[,i], nperm=400 )[,1:5]
            res <- rbind(res, as.data.frame(res1)[1,])
        }
        rownames(res) <- colnames(S)
        res$padj <- p.adjust( res$pval, method="fdr")
        ##res <- res[order(-res$NES),]
        
        fx <- res$NES
        pv <- res$pval
        names(fx) <- names(pv) <- rownames(res)
        top.up   <- names(sort(fx[which(fx>0)],decreasing=TRUE))
        top.down <- names(sort(fx[which(fx<0)]))
        top.q    <- names(pv)[order(pv,-abs(fx))]

        par(mfrow=c(3,3), mar=c(0.2,3.2,3.2,0.2), mgp=c(1.8,0.7,0))
        i=1
        for(i in 1:9) {
            if(i > length(top.q)) {
                frame()
            } else {
                cmp <- top.q[i]
                gsea.enplot(S[,cmp], targets, names=NULL, ##main=gs,
                            main = paste0("#",toupper(keyword),"\n@",cmp),
                            cex.main=0.9, len.main=80)
                qv1 = formatC(res[cmp,"padj"],format="e", digits=2)
                nes1 = formatC(res[cmp,"NES"],format="f", digits=2)
                tt <- c(paste("NES=",nes1),paste("q=",qv1))
                legend("topright", tt, bty="n",cex=0.85)
            }
        }
        
    })

    wordcloud_enrichmentTable.RENDER <- reactive({    
        df <- enrich_getWordEnrichment()
        req(df)
        df <- df[,c("word","pval","padj","ES","NES","size")]
        cat("<wordcloud_enrichmentTable.RENDER> dim(df)=",dim(df),"\n")
        
        numeric.cols <- which(sapply(df, is.numeric))
        numeric.cols
        tbl <- DT::datatable( df, rownames=FALSE,
                             class = 'compact cell-border stripe hover',                  
                             extensions = c('Scroller'),
                             selection=list(mode='single', target='row', selected=1),
                             options=list(
                                 dom = 'lfrtip', 
                                 scrollX = TRUE, scrollY = tabH,
                                 scroller=TRUE, deferRender=TRUE
                             )  ## end of options.list 
                             ) %>%
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
            DT::formatStyle( "NES",
                            background = color_from_middle( df[,"NES"], 'lightblue', '#f5aeae'),
                            backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center') 
        ##tbl <- DT::datatable(df)
        return(tbl)
    })

    wordcloud_leadingEdgeTable.RENDER <- reactive({    

        ngs <- inputData()
        req(ngs, input$fa_contrast)

        df <- enrich_getWordEnrichment()

        sel.row=1
        sel.row <- input$wordcloud_enrichmentTable_rows_selected
        req(df, sel.row)
        if(is.null(sel.row)) return(NULL)
        
        ee <- unlist(df$leadingEdge[sel.row])
        ee <- strsplit(ee, split="//")[[1]]

        ##fx <- ngs$gset.meta$meta[[1]][ee,"meta.fx"]
        ##fx <- ngs$gset.meta$meta[[1]][ee,"fc"][,"fgsea"]  ## real NES
        fx <- ngs$gset.meta$meta[[input$fa_contrast]][ee,"meta.fx"]
        names(fx) <- ee
        df <- data.frame("leading.edge"=ee, fx=fx)
        df <- df[order(-abs(df$fx)),]
        
        numeric.cols <- which(sapply(df, is.numeric))
        numeric.cols

        df$leading.edge <- wrapHyperLink(df$leading.edge, df$leading.edge)  ## add link
        
        tbl <- DT::datatable( df, rownames=FALSE, escape = c(-1,-2),
                             class = 'compact cell-border stripe hover',                  
                             extensions = c('Scroller'),
                             selection=list(mode='single', target='row', selected=1),
                             options=list(
                                 dom = 'lfrtip', 
                                 scrollX = TRUE, scrollY = tabH,
                                 scroller=TRUE, deferRender=TRUE
                             )  ## end of options.list 
                             ) %>%
            formatSignif(numeric.cols,4) %>%
            DT::formatStyle(0, target='row', fontSize='11px', lineHeight='70%') %>% 
            DT::formatStyle( "fx",
                            background = color_from_middle( df[,"fx"], 'lightblue', '#f5aeae'),
                            backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat',
                            backgroundPosition = 'center') 
        ##tbl <- DT::datatable(df)
        return(tbl)
    })

    wordcloud_actmap.RENDER <- reactive({

        cat("<wordcloud_actmap> called\n")
        
        ##df <- enrich_getWordEnrichment()
        ##req(df)
        res <- enrich_getWordFreqResults()   
        score <- sapply(res$gsea, function(x) x$NES)
        rownames(score) <- res$gsea[[1]]$word
        
        ## reduce score matrix
        ##score = head(score[order(-rowSums(abs(score))),],40)
        ##score = score[head(order(-rowSums(score**2)),50),] ## max number of terms
        ##score = score[head(order(-score[,comparison]**2),50),,drop=FALSE] ## max terms
        score = score[head(order(-rowMeans(score[,]**2)),50),,drop=FALSE] ## max terms    
        score = score[,head(order(-colSums(score**2)),25),drop=FALSE] ## max comparisons/FC

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
        
        cex2=1
        colnames(score) = substring(colnames(score),1,30)
        rownames(score) = substring(rownames(score),1,50)
        if(ncol(score)>15) {
            rownames(score) = substring(rownames(score),1,40)
            cex2=0.85
        }
        if(ncol(score)>25) {
            rownames(score) = substring(rownames(score),1,30)
            colnames(score) <- rep("",ncol(score))
            cex2=0.7
        }

        par(mfrow=c(1,1), mar=c(1,1,1,1), oma=c(0,2,0,1))
        require(corrplot)
        score2 <- score
        score2 <- t( t(score2) / apply(abs(score2),2,max)) ## normalize rows???
        score2 <- sign(score2) * abs(score2/max(abs(score2)))**3   ## fudging
        bmar <- 0 + pmax((50 - nrow(score2))*0.25,0)
        corrplot( score2, is.corr=FALSE, cl.pos = "n", col=BLUERED(100),
                 tl.cex=cex2, tl.col="grey20", mar=c(bmar,0,0,0) )
        
    })    

    ##---------------------------------------------------------------
    ##------------- modules for WordCloud ---------------------------
    ##---------------------------------------------------------------
    require(wordcloud2)
    require(rWordCloud)

    enrich_wordtsne_module <- plotModule(
        id="enrich_wordtsne", label="c",
        ##plotlib="ggplot", func=enrich_wordtsne.RENDER,
        plotlib="plotly", func=enrich_wordtsne.PLOTLY, 
        info.text = "<strong>Word t-SNE.</strong> T-SNE of keywords that were found in the title/description of gene sets. Keywords that are often found together in title/descriptions are placed close together in the t-SNE. For each keyword we computed enrichment using GSEA on the mean (absolute) enrichment profiles (averaged over all contrasts). Statistically significant gene sets (q<0.05) are colored in red. The sizes of the nodes are proportional to the normalized enrichment score (NES) of the keyword.",
        options = tagList(
            tipify(radioButtons(ns("enrich_wordtsne_algo"),"Clustering algorithm:",
                                choices=c("tsne","umap"),inline=TRUE),
                   "Choose a clustering algorithm: t-SNE or UMAP.")
        ),
        pdf.width=8, pdf.height=8, res=72, pdf.pointsize=13,
        ##datacsv = enrich_getWordFreq,
        title = "Word t-SNE"
    )
    output <- attachModule(output, enrich_wordtsne_module)

    enrich_wordcloud_opts = tagList(
        tipify(selectInput(ns("enrich_wordcloud_exclude"),"Exclude words:", choices=NULL, multiple=TRUE),
               "Paste a keyword to exclude it from the plot.", placement="top", options = list(container = "body")),
        tipify(selectInput(ns("enrich_wordcloud_colors"),"Colors:", choices=c("Blues","Greys","Accent","Dark2"),
                           multiple=FALSE),
               "Choose a set of colors.", placement="top", options = list(container = "body"))
    )

    enrich_wordcloud_module <- plotModule(
        id="enrich_wordcloud", label="b",
        func = enrich_wordcloud.RENDER,
        plotlib="base", renderFunc="renderPlot", outputFunc="plotOutput",
        ##plotlib="htmlwidget", renderFunc="renderWordcloud2", outputFunc="wordcloud2Output",
        ##plotlib="htmlwidget", renderFunc="renderd3Cloud", outputFunc="d3CloudOutput",
        ##download.fmt = NULL,
        info.text = "<strong>Word cloud.</strong> Word cloud of the most enriched keywords for the data set. Select a keyword in the 'Enrichment table'. In the plot settings, users can exclude certain words from the figure, or choose the color palette. The sizes of the words are relative to the normalized enrichment score (NES) from the GSEA computation. Keyword enrichment is computed by running GSEA on the mean (squared) enrichment profile (averaged over all contrasts). For each keyword, we defined the 'keyword set' as the collection of genesets that contain that keyword in the title/description.",
        options = enrich_wordcloud_opts,
        pdf.width=6, pdf.height=6, res=72,
        title = "Word cloud"
    )
    output <- attachModule(output, enrich_wordcloud_module)

    enrich_keyword_info = "<strong>Keyword enrichment analysis.</strong> Computes enrichment of a selected keyword across all contrasts. Select a keyword by clicking a word in the 'Enrichment table'.

<br><br>Keyword enrichment is computed by running GSEA on the enrichment score profile for all contrasts. We defined the test set as the collection of genesets that contain the keyword in the title/description. Black vertical bars indicate the position of gene sets that contains the *keyword* in the ranked list of enrichment scores. The curve in green corresponds to the 'running statistic' of the keyword enrichment score. The more the green ES curve is shifted to the upper left of the graph, the more the keyword is enriched in the first group. Conversely, a shift of the green ES curve to the lower right, corresponds to keyword enrichment in the second group."

    ##myTextInput('enrich_keyword_keywords','Keyword:',"cell cycle"),

    enrich_keyword_opts = tagList(
        tipify( textInput(ns('enrich_keyword_keywords'),'Keyword:',"cell cycle"),
               "Paste a keyword such as 'apoptosis', 'replication' or 'cell cycle'.",
               placement="top", options = list(container = "body"))
    )

    ## enrich_keyword_opts = textInput('enrich_keyword_keywords','Keyword:',"cell cycle")

    enrich_keyword_module <- plotModule(
        id="enrich_keyword", label="a",
        plotlib="base", func=enrich_keyword.RENDER,
        info.text = enrich_keyword_info,
        ## options = enrich_keyword_opts,
        pdf.width=6, pdf.height=6, res=90,
        title = "Enrichment plots"
    )
    output <- attachModule(output, enrich_keyword_module)

    ##--------buttons for enrichment table
    wordcloud_enrichmentTable_module <- tableModule(
        id = "wordcloud_enrichmentTable", label="e",
        func = wordcloud_enrichmentTable.RENDER,
        info.text="Keyword enrichment table.", 
        title = "Enrichment table"
    )
    output <- attachModule(output, wordcloud_enrichmentTable_module)

    ##--------buttons for leading edge table
    wordcloud_leadingEdgeTable_module <- tableModule(
        id = "wordcloud_leadingEdgeTable", label="f",
        func = wordcloud_leadingEdgeTable.RENDER,
        info.text="Keyword leading edge table.", 
        title = "Leading-edge table"
    )
    output <- attachModule(output, wordcloud_leadingEdgeTable_module)

    ##-------- Activation map plotting module
    wordcloud_actmap.opts = tagList()
    wordcloud_actmap_module <- plotModule(
        id="wordcloud_actmap", wordcloud_actmap.RENDER,
        title = "Activation matrix", label="d",
        info.text = "The <strong>Activation Matrix</strong> visualizes the activation of drug activation enrichment across the conditions. The size of the circles correspond to their relative activation, and are colored according to their upregulation (red) or downregulation (blue) in the contrast profile.",
        options = wordcloud_actmap.opts,
        pdf.width=6, pdf.height=10, res=72
    )
    output <- attachModule(output, wordcloud_actmap_module)

    ##---------------------------------------------------------------
    ##-------------- UI Layout for WordCloud ------------------------
    ##---------------------------------------------------------------

    wordcloud_caption = "<b>Keyword enrichment analysis.</b> <b>(a)</b> Enrichment plots for the top most significant contrasts. Black vertical bars indicate the position of gene sets, in the ranked enrichment scores, that contains the *keyword*. The green curve corresponds to 'running statistics' of the keyword enrichment score. <b>(b)</b> Colored word cloud. The size of the words are relative to the normalized enrichment score (NES) from the GSEA computation. <b>(c)</b> T-SNE plot of keywords extracted from the titles/descriptions of the genesets. <b>(d)</b> Activation matrix showing keyword enrichment across contrasts. <b>(e)</b> Keyword enrichment table for selected contrast. <b>(f)</b> Leading edge terms for selected keyword enrichment."

    output$wordcloud_UI <- renderUI({
        fillCol(
            height = rowH,
            flex = c(1,NA),
            fillRow(
                height = 650,
                flex = c(3.6,1),
                fillCol(
                    flex=c(1.3,0.15,1),
                    height = 650,
                    fillRow(
                        flex = c(1.2,0.05,1,0.05,1),
                        moduleWidget(enrich_keyword_module, ns=ns),
                        br(),
                        moduleWidget(enrich_wordcloud_module, ns=ns),
                        br(),
                        ##moduleWidget(enrich_wordtsne_module), 
                        moduleWidget(enrich_wordtsne_module, outputFunc="plotlyOutput", ns=ns)
                    ),
                    br(),
                    fillRow(
                        flex=c(1,0.08,1),
                        moduleWidget(wordcloud_enrichmentTable_module, outputFunc="dataTableOutput", ns=ns),
                        br(),
                        moduleWidget(wordcloud_leadingEdgeTable_module, outputFunc="dataTableOutput", ns=ns)
                    )
                ),
                moduleWidget(wordcloud_actmap_module, ns=ns)
            ),
            div(HTML(wordcloud_caption),class="caption")
        )
    })

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
        plot( xx0[ii], yy0[ii] + dy, col=klr[ii],
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

    output$fa_fireplot <- renderPlot({
        require(igraph)
        ngs <- inputData()
        req(ngs)
        ##df <- getFilteredKeggTable()
        
        cmp=1
        cmp <- input$fa_contrast

        sel="B-cell related"
        sel="KEGG metabolic pathways"
        sel="Hallmark collection"
        sel=5
        sel <- input$fire_xpcollection    
        if(is.null(sel)) return(NULL)

        gsets <- GSETS[ COLLECTIONS[[sel]] ]
        ##gsets0 <- gsets0[ intersect( names(gsets0), names(zx)) ]
        shownames <- input$fire_shownames
        pgx.firePlot(ngs, cmp, gsets, shownames=shownames)     
        
    }, res=85)

    output$fireplot_UI <- renderUI({
        fillRow(
            height = rowH,
            flex = c(2,1),
            fillCol(
                flex = c(NA,1), 
                inputPanel(
                    selectInput(ns("fire_xpcollection"),NULL,
                                choices=setdiff(names(COLLECTIONS),"<all>") ),
                    checkboxInput(ns('fire_shownames'),'shownames',TRUE),            
                    cellArgs=list(width='100%')
                ),
                plotOutput('fa_fireplot')
            ),
            br()
        )
    })



}
