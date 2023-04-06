##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## =====================================================================================
## ============================= GADGET UI =============================================
## =====================================================================================


MakeContrastGadget <- function(X, pheno, height = 720) {
  gadgetize(MakeContrastUI, MakeContrastServer,
    title = "MakeContrastGadget",
    pheno = pheno, height = height
  )
}

MakeContrastUI <- function(id) {
  ns <- shiny::NS(id)

  shiny::fillCol(
    height = 750,
    flex = c(1, NA, NA, 1),
    shiny::fillRow(
      flex = c(3, 0.06, 1.0),
      shiny::fillCol(
        flex = c(NA, NA, 1.0),
        shiny::h4("Create comparisons"),
        shiny::fillRow(
          flex = c(1, 4),
          shiny::fillCol(
            flex = c(NA, NA, NA, NA, 1),
            tipifyL(
              shiny::selectInput(ns("param"),
                "Phenotype:",
                choices = NULL,
                selected = NULL,
                multiple = TRUE
              ),
              "Select phenotype(s) to create conditions for your groups. Select &ltgene&gt if you want to split by high/low expression of some gene. Select &ltsamples&gt if you want to group manually on sample names. You can select multiple phenotypes to create combinations."
            ),
            shiny::conditionalPanel(
              "input.param == '<gene>'",
              ns = ns,
              shiny::selectizeInput(ns("gene"),
                "Gene:",
                choices = NULL,
                multiple = FALSE
              )
            ),
            shiny::br(),
            tipifyL(
              shiny::textInput(ns("newname"),
                "Comparison name:",
                placeholder = "e.g. MAIN_vs_CONTROL"
              ),
              "Give a name for your contrast as MAIN_vs_CONTROL, with the name of the main group first. You must keep _vs_ in the name to separate the names of the two groups."
            ),
            shiny::br(),
            shiny::actionButton(ns("addcontrast"),
              "add comparison",
              icon = icon("plus"),
              class = "btn-outline-primary"
            ),
            shiny::br()
          ),
          withTooltip(
            shiny::uiOutput(ns("createcomparison"),
              style = "font-size:13px; height: 280px; overflow-y: scroll;"
            ),
            "Create comparisons by dragging conditions into the main or control groups on the right. Then press add comparison to add the contrast to the table.",
            placement = "top", options = list(container = "body")
          )
        )
      ),
      shiny::br(),
      plotWidget(ns("pcaplot"))
    ),
    shiny::h4("Contrast table"),
    shiny::fillRow(
      height = 24,
      flex = c(NA, 0.05, NA, NA, 1),
      withTooltip(
        shiny::actionButton(ns("autocontrast"),
          "add auto-contrasts",
          icon = icon("plus"),
          class = "small-button btn-outline-primary"
        ),
        "If you are feeling lucky, try this to automatically create contrasts.",
        placement = "top", options = list(container = "body")
      ),
      shiny::br(),
      shiny::div(shiny::HTML("<b>Strata:</b>"), style = "padding: 4px 4px;"),
      shiny::selectInput(ns("strata"), NULL, choices = NULL, width = "120px"),
      shiny::br()
    ),
    shiny::div(DT::dataTableOutput(ns("contrastTable")),
      style = "font-size:13px; height: 300px; margin-top: 20px;overflow-y: scroll;"
    )
  )
}

MakeContrastServerRT <- function(id, phenoRT, contrRT, countsRT, height = 720) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      rv <- shiny::reactiveValues(contr = NULL, pheno = NULL)

      shiny::observe({
        rv$contr <- contrRT()
      })

      shiny::observe({
        rv$pheno <- phenoRT()
      })

      shiny::observe({
        shiny::req(phenoRT())
        px <- colnames(phenoRT())
        shiny::updateSelectInput(session, "pcaplot.colvar", choices = px)
        shiny::updateSelectInput(session, "strata", choices = c("<none>", px))
      })

      observeEvent(countsRT(), {
        genes <- sort(rownames(countsRT()))
        updateSelectizeInput(inputId = "gene", choices = genes)

        phenotypes <- c(sort(unique(colnames(phenoRT()))), "<samples>", "<gene>")
        phenotypes <- grep("_vs_", phenotypes, value = TRUE, invert = TRUE) ## no comparisons...
        psel <- c(grep("sample|patient|name|id|^[.]",
          phenotypes,
          value = TRUE, invert = TRUE
        ), phenotypes)[1]
        updateSelectInput(inputId = "param", choices = phenotypes, selected = psel)
      })

      sel.conditions <- shiny::reactive({
        shiny::req(phenoRT(), countsRT())
        df <- phenoRT()

        if ("<samples>" %in% input$param) {
          df$"<samples>" <- rownames(df)
        }
        if ("<gene>" %in% input$param) {
          gene <- input$gene
          if (gene %in% rownames(countsRT())) {
            gx <- log2(1 + countsRT()[gene, ])
            ## df$"<gene>" <- c("low","high")[1 + 1*(gx >= mean(gx,na.rm=TRUE))]
            df$"<gene>" <- gx
          } else {
            return(NULL)
          }
        }

        df <- type.convert(df)
        ii <- which(sapply(type.convert(df), class) %in% c("numeric", "integer"))
        if (length(ii)) {
          for (i in ii) {
            x <- df[, i]
            df[, i] <- c("low", "high")[1 + 1 * (x >= mean(x, na.rm = TRUE))]
          }
        }

        pp <- intersect(input$param, colnames(df))
        ss <- colnames(countsRT())
        cond <- apply(df[ss, pp, drop = FALSE], 1, paste, collapse = "_")
        cond <- gsub("^_|_$", "", cond)
        cond
      })

      output$createcomparison <- shiny::renderUI({
        shiny::req(input$param)
        cond <- sel.conditions()
        if (length(cond) == 0 || is.null(cond)) {
          return(NULL)
        }

        items <- c("<others>", sort(unique(cond)))

        shiny::tagList(
          shiny::tags$head(shiny::tags$style(".default-sortable .rank-list-item {padding: 2px 15px;}")),
          sortable::bucket_list(
            header = NULL,
            sortable::add_rank_list(
              text = "Conditions:",
              labels = items
            ),
            sortable::add_rank_list(
              input_id = ns("group1"),
              text = "Main group:"
            ),
            sortable::add_rank_list(
              input_id = ns("group2"),
              text = "Control group:"
            ),
            group_name = "cmpbucket"
          )
        )
      })

      buttonInput <- function(FUN, len, id, ...) {
        inputs <- character(len)
        for (i in seq_len(len)) {
          inputs[i] <- as.character(FUN(paste0(id, i), ...))
        }
        inputs
      }

      shiny::observeEvent(c(input$group1, input$group2), {
        g1 <- gsub("[-_.,<> ]", ".", input$group1)
        g2 <- gsub("[-_.,<> ]", ".", input$group2)
        g1 <- gsub("[.]+", ".", g1)
        g2 <- gsub("[.]+", ".", g2)
        g1 <- paste(g1, collapse = "")
        g2 <- paste(g2, collapse = "")
        if (is.null(g1) || length(g1) == 0) g1 <- ""
        if (is.null(g2) || length(g2) == 0) g2 <- ""
        if (is.na(g1)) g1 <- ""
        if (is.na(g2)) g2 <- ""
        g1 <- substring(g1, 1, 20)
        g2 <- substring(g2, 1, 20)
        prm.name <- paste(input$param, collapse = ".")
        prm.name <- gsub("[-_.,<> ]", "", prm.name)
        if (any(input$param %in% "<gene>")) {
          prm.name <- sub("gene", input$gene, prm.name)
        }
        tt <- paste0(prm.name, ":", g1, "_vs_", g2)
        if (g1 == "" && g2 == "") tt <- ""
        shiny::updateTextInput(session, "newname", value = tt)
      })

      shiny::observeEvent(input$contrast_delete, {
        ## Observe if a contrast is to be deleted
        ##
        id <- as.numeric(gsub(".*_", "", input$contrast_delete))
        if (length(id) == 0) {
          return(NULL)
        }
        if (!is.null(rv$contr) && NCOL(rv$contr) <= 1) {
          rv$contr <- rv$contr[, 0, drop = FALSE]
        } else {
          rv$contr <- rv$contr[, -id, drop = FALSE]
        }
      })



      shiny::observeEvent(input$addcontrast, {

        cond <- sel.conditions()
        if (length(cond) == 0 || is.null(cond)) {
          return(NULL)
        }

        group1 <- input$group1
        group2 <- input$group2
        in.main <- 1 * (cond %in% group1)
        in.ref1 <- 1 * (cond %in% group2)
        in.ref2 <- ("<others>" %in% group2) & (!cond %in% group1)
        in.ref <- in.ref1 | in.ref2

        ## ctx <- 1*(in.main) - 1*(in.ref)
        ## ct.name <- paste0(input$group1name,"_vs_",input$group2name)
        ct.name <- input$newname
        gr1 <- gsub(".*:|_vs_.*", "", ct.name) ## first is MAIN group!!!
        gr2 <- gsub(".*_vs_|@.*", "", ct.name)
        ctx <- c(NA, gr1, gr2)[1 + 1 * in.main + 2 * in.ref]

        if (sum(in.main) == 0 || sum(in.ref) == 0) {
          shinyalert::shinyalert("ERROR", "Both groups must have samples")
          return(NULL)
        }
        if (ct.name %in% c(NA, "", " ")) {
          shinyalert::shinyalert("ERROR", "You must give a contrast name")
          return(NULL)
        }
        if (1 && gr1 == gr2) {
          shinyalert::shinyalert("ERROR", "Invalid contrast name")
          return(NULL)
        }
        if (!is.null(rv$contr) && ct.name %in% colnames(rv$contr)) {
          shinyalert::shinyalert("ERROR", "Contrast name already exists.")
          return(NULL)
        }
        if (!grepl("_vs_", ct.name)) {
          shinyalert::shinyalert("ERROR", "Contrast must include _vs_ in name")
          return(NULL)
        }

        ## update reactive value
        samples <- colnames(countsRT())

        ctx1 <- matrix(ctx, ncol = 1, dimnames = list(samples, ct.name))
        if (is.null(rv$contr)) {
          rv$contr <- ctx1
        } else {
          rv$contr <- cbind(rv$contr, ctx1)
        }

        ## if(any(input$param %in% c('<gene>','<samples>'))) {
        if (any(input$param %in% c("<gene>"))) {
          if (is.null(rv$pheno) || NCOL(rv$pheno) == 0) {
            rv$pheno <- ctx1
          } else {
            message("[MakeContrastServer:addcontrast] add to cond : dim(ctx1) = ", dim(ctx1))
            if (!ct.name %in% colnames(rv$pheno)) {
              rv$pheno <- cbind(rv$pheno, ctx1)
            }
          }
        }
      })

      shiny::observeEvent(input$autocontrast, {
        shiny::req(phenoRT())
        df <- phenoRT()
        strata.var <- input$strata

        ctx <- NULL
        if (strata.var == "<none>") {
          ct <- playbase::pgx.makeAutoContrasts(
            df,
            mingrp = 3, slen = 20, ref = NULL, fix.degenerate = FALSE
          )
          if (!is.null(ct)) {
            ctx <- ct$exp.matrix
            rownames(ctx) <- rownames(df)
          }
        } else {
          ctx <- playbase::pgx.makeAutoContrastsStratified(
            df,
            strata = strata.var,
            mingrp = 3, slen = 20, ref = NULL, fix.degenerate = FALSE
          )
        }
        if (is.null(ctx)) {
          return(NULL)
        }

        ## update reactive value
        ctx2 <- playbase::contrastAsLabels(ctx)
        if (!is.null(rv$contr)) {
          rv$contr <- cbind(rv$contr, ctx2)
        } else {
          rv$contr <- ctx2
        }
      })

      output$contrastTable <- DT::renderDataTable(
        {
          ct <- rv$contr

          if (is.null(ct) || NCOL(ct) == 0) {
            df <- data.frame(
              delete = 0,
              comparison = "",
              n1 = 0,
              n0 = 0,
              "main.group" = "",
              "control.group" = ""
            )[0, ]
          } else {
            paste.max <- function(x, n = 6) {
              ## x <- unlist(x)
              if (length(x) > n) {
                x <- c(x[1:n], paste("+", length(x) - n, "others"))
              }
              paste(x, collapse = " ")
            }

            ct1 <- playbase::makeContrastsFromLabelMatrix(ct)
            ct1[is.na(ct1)] <- 0

            if (NCOL(ct) == 1) {
              ss1 <- names(which(ct1[, 1] > 0))
              ss2 <- names(which(ct1[, 1] < 0))
              ss1 <- paste.max(ss1, 6)
              ss2 <- paste.max(ss2, 6)
            } else {
              ss0 <- rownames(ct)
              ss1 <- apply(ct1, 2, function(x) paste.max(ss0[which(x > 0)]))
              ss2 <- apply(ct1, 2, function(x) paste.max(ss0[which(x < 0)]))
            }

            deleteButtons <- buttonInput(
              FUN = actionButton,
              len = ncol(ct),
              ## id = 'contrast_delete_',
              id = paste0("contrast_delete_", sample(99999, 1), "_"), ## hack to allow double click
              label = "",
              ## size = "mini",
              width = "50px",
              inline = TRUE,
              icon = shiny::icon("trash-alt"),
              class = "btn-inline btn-outline-danger-hover",
              style = "padding:2px; margin:2px; font-size:95%; color: #B22222;",
              ## onclick = 'Shiny.onInputChange(\"contrast_delete\",this.id)'
              onclick = paste0('Shiny.onInputChange(\"', ns("contrast_delete"), '\",this.id)')
            )

            df <- data.frame(
              delete = deleteButtons,
              comparison = colnames(ct1),
              n1 = colSums(ct1 > 0),
              n0 = colSums(ct1 < 0),
              "main.group" = ss1,
              "control.group" = ss2
            )
          }
          rownames(df) <- NULL

          DT::datatable(
            df,
            rownames = FALSE,
            escape = c(-1),
            selection = "none",
            class = "compact cell-border",
            options = list(
              dom = "t",
              pageLength = 999,
              ## autoWidth = TRUE, ## scrollX=TRUE,
              columnDefs = list(
                list(width = "20px", targets = c(0, 2, 3)),
                list(width = "150px", targets = c(1)),
                list(width = "400px", targets = c(4, 5))
              )
            )
          ) %>%
            DT::formatStyle(0, target = "row", fontSize = "12px", lineHeight = "99%")
        },
        server = FALSE
      )


      pcaplot.RENDER <- shiny::reactive({
        ## X <- ngs$X
        pheno <- phenoRT()
        counts <- countsRT()
        if (is.null(pheno) || is.null(counts)) {
          return(NULL)
        }
        if (NCOL(pheno) == 0 || NCOL(counts) == 0) {
          return(NULL)
        }
        shiny::req(pheno)
        shiny::req(counts)

        method <- input$pcaplot.method
        X <- log2(1 + counts)
        clust <- playbase::pgx.clusterMatrix(X, dims = 2, method = method)
        names(clust)

        cond <- sel.conditions()
        if (length(cond) == 0 || is.null(cond)) {
          return(NULL)
        }
        ## par(mar=c(4,1,1,1))
        playbase::pgx.scatterPlotXY(
          clust$pos2d,
          var = cond, plotlib = "plotly",
          legend = FALSE ## , labels=TRUE
        )
      })

      pcaplot.opts <- shiny::tagList(
        withTooltip(
          shiny::selectInput(ns("pcaplot.method"), "Method:",
            choices = c("pca", "tsne", "umap"),
            width = "100%"
          ), "Choose clustering method.",
          placement = "right", options = list(container = "body")
        )
      )

      shiny::callModule(
        plotModule,
        id = "pcaplot",
        func = pcaplot.RENDER, ## ns=ns,
        plotlib = "plotly",
        options = pcaplot.opts,
        height = c(320, 700), width = c("auto", 800),
        pdf.width = 8, pdf.height = 8,
        title = "PCA/tSNE plot",
        ## info.text = hm_PCAplot_text
        ## caption = pca_caption_static
      )

      return(
        shiny::reactive({
          if (is.null(rv$contr)) {
            return(NULL)
          }
          rv
        })
      ) ## pointing to reactive
    } ## end-of-server
  )
}
