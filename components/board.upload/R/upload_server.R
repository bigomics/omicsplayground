##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2022 BigOmics Analytics Sagl. All rights reserved.
##

UploadBoard <- function(id,
                        pgx_dir,
                        pgx,
                        auth,
                        limits = c(
                          "samples" = 1000, "comparisons" = 20,
                          "genes" = 20000, "genesets" = 10000,
                          "datasets" = 10
                        ),
                        enable_upload = TRUE,
                        enable_save = TRUE,
                        enable_userdir = TRUE) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    loadedDataset <- shiny::reactiveVal(0) ## counts/trigger dataset upload

    phenoRT <- shiny::reactive(uploaded$samples.csv)
    contrRT <- shiny::reactive(uploaded$contrasts.csv)

    rv <- shiny::reactiveValues(contr = NULL, pheno = NULL)

    shiny::observe({
      rv$contr <- contrRT()
    })

    shiny::observe({
      rv$pheno <- phenoRT()
    })

    observe({
      phenotypes <- c(sort(unique(colnames(phenoRT()))), "<samples>", "<gene>")
      phenotypes <- grep("_vs_", phenotypes, value = TRUE, invert = TRUE) ## no comparisons...
      psel <- c(grep("sample|patient|name|id|^[.]", phenotypes,
        value = TRUE,
        invert = TRUE
      ), phenotypes)[1]
      shiny::updateSelectInput(
        session = session,
        inputId = "param",
        choices = phenotypes,
        selected = psel
      )
      genes <- sort(rownames(corrected_counts()))
      shiny::updateSelectInput(
        session = session,
        inputId = "gene",
        choices = genes,
        selected = genes[1]
      )
    })

    output$navheader <- shiny::renderUI({
      fillRow(
        flex = c(NA, 1, NA),
        ## h2(input$nav),
        shiny::div(
          id = "navheader-current-section",
          HTML("Upload data &nbsp;"),
          shiny::actionLink(
            ns("module_info"), "",
            icon = shiny::icon("info-circle"),
            style = "color: #ccc;"
          )
        ),
        shiny::br(),
        shiny::div(pgx$name, id = "navheader-current-dataset")
      )
    })

    shiny::observeEvent(input$module_info, {
      shiny::showModal(shiny::modalDialog(
        title = shiny::HTML("<strong>Upload data</strong>"),
        shiny::HTML(module_infotext),
        easyClose = TRUE, size = "l"
      ))
    })

    module_infotext <- paste0(
      'Under the <b>Upload data</b> panel users can upload their transcriptomics and proteomics data to the platform. The platform requires 3 data files as listed below: a data file containing counts/expression (counts.csv), a sample information file (samples.csv) and a file specifying the statistical comparisons as contrasts (contrasts.csv). It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, row names and column names match for all files. On the left side of the panel, users need to provide a unique name and brief description for the dataset while uploading. N.B. Users can now create contrasts from the platform itself, so the contrasts.csv file is optional.

<br><br>
<ol>
<li>counts.csv: Count/expression file with gene on rows, samples as columns.
<li>samples.csv: Samples file with samples on rows, phenotypes as columns.
<li>contrasts.csv: Contrast file with conditions on rows, contrasts as columns.
</ol>

<br><br><br>
<center><iframe width="560" height="315" src="https://www.youtube.com/embed/elwT6ztt3Fo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe><center>

'
    )

    ## ================================================================================
    ## ====================== NEW DATA UPLOAD =========================================
    ## ================================================================================
    ## reload_pgxdir()

    getPGXDIR <- shiny::reactive({
      ## reload_pgxdir()  ## force reload

      email <- "../me@company.com"
      email <- auth$email()
      email <- gsub(".*\\/", "", email)
      pdir <- pgx_dir ## from module input

      ## USERDIR=FALSE
      if (enable_userdir) {
        pdir <- paste0(pdir, "/", email)
        if (!is.null(email) && !is.na(email) && email != "") pdir <- paste0(pdir, "/")
        if (!dir.exists(pdir)) {
          dbg("[LoadingBoard:getPGXDIR] userdir does not exists. creating pdir = ", pdir)
          dir.create(pdir)
          dbg("[LoadingBoard:getPGXDIR] copy example pgx")
          file.copy(file.path(pgx_dir, "example-data.pgx"), pdir)
        }
      }
      pdir
    })

    if (enable_upload) {
      uploaded_pgx <- UploadModuleServer(
        id = "upload_panel",
        FILES = FILES,
        pgx.dirRT = shiny::reactive(getPGXDIR()),
        height = 720,
        ## limits = c(samples=20, comparisons=20, genes=8000),
        limits = limits
      )

      shiny::observeEvent(uploaded_pgx(), {
        dbg("[observe::uploaded_pgx] uploaded PGX detected!")
        new_pgx <- uploaded_pgx()

        dbg("[observe::uploaded_pgx] initializing PGX object")
        new_pgx <- pgx.initialize(new_pgx)

        ## update Session PGX
        dbg("[UploadBoard@load_react] **** copying current pgx to session.pgx  ****")
        for (i in 1:length(new_pgx)) {
          pgx[[names(new_pgx)[i]]] <- new_pgx[[i]]
        }

        DT::selectRows(proxy = DT::dataTableProxy(ns("pgxtable")), selected = NULL)

        savedata_button <- NULL
        if (enable_save) {
          ## -------------- save PGX file/object ---------------
          pgxname <- sub("[.]pgx$", "", new_pgx$name)
          pgxname <- gsub("^[./-]*", "", pgxname) ## prevent going to parent folder
          pgxname <- paste0(gsub("[ \\/]", "_", pgxname), ".pgx")
          pgxname

          pgxdir <- getPGXDIR()
          fn <- file.path(pgxdir, pgxname)
          fn <- iconv(fn, from = "", to = "ASCII//TRANSLIT")

          ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ## Note: Currently we use 'ngs' as object name but want to go
          ## towards 'pgx' as standard name. Actually saving as RDS
          ## should be better.
          ngs <- new_pgx
          save(ngs, file = fn)

          remove(ngs)
          remove(new_pgx)


          message("[UploadBoard::@savedata] updating PGXINFO")
          pgx.initDatasetFolder(pgxdir, force = FALSE, verbose = TRUE)
          ## reload_pgxdir(reload_pgxdir()+1)
        }

        ## shiny::removeModal()
        msg1 <- "<b>Ready!</b>"
        ## beepr::beep(sample(c(3,4,5,6,8),1))  ## music!!
        beepr::beep(10) ## short beep

        if (enable_save) {
          msg1 <- "<b>Ready!</b><br>Your data is ready and has been saved in your library. You can now start exploring your data."
        } else {
          msg1 <- "<b>Ready!</b><br>Your data is ready. You can now start exploring your data."
        }
        loadedDataset(loadedDataset() + 1) ## notify new data uploaded

        showModal(
          modalDialog(
            HTML(msg1),
            title = NULL,
            size = "s",
            footer = tagList(
              modalButton("Start!")
            )
          )
        )

        shinyjs::runjs("$('.tab-sidebar:eq(1)').trigger('click');")

      })
    }

    # Some 'global' reactive variables used in this file
    uploaded <- shiny::reactiveValues()

    ## Hide/show tabpanels upon available data like a wizard dialog
    shiny::observe({
      has.upload <- Vectorize(function(f) {
        (f %in% names(uploaded) && !is.null(nrow(uploaded[[f]])))
      })
      need2 <- c("counts.csv", "samples.csv")
      need3 <- c("counts.csv", "samples.csv", "contrasts.csv")
      if (all(has.upload(need3))) {
        shiny::showTab("tabs", "Contrasts")
        shiny::showTab("tabs", "Compute")
        if (input$advanced_mode) {
          shiny::showTab("tabs", "Normalize")
          shiny::showTab("tabs", "BatchCorrect")
        }
      } else if (all(has.upload(need2))) {
        if (input$advanced_mode) {
          shiny::showTab("tabs", "Normalize")
          shiny::showTab("tabs", "BatchCorrect")
        }
        shiny::showTab("tabs", "Contrasts")
        shiny::hideTab("tabs", "Compute")
      } else {
        shiny::hideTab("tabs", "Normalize")
        shiny::hideTab("tabs", "BatchCorrect")
        shiny::hideTab("tabs", "Contrasts")
        shiny::hideTab("tabs", "Compute")
      }
    })

    ## =====================================================================
    ## ======================= UI OBSERVERS ================================
    ## =====================================================================

    shiny::observeEvent(input$advanced_mode, {
      if (input$advanced_mode) {
        shiny::showTab("tabs", "Normalize") ## NOT YET!!!
        shiny::showTab("tabs", "BatchCorrect")
      } else {
        shiny::hideTab("tabs", "Normalize")
        shiny::hideTab("tabs", "BatchCorrect")
      }
    })

    ## ========================================================================
    ## ================================= UI ===================================
    ## ========================================================================

    # leaving this for now because not sure about the "suspendWhenHidden" thing (-NC)
    output$contrasts_UI <- shiny::renderUI({
      shiny::fillCol(
        height = height, ## width = 1200,
        MakeContrastUI(ns("makecontrast"))
      )
    })
    shiny::outputOptions(output, "contrasts_UI", suspendWhenHidden = FALSE) ## important!!!

    ## =====================================================================
    ## ================== DATA LOADING OBSERVERS ===========================
    ## =====================================================================

    ## ------------------------------------------------------------------
    ## Observer for uploading data files using fileInput widget.
    ##
    ## Reads in the data files from the file names, checks and
    ## puts in the reactive values object 'uploaded'. Then
    ## uploaded should trigger the computePGX module.
    ## ------------------------------------------------------------------
    shiny::observeEvent(input$upload_files, {
      message("[upload_files] >>> reading uploaded files")
      message("[upload_files] upload_files$name=", input$upload_files$name)
      message("[upload_files] upload_files$datapath=", input$upload_files$datapath)

      ## for(i in 1:length(uploaded)) uploaded[[i]] <- NULL
      uploaded[["pgx"]] <- NULL
      uploaded[["last_uploaded"]] <- NULL

      ## read uploaded files
      pgx.uploaded <- any(grepl("[.]pgx$", input$upload_files$name))
      matlist <- list()

      if (pgx.uploaded) {
        message("[upload_files] PGX upload detected")

        ## If the user uploaded a PGX file, we extract the matrix
        ## dimensions from the given PGX/NGS object. Really?
        ##
        i <- grep("[.]pgx$", input$upload_files$name)
        
        # load(input$upload_files$datapath[i], ngs <- new.env(), verbose = TRUE) ## load NGS/PGX

        ngs <- local(get(load(input$upload_files$datapath[i], verbose = 0))) ## override any name

  
        ## matlist[["counts.csv"]] <- ngs$counts
        ## matlist[["samples.csv"]] <- type.convert(ngs$samples)
        ## matlist[["contrasts.csv"]] <- ngs$model.parameters$exp.matrix
        uploaded[["pgx"]] <- ngs
      } else {
        ## If the user uploaded CSV files, we read in the data
        ## from the files.
        ##
        message("[upload_files] getting matrices from CSV")

        ii <- grep("csv$", input$upload_files$name)
        ii <- grep("sample|count|contrast|expression",
          input$upload_files$name,
          ignore.case = TRUE
        )
        if (length(ii) == 0) {
          return(NULL)
        }

        inputnames <- input$upload_files$name[ii]
        uploadnames <- input$upload_files$datapath[ii]

        if (length(uploadnames) > 0) {
          i <- 1
          for (i in 1:length(uploadnames)) {
            fn1 <- inputnames[i]
            fn2 <- uploadnames[i]
            matname <- NULL
            df <- NULL
            if (grepl("count", fn1, ignore.case = TRUE)) {
              dbg("[upload_files] counts.csv : fn1 = ", fn1)
              ## allows duplicated rownames
              df0 <- read.as_matrix(fn2)
              if (TRUE && any(duplicated(rownames(df0)))) {
                ndup <- sum(duplicated(rownames(df0)))
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Duplicated gene names",
                  text = paste("Your counts matrix has", ndup, "duplicated gene names.\nCounts of those genes will be merged."),
                  type = "warning",
                  btn_labels = "OK",
                  closeOnClickOutside = FALSE,
                )
              }
              dbg(
                "[upload_files] counts.csv : 1 : dim(df0) = ",
                paste(dim(df0), collapse = "x")
              )

              if (nrow(df0) > 1 && NCOL(df0) > 1) {
                df <- as.matrix(df0)
                matname <- "counts.csv"
              }
            } else if (grepl("expression", fn1, ignore.case = TRUE)) {
              dbg("[upload_files] expression.csv : fn1 = ", fn1)
              ## allows duplicated rownames
              df0 <- read.as_matrix(fn2)
              if (TRUE && any(duplicated(rownames(df0)))) {
                ndup <- sum(duplicated(rownames(df0)))
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Duplicated gene names",
                  text = paste("Your counts matrix has", ndup, "duplicated gene names.\nCounts of those genes will be merged."),
                  type = "warning",
                  btn_labels = "OK",
                  closeOnClickOutside = FALSE,
                )
              }
              if (nrow(df0) > 1 && NCOL(df0) > 1) {
                df <- as.matrix(df0)
                message("[UploadModule::upload_files] converting expression to counts...")
                df <- 2**df
                matname <- "counts.csv"
              }
            } else if (grepl("sample", fn1, ignore.case = TRUE)) {
              dbg("[upload_files] samples.csv : fn1 = ", fn1)
              df0 <- read.as_matrix(fn2)
              if (any(duplicated(rownames(df0)))) {
                dup.rows <- rownames(df0)[which(duplicated(rownames(df0)))]
                msg <- paste(
                  "Your samples file has duplicated entries: ",
                  dup.rows, ". This is not allowed, please correct."
                )
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Duplicated sample name",
                  text = msg,
                  type = "error",
                  btn_labels = "OK",
                  ## btn_colors = "red",
                  closeOnClickOutside = FALSE,
                )
              } else if (nrow(df0) > 1 && NCOL(df0) >= 1) {
                df <- as.data.frame(df0)
                matname <- "samples.csv"
              }
            } else if (grepl("contrast", fn1, ignore.case = TRUE)) {
              dbg("[upload_files] contrasts.csv : fn1 = ", fn1)
              df0 <- read.as_matrix(fn2)
              if (any(duplicated(rownames(df0)))) {
                dup.rows <- rownames(df0)[which(duplicated(rownames(df0)))]
                msg <- paste(
                  "Your contrasts file has duplicated entries: ",
                  dup.rows, ". This is not allowed, please correct."
                )
                shinyWidgets::sendSweetAlert(
                  session = session,
                  title = "Duplicated contrast name",
                  text = msg,
                  type = "error",
                  btn_labels = "OK",
                  ## btn_colors = "red",
                  closeOnClickOutside = FALSE,
                )
              } else if (nrow(df0) > 1 && NCOL(df0) >= 1) {
                df <- as.matrix(df0)
                matname <- "contrasts.csv"
              }
            }
            if (!is.null(matname)) {
              matlist[[matname]] <- df
            }
          }
        }
      }

      if ("counts.csv" %in% names(matlist)) {
        ## Convert to gene names (need for biological effects)
        dbg("[upload_files] converting probe names to symbols")
        X0 <- matlist[["counts.csv"]]
        pp <- rownames(X0)
        rownames(X0) <- probe2symbol(pp)
        sel <- !(rownames(X0) %in% c(NA, "", "NA"))
        X0 <- X0[sel, ]
        xx <- tapply(1:nrow(X0), rownames(X0), function(i) colSums(X0[i, , drop = FALSE]))
        X0 <- do.call(rbind, xx)
        matlist[["counts.csv"]] <- X0
      }

      ## put the matrices in the reactive values 'uploaded'
      files.needed <- c("counts.csv", "samples.csv", "contrasts.csv")
      if (length(matlist) > 0) {
        matlist <- matlist[which(names(matlist) %in% files.needed)]
        for (i in 1:length(matlist)) {
          colnames(matlist[[i]]) <- gsub("[\n\t ]", "_", colnames(matlist[[i]]))
          rownames(matlist[[i]]) <- gsub("[\n\t ]", "_", rownames(matlist[[i]]))
          if (names(matlist)[i] %in% c("counts.csv", "contrasts.csv")) {
            matlist[[i]] <- as.matrix(matlist[[i]])
          } else {
            matlist[[i]] <- type.convert(matlist[[i]])
          }
          m1 <- names(matlist)[i]
          message("[upload_files] updating matrix ", m1)
          uploaded[[m1]] <- matlist[[i]]
        }
        uploaded[["last_uploaded"]] <- names(matlist)
      }

      message("[upload_files] done!\n")
    })

    ## ------------------------------------------------------------------
    ## Observer for loading from local exampledata.zip file
    ##
    ## Reads in the data files from zip and puts in the
    ## reactive values object 'uploaded'. Then uploaded should
    ## trigger the computePGX module.
    ## ------------------------------------------------------------------
    shiny::observeEvent(input$load_example, {
      if (input$load_example) {
        zipfile <- file.path(FILES, "exampledata.zip")
        readfromzip1 <- function(file) {
          read.csv(unz(zipfile, file),
            check.names = FALSE, stringsAsFactors = FALSE,
            row.names = 1
          )
        }
        readfromzip2 <- function(file) {
          ## allows for duplicated names
          df0 <- read.csv(unz(zipfile, file), check.names = FALSE, stringsAsFactors = FALSE)
          mat <- as.matrix(df0[, -1])
          rownames(mat) <- as.character(df0[, 1])
          mat
        }
        uploaded$counts.csv <- readfromzip2("exampledata/counts.csv")
        uploaded$samples.csv <- readfromzip1("exampledata/samples.csv")
        uploaded$contrasts.csv <- readfromzip1("exampledata/contrasts.csv")
      } else {
        ## Remove files
        uploaded$counts.csv <- NULL
        uploaded$samples.csv <- NULL
        uploaded$contrasts.csv <- NULL
      }
    })

    ## ------------------------------------------------------------------
    ## Observer for loading CSV from local folder on
    ## host/server using URL. Reads the CSV files from folder
    ## and puts in the reactive values object 'uploaded'.
    ## ------------------------------------------------------------------

    if (ALLOW_URL_QUERYSTRING) {
      shiny::observeEvent(session$clientData$url_search, {
        ## -------------------------------------------------------------
        ## Parse URL query string
        ## -------------------------------------------------------------
        query <- parseQueryString(session$clientData$url_search)
        if (length(query) > 0) {
          dbg("[UploadModule:parseQueryString] names.query =", names(query))
          for (i in 1:length(query)) {
            dbg("[UploadModule:parseQueryString]", names(query)[i], "=>", query[[i]])
          }
        } else {
          dbg("[UploadModule:parseQueryString] no queryString!")
        }

        if (!is.null(query[["csv"]])) {
          qdir <- query[["csv"]]
          dbg("[UploadModule:parseQueryString] *** parseQueryString ***")
          dbg("[UploadModule:parseQueryString] qdir = ", qdir)

          counts_file <- file.path(qdir, "counts.csv")
          samples_file <- file.path(qdir, "samples.csv")
          if (!file.exists(counts_file)) {
            dbg("[SERVER:parseQueryString] ***ERROR*** missing counts.csv in dir = ", qdir)
          }
          if (!file.exists(samples_file)) {
            dbg("[SERVER:parseQueryString] ***ERROR*** missing samples.csv in dir = ", qdir)
          }
          if (!file.exists(counts_file) || !file.exists(samples_file)) {
            return(NULL)
          }

          FUN.readfromdir <- function() {
            dbg("[UploadModule:parseQueryString] *** loading CSV from dir = ", qdir, "***")

            readfromdir1 <- function(file) {
              read.csv(file,
                check.names = FALSE, stringsAsFactors = FALSE,
                row.names = 1
              )
            }
            readfromdir2 <- function(file) {
              ## allows for duplicated names
              df0 <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
              mat <- as.matrix(df0[, -1])
              rownames(mat) <- as.character(df0[, 1])
              mat
            }

            dbg("[UploadModule:parseQueryString] reading samples_csv = ", samples_file)
            uploaded$samples.csv <- readfromdir1(samples_file)

            dbg("[UploadModule:parseQueryString] reading samples_csv = ", samples_file)
            uploaded$counts.csv <- readfromdir2(counts_file)
            uploaded$contrasts.csv <- NULL

            meta_file <- file.path(qdir, "meta.txt")
            uploaded$meta <- NULL
            if (file.exists(meta_file)) {
              dbg("[UploadModule:parseQueryString] reading meta file = ", meta_file)
              ## meta <- read.table(meta_file,sep='\t',header=TRUE,row.names=1)
              meta <- read.table(meta_file, sep = "", header = TRUE, row.names = 1)
              meta <- as.list(array(meta[, 1], dimnames = list(rownames(meta))))
              uploaded$meta <- meta
            }
          }

          shinyalert::shinyalert(
            title = "Load CSV data from folder?",
            text = paste0("folder = ", qdir),
            callbackR = FUN.readfromdir,
            confirmButtonText = "Load!",
            type = "info"
          )

          dbg("[UploadModule:parseQueryString] dim(samples) = ", dim(uploaded$samples.csv))
          dbg("[UploadModule:parseQueryString] dim(counts) = ", dim(uploaded$counts.csv))

          ## focus on this tab
          updateTabsetPanel(session, "tabs", selected = "Upload data")
        }

        if (0 && !is.null(query[["pgx"]])) {
          qdir <- query[["pgx"]]
          dbg("[UploadModule:parseQueryString] pgx =>", qdir)

          pgx_file <- query[["pgx"]]
          pgx_file <- paste0(sub("[.]pgx$", "", pgx_file), ".pgx")
          dbg("[UploadModule:parseQueryString] pgx_file = ", pgx_file)

          if (!file.exists(pgx_file)) {
            dbg("[SERVER:parseQueryString] ***ERROR*** missing pgx_file", pgx_file)
            return(NULL)
          }

          dbg("[UploadModule:parseQueryString] 1:")

          FUN.readPGX <- function() {
            dbg("[UploadModule:parseQueryString] *** loading PGX file = ", pgx_file, "***")

            load(pgx_file) ## load NGS/PGX
            uploaded$pgx <- ngs
            remove(ngs)

            uploaded$meta <- NULL
          }

          dbg("[UploadModule:parseQueryString] 2:")

          shinyalert::shinyalert(
            title = "Load PGX data from folder?",
            text = paste0("folder = ", qdir),
            callbackR = FUN.readPGX,
            confirmButtonText = "Load!",
            type = "info"
          )

          dbg("[UploadModule:parseQueryString] 3:")

          ## focus on this tab
          updateTabsetPanel(session, "tabs", selected = "Upload data")

          dbg("[UploadModule:parseQueryString] 4:")
        }
      })
    }

    ## =====================================================================
    ## ===================== checkTables ===================================
    ## =====================================================================

    checkTables <- shiny::reactive({
      ## check dimensions
      status <- rep("please upload", 3)
      files.needed <- c("counts.csv", "samples.csv", "contrasts.csv")
      names(status) <- files.needed
      files.nrow <- rep(NA, 3)
      files.ncol <- rep(NA, 3)

      for (i in 1:3) {
        fn <- files.needed[i]
        upfile <- uploaded[[fn]]
        if (fn %in% names(uploaded) && !is.null(upfile)) {
          status[i] <- "OK"
          files.nrow[i] <- nrow(upfile)
          files.ncol[i] <- ncol(upfile)
        }
      }

      has.pgx <- ("pgx" %in% names(uploaded))
      if (has.pgx) has.pgx <- has.pgx && !is.null(uploaded[["pgx"]])
      if (has.pgx == TRUE) {
        ## Nothing to check. Always OK.
      } else if (!has.pgx) {
        ## check rownames of samples.csv
        if (status["samples.csv"] == "OK" && status["counts.csv"] == "OK") {
          samples1 <- uploaded[["samples.csv"]]
          counts1 <- uploaded[["counts.csv"]]
          a1 <- mean(rownames(samples1) %in% colnames(counts1))
          a2 <- mean(samples1[, 1] %in% colnames(counts1))

          if (a2 > a1 && NCOL(samples1) > 1) {
            message("[UploadModuleServer] getting sample names from first column\n")
            rownames(samples1) <- samples1[, 1]
            uploaded[["samples.csv"]] <- samples1[, -1, drop = FALSE]
          }
        }

        ## check files: matching dimensions
        if (status["counts.csv"] == "OK" && status["samples.csv"] == "OK") {
          nsamples <- max(ncol(uploaded[["counts.csv"]]), nrow(uploaded[["samples.csv"]]))
          ok.samples <- intersect(
            rownames(uploaded$samples.csv),
            colnames(uploaded$counts.csv)
          )
          n.ok <- length(ok.samples)
          message("[UploadModule::checkTables] n.ok = ", n.ok)
          if (n.ok > 0 && n.ok < nsamples) {
            ## status["counts.csv"]  = "WARNING: some samples with missing annotation)"
          }

          if (n.ok > 0) {
            message("[UploadModule::checkTables] conforming samples/counts...")
            uploaded[["samples.csv"]] <- uploaded$samples.csv[ok.samples, , drop = FALSE]
            uploaded[["counts.csv"]] <- uploaded$counts.csv[, ok.samples, drop = FALSE]
          }

          if (n.ok == 0) {
            status["counts.csv"] <- "ERROR: colnames do not match (with samples)"
            status["samples.csv"] <- "ERROR: rownames do not match (with counts)"
          }

          dbg("[UploadModule::checkTables] dim(samples.csv) = ", dim(uploaded$samples.csv))
          dbg("[UploadModule::checkTables] dim(counts.csv) = ", dim(uploaded$counts.csv))
        }

        if (status["contrasts.csv"] == "OK" && status["samples.csv"] == "OK") {
          samples1 <- uploaded[["samples.csv"]]
          contrasts1 <- uploaded[["contrasts.csv"]]
          group.col <- grep("group", tolower(colnames(samples1)))
          old1 <- (length(group.col) > 0 &&
            nrow(contrasts1) < nrow(samples1) &&
            all(rownames(contrasts1) %in% samples1[, group.col])
          )
          old2 <- all(rownames(contrasts1) == rownames(samples1)) &&
            all(unique(as.vector(contrasts1)) %in% c(-1, 0, 1, NA))

          old.style <- (old1 || old2)
          if (old.style && old1) {
            message("[UploadModule] WARNING: converting old1 style contrast to new format")
            new.contrasts <- samples1[, 0]
            if (NCOL(contrasts1) > 0) {
              new.contrasts <- contrastAsLabels(contrasts1)
              grp <- as.character(samples1[, group.col])
              new.contrasts <- new.contrasts[grp, , drop = FALSE]
              rownames(new.contrasts) <- rownames(samples1)
            }
            dbg("[UploadModule] old.ct1 = ", paste(contrasts1[, 1], collapse = " "))
            dbg("[UploadModule] old.nn = ", paste(rownames(contrasts1), collapse = " "))
            dbg("[UploadModule] new.ct1 = ", paste(new.contrasts[, 1], collapse = " "))
            dbg("[UploadModule] new.nn = ", paste(rownames(new.contrasts), collapse = " "))

            contrasts1 <- new.contrasts
          }
          if (old.style && old2) {
            message("[UploadModule] WARNING: converting old2 style contrast to new format")
            new.contrasts <- samples1[, 0]
            if (NCOL(contrasts1) > 0) {
              new.contrasts <- contrastAsLabels(contrasts1)
              rownames(new.contrasts) <- rownames(samples1)
            }
            contrasts1 <- new.contrasts
          }

          dbg("[UploadModule] 1 : dim.contrasts1 = ", dim(contrasts1))
          dbg("[UploadModule] 1 : dim.samples1   = ", dim(samples1))

          ok.contrast <- length(intersect(rownames(samples1), rownames(contrasts1))) > 0
          if (ok.contrast && NCOL(contrasts1) > 0) {
            ## always clean up
            contrasts1 <- apply(contrasts1, 2, as.character)
            rownames(contrasts1) <- rownames(samples1)
            for (i in 1:ncol(contrasts1)) {
              isz <- (contrasts1[, i] %in% c(NA, "NA", "NA ", "", " ", "  ", "   ", " NA"))
              if (length(isz)) contrasts1[isz, i] <- NA
            }
            uploaded[["contrasts.csv"]] <- contrasts1
            status["contrasts.csv"] <- "OK"
          } else {
            uploaded[["contrasts.csv"]] <- NULL
            status["contrasts.csv"] <- "ERROR: dimension mismatch"
          }
        }

        MAXSAMPLES <- 25
        MAXCONTRASTS <- 5
        MAXSAMPLES <- as.integer(limits["samples"])
        MAXCONTRASTS <- as.integer(limits["comparisons"])

        ## check files: maximum contrasts allowed
        if (status["contrasts.csv"] == "OK") {
          if (ncol(uploaded[["contrasts.csv"]]) > MAXCONTRASTS) {
            status["contrasts.csv"] <- paste("ERROR: max", MAXCONTRASTS, "contrasts allowed")
          }
        }

        ## check files: maximum samples allowed
        if (status["counts.csv"] == "OK" && status["samples.csv"] == "OK") {
          if (ncol(uploaded[["counts.csv"]]) > MAXSAMPLES) {
            status["counts.csv"] <- paste("ERROR: max", MAXSAMPLES, " samples allowed")
          }
          if (nrow(uploaded[["samples.csv"]]) > MAXSAMPLES) {
            status["samples.csv"] <- paste("ERROR: max", MAXSAMPLES, "samples allowed")
          }
        }

        ## check samples.csv: must have group column defined
        if (status["samples.csv"] == "OK" && status["contrasts.csv"] == "OK") {
          samples1 <- uploaded[["samples.csv"]]
          contrasts1 <- uploaded[["contrasts.csv"]]
          if (!all(rownames(contrasts1) %in% rownames(samples1))) {
            status["contrasts.csv"] <- "ERROR: contrasts do not match samples"
          }
        }
      } ## end-if-from-pgx

      e1 <- grepl("ERROR", status["samples.csv"])
      e2 <- grepl("ERROR", status["contrasts.csv"])
      e3 <- grepl("ERROR", status["counts.csv"])
      s1 <- "samples.csv" %in% uploaded$last_uploaded
      s2 <- "contrasts.csv" %in% uploaded$last_uploaded
      s3 <- "counts.csv" %in% uploaded$last_uploaded

      if (e1 || e2 || e3) {
        message("[checkTables] ERROR in samples table : e1 = ", e1)
        message("[checkTables] ERROR in contrasts table : e2 = ", e2)
        message("[checkTables] ERROR in counts table : e2 = ", e3)

        if (e1 && !s1) {
          uploaded[["samples.csv"]] <- NULL
          status["samples.csv"] <- "please upload"
        }
        if (e2 && !s2) {
          uploaded[["contrasts.csv"]] <- NULL
          status["contrasts.csv"] <- "please upload"
        }
        if (e3 && !s3) {
          uploaded[["counts.csv"]] <- NULL
          status["counts.csv"] <- "please upload"
        }
      }


      if (!is.null(uploaded$contrasts.csv) &&
        (is.null(uploaded$counts.csv) ||
          is.null(uploaded$samples.csv))) {
        uploaded[["contrasts.csv"]] <- NULL
        status["contrasts.csv"] <- "please upload"
      }


      ## check files
      description <- c(
        "Count/expression file with gene on rows, samples as columns",
        "Samples file with samples on rows, phenotypes as columns",
        ## "Gene information file with genes on rows, gene info as columns.",
        "Contrast file with conditions on rows, contrasts as columns"
      )
      df <- data.frame(
        filename = files.needed,
        description = description,
        nrow = files.nrow,
        ncol = files.ncol,
        status = status
      )
      rownames(df) <- files.needed

      ## deselect
      ## DT::selectRows(proxy = DT::dataTableProxy("pgxtable"), selected=NULL)
      return(df)
    })

    output$downloadExampleData <- shiny::downloadHandler(
      filename = "exampledata.zip",
      content = function(file) {
        zip <- file.path(FILES, "exampledata.zip")
        file.copy(zip, file)
      }
    )

    output$upload_info <- shiny::renderUI({
      upload_info <- "<h4>User file upload</h4><p>Please prepare the data files in CSV format as listed below. It is important to name the files exactly as shown. The file format must be comma-separated-values (CSV) text. Be sure the dimensions, rownames and column names match for all files. You can download a zip file with example files here: EXAMPLEZIP. You can upload a maximum of <u>LIMITS</u>."
      DLlink <- shiny::downloadLink(ns("downloadExampleData"), "exampledata.zip")
      upload_info <- sub("EXAMPLEZIP", DLlink, upload_info)

      limits0 <- paste(
        limits["datasets"], "datasets (with each up to",
        limits["samples"], "samples and",
        limits["comparisons"], "comparisons)"
      )
      upload_info <- sub("LIMITS", limits0, upload_info)
      shiny::HTML(upload_info)
    })

    ## =====================================================================
    ## ========================= SUBMODULES/SERVERS ========================
    ## =====================================================================

    ## correctedX <- shiny::reactive({
    normalized_counts <- NormalizeCountsServerRT(
      id = "normalize",
      counts = shiny::reactive(uploaded$counts.csv),
      height = height
    )

    ## correctedX <- shiny::reactive({
    correctedX <- BatchCorrectServer(
      id = "batchcorrect",
      X = shiny::reactive(uploaded$counts.csv),
      ## X = normalized_counts,  ## NOT YET!!!!
      is.count = TRUE,
      pheno = shiny::reactive(uploaded$samples.csv),
      height = height
    )

    corrected_counts <- shiny::reactive({
      counts <- NULL
      advanced_mode <- (length(input$advanced_mode) > 0 &&
        input$advanced_mode[1] == 1)
      if (advanced_mode) {
        out <- correctedX()
        counts <- pmax(2**out$X - 1, 0)
      } else {
        counts <- uploaded$counts.csv
      }
      counts
    })

    ## mkContrast <- shiny::reactive({
    modified_ct <- MakeContrastServerRT(
      id = "makecontrast",
      phenoRT = shiny::reactive(uploaded$samples.csv),
      contrRT = shiny::reactive(uploaded$contrasts.csv),
      ## countsRT = shiny::reactive(uploaded$counts.csv),
      countsRT = corrected_counts,
      height = height
    )

    shiny::observeEvent(modified_ct(), {
      ## Monitor for changes in the contrast matrix and if
      ## so replace the uploaded reactive values.
      ##
      modct <- modified_ct()
      uploaded$contrasts.csv <- modct$contr
      uploaded$samples.csv <- modct$pheno
    })

    upload_ok <- shiny::reactive({
      check <- checkTables()
      all(check[, "status"] == "OK")
      all(grepl("ERROR", check[, "status"]) == FALSE)
    })

    batch_vectors <- shiny::reactive({
      correctedX()$B
    })

    ## computed_pgx <- ComputePgxServer(
    computed_pgx <- ComputePgxServer(
      id = "compute",
      ## countsRT = shiny::reactive(uploaded$counts.csv),
      countsRT = corrected_counts,
      samplesRT = shiny::reactive(uploaded$samples.csv),
      contrastsRT = shiny::reactive(uploaded$contrasts.csv),
      batchRT = batch_vectors,
      metaRT = shiny::reactive(uploaded$meta),
      enable_button = upload_ok,
      alertready = FALSE,
      FILES = FILES,
      pgx.dirRT = shiny::reactive(getPGXDIR()),
      max.genes = as.integer(limits["genes"]),
      max.genesets = as.integer(limits["genesets"]),
      max.datasets = as.integer(limits["datasets"]),
      height = height
    )

    uploaded_pgx <- shiny::reactive({
      if (!is.null(uploaded$pgx)) {
        pgx <- uploaded$pgx
        ## pgx <- pgx.initialize(pgx)
      } else {
        pgx <- computed_pgx()
      }
      return(pgx)
    })

    ## =====================================================================
    ## ===================== PLOTS AND TABLES ==============================
    ## =====================================================================

    output$countStats <- shiny::renderPlot({

      check <- checkTables()
      status.ok <- check["counts.csv", "status"]
      dbg("[countStats] status.ok = ", status.ok)

      if (status.ok != "OK") {
        frame()
        status.ds <- check["counts.csv", "description"]
        msg <- paste(
          toupper(status.ok), "\n", "(Required) Upload 'counts.csv'",
          tolower(status.ds)
        )
        graphics::text(0.5, 0.5, paste(strwrap(msg, 30), collapse = "\n"), col = "grey25")
        graphics::box(lty = 2, col = "grey60")
        return(NULL)
      }

      counts <- uploaded[["counts.csv"]]
      xx <- log2(1 + counts)
      if (nrow(xx) > 1000) xx <- xx[sample(1:nrow(xx), 1000), , drop = FALSE]
      ## dc <- reshape::melt(xx)
      suppressWarnings(dc <- data.table::melt(xx))
      dc$value[dc$value == 0] <- NA
      tt2 <- paste(nrow(counts), "genes x", ncol(counts), "samples")
      ggplot2::ggplot(dc, ggplot2::aes(x = value, color = Var2)) +
        ggplot2::geom_density() +
        ggplot2::xlab("log2(1+counts)") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle("COUNTS", subtitle = tt2)
    })

    output$phenoStats <- shiny::renderPlot({
      dbg("[phenoStats] renderPlot called \n")
      ## req(uploaded$samples.csv)

      check <- checkTables()
      status.ok <- check["samples.csv", "status"]
      if (status.ok != "OK") {
        frame()
        status.ds <- check["samples.csv", "description"]
        msg <- paste(
          toupper(status.ok), "\n", "(Required) Upload 'samples.csv'",
          tolower(status.ds)
        )
        graphics::text(0.5, 0.5, paste(strwrap(msg, 30), collapse = "\n"), col = "grey25")
        graphics::box(lty = 2, col = "grey60")
        return(NULL)
      }

      pheno <- uploaded[["samples.csv"]]
      px <- head(colnames(pheno), 20) ## show maximum??

      df <- type.convert(pheno[, px, drop = FALSE])
      vt <- df %>% inspectdf::inspect_types()
      vt

      ## discretized continuous variable into 10 bins
      ii <- unlist(vt$col_name[c("numeric", "integer")])
      ii
      if (!is.null(ii) && length(ii)) {
        cat("[UploadModule::phenoStats] discretizing variables:", ii, "\n")
        df[, ii] <- apply(df[, ii, drop = FALSE], 2, function(x) {
          if (any(is.infinite(x))) x[which(is.infinite(x))] <- NA
          cut(x, breaks = 10)
        })
      }

      p1 <- df %>%
        inspectdf::inspect_cat() %>%
        inspectdf::show_plot()
      tt2 <- paste(nrow(pheno), "samples x", ncol(pheno), "phenotypes")
      ## tt2 <- paste(ncol(pheno),"phenotypes")
      p1 <- p1 + ggplot2::ggtitle("PHENOTYPES", subtitle = tt2) +
        ggplot2::theme(
          ## axis.text.x = ggplot2::element_text(size=8, vjust=+5),
          axis.text.y = ggplot2::element_text(
            size = 12,
            margin = ggplot2::margin(0, 0, 0, 25),
            hjust = 1
          )
        )

      p1
    })

    output$contrastStats <- shiny::renderPlot({
      ## req(uploaded$contrasts.csv)
      ct <- uploaded$contrasts.csv
      has.contrasts <- !is.null(ct) && NCOL(ct) > 0
      check <- checkTables()
      status.ok <- check["contrasts.csv", "status"]

      if (status.ok != "OK" || !has.contrasts) {
        frame()
        status.ds <- check["contrasts.csv", "description"]
        msg <- paste(
          toupper(status.ok), "\n", "(Optional) Upload 'contrasts.csv'",
          tolower(status.ds)
        )
        ## text(0.5,0.5,"Please upload contrast file 'contrast.csv' with conditions on rows, contrasts as columns")
        graphics::text(0.5, 0.5, paste(strwrap(msg, 30), collapse = "\n"), col = "grey25")
        graphics::box(lty = 2, col = "grey60")
        return(NULL)
      }

      dbg("[output$contrastStats] 2 : ")

      contrasts <- uploaded$contrasts.csv

      dbg("[output$contrastStats] 3 : ")

      ## contrasts <- sign(contrasts)
      ## df <- contrastAsLabels(contrasts)
      df <- contrasts
      px <- head(colnames(df), 20) ## maximum to show??
      df <- data.frame(df[, px, drop = FALSE], check.names = FALSE)
      tt2 <- paste(nrow(contrasts), "samples x", ncol(contrasts), "contrasts")
      ## tt2 <- paste(ncol(contrasts),"contrasts")
      dbg("[output$contrastStats] 4a : dim.df=", dim(df))

      p1 <- df %>%
        inspectdf::inspect_cat() %>%
        inspectdf::show_plot()
      dbg("[output$contrastStats] 4b : ")

      p1 <- p1 + ggplot2::ggtitle("CONTRASTS", subtitle = tt2) +
        ggplot2::theme(
          ## axis.text.x = ggplot2::element_text(size=8, vjust=+5),
          axis.text.y = ggplot2::element_text(
            size = 12,
            margin = ggplot2::margin(0, 0, 0, 25),
            hjust = 1
          )
        )

      dbg("[output$contrastStats] 5 : ")

      p1
    })

    sel.conditions <- shiny::reactive({
      shiny::req(phenoRT(), corrected_counts())
      df <- phenoRT()

      if ("<samples>" %in% input$param) {
        df$"<samples>" <- rownames(df)
      }
      if ("<gene>" %in% input$param) {
        gene <- input$gene
        if (gene %in% rownames(corrected_counts())) {
          gx <- log2(1 + corrected_counts()[gene, ])
          ## df$"<gene>" <- c("low","high")[1 + 1*(gx >= mean(gx,na.rm=TRUE))]
          df$"<gene>" <- gx
        } else {
          return(NULL)
        }
      }

      df <- type.convert(df)
      ii <- which(sapply(type.convert(df), class) %in% c("numeric", "integer"))
      ii
      if (length(ii)) {
        for (i in ii) {
          x <- df[, i]
          df[, i] <- c("low", "high")[1 + 1 * (x >= mean(x, na.rm = TRUE))]
        }
      }

      pp <- intersect(input$param, colnames(df))
      ss <- colnames(corrected_counts())
      cond <- apply(df[ss, pp, drop = FALSE], 1, paste, collapse = "_")
      cond <- gsub("^_|_$", "", cond)
      cond
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

      message("[MakeContrastServer:addcontrast] 1 : ")

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

      message("[MakeContrastServer:addcontrast] update reactive values : 1")

      ## update reactive value
      samples <- colnames(corrected_counts())

      message("[MakeContrastServer:addcontrast] 1 : samples = ", samples)
      message("[MakeContrastServer:addcontrast] 1 : ct.name = ", ct.name)
      message("[MakeContrastServer:addcontrast] 1 : len.ctx = ", length(ctx))

      ctx1 <- matrix(ctx, ncol = 1, dimnames = list(samples, ct.name))
      if (is.null(rv$contr)) {
        rv$contr <- ctx1
      } else {
        rv$contr <- cbind(rv$contr, ctx1)
      }

      message("[MakeContrastServer:addcontrast] update reactive values : 2")
      message("[MakeContrastServer:addcontrast] ct.name in pheno = ", ct.name %in% colnames(rv$pheno))

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
          ## header = shiny::h4("Create comparison:"),
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

          ct1 <- makeContrastsFromLabelMatrix(ct)
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

    shiny::observeEvent(input$contrast_delete, {
      ## Observe if a contrast is to be deleted
      ##
      id <- as.numeric(gsub(".*_", "", input$contrast_delete))
      message("[contrast_delete] clicked on delete contrast", id)
      if (length(id) == 0) {
        return(NULL)
      }
      ## updateActionButton(session, paste0("contrast_delete_",id),label="XXX")
      if (!is.null(rv$contr) && NCOL(rv$contr) <= 1) {
        rv$contr <- rv$contr[, 0, drop = FALSE]
      } else {
        rv$contr <- rv$contr[, -id, drop = FALSE]
      }
    })

    output$checkTablesOutput <- DT::renderDataTable({
      ## Render the upload status table
      ##
      if (!input$advanced_mode) {
        return(NULL)
      }
      df <- checkTables()
      dt <- DT::datatable(
        df,
        rownames = FALSE,
        selection = "none",
        class = "compact cell-border",
        options = list(
          dom = "t"
        )
      ) %>%
        DT::formatStyle(0, target = "row", fontSize = "12px", lineHeight = "100%")
    })

    upload_plot_pcaplot_server(
      "pcaplot",
      phenoRT = phenoRT,
      countsRT = corrected_counts,
      sel.conditions = sel.conditions,
      watermark = FALSE
    )

    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    res <- list(
      loaded = loadedDataset
    )
    return(res)
  })
}
