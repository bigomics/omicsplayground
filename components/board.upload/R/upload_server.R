##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

UploadBoard <- function(id,
                        pgx_dir,
                        pgx,
                        auth,
                        reload_pgxdir,
                        load_uploaded_data,
                        recompute_pgx,
                        recompute_info) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    # Some 'global' reactive variables used in this file
    uploaded <- shiny::reactiveValues()

    output$navheader <- shiny::renderUI({
      fillRow(
        flex = c(NA, 1, NA),
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
        title = shiny::HTML("<strong>How to upload new data</strong>"),
        shiny::HTML(module_infotext),
        easyClose = TRUE,
        size = "xl"
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
<center><iframe width="560" height="315" src="https://www.youtube.com/embed/elwT6ztt3Fo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe><center>'
    )

    module_infotext <- HTML('<center><iframe width="1120" height="630" src="https://www.youtube.com/embed/elwT6ztt3Fo" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe><center>')


    ## ================================================================================
    ## ====================== NEW DATA UPLOAD =========================================
    ## ================================================================================

    ## keeps track of how pgx was obtained: uploaded or computed. NEED
    ## RETHINK: is this robust to multiple users on same R process?
    uploaded_method <- NA

    shiny::observeEvent(uploaded_pgx(), {
      dbg("[UploadBoard:observe:uploaded_pgx] uploaded PGX detected!")

      new_pgx <- uploaded_pgx()

      ## NEED RETHINK: if "uploaded" we unneccessarily saving the pgx
      ## object again.  We should skip saving and pass the filename to
      ## pgxfile to be sure the filename is correct.

      ## new_pgx <- playbase::pgx.initialize(new_pgx)  ## already done later
      ## -------------- save PGX file/object ---------------
      pgxfile <- sub("[.]pgx$", "", new_pgx$name)
      pgxfile <- gsub("^[./-]*", "", pgxfile) ## prevent going to parent folder
      pgxfile <- paste0(gsub("[ \\/]", "_", pgxfile), ".pgx")
      pgxdir <- auth$user_dir
      fn <- file.path(pgxdir, pgxfile)
      fn <- iconv(fn, from = "", to = "ASCII//TRANSLIT")
      ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ## switch 'pgx' as standard name. Actually saving as RDS
      ## would have been better...
      dbg("[UploadBoard:observe::uploaded_pgx] saving pgx as = ", fn)
      playbase::pgx.save(new_pgx, file = fn)

      shiny::withProgress(message = "Scanning dataset library...", value = 0.33, {
        playbase::pgxinfo.updateDatasetFolder(
          pgxdir,
          new.pgx = pgxfile,
          update.sigdb = FALSE
        )
      })

      ## trigger reload of pgx table
      reload_pgxdir(reload_pgxdir() + 1)

      beepr::beep(10) ## short beep

      load_my_dataset <- function() {
        if (input$confirmload) {
          load_uploaded_data(pgxfile)
        }
      }

      if (uploaded_method == "computed") {
        shinyalert::shinyalert(
          title = paste("Your dataset is ready!"),
          text = paste("Your dataset", new_pgx$name, "is ready for visualization. Happy discoveries!"),
          confirmButtonText = "Load my new data!",
          showCancelButton = TRUE,
          cancelButtonText = "Stay here.",
          inputId = "confirmload",
          closeOnEsc = FALSE,
          immediate = TRUE,
          callbackR = load_my_dataset
        )
      } else {
        shinyalert::shinyalert(
          title = paste("Dataset is loaded!"),
          text = paste("Your uploaded dataset", new_pgx$name, "is ready for visualization. Happy discoveries!"),
          confirmButtonText = "Show dataset!",
          showCancelButton = FALSE,
          inputId = "confirmload",
          closeOnEsc = FALSE,
          immediate = TRUE,
          callbackR = load_my_dataset
        )
      }
    })

    # create an observer that will hide tabs Upload Files if selected organism if null and show if the button proceed_to_upload is clicked
    observeEvent(input$proceed_to_upload,{
      # show tab Upload Files
      shiny::showTab("tabs", "Upload Files")
      # for some reason, we need to click in another tab to make the upload tab run
      shinyjs::runjs('document.querySelector("a[data-value=\'Upload Files\']").click();')
      shinyjs::runjs('document.querySelector("a[data-value=\'Select Organism\']").click();')
      shinyjs::runjs('document.querySelector("a[data-value=\'Upload Files\']").click();')
    })

    
    # upload tab should be hidden when when no selected organism is selected (when componene is loaded)
    shiny::hideTab("tabs", "Upload Files")


    ## Hide/show tabpanels upon available data like a wizard dialog
    shiny::observe({
      has.upload <- Vectorize(function(f) {
        (f %in% names(uploaded) && !is.null(nrow(uploaded[[f]])))
      })
      need2 <- c("counts.csv", "samples.csv")
      need3 <- c("counts.csv", "samples.csv", "contrasts.csv")
      if (all(has.upload(need3))) {
        shiny::showTab("tabs", "Upload Files")
        shiny::showTab("tabs", "Comparisons")
        shiny::showTab("tabs", "Compute")
        if (input$advanced_mode) {
          shiny::showTab("tabs", "BatchCorrect")
        }
      } else if (all(has.upload(need2))) {
        if (input$advanced_mode) {
          shiny::showTab("tabs", "BatchCorrect")
        }
        shiny::showTab("tabs", "Comparisons")
        shiny::hideTab("tabs", "Compute")
      } else {
        #shiny::hideTab("tabs", "Upload Files")
        shiny::hideTab("tabs", "BatchCorrect")
        shiny::hideTab("tabs", "Comparisons")
        shiny::hideTab("tabs", "Compute")
      }
    })

    ## =====================================================================
    ## ======================= UI OBSERVERS ================================
    ## =====================================================================

    shiny::observeEvent(input$advanced_mode, {
      if (input$advanced_mode) {
        shiny::showTab("tabs", "BatchCorrect")
      } else {
        shiny::hideTab("tabs", "BatchCorrect")
      }
    })



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

    # this directory is used to save pgx files, logs, inputs, etc..
    raw_dir <- reactiveVal(NULL)

    create_raw_dir <- function(auth) {
      auth_id <- ifelse(!auth$email %in% c("", NA), auth$email, auth$username)
      prefix <- paste0("raw_", auth_id, "_")
      raw_dir <- tempfile(pattern = prefix, tmpdir = file.path(PGX.DIR, "USER_INPUT"))
      dir.create(raw_dir, recursive = TRUE)
      dbg("[UploadBoard:raw_dir<-eventReactive] creating raw_dir", raw_dir)
      raw_dir
    }

    shiny::observeEvent(input$upload_files, {
      if (is.null(raw_dir())) {
        raw_dir(create_raw_dir(auth))
      }

      upload_table <- input$upload_files

      if(class(upload_table) != "data.frame" && upload_table == "hello_example"){
        upload_table <- data.frame(
          name = c("counts.csv", "samples.csv", "contrasts.csv"),
          type = c("text/csv", "text/csv", "text/csv"),
          datapath = c("examplecounts","examplesamples","examplecontrasts")
        )
      }

      message("[upload_files] >>> reading uploaded files")
      message("[upload_files] upload_files$name=", upload_table$name)
      message("[upload_files] upload_files$datapath=", upload_table$datapath)

      uploaded[["counts.csv"]] <- NULL
      uploaded[["samples.csv"]] <- NULL
      uploaded[["contrasts.csv"]] <- NULL
      uploaded[["pgx"]] <- NULL
      uploaded[["last_uploaded"]] <- NULL
      uploaded[["checklist"]] <- NULL

      ## read uploaded files
      pgx.uploaded <- any(grepl("[.]pgx$", upload_table$name))
      matlist <- list()
      checklist <- list()

      if (pgx.uploaded) {
        ## If the user uploaded a PGX file, we extract the matrix
        ## dimensions from the given PGX/NGS object. Really?
        i <- grep("[.]pgx$", upload_table$name)
        pgxfile <- upload_table$datapath[i]
        uploaded[["pgx"]] <- local(get(load(pgxfile, verbose = 0))) ## override any name
      } else {
        ## If the user uploaded CSV files, we read in the data
        ## from the files.

        ii <- grep("csv$", upload_table$name)
        ii <- grep("sample|count|contrast|expression|comparison",
          upload_table$name,
          ignore.case = TRUE
        )
        if (length(ii) == 0) {
          return(NULL)
          }

        inputnames <- upload_table$name[ii]
        uploadnames <- upload_table$datapath[ii]
        message("[upload_files] uploaded files: ", inputnames)

        ## remove any old gui_contrasts.csv
        user_ctfile <- file.path(raw_dir(), "user_contrasts.csv")
        if (file.exists(user_ctfile)) unlink(user_ctfile)

        error_list <- playbase::PGX_CHECKS

        if (length(uploadnames) > 0) {
          for (i in 1:length(uploadnames)) {
            #i = 1
            fn1 <- inputnames[i]
            fn2 <- uploadnames[i]
            matname <- NULL
            df <- NULL
            IS_COUNT <- grepl("count", fn1, ignore.case = TRUE)
            IS_EXPRESSION <- grepl("expression", fn1, ignore.case = TRUE)
            IS_SAMPLE <- grepl("sample", fn1, ignore.case = TRUE)
            IS_CONTRAST <- grepl("contrast|comparison", fn1, ignore.case = TRUE)

            if (IS_COUNT || IS_EXPRESSION) {
              ## allows duplicated rownames
              if(fn2 == "examplecounts"){
                df0 <- playbase::COUNTS
                # save a warning file telling this folder is example data
                writeLines("", file.path(raw_dir(), "EXAMPLE_DATA"))
              } else {
                df0 <- playbase::read.as_matrix(fn2)
              }
              
              # save input as raw file in raw_dir
              file.copy(fn2, file.path(raw_dir(), "raw_counts.csv"))

              COUNTS_check <- playbase::pgx.checkINPUT(df0, "COUNTS")
              check <- COUNTS_check$check

              if (COUNTS_check$PASS && IS_COUNT) {
                df <- as.matrix(COUNTS_check$df)
                matname <- "counts.csv"
              }

              if (COUNTS_check$PASS && IS_EXPRESSION) {
                df <- as.matrix(COUNTS_check$df)
                message("[UploadModule::upload_files] converting expression to counts...")
                df <- 2**df
                matname <- "counts.csv"
              }

              # store check and data regardless of it errors
              checklist[["counts.csv"]]$check <- check
              checklist[["counts.csv"]]$file <- COUNTS_check$df
            }

            if (IS_SAMPLE) {
              if(fn2 == "examplesamples"){
                df0 <- playbase::SAMPLES
              } else {
                df0 <- playbase::read.as_matrix(fn2)
              }
              # save input as raw file in raw_dir
              file.copy(fn2, file.path(raw_dir(), "raw_samples.csv"))

              SAMPLES_check <- playbase::pgx.checkINPUT(df0, "SAMPLES")
              check <- SAMPLES_check$check

              if (SAMPLES_check$PASS && IS_SAMPLE) {
                df <- as.data.frame(SAMPLES_check$df)
                matname <- "samples.csv"
              }
              # store check and data regardless of it errors
              checklist[["samples.csv"]]$check <- check
              checklist[["samples.csv"]]$file <- SAMPLES_check$df
            }

            if (IS_CONTRAST) {
              if(fn2 == "examplecontrasts"){
                df0 <- playbase::CONTRASTS
              } else {
                df0 <- playbase::read.as_matrix(fn2)
              }
              # save input as raw file in raw_dir
              file.copy(fn2, file.path(raw_dir(), "raw_contrasts.csv"))

              CONTRASTS_check <- playbase::pgx.checkINPUT(df0, "CONTRASTS")
              check <- CONTRASTS_check$check

              if (CONTRASTS_check$PASS) {
                df <- as.matrix(CONTRASTS_check$df)
                matname <- "contrasts.csv"
              }
              # store check and data regardless of it errors
              checklist[["contrasts.csv"]]$check <- check
              checklist[["contrasts.csv"]]$file <- CONTRASTS_check$df
            }

            if (!is.null(matname)) {
              matlist[[matname]] <- df
            }
          }
        }
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
        uploaded[["checklist"]] <- checklist
      }

      message("[upload_files] done!\n")
    })

    # In case the user is reanalysing the data, get the info from pgx
    observeEvent(recompute_pgx(), {
      pgx <- recompute_pgx()
      uploaded$samples.csv <- pgx$samples
      uploaded$contrasts.csv <- pgx$contrast
      uploaded$counts.csv <- pgx$counts
      corrected_counts <- pgx$counts
      recompute_info(list("name" = pgx$name, "description" = pgx$description))
    })


    ## ------------------------------------------------------------------
    ## Observer for loading example data
    ##
    ## Reads in the data files from zip and puts in the
    ## reactive values object 'uploaded'. Then uploaded should
    ## trigger the computePGX module.
    ## ------------------------------------------------------------------
    shiny::observeEvent(input$load_example, {
        shiny::updateSelectInput(session, "selected_organism", choices = "Human")
        shiny::showTab("tabs", "Upload Files")

        shinyjs::runjs('document.querySelector("a[data-value=\'Upload Files\']").click();')

        # go to upload tab
                
        shinyjs::runjs('document.queryselector("a[data-value=\'Comparisons\']").click();')

        shinyjs::runjs("Shiny.setInputValue('upload-upload_files', 'hello_example',{priority: 'event'});")

      })

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

      error_list <- playbase::PGX_CHECKS

      has.pgx <- ("pgx" %in% names(uploaded))
      if (has.pgx) has.pgx <- has.pgx && !is.null(uploaded[["pgx"]])
      if (has.pgx == TRUE) {
        ## Nothing to check. Always OK.
      } else if (!has.pgx) {
        ## check rownames of samples.csv
        if (status["samples.csv"] == "OK" && status["counts.csv"] == "OK") {
          FILES_check <- playbase::pgx.crosscheckINPUT(
            SAMPLES = uploaded[["samples.csv"]],
            COUNTS = uploaded[["counts.csv"]]
          )
          uploaded[["checklist"]][["samples_counts"]] <- FILES_check$check


          uploaded[["samples.csv"]] <- FILES_check$SAMPLES
          uploaded[["counts.csv"]] <- FILES_check$COUNTS
          samples1 <- FILES_check$SAMPLES
          counts1 <- FILES_check$COUNTS
          a1 <- mean(rownames(samples1) %in% colnames(counts1))
          a2 <- mean(samples1[, 1] %in% colnames(counts1))

          if (a2 > a1 && NCOL(samples1) > 1) {
            message("[UploadModuleServer] getting sample names from first column\n")
            rownames(samples1) <- samples1[, 1]
            uploaded[["samples.csv"]] <- samples1[, -1, drop = FALSE]
          }

          if (FILES_check$PASS == FALSE) {
            status["samples.csv"] <- "Error, please check your samples files."
            status["counts.csv"] <- "Error, please check your counts files."
            uploaded[["counts.csv"]] <- NULL
            uploaded[["samples.csv"]] <- NULL
          }

          if (FILES_check$PASS == TRUE) {
            status["samples.csv"] <- "OK"
            status["counts.csv"] <- "OK"
          }
        }
      }

      if (status["contrasts.csv"] == "OK" && status["samples.csv"] == "OK") {
        FILES_check <- playbase::pgx.crosscheckINPUT(
          SAMPLES = uploaded[["samples.csv"]],
          CONTRASTS = uploaded[["contrasts.csv"]]
        )
        uploaded[["checklist"]][["samples_contrasts"]] <- FILES_check$check

        uploaded[["samples.csv"]] <- FILES_check$SAMPLES
        uploaded[["contrasts.csv"]] <- FILES_check$CONTRASTS


        if (FILES_check$PASS == FALSE) {
          status["samples.csv"] <- "Error, please check your samples files."
          status["contrasts.csv"] <- "Error, please check your contrasts files."
          uploaded[["samples.csv"]] <- NULL
          uploaded[["contrasts.csv"]] <- NULL
        }
      }

      MAXSAMPLES <- 25
      MAXCONTRASTS <- 5
      MAXSAMPLES <- as.integer(auth$options$MAX_SAMPLES)
      MAXCONTRASTS <- as.integer(auth$options$MAX_COMPARISONS)

      ## check files: maximum contrasts allowed
      if (status["contrasts.csv"] == "OK") {
        if (ncol(uploaded[["contrasts.csv"]]) > MAXCONTRASTS) {
          status["contrasts.csv"] <- paste("ERROR: max", MAXCONTRASTS, "contrasts allowed")
        }
      }

      ## check files: maximum samples allowed
      if (status["counts.csv"] == "OK") {
        if (ncol(uploaded[["counts.csv"]]) > MAXSAMPLES) {
          status["counts.csv"] <- paste("ERROR: max", MAXSAMPLES, " samples allowed")
        }
      }
      if (status["samples.csv"] == "OK") {
        if (nrow(uploaded[["samples.csv"]]) > MAXSAMPLES) {
          status["samples.csv"] <- paste("ERROR: max", MAXSAMPLES, "samples allowed")
        }
      }

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
      return(df)
    })

    upload_module_preview_server("upload_preview", uploaded)

    output$downloadExampleData <- shiny::downloadHandler(
      filename = "exampledata.zip",
      content = function(file) {
        # save samples, counts and contrasts locally
        samples <- playbase::SAMPLES
        counts <- playbase::COUNTS
        contrasts <- playbase::CONTRASTS

        # Write each data.frame to a temporary CSV file
        tempdir <- tempdir()
        samples_csv <- file.path(tempdir, "samples.csv")
        counts_csv <- file.path(tempdir, "counts.csv")
        contrasts_csv <- file.path(tempdir, "contrasts.csv")
        write.csv(samples, samples_csv, row.names = TRUE)
        write.csv(counts, counts_csv, row.names = TRUE)
        write.csv(contrasts, contrasts_csv, row.names = TRUE)

        # Create a zip file containing the CSV files
        zipfile <- file.path(tempdir,"data.zip")
        zip(zipfile, files = c(samples_csv, counts_csv, contrasts_csv), flags = "-r9Xj")

        # clean up
        withr::defer(unlink(tempdir, recursive = TRUE), env = globalenv())
        
        # Return the zip file
        file.copy(zipfile, file)
      }
    )

    output$upload_info <- shiny::renderUI({
      
      upload_info <- glue::glue("<div><h4>Organism: {input$selected_organism}<h4></div><br><h4>How to upload your files:</h4><p>Please prepare the data files in CSV format as shown in the example data. The file format must be comma-separated-values (.CSV). Be sure the dimensions, rownames and column names match for all files. You can upload a maximum of <u>LIMITS</u>. <a target='_blank' href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep.html'>Click here to read more about data preparation.</a>.</p>")

      #DLlink <- shiny::downloadLink(ns("downloadExampleData"), "exampledata.zip")
      #upload_info <- sub("EXAMPLEZIP", upload_info)

      limits.text <- paste(
        auth$options$MAX_DATASETS, "datasets (with each up to",
        auth$options$MAX_SAMPLES, "samples and",
        auth$options$MAX_COMPARISONS, "comparisons)"
      )
      upload_info <- sub("LIMITS", limits.text, upload_info)
      shiny::HTML(upload_info)
    })

    ## =====================================================================
    ## ========================= SUBMODULES/SERVERS ========================
    ## =====================================================================

    correctedX <- upload_module_batchcorrect_server(
      id = "batchcorrect",
      X = shiny::reactive(uploaded$counts.csv),
      is.count = TRUE,
      pheno = shiny::reactive(uploaded$samples.csv),
      height = height
    )

    corrected_counts <- shiny::reactive({
      counts <- NULL
      if (input$advanced_mode) {
        out <- correctedX()
        counts <- pmax(2**out$X - 1, 0)
      } else {
        counts <- uploaded$counts.csv
      }
      counts
    })

    modified_ct <- upload_module_makecontrast_server(
      id = "makecontrast",
      phenoRT = shiny::reactive(uploaded$samples.csv),
      contrRT = shiny::reactive(uploaded$contrasts.csv),
      countsRT = corrected_counts,
      height = height
    )

    shiny::observeEvent(modified_ct(), {
      ## Monitor for changes in the contrast matrix and if
      ## so replace the uploaded reactive values.
      modct <- modified_ct()
      uploaded$contrasts.csv <- modct$contr
      if (!is.null(raw_dir()) && dir.exists(raw_dir())) {
        write.csv(modct$contr, file.path(raw_dir(), "user_contrasts.csv"), row.names = TRUE)
      }
    })

    upload_ok <- shiny::reactive({
      check <- checkTables()
      all(check[, "status"] == "OK")
      all(grepl("ERROR", check[, "status"]) == FALSE)
    })

    batch_vectors <- shiny::reactive({
      correctedX()$B
    })

    computed_pgx <- upload_module_computepgx_server(
      id = "compute",
      countsRT = corrected_counts,
      samplesRT = shiny::reactive(uploaded$samples.csv),
      contrastsRT = shiny::reactive(uploaded$contrasts.csv),
      raw_dir = raw_dir,
      batchRT = batch_vectors,
      metaRT = shiny::reactive(uploaded$meta),
      selected_organism = shiny::reactive(input$selected_organism),
      enable_button = upload_ok,
      alertready = FALSE,
      lib.dir = FILES,
      auth = auth,
      create_raw_dir = create_raw_dir,
      height = height,
      recompute_info
    )

    uploaded_pgx <- shiny::reactive({
      if (!is.null(uploaded$pgx)) {
        pgx <- uploaded$pgx
        uploaded_method <<- "uploaded"
      } else {
        pgx <- computed_pgx()
        uploaded_method <<- "computed"
      }
      return(pgx)
    })

    ## =====================================================================
    ## ===================== PLOTS AND TABLES ==============================
    ## =====================================================================

    upload_plot_countstats_server(
      "countStats",
      checkTables,
      uploaded
    )

    upload_plot_phenostats_server(
      "phenoStats",
      checkTables,
      uploaded
    )

    upload_plot_contraststats_server(
      "contrastStats",
      checkTables,
      uploaded
    )

    buttonInput <- function(FUN, len, id, ...) {
      inputs <- character(len)
      for (i in seq_len(len)) {
        inputs[i] <- as.character(FUN(paste0(id, i), ...))
      }
      inputs
    }

    output$checkTablesOutput <- DT::renderDataTable({
      ## Render the upload status table
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
      phenoRT = shiny::reactive(uploaded$samples.csv),
      countsRT = corrected_counts,
      sel.conditions = sel.conditions,
      watermark = WATERMARK
    )

    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    # board does not return anything
  })
}
