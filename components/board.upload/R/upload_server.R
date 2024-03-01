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
                        recompute_info,
                        inactivityCounter,
                        new_upload) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns ## NAMESPACE

    # Some 'global' reactive variables used in this file
    uploaded <- shiny::reactiveValues()
    ##    checked   <- shiny::reactiveValues()
    checklist <- shiny::reactiveValues()

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

    shiny::observeEvent( input$upload_files_btn, {
      shinyjs::click(id = "upload_files")
    }, ignoreNULL = TRUE )

    
    shiny::observeEvent({
      list( uploaded, input$tabs )
    }, {
      dbg("[upload_server:observeEvent(names.uploaded)] names(uploaded) = ",names(uploaded))
      dbg("[upload_server:observeEvent(names.uploaded)] input$tabs = ",input$tabs)      
      if( is.null(uploaded$counts.csv) && is.null(uploaded$samples.csv) ) {
        dbg("[upload_server:observeEvent(names.uploaded)] *** UPLOADED EMPTY *** ")
      }
    }, ignoreNULL = FALSE)
    
    
    shiny::observeEvent(uploaded_pgx(), {
      new_pgx <- uploaded_pgx()

      dbg("[upload_server:observeEvent(uploaded_pgx()] names(new_pgx) = ", names(new_pgx))      
      dbg("[upload_server:observeEvent(uploaded_pgx()] dim(new_pgx$X) = ", dim(new_pgx$X))
      dbg("[upload_server:observeEvent(uploaded_pgx()] new_pgx$name = ", new_pgx$name)
      dbg("[upload_server:observeEvent(uploaded_pgx()] uploaded_method = ", uploaded_method)
      
      ## NEED RETHINK: if "uploaded" we unneccessarily saving the pgx
      ## object again.  We should skip saving and pass the filename to
      ## pgxfile to be sure the filename is correct.

      ## new_pgx <- playbase::pgx.initialize(new_pgx)  ## already done later
      ## -------------- save PGX file/object ---------------
      # Old pgx does not have name slot, overwrite it with file name
      if (is.null(new_pgx$name)) {
        new_pgx$name <- sub("[.]pgx$", "", input$upload_files$name)
      }
      pgxfile <- sub("[.]pgx$", "", new_pgx$name)
      pgxfile <- gsub("^[./-]*", "", pgxfile) ## prevent going to parent folder
      pgxfile <- paste0(gsub("[ \\/]", "_", pgxfile), ".pgx")
      pgxdir  <- auth$user_dir
      fn <- file.path(pgxdir, pgxfile)
      fn <- iconv(fn, from = "", to = "ASCII//TRANSLIT")
      ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ## switch 'pgx' as standard name. Actually saving as RDS
      ## would have been better...
      dbg("[UploadBoard:observe:uploaded_pgx] saving pgx as = ", fn)
      playbase::pgx.save(new_pgx, file = fn)

      shiny::withProgress(
        message = "Scanning dataset library...", value = 0.33, {
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

      ## clean up reactiveValues
      isolate({
        lapply(names(uploaded), function(i) uploaded[[i]] <- NULL)
        lapply(names(checklist), function(i) checklist[[i]] <- NULL)
      })

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

    # hide upload tab at server start
##    shinyjs::runjs('document.querySelector(\'[data-value="Upload"]\').style.display = "none";')

    ## Hide/show tabpanels upon available data like a wizard dialog
    shiny::observe({

      return(NULL)  ### TEMPORARY FOR DEVELOPMENT
      
      has.counts <- !is.null(checked_counts()$matrix)
      has.samples <- !is.null(checked_samples()$matrix)
      has.contrasts <- !is.null(checked_contrasts()$matrix)
      # check that modified contrast is not NULL and has col dim >0
      has.contrasts <- !is.null(modified_ct()) && ncol(modified_ct()) > 0

      need2 <- has.counts && has.samples
      need3 <- need2 && has.contrasts

      if (need3) {
        # show compute if contrast is done
        shiny::showTab("tabs", "Compute")
        shiny::showTab("tabs", "Comparisons")
        shiny::showTab("tabs", "QC/BC")
        if (input$expert_mode) {
          shiny::showTab("tabs", "BatchEffects")
        }
      } else if (need2) {
        shiny::hideTab("tabs", "Compute")
        shiny::showTab("tabs", "Comparisons")
        shiny::showTab("tabs", "QC/BC")
        if (input$expert_mode) {
          shiny::showTab("tabs", "BatchEffects")
        }
      } else {
        shiny::hideTab("tabs", "Compute")
        shiny::hideTab("tabs", "BatchEffects")
        shiny::hideTab("tabs", "QC/BC")
        shiny::hideTab("tabs", "Comparisons")
      }
    })

    ## =====================================================================
    ## ======================= UI OBSERVERS ================================
    ## =====================================================================

    shiny::observeEvent( input$expert_mode, {
      if (input$expert_mode) {
        shiny::showTab("tabs", "BatchEffects")
      } else {
        shiny::hideTab("tabs", "BatchEffects")
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
    last_hash <- 1234

    create_raw_dir <- function(auth) {
      auth_id <- ifelse(!auth$email %in% c("", NA), auth$email, auth$username)
      prefix <- paste0("raw_", auth_id, "_")
      raw_dir <- tempfile(pattern = prefix, tmpdir = file.path(PGX.DIR, "USER_INPUT"))
      dir.create(raw_dir, recursive = TRUE)
      dbg("[UploadBoard:raw_dir<-eventReactive] creating raw_dir", raw_dir)
      raw_dir
    }

    # In case the user is reanalysing the data, get the info from pgx
    observeEvent(recompute_pgx(), {
      pgx <- recompute_pgx()
      uploaded$samples.csv <- pgx$samples
      uploaded$contrasts.csv <- pgx$contrast
      uploaded$counts.csv <- pgx$counts
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
      # show tab Upload
      shinyjs::runjs('document.querySelector(\'[data-value="Upload"]\').style.display = "";')
      # check on upload tab
      shinyjs::runjs('document.querySelector("a[data-value=\'Upload\']").click();')

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
        # this was re-done in multi-species, it will be much better. Temporary solution for legacy code. MMM
        checklist[["contrasts.csv"]]$checks <- list()
        checklist[["samples.csv"]]$checks <- list()
        checklist[["counts.csv"]]$checks <- list()
        checklist[["contrasts.csv"]]$PASS <- TRUE
        checklist[["samples_counts"]] <- NULL
        checklist[["samples_contrasts"]] <- NULL
        uploaded[["last_uploaded"]] <- c("counts.csv", "samples.csv", "contrasts.csv")

        shinyalert::shinyalert(
          title = "Example dataset",
          text = 'This example dataset is a time course experiment measuring the protein abundances in T cells comparing activated vs. resting cells at different time points (Geiger et al., Cell 2016)',
          html = TRUE,
        )
        
      } else {
        ## clear files
        lapply(names(uploaded), function(i) uploaded[[i]] <- NULL)
        lapply(names(checklist), function(i) checklist[[i]] <- NULL)
      }
    })


    ## =====================================================================
    ## ===================== checkTables ===================================
    ## =====================================================================

    ## --------------------------------------------------------
    ## Check COUNTS matrix
    ## --------------------------------------------------------
    checked_counts <- shiny::eventReactive(
      {
        list(uploaded$counts.csv)
      },
      {
        ## get uploaded counts
        df0 <- uploaded$counts.csv
        if (is.null(df0)) {
          return(list(status = "Missing counts.csv", matrix = NULL))
        }

        ## --------------------------------------------------------
        ## Single matrix counts check
        ## --------------------------------------------------------
        res <- playbase::pgx.checkINPUT(df0, "COUNTS")
        # store check and data regardless of it errors
        checklist[["counts.csv"]]$checks <- res$checks
        if (res$PASS) {
          checked <- res$df
          status <- "OK"
        } else {
          checked <- NULL
          status <- "ERROR: incorrect counts matrix"
          uploaded$counts.csv <- NULL
        }

        ## --------------------------------------------------------
        ## check files: maximum samples allowed
        ## --------------------------------------------------------
        MAXSAMPLES <- as.integer(auth$options$MAX_SAMPLES)
        if (!is.null(checked)) {
          if (ncol(checked) > MAXSAMPLES) {
            status <- paste("ERROR: max", MAXSAMPLES, " samples allowed")
            checked <- NULL
            # remove only counts.csv from last_uploaded
            uploaded[["last_uploaded"]] <- setdiff(uploaded[["last_uploaded"]], "counts.csv")
            ## uploaded[["counts.csv"]] <- NULL
            # pop up telling user max contrasts reached
            shinyalert::shinyalert(
              title = "Maximum counts reached",
              text = paste(
                "You have reached the maximum number of counts allowed. Please",
                "upload a new COUNTS file with a maximum of", MAXSAMPLES, "samples."
              ),
              type = "error"
            )
          }
        }
        if (is.null(checked)) {
          uploaded[["last_uploaded"]] <- setdiff(uploaded[["last_uploaded"]], "counts.csv")
        }

        list(status = status, matrix = checked)
      }
    )


    ## --------------------------------------------------------
    ## Check SAMPLES matrix
    ## --------------------------------------------------------
    checked_samples <- shiny::eventReactive(
      {
        list(uploaded$counts.csv, uploaded$samples.csv)
      },
      {
        ## get uploaded counts
        df0 <- uploaded$samples.csv
        if (is.null(df0)) {
          return(list(status = "Missing samples.csv", matrix = NULL))
        }

        ## Single matrix counts check
        res <- playbase::pgx.checkINPUT(df0, "SAMPLES")
        # store check and data regardless of it errors
        checklist[["samples.csv"]]$checks <- res$checks
        if (res$PASS) {
          checked <- res$df
          status <- "OK"
        } else {
          checked <- NULL
          status <- "ERROR: incorrect samples matrix"
        }

        MAXSAMPLES <- as.integer(auth$options$MAX_SAMPLES)
        if (!is.null(checked)) {
          if (nrow(checked) > MAXSAMPLES) {
            status <- paste("ERROR: max", MAXSAMPLES, "samples allowed")
            ## uploaded[["samples.csv"]] <- NULL
            checked <- NULL
            # pop up telling user max samples reached
            shinyalert::shinyalert(
              title = "Maximum samples reached",
              text = paste("You have reached the maximum number of samples allowed. Please upload a new SAMPLES file with a maximum of", MAXSAMPLES, "samples."),
              type = "error"
            )
          }
        }

        ## -------------- cross-check with counts ------------------
        cc <- checked_counts()
        if (!is.null(checked) && !is.null(cc$matrix)) {
          cross_check <- playbase::pgx.crosscheckINPUT(
            SAMPLES = checked,
            COUNTS = cc$matrix
          )
          checklist[["samples_counts"]]$checks <- cross_check$checks

          if (cross_check$PASS) {
            checked <- res$df
            status <- "OK"
          } else {
            checked <- NULL
            status <- "ERROR: samples matrix does not match counts"
          }
        }

        if (is.null(checked)) {
          uploaded[["last_uploaded"]] <- setdiff(uploaded[["last_uploaded"]], "samples.csv")
          ## uploaded[["samples.csv"]] <- NULL
          uploaded[["contrasts.csv"]] <- NULL
        }

        list(status = status, matrix = checked)
      }
    )


    ## --------------------------------------------------------
    ## Check contrast matrix
    ## --------------------------------------------------------
    checked_contrasts <- shiny::eventReactive(
      {
        list(uploaded$contrasts.csv, uploaded$samples.csv)
      },
      {
        ## get uploaded counts
        df0 <- uploaded$contrasts.csv
        if (is.null(df0)) {
          return(list(status = "Missing contrasts.csv", matrix = NULL))
        }

        ## insanity check
        if (NCOL(df0) == 0) {
          checked <- NULL
          status <- "ERROR: no contrasts. please check your input file."
        }

        ## --------- Single matrix counts check----------
        res <- playbase::pgx.checkINPUT(df0, "CONTRASTS")
        # store check and data regardless of it errors
        checklist[["contrasts.csv"]]$checks <- res$checks
        if (res$PASS) {
          checked <- res$df
          status <- "OK"
        } else {
          checked <- NULL
          status <- "ERROR: invalid contrast. please check your input file."
        }


        ## Check if samples.csv exists before uploading contrast.csv
        cc <- checked_samples()
        if (!is.null(checked) && is.null(cc$matrix)) {
          status <- "ERROR: please upload samples file first."
          checked <- NULL
          # pop up telling the user to upload samples.csv first
          shinyalert::shinyalert(
            title = "Samples.csv file missing",
            text = "Please upload the samples.csv file first.",
            type = "error"
          )
        }

        ## -------------- max contrast check ------------------
        MAXCONTRASTS <- as.integer(auth$options$MAX_COMPARISONS)
        if (!is.null(checked)) {
          if (ncol(checked) > MAXCONTRASTS) {
            status <- paste("ERROR: max", MAXCONTRASTS, "contrasts allowed")
            checked <- NULL
            # pop up telling user max contrasts reached
            shinyalert::shinyalert(
              title = "Maximum contrasts reached",
              text = paste("You have reached the maximum number of contrasts allowed. Please upload a new contrasts file with a maximum of", MAXCONTRASTS, "contrasts."),
              type = "error"
            )
          }
        }

        ## -------------- cross-check with samples ------------------
        cc <- checked_samples()
        if (!is.null(checked) && !is.null(cc$matrix)) {
          cross_check <- playbase::pgx.crosscheckINPUT(
            SAMPLES = cc$matrix,
            CONTRASTS = checked
          )
          checklist[["samples_contrasts"]]$checks <- cross_check$checks
          if (cross_check$PASS) {
            checked <- res$df
            status <- "OK"
          } else {
            checked <- NULL
            status <- "ERROR: contrasts do not match samples."
          }
        }

        if (!is.null(checked)) {
          checked <- playbase::contrasts.convertToLabelMatrix(
            contrasts = checked, samples = cc$matrix
          )
        }
        if (is.null(checked)) {
          uploaded[["last_uploaded"]] <- setdiff(uploaded[["last_uploaded"]], "contrasts.csv")
        }

        list(status = status, matrix = checked)
      }
    )


    ## --------------------------------------------------------
    ## Gather all checks in table
    ## --------------------------------------------------------

    checkTables <- shiny::reactive({
      # check if status exists
      status <- rep("please upload", 3)
      files.needed <- c("counts.csv", "samples.csv", "contrasts.csv")
      names(status) <- files.needed
      files.nrow <- rep(NA, 3)
      files.ncol <- rep(NA, 3)

      ## check if all files in uploaded
      for (i in 1:3) {
        fn <- files.needed[i]
        upfile <- uploaded[[fn]]
        if (fn %in% names(uploaded) && !is.null(upfile)) {
          status[i] <- "OK"
          files.nrow[i] <- nrow(upfile)
          files.ncol[i] <- ncol(upfile)
        }
      }
      ##      ERROR_CODES <- playbase::PGX_CHECKS
      has.pgx <- ("pgx" %in% names(uploaded))
      has.csv <- any(grepl("csv", names(uploaded)))
      if (has.pgx) has.pgx <- has.pgx && !is.null(uploaded[["pgx"]])
      if (has.pgx == TRUE) {
        ## Nothing to check. Always OK.
        status <- c("counts.csv" = "OK", "samples.csv" = "OK", "contrasts.csv" = "OK")
      } else if (!has.pgx) {
        c1 <- checked_counts()
        c2 <- checked_samples()
        c3 <- checked_contrasts()

        status <- c(
          "counts.csv" = c1$status,
          "samples.csv" = c2$status,
          "contrasts.csv" = c3$status
        )
      }

      ## --------------------------------------------------------
      ## Build summary table
      ## --------------------------------------------------------
      description <- c(
        "Count/expression file with gene on rows, samples as columns",
        "Samples file with samples on rows, phenotypes as columns",
        ## "Gene information file with genes on rows, gene info as columns.",
        "Contrast file with conditions on rows, contrasts as columns"
      )
      description <- c(
        "genes x samples",
        "samples x phenotypes",
        ## "Gene information file with genes on rows, gene info as columns.",
        "conditions x comparisons"
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
        zipfile <- file.path(tempdir, "data.zip")
        zip(zipfile, files = c(samples_csv, counts_csv, contrasts_csv), flags = "-r9Xj")

        # clean up
        withr::defer(unlink(tempdir, recursive = TRUE), env = globalenv())

        # Return the zip file
        file.copy(zipfile, file)
      }
    )

    output$upload_info <- shiny::renderUI({
      upload_info <- "Please prepare the data files in CSV format with the names 'counts.csv', 'samples.csv' and 'contrasts.csv'. Be sure the dimensions, rownames and column names match for all files. You can upload a maximum of _LIMITS_. Click <u><a target='_blank' href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep.html'>here</a></u> to read more about data preparation.</p>"
      limits.text <- paste(
        auth$options$MAX_DATASETS, "datasets (with each up to",
        auth$options$MAX_SAMPLES, "samples and",
        auth$options$MAX_COMPARISONS, "comparisons)"
      )
      upload_info <- sub("_LIMITS_", limits.text, upload_info, fixed=TRUE)
      shiny::HTML(upload_info)
    })

    ## =====================================================================
    ## ========================= SUBMODULES/SERVERS ========================
    ## =====================================================================

    modified_ct <- upload_module_makecontrast_server(
      id = "makecontrast",
      phenoRT = reactive(checked_samples()$matrix),
      contrRT = reactive(checked_contrasts()$matrix),
      countsRT = reactive(checked_counts()$matrix)
    )

    shiny::observeEvent( modified_ct(), {
      ## Monitor for changes in the contrast matrix and if
      ## so replace the uploaded reactive values.
      modct <- modified_ct()
      if (!is.null(raw_dir()) && dir.exists(raw_dir())) {
        write.csv(modct, file.path(raw_dir(), "user_contrasts.csv"), row.names = TRUE)
      }
    })
    
    corrected1 <- upload_module_outliers_server(
      id = "checkqc",
      r_X = shiny::reactive(checked_counts()$matrix),
      r_samples = shiny::reactive(checked_samples()$matrix),
      r_contrasts = modified_ct,
      is.count = TRUE,
      height = height
    )

    upload_module_batchcorrect_server(
      id = "batchcorrect",
      r_X = shiny::reactive(checked_counts()$matrix),
      r_samples = shiny::reactive(checked_samples()$matrix),
      r_contrasts = modified_ct,
      r_results = corrected1$results,
      is.count = TRUE
    )

    upload_ok <- shiny::reactive({
      check <- checkTables()
      all(check[, "status"] == "OK")
      all(grepl("ERROR", check[, "status"]) == FALSE)
    })

    computed_pgx <- upload_module_computepgx_server(
      id = "compute",
      countsRT = corrected1$correctedCounts,
      samplesRT = shiny::reactive(checked_samples()$matrix),
      contrastsRT = modified_ct,
      raw_dir = raw_dir,
      metaRT = shiny::reactive(uploaded$meta),
      selected_organism = shiny::reactive(input$selected_organism),
      enable_button = upload_ok,
      alertready = FALSE,
      lib.dir = FILES,
      auth = auth,
      create_raw_dir = create_raw_dir,
      height = "100%",
      recompute_info = recompute_info,
      inactivityCounter = inactivityCounter
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

    # wizard lock/unlock logic

    
    # lock/unlock wizard for samples.csv
    observeEvent(
      list(uploaded$counts.csv, checked_counts), {
        if (is.null(checked_counts()$status) || checked_counts()$status != "OK"){
          wizardR::lock("upload-wizard")
        } else if (!is.null(checked_counts()$status) && checked_counts()$status == "OK"){
          wizardR::unlock("upload-wizard")
        }
    })

    # lock/unlock wizard for samples.csv
    observeEvent(
      list(uploaded$samples.csv, checked_samples), {
        if (is.null(checked_samples()$status) || checked_samples()$status != "OK"){
          wizardR::lock("upload-wizard")
        } else if (!is.null(checked_samples()$status) && checked_samples()$status == "OK"){
          wizardR::unlock("upload-wizard")
        }
    })

    
    
    ## =====================================================================
    ## ===================== PLOTS AND TABLES ==============================
    ## =====================================================================

    upload_plot_countstats_server(
      "countStats",
      checkTables,
      countsRT = reactive(checked_counts()$matrix)
    )

    upload_plot_phenostats_server(
      "phenoStats",
      checkTables,
      samplesRT = reactive(checked_samples()$matrix)
    )

    upload_plot_contraststats_server(
      "contrastStats",
      checkTables,
      contrastsRT = reactive(checked_contrasts()$matrix),
      samplesRT = reactive(checked_samples()$matrix)
    )

    upload_table_preview_counts_server(
      "counts_preview",
      uploaded,
      checklist = checklist,
      scrollY = "calc(50vh - 140px)",
      width =  c("auto", "100%"),
      height = c("100%", TABLE_HEIGHT_MODAL),
      title = "Uploaded Counts",
      info.text = "This is the uploaded counts data.",
      caption = "This is the uploaded counts data."
    )

    upload_table_preview_samples_server(
      "samples_preview",
      uploaded,
      checklist = checklist,
      scrollY = "calc(50vh - 140px)",
      width =  c("auto", "100%"),
      height = c("100%", TABLE_HEIGHT_MODAL),
      title = "Uploaded Samples",
      info.text = "This is the uploaded samples data.",
      caption = "This is the uploaded samples data."
    )

    upload_table_preview_contrasts_server(
      "contrasts_preview",
      uploaded,
      checklist,
      scrollY = "calc(50vh - 140px)",
      height = c("100%", TABLE_HEIGHT_MODAL),
      width = c("auto", "100%"),
      title = "Uploaded Contrasts",
      info.text = "This is the uploaded comparison data.",
      caption = "This is the uploaded comparison data."
    )

    # observe show_modal and start modal
    shiny::observeEvent(new_upload(), {
        shiny::req(auth$options)
        enable_upload <- auth$options$ENABLE_UPLOAD

        # skip upload trigger at first startup
        if (new_upload() == 0) {
          return(NULL)
        }
        
        if (enable_upload) {
          wizardR::wizard_show(ns("upload-wizard"))
        } else {
          shinyalert::shinyalert(
            title = "Upload disabled",
            text = "Sorry, upload of new data is disabled for this account.",
            type = "warning",
            #
            closeOnClickOutside = FALSE
          )
        }
    })

    shiny::observeEvent(input$`upload-upload-wizard`,{
      print("upload-wizard fired")
    })
    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
  })
}
