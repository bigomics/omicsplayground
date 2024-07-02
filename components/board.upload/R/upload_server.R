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
    # this directory is used to save pgx files, logs, inputs, etc..
    raw_dir <- reactiveVal(NULL)
    upload_organism <- reactiveVal(NULL)
    upload_name <- reactiveVal(NULL)
    upload_description <- reactiveVal(NULL)
    upload_datatype <- reactiveVal(NULL)
    upload_gset_methods <- reactiveVal(NULL)
    upload_gx_methods <- reactiveVal(NULL)
    process_counter <- reactiveVal(0)
    show_comparison_builder <- shiny::reactiveVal(TRUE)
    selected_contrast_input <- shiny::reactiveVal(TRUE)
    reset_upload_text_input <- shiny::reactiveVal(0)

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

    shiny::observeEvent(input$upload_files_btn,
      {
        shinyjs::click(id = "upload_files")
      },
      ignoreNULL = TRUE
    )


    shiny::observeEvent(
      {
        list(uploaded, input$tabs)
      },
      {
        dbg("[upload_server:observeEvent(names.uploaded)] names(uploaded) = ", names(uploaded))
        dbg("[upload_server:observeEvent(names.uploaded)] input$tabs = ", input$tabs)
        if (is.null(uploaded$counts.csv) && is.null(uploaded$samples.csv)) {
          dbg("[upload_server:observeEvent(names.uploaded)] *** UPLOADED EMPTY *** ")
        }
      },
      ignoreNULL = FALSE
    )


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
      pgxdir <- auth$user_dir
      fn <- file.path(pgxdir, pgxfile)
      fn <- iconv(fn, from = "", to = "ASCII//TRANSLIT")
      ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ## switch 'pgx' as standard name. Actually saving as RDS
      ## would have been better...
      dbg("[UploadBoard:observe:uploaded_pgx] saving pgx as = ", fn)
      playbase::pgx.save(new_pgx, file = fn)

      shiny::withProgress(
        message = "Scanning dataset library...",
        value = 0.33,
        {
          playbase::pgxinfo.updateDatasetFolder(
            pgxdir,
            new.pgx = pgxfile,
            update.sigdb = FALSE
          )
        }
      )

      ## trigger reload of pgx table
      reload_pgxdir(reload_pgxdir() + 1)

      beepr::beep(10) ## short beep

      load_my_dataset <- function() {
        if (input$confirmload) {
          load_uploaded_data(pgxfile)
        }
      }

      # reset new_upload to 0, so upload will not trigger when computation is done
      new_upload(0)

      if (uploaded_method == "computed") {
        shinyalert::shinyalert(
          title = paste("Your dataset is ready!"),
          text = paste("Your dataset", new_pgx$name, "is ready for visualization. Happy discoveries!"),
          confirmButtonText = "Analyze my new data!",
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


    ## =====================================================================
    ## ======================= UI OBSERVERS ================================
    ## =====================================================================

    shiny::observeEvent(input$expert_mode, {
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


    last_hash <- 1234

    create_raw_dir <- function(auth) {
      auth_id <- ifelse(!auth$email %in% c("", NA), auth$email, auth$username)
      prefix <- paste0("raw_", auth_id, "_")
      raw_dir <- tempfile(pattern = prefix, tmpdir = file.path(PGX.DIR, "USER_INPUT"))
      dir.create(raw_dir, recursive = TRUE)
      raw_dir
    }

    ## =====================================================================
    ## ===================== checkTables ===================================
    ## =====================================================================

    ## --------------------------------------------------------
    ## Check COUNTS matrix
    ## --------------------------------------------------------

    checked_for_log <- reactiveVal(FALSE)

    checked_counts <- shiny::eventReactive(
      {
        list(uploaded$counts.csv)
      },
      {
        ## get uploaded counts

        checked <- NULL
        df0 <- uploaded$counts.csv
        if (is.null(df0)) {
          return(list(status = "Missing counts.csv", matrix = NULL))
        }

        ## --------------------------------------------------------
        ## Single matrix counts check
        ## --------------------------------------------------------
        browser()
        res <- playbase::pgx.checkINPUT(df0, "COUNTS", organism = upload_organism())
        write_check_output(res$checks, "COUNTS", raw_dir())

        # check if error 29 exists (log2 transform detected), give action to user revert to intensities or skip correction
        #   counts_log_correction <- function(isConfirmed) {
        #   if (isConfirmed) {
        #     res$df <- 2**res$df
        #     if(min(res$df,na.rm=TRUE) > 0) res$df <- res$df - 1
        #     checked <<- res$df
        #     checked_for_log(TRUE)
        #   }
        # }

        if ("e29" %in% names(res$checks)) {
          # TODO find a way to give user option to correct or not (shinyalert callback breaks req and eventReactive)
          # shinyalert::shinyalert(
          #   title = paste("log-transformed counts?"),
          #   text = paste("Your counts data seems to be log-transformed. Would you like to revert to intensities?"),
          #   confirmButtonText = "Convert to intensities.",
          #   showCancelButton = TRUE,
          #   cancelButtonText = "My counts are not log transformed.",
          #   inputId = "logCorrectCounts",
          #   closeOnEsc = FALSE,
          #   immediate = TRUE,
          #   callbackR = counts_log_correction
          # )

          # inform user that we are applying the correction
          shinyalert::shinyalert(
            title = "Converting log-transformed counts",
            text = "Your counts data seems to be log-transformed. We are converting it to intensities.",
            type = "info"
          )


          # this should be run only when user confirms to convert to intensities in shinyalert (counts_log_correction function)
          res$df <- 2**res$df
          if (min(res$df, na.rm = TRUE) > 0) res$df <- res$df - 1
          checked <- res$df
          checked_for_log(TRUE)
        } else {
          # no correction needed
          checked <- res$df
          checked_for_log(TRUE)
        }

        # TODO if you use the req, eventReactive will return at shiny alert execution, and data will not be corrected
        # req(checked_for_log(), !is.null(checked))

        checklist[["counts.csv"]]$checks <- res$checks

        if (res$PASS) {
          status <- "OK"
        } else {
          checked <- NULL
          status <- "ERROR: incorrect counts matrix"
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
    checked_samples_counts <- shiny::eventReactive(
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

        write_check_output(res$checks, "SAMPLES", raw_dir())
        # store check and data regardless of it errors
        checklist[["samples.csv"]]$checks <- res$checks
        checked <- res$df
        if (res$PASS) {
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

          write_check_output(cross_check$checks, "SAMPLES_COUNTS", raw_dir())

          checklist[["samples_counts"]]$checks <- cross_check$checks

          # initialize results
          res_samples <- NULL
          res_counts <- NULL

          if (cross_check$PASS) {
            res_samples <- cross_check$SAMPLES
            res_counts <- cross_check$COUNTS
            status <- "OK"
          } else {
            checked <- NULL
            status <- "ERROR: samples matrix does not match counts"
          }
        }

        if (is.null(checked)) {
          uploaded[["last_uploaded"]] <<- setdiff(uploaded[["last_uploaded"]], "samples.csv")
          ## uploaded[["samples.csv"]] <<- NULL
          uploaded[["contrasts.csv"]] <<- NULL
        }

        list(status = status, SAMPLES = res_samples, COUNTS = res_counts)
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

        ## --------- Single matrix counts check----------
        res <- playbase::pgx.checkINPUT(df0, "CONTRASTS")
        # store check and data regardless of it errors
        checklist[["contrasts.csv"]]$checks <- res$checks
        write_check_output(res$checks, "CONTRASTS", raw_dir())
        checked <- res$df
        if (res$PASS) {
          status <- "OK"
        } else {
          checked <- NULL
          status <- "ERROR: invalid contrast. please check your input file."
        }

        ## Check if samples.csv exists before uploading contrast.csv
        cc <- checked_samples_counts()

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
        if (!is.null(checked) && !is.null(cc$SAMPLES)) {
          cross_check <- playbase::pgx.crosscheckINPUT(
            SAMPLES = cc$SAMPLES,
            CONTRASTS = checked
          )

          write_check_output(cross_check$checks, "SAMPLES_CONTRASTS", raw_dir())
          checklist[["samples_contrasts"]]$checks <- cross_check$checks
          checked <- res$df
          if (cross_check$PASS) {
            # checked <- res$df
            status <- "OK"
          } else {
            checked <- NULL
            status <- "ERROR: contrasts do not match samples."
          }
        }

        if (!is.null(checked)) {
          checked <- playbase::contrasts.convertToLabelMatrix(
            contrasts = checked, samples = cc$SAMPLES
          )
        }
        if (is.null(checked)) {
          uploaded[["last_uploaded"]] <<- setdiff(uploaded[["last_uploaded"]], "contrasts.csv")
        }

        list(status = status, matrix = checked)
      }
    )

    output$input_recap <- renderUI({
      shiny::fluidRow(
        shiny::column(3, tags$h3(shiny::HTML(paste("<b>Organism:</b> ", upload_organism(), sep = "<br>")))),
        shiny::column(3, tags$h3(shiny::HTML(paste("<b>Name:</b> ", upload_name(), sep = "<br>")))),
        shiny::column(3, tags$h3(shiny::HTML(paste("<b>Description:</b> ", upload_description(), sep = "<br>")))),
        shiny::column(3, tags$h3(shiny::HTML(paste("<b>Data type:</b> ", upload_datatype(), sep = "<br>"))))
      )
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
      upload_info <- sub("_LIMITS_", limits.text, upload_info, fixed = TRUE)
      shiny::HTML(upload_info)
    })

    ## =====================================================================
    ## ========================= SUBMODULES/SERVERS ========================
    ## =====================================================================

    shiny::observeEvent(modified_ct(), {
      ## Monitor for changes in the contrast matrix and if
      ## so replace the uploaded reactive values.
      modct <- modified_ct()
      if (!is.null(raw_dir()) && dir.exists(raw_dir())) {
        write.csv(modct, file.path(raw_dir(), "user_contrasts.csv"), row.names = TRUE)
      }
    })

    corrected1 <- upload_module_outliers_server(
      id = "checkqc",
      r_X = shiny::reactive(checked_samples_counts()$COUNTS),
      r_samples = shiny::reactive(checked_samples_counts()$SAMPLES),
      r_contrasts = modified_ct,
      is.count = TRUE,
      height = height
    )

    correctedX <- upload_module_batchcorrect_server(
      id = "batchcorrect",
      r_X = shiny::reactive(checked_samples_counts()$COUNTS),
      r_samples = shiny::reactive(checked_samples_counts()$SAMPLES),
      r_contrasts = modified_ct,
      r_results = modified_ct,
      is.count = TRUE
    )


    computed_pgx <- upload_module_computepgx_server(
      id = "compute",
      countsRT = shiny::reactive(checked_samples_counts()$COUNTS), # TODO add return from new-bc module: corrected1$correctedCounts,
      samplesRT = shiny::reactive(checked_samples_counts()$SAMPLES),
      contrastsRT = modified_ct,
      raw_dir = raw_dir,
      metaRT = shiny::reactive(uploaded$meta),
      upload_organism = upload_organism,
      alertready = FALSE,
      lib.dir = FILES,
      auth = auth,
      create_raw_dir = create_raw_dir,
      height = "100%",
      recompute_info = recompute_info,
      inactivityCounter = inactivityCounter,
      upload_wizard = reactive(input$upload_wizard),
      upload_datatype = upload_datatype,
      upload_name = upload_name,
      upload_description = upload_description,
      upload_gx_methods = upload_gx_methods,
      upload_gset_methods = upload_gset_methods,
      process_counter = process_counter,
      reset_upload_text_input = reset_upload_text_input
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

    # warn user when locked button is clicked (UX)
    observeEvent(
      input$upload_wizard_locked,
      {
        # summaryze all error logs
        summary_checks <- list(
          checklist$samples.csv$checks,
          checklist$counts.csv$checks,
          checklist$contrasts.csv$checks,
          checklist$samples_counts$checks,
          checklist$samples_contrasts$checks
        )

        summary_check_content <- length(unlist(summary_checks, recursive = FALSE))

        result_alert <- NULL

        if (summary_check_content > 0) {
          # chekc which checks have error results
          find_content <- !sapply(
            summary_checks,
            function(x) is.null(x) || length(x) == 0
          )

          summary_checks <- summary_checks[find_content]


          # get the names of each list within summary checks
          get_all_codes <- sapply(summary_checks, function(x) names(x))


          # check if any any code is error code
          error_list <- playbase::PGX_CHECKS
          error_list <- error_list[error_list$warning_type == "error", ]

          if (any(get_all_codes %in% error_list$error)) {
            # which warning is error
            which_error <- which(get_all_codes %in% error_list$error)

            result_alert <- check_to_html(
              unlist(summary_checks[which_error], recursive = FALSE),
              pass_msg = "All counts checks passed",
              null_msg = "Fix any errors with your input first."
            )

            # return shiny alert with result_alert
            shinyalert::shinyalert(
              title = "Please review your input:",
              text = result_alert,
              type = "error",
              html = TRUE
            )
          }
        }

        if (is.null(result_alert)) {
          if (input$upload_wizard == "step_samples") {
            shinyalert::shinyalert(
              title = "Upload your samples!",
              text = "Please finish the current step before proceeding.",
              type = "warning"
            )
          } else if (input$upload_wizard == "step_counts") {
            shinyalert::shinyalert(
              title = "Upload your counts!",
              text = "Please finish the current step before proceeding.",
              type = "warning"
            )
          } else if (input$upload_wizard == "step_comparisons") {
            shinyalert::shinyalert(
              title = "Create at least one comparison!",
              text = "Upload or build your comparisons.",
              type = "warning"
            )
          } else if (input$upload_wizard == "step_compute") {
            shinyalert::shinyalert(
              title = "Start the computation!",
              text = "Please finish the current step before proceeding.",
              type = "warning"
            )
          }
        }
      }
    )

    # wizard lock/unlock logic

    # lock/unlock wizard for counts.csv
    observeEvent(
      list(uploaded$counts.csv, checked_counts, input$upload_wizard),
      {
        req(input$upload_wizard == "step_counts")

        if (is.null(checked_counts()$status) || checked_counts()$status != "OK") {
          wizardR::lock("upload_wizard")
        } else if (!is.null(checked_counts()$status) && checked_counts()$status == "OK") {
          wizardR::unlock("upload_wizard")
        }
      }
    )

    # lock/unlock wizard for samples.csv
    observeEvent(
      list(uploaded$samples.csv, checked_samples_counts(), input$upload_wizard),
      {
        req(input$upload_wizard == "step_samples")
        if (is.null(checked_samples_counts()$status) || checked_samples_counts()$status != "OK") {
          wizardR::lock("upload_wizard")
        } else if (!is.null(checked_samples_counts()$status) && checked_samples_counts()$status == "OK") {
          wizardR::unlock("upload_wizard")
        }
      }
    )

    # lock wizard at Comparison step
    observeEvent(
      list(input$upload_wizard, modified_ct()),
      {
        req(input$upload_wizard == "step_comparisons")
        if (is.null(modified_ct()) || ncol(modified_ct()) == 0 || is.null(checked_contrasts()) || is.null(checked_samples_counts()) || is.null(checked_counts())) {
          wizardR::lock("upload_wizard")
        } else {
          wizardR::unlock("upload_wizard")
        }
      }
    )

    # lock wizard it compute step
    observeEvent(
      list(
        input$upload_wizard,
        upload_name(),
        upload_datatype(),
        upload_description(),
        upload_organism(),
        upload_gset_methods(),
        upload_gx_methods()
      ),
      {
        req(input$upload_wizard == "step_compute")

        pgx_files <- playbase::pgxinfo.read(auth$user_dir, file = "datasets-info.csv")

        if (!is.null(upload_name()) && upload_name() %in% pgx_files$dataset) {
          shinyalert::shinyalert(
            title = "Invalid name",
            text = "This dataset name already exists.",
            type = "error"
          )
          upload_name(NULL)
        }

        if (is.null(upload_gx_methods())) {
          shinyalert::shinyalert(
            title = "ERROR",
            text = "You must select at least one gene test method",
            type = "error"
          )
        }
        if (is.null(upload_gset_methods())) {
          shinyalert::shinyalert(
            title = "ERROR",
            text = "You must select at least one geneset (enrichment) test method",
            type = "error"
          )
        }

        if (!is.null(upload_name()) && !isValidFileName(upload_name())) {
          message("[ComputePgxServer:input$compute] WARNING:: Invalid name")
          shinyalert::shinyalert(
            title = "Invalid name",
            text = "Please remove any slashes (/) from the name",
            type = "error"
          )
          upload_name(NULL)
        }

        if (
          is.null(upload_name()) ||
            upload_name() == "" ||
            upload_description() == "" ||
            is.null(upload_description()) ||
            is.null(upload_gx_methods()) ||
            is.null(upload_gset_methods())
        ) {
          wizardR::lock("upload_wizard")
        } else {
          wizardR::unlock("upload_wizard")
        }
      }
    )


    ## =====================================================================
    ## ===================== PLOTS AND TABLES ==============================
    ## =====================================================================

    upload_table_preview_counts_server(
      "counts_preview",
      uploaded,
      checklist = checklist,
      scrollY = "calc(50vh - 140px)",
      width = c("auto", "100%"),
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
      width = c("auto", "100%"),
      height = c("100%", TABLE_HEIGHT_MODAL),
      title = "Uploaded Samples",
      info.text = "This is the uploaded samples data.",
      caption = "This is the uploaded samples data."
    )

    modified_ct <- upload_table_preview_contrasts_server(
      "contrasts_preview",
      uploaded,
      checklist,
      scrollY = "calc(50vh - 140px)",
      height = c("100%", TABLE_HEIGHT_MODAL),
      width = c("auto", "100%"),
      title = "Uploaded Contrasts",
      info.text = "This is the uploaded comparison data.",
      caption = "This is the uploaded comparison data.",
      checked_samples = checked_samples_counts,
      checked_counts = checked_samples_counts,
      checked_contrasts = checked_contrasts,
      show_comparison_builder = show_comparison_builder,
      selected_contrast_input = selected_contrast_input,
      upload_wizard = shiny::reactive(input$upload_wizard)
    )

    # observe show_modal and start modal
    shiny::observeEvent(
      list(new_upload()),
      {
        shiny::req(auth$options)
        enable_upload <- auth$options$ENABLE_UPLOAD

        isolate({
          lapply(names(uploaded), function(i) uploaded[[i]] <- NULL)
          lapply(names(checklist), function(i) checklist[[i]] <- NULL)
          upload_datatype(NULL)
          upload_name(NULL)
          upload_description(NULL)
          # upload_organism(NULL)
          show_comparison_builder(TRUE)
          selected_contrast_input(FALSE)
        })

        reset_upload_text_input(reset_upload_text_input() + 1)

        wizardR::reset("upload_wizard")

        # skip upload trigger at first startup
        if (new_upload() == 0) {
          return(NULL)
        }

        if (enable_upload) {
          MAX_DS_PROCESS <- 1
          if (process_counter() < MAX_DS_PROCESS) {
            wizardR::wizard_show(ns("upload_wizard"))
            if (!is.null(recompute_pgx())) {
              pgx <- recompute_pgx()
              uploaded$samples.csv <- pgx$samples
              uploaded$contrasts.csv <- pgx$contrast
              uploaded$counts.csv <- pgx$counts
              recompute_info(
                list(
                  "name" = pgx$name,
                  "description" = pgx$description
                )
              )
            }

            # if recomputing pgx, add data to wizard
          } else {
            shinyalert::shinyalert(
              title = "Computation in progress",
              text = "Sorry, only one computation is allowed at a time. Please wait for the current computation to finish.",
              type = "warning",
              closeOnClickOutside = FALSE
            )
          }
        } else {
          shinyalert::shinyalert(
            title = "Upload disabled",
            text = "Sorry, upload of new data is disabled for this account.",
            type = "warning",
            closeOnClickOutside = FALSE
          )
        }
      }
    )

    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
  })
}
