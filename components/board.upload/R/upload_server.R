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
                        inactivityCounter) {
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

    shiny::observeEvent(uploaded_pgx(), {
      new_pgx <- uploaded_pgx()

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
    shinyjs::runjs('document.querySelector(\'[data-value="Upload"]\').style.display = "none";')

    ## Hide/show tabpanels upon available data like a wizard dialog
    shiny::observe({
      has.counts <- !is.null(checked_counts()$matrix)
      has.samples <- !is.null(checked_samples()$matrix)
      has.contrasts <- !is.null(checked_contrasts()$matrix)
      # check that modified contrast is not NULL and has col dim >0
      has.contrasts <- !is.null(modified_ct()$contr) && ncol(modified_ct()$contr) > 0

      need2 <- has.counts && has.samples
      need3 <- need2 && has.contrasts

      if (need3) {
        # show compute if contrast is done
        shiny::showTab("tabs", "Compute")
        shiny::showTab("tabs", "Comparisons")
        if (input$advanced_mode) {
          shiny::showTab("tabs", "BatchCorrect")
        }
      } else if (need2) {
        shiny::hideTab("tabs", "Compute")
        shiny::showTab("tabs", "Comparisons")
        if (input$advanced_mode) {
          shiny::showTab("tabs", "BatchCorrect")
        }
      } else {
        shiny::hideTab("tabs", "Compute")
        shiny::hideTab("tabs", "BatchCorrect")
        shiny::hideTab("tabs", "Comparisons")
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
    last_hash <- 1234

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
      if (class(upload_table) != "data.frame" && upload_table == "hello_example") {
        upload_table <- data.frame(
          name = c("counts.csv", "samples.csv", "contrasts.csv"),
          type = c("text/csv", "text/csv", "text/csv"),
          datapath = c("examplecounts", "examplesamples", "examplecontrasts")
        )
      }

      uploaded[["counts.csv"]] <- NULL
      uploaded[["samples.csv"]] <- NULL
      uploaded[["contrasts.csv"]] <- NULL
      uploaded[["pgx"]] <- NULL
      uploaded[["last_uploaded"]] <- NULL
      uploaded[["checklist"]] <- NULL
      checklist[["samples_counts"]] <- NULL
      checklist[["samples_contrasts"]] <- NULL

      ## read uploaded files
      pgx.uploaded <- any(grepl("[.]pgx$", upload_table$name))
      matlist <- list()

      if (pgx.uploaded) {
        ## If the user uploaded a PGX file, we extract the matrix
        ## dimensions from the given PGX/NGS object. Really?
        i <- grep("[.]pgx$", upload_table$name)
        pgxfile <- upload_table$datapath[i]
        uploaded[["pgx"]] <- local(get(load(pgxfile, verbose = 0))) ## override any name
        return(NULL)
      } else {
        ## If the user uploaded CSV files, we read in the data
        ## from the files.
        ii <- grep("csv$", input$upload_files$name)
        ii <- grep("sample|count|contrast|expression|comparison",
          upload_table$name,
          ignore.case = TRUE
        )
        if (length(ii) == 0) {
          return(NULL)
        }

        inputnames <- upload_table$name[ii]
        uploadnames <- upload_table$datapath[ii]

        ## remove any old gui_contrasts.csv
        user_ctfile <- file.path(raw_dir(), "user_contrasts.csv")
        if (file.exists(user_ctfile)) unlink(user_ctfile)

        ## ERROR_CODES <- playbase::PGX_CHECKS  ##??

        ## ---------------------------------------------------------------------
        ## This goes one-by-one over the uploaded matrices and
        ## performs a single matrix check. The cross-checks are done
        ## later.
        ## ---------------------------------------------------------------------
        if (length(uploadnames) > 0) {
          for (i in 1:length(uploadnames)) {
            # i = 1
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
              if (fn2 == "examplecounts") {
                df0 <- playbase::COUNTS
                # save a warning file telling this folder is example data
                writeLines("", file.path(raw_dir(), "EXAMPLE_DATA"))
              } else {
                df0 <- playbase::read.as_matrix(fn2)
              }

              # save input as raw file in raw_dir
              file.copy(fn2, file.path(raw_dir(), "raw_counts.csv"))

              if (IS_COUNT) {
                df <- df0
                matname <- "counts.csv"
              }

              if (IS_EXPRESSION) {
                message("[UploadModule::upload_files] converting expression to counts...")
                df <- 2**df0 - 1
                matname <- "counts.csv"
              }
            }

            if (IS_SAMPLE) {
              if (fn2 == "examplesamples") {
                df0 <- playbase::SAMPLES
              } else {
                df0 <- playbase::read.as_matrix(fn2)
              }
              # save input as raw file in raw_dir
              file.copy(fn2, file.path(raw_dir(), "raw_samples.csv"))
              df <- df0
              matname <- "samples.csv"
            }

            if (IS_CONTRAST) {
              if (fn2 == "examplecontrasts") {
                df0 <- playbase::CONTRASTS
              } else {
                df0 <- playbase::read.as_matrix(fn2)
              }
              # save input as raw file in raw_dir
              file.copy(fn2, file.path(raw_dir(), "raw_contrasts.csv"))
              df <- df0
              matname <- "contrasts.csv"
            }
            if (!is.null(matname)) {
              matlist[[matname]] <- df
            }
          }
        }
      }

      if (is.null(uploaded$counts.csv) && !"counts.csv" %in% names(matlist)) {
        shinyalert::shinyalert(
          title = "Please upload counts.csv matrix first",
          text = NULL,
          type = "info"
        )
        ## cancel upload!!
        matlist <- NULL
      }

      ## check order
      no.samples <- !("samples.csv" %in% names(matlist) || "samples.csv" %in% names(uploaded))
      no.counts <- !("counts.csv" %in% names(matlist) || "counts.csv" %in% names(uploaded))
      if ("contrasts.csv" %in% names(matlist) && (no.samples || no.counts)) {
        shinyalert::shinyalert(
          title = "Please upload counts.csv and samples.csv matrices first",
          text = NULL,
          type = "info"
        )
        ## cancel upload!!
        matlist <- NULL
      }

      ## check hash of new counts file
      if ("counts.csv" %in% names(matlist)) {
        new_hash <- rlang::hash(matlist[["counts.csv"]])
        if (new_hash != last_hash) {
          uploaded[["samples.csv"]] <- NULL
          uploaded[["contrasts.csv"]] <- NULL
          uploaded[["last_uploaded"]] <- NULL
          last_hash <<- new_hash
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
            matlist[[i]] <- type.convert(matlist[[i]], as.is = TRUE)
          }
          m1 <- names(matlist)[i]
          message("[upload_files] updating matrix ", m1)
          uploaded[[m1]] <- matlist[[i]]
        }
        uploaded[["last_uploaded"]] <- names(matlist)
      }
      message("[upload_files] done!\n")
    })


    # In case the user is reanalysing the data, get the info from pgx
    observeEvent(recompute_pgx(), {
      pgx <- recompute_pgx()
      uploaded$samples.csv <- pgx$samples
      uploaded$contrasts.csv <- pgx$contrast
      uploaded$counts.csv <- pgx$counts
      ##      corrected_counts <- pgx$counts  ## ?? IK
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
        write_check_output(res$checks, "COUNTS", raw_dir())
        # store check and data regardless of it errors
        checklist[["counts.csv"]]$checks <- res$checks
        if (res$PASS) {
          checked <- res$df
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

        write_check_output(res$checks, "SAMPLES", raw_dir())
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
          write_check_output(cross_check$checks, "SAMPLES_COUNTS", raw_dir())

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
        write_check_output(res$checks, "CONTRASTS", raw_dir())

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

          write_check_output(cross_check$checks, "SAMPLES_CONTRASTS", raw_dir())
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
      upload_info <- glue::glue("<div><h4>How to upload your files:</h4><p>Please prepare the data files in CSV format as shown in the example data. The file format must be comma-separated-values (.CSV). Be sure the dimensions, rownames and column names match for all files. You can upload a maximum of <u>LIMITS</u>. <a target='_blank' href='https://omicsplayground.readthedocs.io/en/latest/dataprep/dataprep.html'>Click here to read more about data preparation.</a>.</p>")
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

    ## Preview Server
    upload_module_preview_server("upload_preview", uploaded, checklist, checkTables)

    ## correctedX <- shiny::reactive({
    correctedX <- upload_module_batchcorrect_server(
      id = "batchcorrect",
      X = shiny::reactive(checked_counts()$matrix),
      is.count = TRUE,
      pheno = shiny::reactive(checked_samples()$matrix),
      height = height
    )

    corrected_counts <- shiny::reactive({
      counts <- NULL
      if (input$advanced_mode) {
        out <- correctedX()
        counts <- pmax(2**out$X - 1, 0)
      } else {
        counts <- checked_counts()$matrix
      }
      counts
    })

    modified_ct <- upload_module_makecontrast_server(
      id = "makecontrast",
      phenoRT = reactive(checked_samples()$matrix),
      contrRT = reactive(checked_contrasts()$matrix),
      countsRT = corrected_counts,
      height = height
    )

    shiny::observeEvent(modified_ct(), {
      ## Monitor for changes in the contrast matrix and if
      ## so replace the uploaded reactive values.
      modct <- modified_ct()
      uploaded[["contrasts.csv"]] <- modct$contr   ## trigger check
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
      samplesRT = shiny::reactive(checked_samples()$matrix),
      contrastsRT = shiny::reactive(modified_ct()$contr),
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



    # create an observer that will hide tabs Upload if selected organism if null and show if the button proceed_to_upload is clicked
    observeEvent(input$proceed_to_upload, {
      # show tab Upload
      shinyjs::runjs('document.querySelector(\'[data-value="Upload"]\').style.display = "";')
      # check on upload tab
      shinyjs::runjs('document.querySelector("a[data-value=\'Upload\']").click();')
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

    ## ------------------------------------------------
    ## Board return object
    ## ------------------------------------------------
    # board does not return anything
  })
}
