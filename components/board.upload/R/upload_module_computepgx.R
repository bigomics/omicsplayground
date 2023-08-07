##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

ComputePgxGadget <- function(counts, samples, contrasts, height = 720) {
  gadgetize(
    ComputePgxUI, ComputePgxServer,
    title = "ComputePGX",
    countsRT = shiny::reactive(counts),
    samplesRT = shiny::reactive(samples),
    contrastsRT = shiny::reactive(contrasts)
  )
}

upload_module_computepgx_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("UI"))
}

upload_module_computepgx_server <- function(
    id,
    countsRT,
    samplesRT,
    contrastsRT,
    raw_dir,
    batchRT,
    metaRT,
    lib.dir,
    auth,
    create_raw_dir,
    enable_button = TRUE,
    alertready = TRUE,
    height = 720) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      ## statistical method for GENE level testing
      GENETEST.METHODS <- c(
        "ttest", "ttest.welch", "voom.limma", "trend.limma", "notrend.limma",
        "deseq2.wald", "deseq2.lrt", "edger.qlf", "edger.lrt"
      )
      GENETEST.SELECTED <- c("trend.limma", "deseq2.wald", "edger.qlf")

      ## statistical method for GENESET level testing
      GENESET.METHODS <- c(
        "fisher", "ssgsea", "gsva", "spearman", "camera", "fry",
        ## "plage","enricher","gsea.permPH","gsea.permGS","gseaPR",
        "fgsea"
      )
      GENESET.SELECTED <- c("fisher", "gsva", "fgsea")

      ## batch correction and extrs methods
      EXTRA.METHODS <- c("deconv", "drugs", "wordcloud", "connectivity", "wgcna")
      EXTRA.NAMES <- c(
        "celltype deconvolution", "drugs connectivity",
        "wordcloud", "experiment similarity", "WGCNA"
      )
      EXTRA.SELECTED <- c("deconv", "drugs", "wordcloud", "connectivity", "wgcna")

      ONESAMPLE.GENE_METHODS <- c('ttest','ttest.welch')
      ONESAMPLE.GENESET_METHODS <- sort(c('spearman','gsva','fgsea','ssgsea','fisher'))
 
      DEV.METHODS <- c("noLM.prune")
      DEV.NAMES <- c("noLM + prune")
      DEV.SELECTED <- c()

      path_gmt <- "https://omicsplayground.readthedocs.io/en/latest/"

      output$UI <- shiny::renderUI({
        shiny::fillCol(
          height = height,
          flex = c(0.2, NA, 0.05, 1.3),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(
              12,
              align = "center", offset = 0,
              shiny::tags$table(
                style = "width:100%;vertical-align:top;padding:4px;",
                shiny::tags$tr(
                  shiny::tags$td("", width = "300"),
                  shiny::tags$td("Name", width = "100"),
                  shiny::tags$td(
                    shiny::textInput(
                      ns("upload_name"), NULL, ## "Dataset:",
                      placeholder = "Name of your dataset"
                    ),
                    width = "600"
                  ),
                  shiny::tags$td("", width = "120")
                ),
                shiny::tags$tr(
                  shiny::tags$td(""),
                  shiny::tags$td("Datatype"),
                  shiny::tags$td(shiny::selectInput(
                    ns("upload_datatype"), NULL,
                    choices = c(
                      "RNA-seq", "scRNA-seq", "proteomics",
                      "mRNA microarray", "other"
                    )
                  )),
                  shiny::tags$td("")
                ),
                shiny::tags$tr(
                  shiny::tags$td(""),
                  shiny::tags$td("Description"),
                  shiny::tags$td(shiny::div(
                    shiny::textAreaInput(
                      ns("upload_description"), NULL, ## "Description:",
                      placeholder = "Give a short description of your dataset",
                      height = 100, resize = "none"
                    ),
                    style = "margin-left: 0px;"
                  )),
                  shiny::tags$td("")
                )
              ),
              shiny::br(),
              shiny::div(
                shiny::actionButton(ns("compute"), "Compute!",
                  icon = icon("running"),
                  class = "btn-outline-primary"
                ),
                shiny::br(), br(),
                shiny::actionLink(ns("options"), "Computation options",
                  icon = icon("cog", lib = "glyphicon")
                ),
                style = "padding-right: 80px;"
              )
            )
          ),
          shiny::br(),
          shiny::conditionalPanel(
            "input.options%2 == 1",
            ns = ns,
            shiny::fillRow(
              shiny::wellPanel(
                style = "width: 95%;",
                shiny::checkboxGroupInput(
                  ns("filter_methods"),
                  shiny::HTML("<h4>Feature filtering:</h4><br/>"),
                  choiceValues =
                    c(
                      "only.hugo",
                      "only.proteincoding",
                      "remove.unknown",
                      "remove.notexpressed"
                      ## "excl.immuno"
                      ## "excl.xy"
                    ),
                  choiceNames =
                    c(
                      "convert to HUGO",
                      "protein-coding only",
                      "remove Rik/ORF/LOC genes",
                      "remove not-expressed"
                      ## "Exclude immunogenes",
                    ),
                  selected = c(
                    "only.hugo",
                    "only.proteincoding",
                    "remove.unknown",
                    "remove.notexpressed"
                  )
                )
              ),
              shiny::wellPanel(
                style = "width: 95%;",
                shiny::checkboxGroupInput(
                  ns("gene_methods"),
                  shiny::HTML("<h4>Gene tests:</h4><br/>"),
                  GENETEST.METHODS,
                  selected = GENETEST.SELECTED
                )
              ),
              shiny::wellPanel(
                style = "width: 95%;",
                shiny::checkboxGroupInput(
                  ns("gset_methods"),
                  shiny::HTML("<h4>Enrichment methods:</h4><br/>"),
                  GENESET.METHODS,
                  selected = GENESET.SELECTED
                ),
              ),
              shiny::wellPanel(
                style = "width: 95%;",
                shiny::checkboxGroupInput(
                  ns("extra_methods"),
                  shiny::HTML("<h4>Extra analysis:</h4><br/>"),
                  choiceValues = EXTRA.METHODS,
                  choiceNames = EXTRA.NAMES,
                  selected = EXTRA.SELECTED
                )
              ),
              shiny::wellPanel(
                fileInput2(
                  ns("upload_custom_genesets"),
                  shiny::tagList(
                    shiny::tags$h4("Custom genesets (.gmt) file:"),
                    shiny::tags$h6(
                      "A GMT file as described",
                      tags$a(
                        "here.",
                        href = path_gmt,
                        target = "_blank",
                        style = "text-decoration: underline;"
                      )
                    )
                  ),
                  multiple = FALSE,
                  accept = c(".txt", ".gmt")
                ),
                shiny::checkboxGroupInput(
                  ns("dev_options"),
                  shiny::HTML("<h4>Developer options:</h4><br/>"),
                  choiceValues = DEV.METHODS,
                  choiceNames = DEV.NAMES,
                  selected = DEV.SELECTED
                )
              )
            ) ## end of fillRow
          ) ## end of conditional panel
        ) ## end of fill Col
      })
      shiny::outputOptions(output,
                           "UI",
                           suspendWhenHidden = FALSE) ## important!!!

      shiny::observeEvent(enable_button(), {
        if (!enable_button()) {
          shinyjs::disable(ns("compute"))
        } else {
          shinyjs::enable(ns("compute"))
        }
      })

      shiny::observeEvent(metaRT(), {
        meta <- metaRT()

        if (!is.null(meta[["name"]])) {
          shiny::updateTextInput(session,
                                 "upload_name",
                                 value = meta[["name"]])
        }
        if (!is.null(meta[["description"]])) {
          shiny::updateTextAreaInput(session,
                                     "upload_description",
                                     value = meta[["description"]])
        }
      })
        
      shiny::observeEvent(contrastsRT(), {
        contrasts <- as.data.frame(contrastsRT())
        has_one <- apply(contrasts, 2, function(x) any(table(x) == 1))

        if (any(has_one)) {
          shinyalert::shinyalert(
            title = "WARNING",
            text = "There are cases where there is only one samples in a group. 
                    Some of the gene tests and enrichment 
                    methods are disabled.",
            type = "warning"
          )
          shiny::updateCheckboxGroupInput(
                                          session,
                                          "gene_methods",
                                          choices = ONESAMPLE.GENE_METHODS,
                                          sel = "ttest")
          shiny::updateCheckboxGroupInput(session,
                                          "gset_methods",
                                          choices = ONESAMPLE.GENESET_METHODS,
                                          sel = c("fisher", "fgsea", "gsva"))
          }
      })

      ## ------------------------------------------------------------------
      ## After confirmation is received, start computing the PGX
      ## object from the uploaded files
      ## ------------------------------------------------------------------

      # Define a reactive value to store the process object
      process_obj <- reactiveVal(NULL)
      computedPGX <- shiny::reactiveVal(NULL)
      process_counter <- reactiveVal(0)
      reactive_timer <- reactiveTimer(20000) # Triggers every 10000 milliseconds (20 second)
      custom.geneset <- reactiveValues(gmt = NULL, info = NULL)
      store_error_from_process <- reactiveValues(user_email = NULL, pgx_name = NULL, pgx_path = NULL, error = NULL)

      shiny::observeEvent(input$upload_custom_genesets, {
        filePath <- input$upload_custom_genesets$datapath

        fileName <- input$upload_custom_genesets$name

        if (endsWith(filePath, ".txt") || endsWith(filePath, ".gmt")) {
          custom.geneset$gmt <- playbase::read.gmt(filePath)
          # perform some basic checks
          gmt.length <- length(custom.geneset$gmt)
          gmt.is.list <- is.list(custom.geneset$gmt)

          # clean genesets

          # eventually we can add multiple files, but for now we only support one
          # just add more gmts to the list below
          custom.geneset$gmt <- list(CUSTOM = custom.geneset$gmt)

          # convert gmt to OPG standard
          custom.geneset$gmt <- playbase::clean_gmt(custom.geneset$gmt, "CUSTOM")

          # compute custom geneset stats
          custom.geneset$gmt <- custom.geneset$gmt[!duplicated(names(custom.geneset$gmt))]
          custom.geneset$info$GSET_SIZE <- sapply(custom.geneset$gmt, length)

          # tell user that custom genesets are "ok"
          # we could perform an additional check to verify that items in lists are genes
          if (gmt.length > 0 && gmt.is.list) {
            shinyalert::shinyalert(
              title = "Custom genesets uploaded!",
              text = "Your genesets will be incorporated in the analysis.",
              type = "success",
              closeOnClickOutside = TRUE
            )
          }
        }

        # error message if custom genesets not detected
        if (is.null(custom.geneset$gmt)) {
          shinyalert::shinyalert(
            title = "Invalid custom genesets",
            text = "Please update a .txt file. See guidelines here <PLACEHOLDER>.",
            type = "error",
            closeOnClickOutside = TRUE
          )
          custom.geneset <- list(gmt = NULL, info = NULL)
          return(NULL)
        }
      })

      shiny::observeEvent(input$compute, {
        ## -----------------------------------------------------------
        ## Check validity
        ## -----------------------------------------------------------
        if (!enable_button()) {
          message("[ComputePgxServer:@compute] WARNING:: *** NOT ENABLED ***")
          return(NULL)
        }

        max.datasets <- as.integer(auth$options$MAX_DATASETS)
        pgxdir <- auth$user_dir
        numpgx <- length(dir(pgxdir, pattern = "*.pgx$"))
        if (!auth$options$ENABLE_DELETE) {
          numpgx <- length(dir(pgxdir, pattern = "*.pgx$|*.pgx_$")) ## count deleted...
        }
        if (numpgx >= max.datasets) {
          shinyalert_storage_full() ## from ui-alerts.R
          return(NULL)
        }

        ## check for name and description
        has.name <- input$upload_name != ""
        has.description <- input$upload_description != ""
        if (!has.name || !has.description) {
          shinyalert::shinyalert(
            title = "ERROR",
            text = "You must give a dataset name and description",
            type = "error"
          )
          return(NULL)
        }

        ## -----------------------------------------------------------
        ## Retrieve the most recent matrices from reactive values
        ## -----------------------------------------------------------
        counts <- countsRT()
        samples <- samplesRT()
        samples <- data.frame(samples, stringsAsFactors = FALSE, check.names = FALSE)
        contrasts <- as.matrix(contrastsRT())

        ## !!!!!!!!!!!!!! This is blocking the computation !!!!!!!!!!!
        ## batch  <- batchRT()

        ## -----------------------------------------------------------
        ## Set statistical methods and run parameters
        ## -----------------------------------------------------------
        max.genes <- as.integer(auth$options$MAX_GENES)
        max.genesets <- as.integer(auth$options$MAX_GENESETS)

        dbg("[upload_module_computepgx_server] max.genes = ", max.genes)
        dbg("[upload_module_computepgx_server] max.genesets = ", max.genesets)

        ## get selected methods from input
        gx.methods <- input$gene_methods
        gset.methods <- input$gset_methods
        extra.methods <- input$extra_methods

        if (length(gx.methods) == 0) {
          shinyalert::shinyalert(
            title = "ERROR",
            text = "You must select at least one gene test method",
            type = "error"
          )
          return(NULL)
        }
        if (length(gset.methods) == 0) {
          shinyalert::shinyalert(
            title = "ERROR",
            text = "You must select at least one geneset test method",
            type = "error"
          )
          return(NULL)
        }

        ## at least do meta.go, infer
        extra.methods <- unique(c("meta.go", "infer", extra.methods))

        ## ----------------------------------------------------------------------
        ## Start computation
        ## ----------------------------------------------------------------------

        flt <- ""
        use.design <- TRUE
        prune.samples <- FALSE
        flt <- input$filter_methods
        only.hugo <- ("only.hugo" %in% flt)
        do.protein <- ("proteingenes" %in% flt)
        remove.unknown <- ("remove.unknown" %in% flt)
        excl.immuno <- ("excl.immuno" %in% flt)
        excl.xy <- ("excl.xy" %in% flt)
        only.proteincoding <- ("only.proteincoding" %in% flt)
        filter.genes <- ("remove.notexpressed" %in% flt)
        use.design <- !("noLM.prune" %in% input$dev_options)
        prune.samples <- ("noLM.prune" %in% input$dev_options)
        this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

        # if no raw_dir (happens when we auto-load example data via
        # button), or user click compute a second time
        if (is.null(raw_dir())) {
          raw_dir(create_raw_dir(auth))
          ## raw_dir(tempfile(pattern = "pgx_"))
          ## dir.create(raw_dir())
          ##
          ## NOTE: should we save any counts/samples/contrast matrices
          ## here?? They might be altered here but at least we have
          ## them saved.
          dbg("[compute PGX process] : tempFile", raw_dir())
        }

        dataset_name <- gsub("[ ]", "_", trimws(input$upload_name))
        creator <- auth$email
        libx.dir <- paste0(sub("/$", "", lib.dir), "x") ## set to .../libx

        pgx_save_folder <- auth$user_dir

        # get rid of reactive container
        custom.geneset <- list(gmt = custom.geneset$gmt, info = custom.geneset$info)
        # Define create_pgx function arguments
        params <- list(
          samples = samples,
          counts = counts,
          contrasts = contrasts,
          batch.correct = FALSE,
          prune.samples = TRUE,
          filter.genes = filter.genes,
          only.known = !remove.unknown,
          only.proteincoding = only.proteincoding,
          only.hugo = only.hugo,
          convert.hugo = only.hugo,
          do.cluster = TRUE,
          cluster.contrasts = FALSE,
          max.genes = max.genes,
          max.genesets = max.genesets,
          custom.geneset = custom.geneset,
          gx.methods = gx.methods,
          gset.methods = gset.methods,
          extra.methods = extra.methods,
          use.design = use.design, ## no.design+prune are combined
          prune.samples = prune.samples,
          do.cluster = TRUE,
          libx.dir = libx.dir, # needs to be replaced with libx.dir
          name = dataset_name,
          datatype = input$upload_datatype,
          description = input$upload_description,
          creator = creator,
          date = this.date,
          pgx.save.folder = pgx_save_folder
        )

        path_to_params <- file.path(raw_dir(), "params.RData")
        saveRDS(params, file = path_to_params)

        # Normalize paths
        script_path <- normalizePath(file.path(get_opg_root(), "bin", "pgxcreate_op.R"))
        tmpdir <- normalizePath(raw_dir())

        # Start the process and store it in the reactive value
        shinyalert::shinyalert(
          title = "Crunching your data!",
          text = "Your dataset will be computed in the background. 
          You can continue to play with a different dataset in the meantime. 
          When it is ready, it will appear in your dataset library.",
          type = "info",
          timer = 60000
        )
        bigdash.selectTab(
          session,
          selected = "load-tab"
        )

        process_counter(process_counter() + 1)
        dbg("[compute PGX process] : starting processx nr: ", process_counter())
        dbg("[compute PGX process] : process tmpdir = ", tmpdir)
        dbg("[compute PGX process] : see error.log => tail -f", paste0(tmpdir, "/processx-error.log"))
        ## append to process list
        process_obj(
          append(
            process_obj(),
            list(
              ## new background computation job
              list(
                process = processx::process$new(
                  "Rscript",
                  args = c(script_path, tmpdir),
                  supervise = TRUE,
                  stderr = "|",
                  stdout = "|"
                ),
                number = process_counter(),
                dataset_name = gsub("[ ]", "_", input$upload_name),
                raw_dir = raw_dir(),
                stderr = c(),
                stdout = c()
              )
            )
          )
        )
      })

      check_process_status <- reactive({
        if (process_counter() == 0) {
          return(NULL)
        }

        reactive_timer()
        active_processes <- process_obj()
        completed_indices <- c()

        for (i in seq_along(active_processes)) {
          active_obj <- active_processes[[i]]
          current_process <- active_processes[[i]]$process
          raw_dir <- active_processes[[i]]$raw_dir

          process_status <- current_process$get_exit_status()
          process_alive <- current_process$is_alive()
          nr <- active_obj$number

          ## [https://processx.r-lib.org/] Always
          ## make sure that you read out the standard
          ## output and/or error of the pipes, otherwise
          ## the background process will stop running!
          stderr_output <- current_process$read_error_lines()
          active_obj$stderr <- c(active_obj$stderr, stderr_output)
          stdout_output <- current_process$read_output_lines()
          active_obj$stdout <- c(active_obj$stdout, stdout_output)

          errlog <- file.path(raw_dir, "processx-error.log")
          outlog <- file.path(raw_dir, "processx-output.log")
          cat(paste(stderr_output, collapse = "\n"), file = errlog, append = TRUE)
          cat(paste(stdout_output, collapse = "\n"), file = outlog, append = TRUE)

          if (!is.null(process_status)) {
            if (process_status == 0) {
              # Process completed successfully
              dbg("[compute PGX process] : process completed")
              on_process_completed(raw_dir = raw_dir, nr = nr)
              ds_name_bold <- paste0("<b>", active_processes[[i]]$dataset_name, "</b>")

              if (!auth$email == "") {
                gmail_creds <- file.path(ETC, "gmail_creds")
                sendSuccessMessageToUser(
                  user_email = auth$email,
                  pgx_name = ds_name_bold,
                  path_to_creds = gmail_creds
                )
              }
              raw_dir(NULL)
            } else {
              on_process_error(nr = nr)

              log_pgx_compute <- ""

              if (length(active_obj$stderr) > 0) {
                ## Copy the error to the stderr of main app
                message("Standard error from processx:")
                err <- paste0("[processx.", nr, ":stderr] ", active_obj$stderr)
                # save err to log_pgx_compute, separated by new lines
                log_pgx_compute <- paste0(log_pgx_compute, "Error:", "<br>")
                # append err to log_pgx_compute
                err <- paste0(err, cat = "<br>")
                log_pgx_compute <- c(log_pgx_compute, err, "<br>")
              }
              if (length(active_obj$stdout) > 0) {
                ## Copy the error to the stderr of main app
                cat("Standard output from processx:")
                out <- paste0("[processx.", nr, ":stdout] ", active_obj$stdout)
                out <- paste0(out, "<br>")
                log_pgx_compute <- c(log_pgx_compute, "Output:", "<br>")
                log_pgx_compute <- c(log_pgx_compute, out, "<br>")
              }

              ds_name_bold <- paste0("<b>", active_processes[[i]]$dataset_name, "</b>")
              title <- shiny::HTML(paste("The dataset", ds_name_bold, "could not be computed."))

              # pass error data to reactive
              store_error_from_process$error <- log_pgx_compute
              store_error_from_process$pgx_name <- ds_name_bold
              store_error_from_process$user_email <- auth$email
              store_error_from_process$pgx_path <- raw_dir

              # if auth$email is empty, then the user is not logged in
              if (auth$email == "") {
                error_popup(
                  title = "Error:",
                  header = title,
                  message = "No email detected! Contact CS not possible!",
                  error = shiny::HTML(log_pgx_compute),
                  btn_id = "send_data_to_support__",
                  onclick = NULL
                )
              } else {
                error_popup(
                  title = "Error:",
                  header = title,
                  message = "Would you like to get support from our customer service?",
                  error = shiny::HTML(log_pgx_compute),
                  btn_id = "send_data_to_support__",
                  onclick = paste0('Shiny.onInputChange(\"', ns("send_data_to_support"), '\", this.id, {priority: "event"})')
                )
                # send error message to user
                gmail_creds <- file.path(ETC, "gmail_creds")

                sendErrorMessageToUser(
                  user_email = store_error_from_process$user_email,
                  pgx_name = store_error_from_process$pgx_name,
                  error = paste0(store_error_from_process$error, collapse = ""),
                  path_to_creds = gmail_creds
                )
              }
              raw_dir(NULL)
            }
            completed_indices <- c(completed_indices, i)

            ## read latest output/error and copy to main stderr/stdout
            if (length(active_obj$stderr) > 0) {
              ## Copy the error to the stderr of main app
              message("Standard error from processx:")
              err <- paste0("[processx.", nr, ":stderr] ", active_obj$stderr)
              writeLines(err, con = stderr())
            }
            if (length(active_obj$stdout) > 0) {
              ## Copy the error to the stderr of main app
              cat("Standard output from processx:")
              out <- paste0("[processx.", nr, ":stdout] ", active_obj$stdout)
              writeLines(out, con = stdout())
            }
          } else {
            # Process is still running, do nothing
            dbg("[compute PGX process] : process still running")
          }
        }

        # Remove completed processes from the list
        if (length(completed_indices) > 0) {
          active_processes <- active_processes[-completed_indices]
          process_obj(active_processes)
        }
        return(NULL)
      })

      # Function to execute when the process is completed successfully
      on_process_completed <- function(raw_dir, nr) {
        dbg("[computePGX:on_process_completed] process", nr, "completed!")
        process_counter(process_counter() - 1) # stop the timer

        path_to_params <- file.path(raw_dir, "params.RData")
        params <- readRDS(path_to_params)
        pgx_save_folder_px <- params$pgx.save.folder

        # check if user folder matches pgx processx folder, it not stop here
        if (pgx_save_folder_px != auth$user_dir) {
          info("[computePGX:on_process_completed] : ERROR: pgx_save_folder != auth$user_dir)")
          dbg("[computePGX:on_process_completed] : pgx_save_folder = ", pgx_save_folder_px)
          dbg("[computePGX:on_process_completed] : auth$user_dir = ", auth$user_dir)
          return()
        }

        ## dataset_name <- paste0(gsub("[ ]", "_", input$upload_name),".pgx")
        dataset_name <- paste0(params$name, ".pgx")
        result_pgx <- file.path(pgx_save_folder_px, dataset_name)
        dbg("[computePGX:on_process_completed] : result_pgx = ", result_pgx)

        if (file.exists(result_pgx)) {
          pgx <- playbase::pgx.load(result_pgx) ## always pgx
          computedPGX(pgx)
        } else {
          info("[computePGX:on_process_completed] : ERROR: Result file not found")
        }
        ## remove temp dir only if "user_input/raw_" is present in raw_dir
        if (grepl("raw_", raw_dir)) {
          unlink(raw_dir, recursive = TRUE)
        }
      }

      on_process_error <- function(nr) {
        info("[computePGX:on_process_error] ERROR: process", nr, "completed with an error!")
        process_counter(process_counter() - 1) # stop the timer
      }

      ## what does this do???
      observe(check_process_status())

      observe({
        if (process_counter() > 0) {
          shiny::insertUI(
            selector = "#current_dataset",
            where = "beforeEnd",
            ui = loading_spinner("Computation in progress...")
          )
        } else if (process_counter() == 0) {
          shiny::removeUI(selector = "#current_dataset > #spinner-container")
        }

        MAX_DS_PROCESS <- 1
        if (process_counter() < MAX_DS_PROCESS) {
          shinyjs::enable("compute")
        } else {
          shinyjs::disable("compute")
        }
      })

      # observer to listed to click on send_data_to_support button
      observeEvent(input$send_data_to_support, {
        # write a message to console with shinyjs
        shinyjs::runjs("console.log('send_data_to_support button clicked')")
        message("send_data_to_support button clicked")

        gmail_creds <- file.path(ETC, "gmail_creds")

        sendErrorMessageToCustomerSuport(
          user_email = store_error_from_process$user_email,
          pgx_name = store_error_from_process$pgx_name,
          pgx_path = store_error_from_process$pgx_path,
          error = paste0(store_error_from_process$error, collapse = ""),
          path_to_creds = gmail_creds
        )

        # close modal

        shinyjs::runjs("document.getElementById('sendLogModal').style.display = 'none';")

        # alert user that a message was sent to CS

        shinyalert::shinyalert(
          title = "Message sent",
          text = "We are sorry you had a problem. We will get back to you as soon as possible.",
          type = "success"
        )
      })

      return(computedPGX)
    } ## end-of-server
  )
}
