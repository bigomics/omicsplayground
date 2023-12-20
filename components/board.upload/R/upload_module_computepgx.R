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
    selected_organism,
    auth,
    create_raw_dir,
    enable_button = shiny::reactive(TRUE),
    alertready = TRUE,
    height = 720,
    recompute_info,
    inactivityCounter) {
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

      ONESAMPLE.GENE_METHODS <- c("ttest", "ttest.welch")
      ONESAMPLE.GENESET_METHODS <- sort(c("spearman", "gsva", "fgsea", "ssgsea", "fisher"))

      DEV.METHODS <- c("noLM.prune")
      DEV.NAMES <- c("noLM + prune")
      DEV.SELECTED <- c()

      readthedocs_url <- "https://omicsplayground.readthedocs.io/en/latest/dataprep/geneset.html"

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
                  shiny::tags$td("Organism"),
                  shiny::tags$td(shiny::tags$h6(selected_organism())),
                  shiny::tags$td("")
                ),
                shiny::tags$tr(
                  shiny::tags$td(""),
                  shiny::tags$td("Description"),
                  shiny::tags$td(shiny::div(
                    shiny::textAreaInput(
                      ns("upload_description"), NULL,
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
                      "remove.notexpressed",
                      "skip.normalization"
                    ),
                  choiceNames =
                    c(
                      "Transform features to gene symbols",
                      "protein-coding only",
                      "remove Rik/ORF/LOC genes",
                      "remove not-expressed",
                      "skip normalization"
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
                ),
                fileInput2(
                  ns("upload_annot_table"),
                  shiny::tags$h4("Counts annotation table (optional):"),
                  multiple = FALSE,
                  accept = c(".csv")
                )
              ),
              shiny::wellPanel(
                fileInput2(
                  ns("upload_gmt"),
                  shiny::tagList(
                    shiny::tags$h4("Custom genesets (.gmt) file:"),
                    shiny::tags$h6(
                      "A GMT file as described",
                      tags$a(
                        "here.",
                        href = readthedocs_url,
                        target = "_blank",
                        style = "text-decoration: underline;"
                      ),
                      "or download an",
                      downloadLink(ns("download_gmt"), shiny::HTML("<u>example GMT</u>")),
                      " (gene targets of the EGFR transcription factor)."
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
        suspendWhenHidden = FALSE
      ) ## important!!!

      shiny::observeEvent(enable_button(), {
        if (!enable_button()) {
          shinyjs::disable(ns("compute"))
        } else {
          shinyjs::enable(ns("compute"))
        }
      })
      # Input name and description
      shiny::observeEvent(list(metaRT(), recompute_info()), {
        meta <- metaRT()
        pgx_info <- recompute_info()
        # If the user recomputes, recycle old names/description
        if (is.null(pgx_info)) {
          if (!is.null(meta[["name"]])) {
            shiny::updateTextInput(session,
              "upload_name",
              value = meta[["name"]]
            )
          }
          if (!is.null(meta[["description"]])) {
            shiny::updateTextAreaInput(session,
              "upload_description",
              value = meta[["description"]]
            )
          }
        } else {
          shiny::updateTextInput(session,
            "upload_name",
            value = gsub(".pgx$", "", pgx_info[["name"]])
          )
          shiny::updateTextAreaInput(session,
            "upload_description",
            value = pgx_info[["description"]]
          )
        }
      })

      shiny::observeEvent(contrastsRT(), {
        contrasts <- as.data.frame(contrastsRT())
        has_one <- apply(contrasts, 2, function(x) all(table(x) == 1))

        if (any(has_one)) {
          shinyalert::shinyalert(
            title = "WARNING",
            text = stringr::str_squish("There are cases where there is only one samples
                    in a group. Some of the gene tests and enrichment
                    methods are disabled."),
            type = "warning"
          )
          shiny::updateCheckboxGroupInput(
            session,
            "gene_methods",
            choices = ONESAMPLE.GENE_METHODS,
            sel = "ttest"
          )
          shiny::updateCheckboxGroupInput(session,
            "gset_methods",
            choices = ONESAMPLE.GENESET_METHODS,
            sel = c("fisher", "fgsea", "gsva")
          )
        }
      })

      ## ------------------------------------------------------------------
      ## After confirmation is received, start computing the PGX
      ## object from the uploaded files
      ## ------------------------------------------------------------------

      # Define a reactive value to store the process object
      PROCESS_LIST <- list()
      computedPGX <- shiny::reactiveVal(NULL)
      process_counter <- reactiveVal(0)
      custom_geneset <- list(gmt = NULL, info = NULL)
      annot_table <- NULL
      processx_error <- list(user_email = NULL, pgx_name = NULL, pgx_path = NULL, error = NULL)

      ## react on custom GMT upload
      shiny::observeEvent(input$upload_gmt, {
        filePath <- input$upload_gmt$datapath
        fileName <- input$upload_gmt$name
        GSET_CHECK <- FALSE

        if (endsWith(filePath, ".txt") || endsWith(filePath, ".gmt")) {
          gmt <- playbase::read.gmt(filePath)

          # Clean genesets. eventually we can add multiple gmt files, but
          # for now we only support one
          gmt <- list(CUSTOM = gmt)
          gmt <- playbase::clean_gmt(gmt, "CUSTOM")

          # compute custom geneset stats
          gmt <- gmt[!duplicated(names(gmt))]
          gset_size <- sapply(gmt, length)
          gmt <- gmt[gset_size >= 3] ## how many?


          # an additional check to verify that items in lists are
          # genes
          if (length(gmt) > 0 && is.list(gmt)) {
            ## put in parent variable
            info <- list(GSET_SIZE = gset_size)
            custom_geneset <<- list(gmt = gmt, info = info)

            # tell user that custom genesets are "ok"
            shinyalert::shinyalert(
              title = "Custom genesets uploaded!",
              text = "Your genesets will be incorporated in the analysis.",
              type = "success",
              closeOnClickOutside = TRUE
            )
            GSET_CHECK <- TRUE
          }
        }

        # error message if custom genesets not detected or GSET_CHECK is FALSE add
        # https://omicsplayground.readthedocs.io/en/latest/dataprep/geneset.html
        # to guidelins, target blank as html
        if (is.null(custom_geneset$gmt) || !GSET_CHECK) {
          shinyalert::shinyalert(
            title = "Invalid custom genesets",
            text = "Please update a tsv file. See guidelines <a href='https://omicsplayground.readthedocs.io/en/latest/dataprep/geneset.html' target='_blank'>here</a>.",
            type = "error",
            html = TRUE,
            closeOnClickOutside = TRUE
          )
          custom_geneset <<- list(gmt = NULL, info = NULL)
          return(NULL)
        }
      })

      # react on upload_annot_table
      shiny::observeEvent(input$upload_annot_table, {
        # trigger a popup
        
        annot_able <<- playbase::read.as_matrix(input$upload_annot_table$datapath)

        shinyalert::shinyalert(
          title = "Annotation table uploaded!",
          text = "Your annotation table will be incorporated in the analysis.",
          type = "success",
          closeOnClickOutside = TRUE
        )
      })

      shiny::observeEvent(input$compute, {
        ## -----------------------------------------------------------
        ## Check validity
        ## -----------------------------------------------------------
        if (!enable_button()) {
          message("[ComputePgxServer:input$compute] WARNING:: *** DISABLED ***")
          shinyalert::shinyalert(
            title = "ERROR",
            text = "Compute is disabled",
            type = "error"
          )
          return(NULL)
        }
        if (!isValidFileName(input$upload_name)) {
          message("[ComputePgxServer:input$compute] WARNING:: Invalid name")
          shinyalert::shinyalert(
            title = "Invalid name",
            text = "Please remove any slashes (/) from the name",
            type = "error"
          )
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
        do.normalization <- !("skip.normalization" %in% flt)
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
        }

        dataset_name <- gsub("[ ]", "_", trimws(input$upload_name))
        creator <- auth$email
        libx.dir <- paste0(sub("/$", "", lib.dir), "x") ## set to .../libx

        pgx_save_folder <- auth$user_dir

        # Define create_pgx function arguments

        params <- list(
          organism = selected_organism(),
          samples = samples,
          counts = counts,
          contrasts = contrasts,
          batch.correct = FALSE,
          normalize = do.normalization,
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
          custom.geneset = custom_geneset,
          annot_table = annot_table,
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
          text = stringr::str_squish("Your dataset will be computed in the background.
            You can continue to play with a different dataset in the meantime.
            When it is ready, it will appear in your dataset library. Most datasets
            take between 30 - 60 minutes to complete."),
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


        new.job <- list(
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

        ## append to process list
        PROCESS_LIST <<- c(PROCESS_LIST, list(new.job))
      }) ## end observe input$compute

      check_process_status <- reactive({
        if (process_counter() == 0) {
          return(NULL)
        }

        # Re-execute this reactive expression after 30 seconds
        shiny::invalidateLater(30 * 1000, session)

        # When computing PGX, reset inactivity counter
        # to avoid session closure while computing
        inactivityCounter(0)

        completed_indices <- c()
        for (i in seq_along(PROCESS_LIST)) {
          active_obj <- PROCESS_LIST[[i]]
          current_process <- PROCESS_LIST[[i]]$process
          raw_dir <- PROCESS_LIST[[i]]$raw_dir

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
              ds_name <- paste0("<b>", PROCESS_LIST[[i]]$dataset_name, "</b>")
              if (!auth$email == "") {
                gmail_creds <- file.path(ETC, "gmail_creds")
                sendSuccessMessageToUser(
                  user_email = auth$email,
                  pgx_name = ds_name,
                  path_to_creds = gmail_creds
                )
              }
              raw_dir(NULL)
            } else {
              on_process_error(nr = nr)

              errors <- ""
              if (length(active_obj$stderr) > 0) {
                ## Copy the error to the stderr of main app
                message("Standard error from processx:")
                err <- paste0("[processx.", nr, ":stderr] ", active_obj$stderr)
                # save err to errors, separated by new lines
                errors <- paste0(errors, "Error:", "<br>")
                # append err to errors
                err <- paste0(err, cat = "<br>")
                errors <- c(errors, err, "<br>")
              }
              if (length(active_obj$stdout) > 0) {
                ## Copy the error to the stderr of main app
                cat("Standard output from processx:")
                out <- paste0("[processx.", nr, ":stdout] ", active_obj$stdout)
                out <- paste0(out, "<br>")
                errors <- c(errors, "Output:", "<br>")
                errors <- c(errors, out, "<br>")
              }

              ds_name <- paste0("<b>", PROCESS_LIST[[i]]$dataset_name, "</b>")
              title <- shiny::HTML(paste("The dataset", ds_name, "could not be computed."))

              # pass error data to parent variable
              processx_error <<- list(
                error      = errors,
                pgx_name   = ds_name,
                user_email = auth$email,
                pgx_path   = raw_dir
              )


              # if auth$email is empty, then the user is not logged in
              if (auth$email == "") {
                error_popup(
                  title   = "Error:",
                  header  = title,
                  message = "No email detected! Contact CS not possible!",
                  error   = shiny::HTML(errors),
                  btn_id  = "send_data_to_support__",
                  onclick = NULL
                )
              } else {
                error_popup(
                  title   = "Error:",
                  header  = title,
                  message = "Would you like to get support from our customer service?",
                  error   = shiny::HTML(errors),
                  btn_id  = "send_data_to_support__",
                  onclick = paste0('Shiny.onInputChange(\"', ns("send_data_to_support"), '\", this.id, {priority: "event"})')
                )
                # send error message to user
                gmail_creds <- file.path(ETC, "gmail_creds")

                sendErrorMessageToUser(
                  user_email = processx_error$user_email,
                  pgx_name = processx_error$pgx_name,
                  error = paste0(processx_error$error, collapse = ""),
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
          PROCESS_LIST <<- PROCESS_LIST[-completed_indices]
        }
        return(NULL)
      })

      # Function to execute when the process is completed successfully
      on_process_completed <- function(raw_dir, nr) {
        dbg("[computePGX:on_process_completed] process", nr, "completed!")
        process_counter(process_counter() - 1)

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
        dbg("[computePGX:on_process_completed] : saved result_pgx in ", result_pgx)

        if (file.exists(result_pgx)) {
          pgx <- playbase::pgx.load(result_pgx) ## always pgx
          computedPGX(pgx)
        } else {
          info("[computePGX:on_process_completed] : ERROR: Result file not found")
        }
        ## remove temp dir only if "user_input/raw_" is present in raw_dir
        if (grepl("raw_", raw_dir)) {
          # check if no ERROR_ files exist in raw_dir
          if (length(list.files(raw_dir, pattern = "ERROR_")) == 0) {
            unlink(raw_dir, recursive = TRUE)
          }
        }
      }

      on_process_error <- function(nr) {
        info("[computePGX:on_process_error] ERROR: process", nr, "completed with an error!")
        process_counter(process_counter() - 1) # stop the timer
      }

      ## This starts the check_process_status() reactive that gets
      ## invalidated/run every some minutes..
      observe(check_process_status())

      ## upon every new process check if we need to show the spinner
      ## or not, and if we are allowed to show the compute button.
      observeEvent(process_counter(), {
        if (process_counter() > 0) {
          shiny::insertUI(
            selector = "#current_dataset",
            where = "beforeBegin",
            ui = loading_spinner("Computation in progress..."),
            session = session
          )
        } else if (process_counter() == 0) {
          # remove UI with JS, had problems with shiny::removeUI
          shinyjs::runjs("document.querySelector('.current-dataset #spinner-container').remove();")
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
          user_email = processx_error$user_email,
          pgx_name = processx_error$pgx_name,
          pgx_path = processx_error$pgx_path,
          error = paste0(processx_error$error, collapse = ""),
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

      output$download_gmt <- downloadHandler(
        filename = function() {
          # Set the filename for the downloaded file
          "EGFR_TARGET_GENES.v2023.1.Hs.gmt"
        },
        content = function(file) {
          gmt_path <- file.path(FILES, "/gmt/EGFR_TARGET_GENES.v2023.1.Hs.gmt")
          gmt <- readBin(gmt_path, what = raw(), n = file.info(gmt_path)$size)
          writeBin(gmt, file)
        }
      )

      return(computedPGX)
    } ## end-of-server
  )
}
