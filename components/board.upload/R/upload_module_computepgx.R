##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_computepgx_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("UI"), fill = TRUE)
}

upload_module_computepgx_server <- function(
    id,
    countsRT,
    countsX,
    impX,
    norm_method,
    samplesRT,
    contrastsRT,
    annotRT = reactive(NULL),
    raw_dir,
    metaRT,
    lib.dir,
    auth,
    create_raw_dir,
    alertready = TRUE,
    height = 720,
    compute_settings,
    inactivityCounter,
    upload_wizard,
    upload_name,
    upload_description,
    upload_datatype,
    upload_organism,
    upload_gx_methods,
    upload_gset_methods,
    process_counter,
    reset_upload_text_input,
    probetype) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      ## statistical method for GENE level testing
      GENETEST.METHODS <- shiny::eventReactive(
        {
          upload_datatype()
        },
        {
          if (grepl("proteomics", upload_datatype(), ignore.case = TRUE)) {
            mm <- c("ttest", "ttest.welch", "trend.limma", "notrend.limma")
          } else {
            mm <- c(
              "ttest", "ttest.welch", "voom.limma", "trend.limma", "notrend.limma",
              "deseq2.wald", "deseq2.lrt", "edger.qlf", "edger.lrt"
            )
          }
          return(mm)
        }
      )

      GENETEST.SELECTED <- shiny::eventReactive(
        {
          upload_datatype()
        },
        {
          if (grepl("proteomics", upload_datatype(), ignore.case = TRUE)) {
            mm <- c("ttest", "ttest.welch", "trend.limma", "notrend.limma")
          } else {
            mm <- c("trend.limma", "voom.limma", "deseq2.wald", "edger.qlf")
          }
          return(mm)
        }
      )

      ## statistical method for GENESET level testing
      GENESET.METHODS <- c(
        "fisher", "ssgsea", "gsva", "spearman", "camera", "fry",
        ## "plage","enricher","gsea.permPH","gsea.permGS","gseaPR",
        "fgsea"
      )
      GENESET.SELECTED <- c("fisher", "gsva", "ssgsea", "fgsea")

      ## batch correction and extrs methods
      EXTRA.METHODS <- c("deconv", "drugs", "wordcloud", "connectivity", "wgcna")
      EXTRA.NAMES <- c(
        "celltype deconvolution", "drugs connectivity",
        "wordcloud", "experiment similarity", "WGCNA"
      )
      EXTRA.SELECTED <- c("deconv", "drugs", "wordcloud", "connectivity", "wgcna")

      ONESAMPLE.GENE_METHODS <- c("ttest", "ttest.welch")
      ONESAMPLE.GENESET_METHODS <- sort(c("fgsea", "fisher"))
      DEV.METHODS <- c("noLM.prune")
      DEV.NAMES <- c("noLM + prune")
      DEV.SELECTED <- c()

      ## Probe filtering defaults
      PROBE_FILTER_SELECTED <- DEFAULTS$computation_options$probe_filtering

      readthedocs_url <- "https://omicsplayground.readthedocs.io/en/latest/dataprep/geneset.html"

      output$UI <- shiny::renderUI({
        upload_annot_table_ui <- NULL
        if (auth$options$ENABLE_ANNOT) {
          upload_annot_table_ui <- fileInput2(
            ns("upload_annot_table"),
            shiny::tags$h4("Probe annotation (alpha):"),
            multiple = FALSE,
            accept = c(".csv")
          )
        }
        div(
          style = "overflow: auto;",
          bslib::as_fill_carrier(),
          bslib::layout_columns(
            width = "100%",
            col_widths = c(-5, 2, -5),
            fill = FALSE,
            div(
              ## style = "display: flex; flex-direction: column; align-items: center; gap: 20px; width: 100%;",
              style = "display: flex; flex-direction: column; gap: 20px; width: 100%;",
              shiny::div(
                style = "margin-left: 0px; text-align: left; width: 100%;",
                shiny::uiOutput(ns("input_recap2"))
              ),
              div(
                p("Dataset name:", style = "text-align: left;  margin: 0 0 2px 0; ;  font-weight: bold;"),
                shiny::textInput(
                  ns("selected_name"), NULL,
                  placeholder = "Name of your dataset"
                )
              ),
              div(
                p("Description:", style = "text-align: left;   margin: 0 0 2px 0;; font-weight: bold;"),
                shiny::textAreaInput(
                  ns("selected_description"), NULL,
                  placeholder = "Give a short description of your dataset",
                  height = 60, resize = "none"
                )
              )
            )
          ), ## end layout_col
          if (!is.null(probetype()) && probetype() == "running") {
            shiny::div(
              style = "display: flex; justify-content: center; align-items: center;",
              shiny::tags$h4(
                "Probe type detection still running, please wait...",
                style = "color: red;"
              )
            )
          },
          shiny::div(
            shiny::uiOutput(ns("probetype_result"))
          ),
          shiny::div(
            shiny::actionLink(ns("options"), "Computation options",
              icon = icon("cog", lib = "glyphicon")
            ),
            style = "display: flex; justify-content: center; margin: 15px 0;"
          ),
          shiny::conditionalPanel(
            "input.options%2 == 1",
            ns = ns,
            bslib::layout_columns(
              width = 12,
              bslib::card(
                shiny::checkboxGroupInput(
                  ns("filter_methods"),
                  ## shiny::HTML("<h4>Probe filtering:</h4>"),
                  shiny::HTML("<h4>Feature filtering:</h4>"),
                  choiceValues =
                    c(
                      "append.symbol",
                      "remove.notexpressed",
                      "remove.unknown",
                      "only.proteincoding"
                    ),
                  choiceNames =
                    c(
                      "Append symbol to feature ID",
                      "Remove not-expressed",
                      "Remove features without symbol",
                      "Remove Rik/ORF/LOC genes"
                    ),
                  selected = PROBE_FILTER_SELECTED
                )
              ),
              bslib::card(
                shiny::checkboxGroupInput(
                  ns("gene_methods"),
                  shiny::HTML("<h4>Gene tests:</h4>"),
                  GENETEST.METHODS(),
                  selected = GENETEST.SELECTED()
                )
              ),
              bslib::card(
                shiny::checkboxGroupInput(
                  ns("gset_methods"),
                  shiny::HTML("
                    <div style='display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; width: 100%;'>
                      <h4>Enrichment methods:</h4>
                      <a href='https://omicsplayground.readthedocs.io/en/latest/methods/' target='_blank' class='info-link' style='margin-left: 15px;'>
                        <i class='fa-solid fa-circle-info info-icon' style='color: blue; font-size: 20px;'></i>
                      </a>
                    </div>
                  "),
                  # <a href='https://example.com' target='_blank' id='infoButton' style='flex-shrink: 0; padding: 10px 20px; background-color: blue; color: white; text-decoration: none; border-radius: 4px;'>Info</a>
                  GENESET.METHODS,
                  selected = GENESET.SELECTED
                ),
              ),
              bslib::card(
                shiny::checkboxGroupInput(
                  ns("extra_methods"),
                  shiny::HTML("<h4>Extra analysis:</h4>"),
                  choiceValues = EXTRA.METHODS,
                  choiceNames = EXTRA.NAMES,
                  selected = EXTRA.SELECTED
                ),
                shiny::checkboxGroupInput(
                  ns("dev_options"),
                  shiny::HTML("<br><h4>Developer options:</h4>"),
                  choiceValues = DEV.METHODS,
                  choiceNames = DEV.NAMES,
                  selected = DEV.SELECTED
                )
              ),
              bslib::card(
                fileInput2(
                  ns("upload_gmt"),
                  shiny::tagList(
                    shiny::tags$h4("Custom genesets:"),
                    shiny::p(
                      "Upload a custom GMT file (.gmt) as described",
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
                upload_annot_table_ui
              )
            ), ## end of fillRow
            tags$style(HTML("#upload-compute-gset_methods-label { width: -webkit-fill-available; }")),
            tags$style(HTML("
            .flex-button {
              display:block;
              width:25px;
              height:25px;
              line-height:25px;
              border-radius: 50%;
              color:#fff;
              text-align:center;
              background: #555777;
              box-shadow: 0 0 3px gray;
              font-size:20px;
              font-weight:bold;
            }
            .flex-button:hover {
              color:#fff;
              background: black;
            }
            "))
          ) ## end of conditional panel
        )
      })

      # Input validators
      iv <- shinyvalidate::InputValidator$new()
      iv$enable()
      shiny::observeEvent(input$selected_name, {
        iv$add_rule("selected_name", shinyvalidate::sv_required())
      })
      shiny::observeEvent(input$selected_description, {
        iv$add_rule("selected_description", shinyvalidate::sv_required())
      })

      shiny::outputOptions(output,
        "UI",
        suspendWhenHidden = FALSE
      ) ## important!!!  Really???

      shiny::observeEvent(
        {
          list(countsRT(), samplesRT(), contrastsRT())
        },
        {
          ## invalidate any previously computed pgx
          computedPGX(NULL)
        }
      )

      # reset dataset name and description
      observeEvent(reset_upload_text_input(), {
        shiny::updateTextInput(session, "selected_name", value = "")
        shiny::updateTextAreaInput(session, "selected_description", value = "")
      })


      # change upload_name to selected_name
      observeEvent(input$selected_name, {
        upload_name(input$selected_name)
      })

      observeEvent(input$selected_description, {
        upload_description(input$selected_description)
      })

      # save input$gene_methods to upload_gx_methods
      observeEvent(input$gene_methods,
        {
          upload_gx_methods(input$gene_methods)
        },
        ignoreNULL = FALSE
      )

      # save input$gset_methods to upload_gset_methods
      observeEvent(input$gset_methods,
        {
          upload_gset_methods(input$gset_methods)
        },
        ignoreNULL = FALSE
      )

      output$input_recap2 <- renderUI({
        tagList(
          shiny::HTML("<b>Data type:</b>", upload_datatype()),
          shiny::HTML("<br><b>Organism:</b>", upload_organism()),
          shiny::HTML("<br><b>Probe type:</b>", probetype())
          ## shiny::HTML("<b>Name:</b><br>", upload_name()),
          ## shiny::HTML("<b>Description:</b><br>", upload_description())
        )
      })

      # handle ah task result
      output$probetype_result <- shiny::renderUI({
        p <- probetype()
        if (is.null(p) || p == "error") {
          shiny::div(
            style = "display: flex; justify-content: center; align-items: center; color: red;",
            shiny::tags$p("Probes not recognized, please check organism or your probe names.")
          )
        }
      })

      # Input name and description. NEED CHECK!!! seems not to
      # work. 18.11.24IK.
      shiny::observeEvent(list(metaRT(), compute_settings), {
        meta <- metaRT()
        pgx_info <- compute_settings
        if (is.null(meta) && is.null(pgx_info)) {
          return(NULL)
        }

        if (!is.null(pgx_info) && length(pgx_info) > 0) {
          meta <- pgx_info
        }

        # If the user recomputes, recycle old names/description
        if (!is.null(meta$name)) {
          shiny::updateTextInput(
            session,
            "selected_name",
            value = gsub(".pgx$", "", meta$name)
          )
        }
        if (!is.null(meta$description)) {
          shiny::updateTextAreaInput(
            session,
            "selected_description",
            value = meta$description
          )
        }
      })

      shiny::observeEvent(contrastsRT(), {
        contrasts <- as.data.frame(contrastsRT())
        has_one <- apply(contrasts, 2, function(x) any(table(x) == 1))
        if (any(has_one)) {
          shinyalert::shinyalert(
            title = "WARNING",
            text = stringr::str_squish("There are cases where there is only one samples in a group. Some of the gene tests and enrichment methods are disabled. Please note: a good experiment should have at least 3 replicates per condition."),
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
            sel = c("fisher", "fgsea")
          )
        } else {
          shiny::updateCheckboxGroupInput(
            session,
            "gene_methods",
            choices = GENETEST.METHODS(),
            sel = GENETEST.SELECTED()
          )
          shiny::updateCheckboxGroupInput(session,
            "gset_methods",
            choices = GENESET.METHODS,
            sel = GENESET.SELECTED
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
      custom_geneset <- list(gmt = NULL, info = NULL)
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

      shiny::observeEvent(upload_wizard(), {
        if (!is.null(upload_wizard()) && upload_wizard() != "wizard_finished") {
          return(NULL)
        }

        ## bail out if probetype task is not finished or has error
        p <- probetype()
        dbg("[computepgx_server:upload_wizard] start compute PGX!!!")
        dbg("[computepgx_server:upload_wizard] probetype = ", p)

        if (is.null(p) || grepl("error", tolower(p)) || p == "") {
          dbg("[computepgx_server:upload_wizard] ERROR probetype failed")
          shinyalert::shinyalert("ERROR", "probetype detection failed",
            type = "error"
          )
          return(NULL)
        }
        shiny::req(!(p %in% c("error", "running", ""))) ## wait for process??

        ## -----------------------------------------------------------
        ## Retrieve the most recent matrices from reactive values
        ## -----------------------------------------------------------

        counts <- countsRT()
        countsX <- countsX()
        impX <- impX()
        samples <- samplesRT()
        samples <- data.frame(samples, stringsAsFactors = FALSE, check.names = FALSE)
        contrasts <- as.matrix(contrastsRT())
        annot_table <- annotRT()

        nmissing.counts <- sum(is.na(counts))
        nmissing.countsX <- sum(is.na(countsX))
        if (nmissing.counts > 0 || nmissing.countsX > 0) {
          shinyalert::shinyalert(
            title = "WARNING",
            text = stringr::str_squish("Missing values are present in your data. You chose not to impute. The following differential gene expression (DGE) tests are currently unsupported with missing values: limma/voom, DESeq2 and edgeR LRT, QL F-test, Wald test. Please adjust your DGE test selection accordingly."),
            type = "warning",
            timer = 60000
          )
        }

        ## -----------------------------------------------------------
        ## Set statistical methods and run parameters
        ## -----------------------------------------------------------
        max.genes <- as.integer(auth$options$MAX_GENES)
        max.genesets <- as.integer(auth$options$MAX_GENESETS)

        ## get selected methods from input
        gx.methods <- input$gene_methods
        gset.methods <- input$gset_methods
        extra.methods <- input$extra_methods

        ## at least do meta.go, infer
        extra.methods <- unique(c("meta.go", "infer", extra.methods))

        ## ----------------------------------------------------------------------
        ## Start computation
        ## ----------------------------------------------------------------------

        flt <- ""
        use.design <- TRUE
        prune.samples <- FALSE
        flt <- input$filter_methods
        append.symbol <- ("append.symbol" %in% flt)
        do.protein <- ("proteingenes" %in% flt)
        remove.unknown <- ("remove.unknown" %in% flt)
        ## do.normalization <- !("skip.normalization" %in% flt)
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

        dataset_name <- gsub("[ ]", "_", trimws(upload_name()))
        creator <- auth$email
        libx.dir <- paste0(sub("/$", "", lib.dir), "x") ## set to .../libx

        pgx_save_folder <- auth$user_dir

        ## Define create_pgx function arguments
        params <- list(
          organism = upload_organism(),
          samples = samples,
          counts = counts,
          countsX = countsX,
          impX = impX,
          contrasts = contrasts,
          probe_type = probetype(),
          # Extra tables
          annot_table = annot_table,
          custom.geneset = custom_geneset,
          # Options
          batch.correct = FALSE,
          norm_method = norm_method(),
          settings = list(
            ## compute settings only for info
            imputation_method = compute_settings$imputation_method,
            bc_method = compute_settings$bc_method,
            remove_outliers = compute_settings$remove_outliers,
            norm_method = norm_method()
          ),
          ## normalize = do.normalization,
          prune.samples = TRUE,
          filter.genes = filter.genes,
          only.known = remove.unknown,
          only.proteincoding = only.proteincoding,
          only.hugo = append.symbol, ## DEPRECATED
          convert.hugo = append.symbol, ## should be renamed
          do.cluster = TRUE,
          cluster.contrasts = FALSE,
          max.genes = max.genes,
          max.genesets = max.genesets,
          gx.methods = gx.methods,
          gset.methods = gset.methods,
          extra.methods = extra.methods,
          use.design = use.design, ## no.design+prune are combined
          prune.samples = prune.samples,
          do.cluster = TRUE,
          libx.dir = libx.dir, # needs to be replaced with libx.dir
          name = dataset_name,
          datatype = upload_datatype(),
          description = input$selected_description,
          creator = creator,
          date = this.date,
          pgx.save.folder = pgx_save_folder,
          ETC = ETC,
          email = auth$email,
          sendSuccessMessageToUser = sendSuccessMessageToUser
        )

        path_to_params <- file.path(raw_dir(), "params.RData")
        saveRDS(params, file = path_to_params)

        # Normalize paths
        script_path <- normalizePath(file.path(get_opg_root(), "bin", "pgxcreate_op.R"))
        tmpdir <- normalizePath(raw_dir())

        # Remove global variables
        try(rm(annot_table))
        try(rm(custom_geneset))

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
          dataset_name = gsub("[ ]", "_", input$selected_name),
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
                ds_name <- paste0("<b>", PROCESS_LIST[[i]]$dataset_name, "</b>")
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
        info("[computePGX:on_process_completed] raw_dir = ", raw_dir)
        if (grepl("raw_", raw_dir)) {
          # check if no ERROR_ files exist in raw_dir
          if (length(list.files(raw_dir, pattern = "ERROR_")) == 0) {
            info("[computePGX:on_process_completed] : SUCCESS: Removing folder ", raw_dir)
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
          shinyjs::runjs("document.querySelector('.current-dataset #spinner-container')?.remove();")
        }
      })

      # observer to listed to click on send_data_to_support button
      observeEvent(input$send_data_to_support, {
        # write a message to console with shinyjs
        shinyjs::runjs("console.log('send_data_to_support button clicked')")
        message("send_data_to_support button clicked")

        credential <- file.path(ETC, "hubspot_creds")

        sendErrorMessageToCustomerSuport(
          user_email = processx_error$user_email,
          pgx_name = processx_error$pgx_name,
          pgx_path = processx_error$pgx_path,
          error = paste0(processx_error$error, collapse = ""),
          path_to_creds = credential
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
