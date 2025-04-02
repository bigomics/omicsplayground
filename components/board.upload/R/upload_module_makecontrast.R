##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

upload_module_makecontrast_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = 12,
    height = "100%",
    row_heights = c(3, 3),
    bslib::card(
      full_screen = FALSE,
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(3, 9),
          ## left 3-column
          shiny::div(
            shiny::HTML("<b>1. Choose phenotype:</b>"),
            withTooltip(
              shiny::selectizeInput(
                inputId = ns("param"),
                NULL,
                width = "100%",
                choices = "<samples>",
                selected = "<samples>",
                multiple = TRUE,
                options = list(maxItems = 3)
              ),
              "Select phenotype(s) to create conditions for your groups. Select &ltsamples&gt if you want to group manually on sample names. You can select multiple phenotypes to create combinations.",
              placement = "left", options = list(container = "body")
            ),
            shiny::div(
              style = "padding-top: 5px;",
              shiny::HTML("<b>3. Comparison name:</b>"),
              withTooltip(
                shiny::textInput(
                  ns("newname"),
                  NULL,
                  width = "100%",
                  placeholder = "e.g. MAIN_vs_CONTROL"
                ),
                "Give a name for your comparison as MAIN_vs_CONTROL, with the name of the main group first. You must keep _vs_ in the name to separate the names of the two groups."
              ),
              div(
                style = "padding-top: -2px;",
                shiny::checkboxInput(
                  ns("add_prefix"), "Add phenotype prefix",
                  TRUE
                )
              )
            ),
            shiny::div(
              style = "padding-top: 3px;",
              shiny::HTML("<b>4. Add my comparison:</b>")
            ),
            shiny::div(
              style = "margin-top: 0px;",
              withTooltip(
                shiny::actionButton(
                  ns("addcontrast"),
                  "add to list",
                  icon = icon("plus"),
                  class = "btn-sm btn-outline-primary",
                  width = "48%"
                ),
                "Add this comparison to the list.",
                placement = "top", options = list(container = "body")
              )
            )
          ),

          ## right 9-column
          shiny::div(
            style = "overflow: auto; margin-left: 30px;",
            shiny::HTML("<b>2. Create comparison:</b>"),
            withTooltip(
              shiny::uiOutput(
                ns("createcomparison"),
                style = "font-size:13px; margin-top: 2px;"
              ),
              "Create comparisons by dragging conditions into the main or control groups on the right. Then press add comparison to add them to the table."
            )
          )
        ) ## end of layout-columns c(3,9)
      ) # end of card-body
    ),
    #-----------------------------------------------------
    bslib::layout_columns(
      col_widths = c(8, 4),
      bslib::card(
        full_screen = FALSE,
        bslib::card_body(
          style = "padding: 10px 20px;",
          DT::dataTableOutput(ns("contrastTable"))
        )
      ),
      bslib::card(
        full_screen = FALSE,
        bslib::card_body(
          style = "padding: 10px 20px;",
          plotOutput(ns("contrastPlot"))
        )
      )
    )
  )
}

upload_module_makecontrast_server <- function(
    id,
    phenoRT,
    contrRT,
    countsRT,
    upload_wizard,
    show_comparison_builder,
    autocontrast) {
  shiny::moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      ## rv <- shiny::reactiveValues(contr = NULL)
      rv_contr <- shiny::reactiveVal(NULL)

      rv <- reactiveValues(
        condition_start = NULL,
        condition_group1 = NULL,
        condition_group2 = NULL
      )

      shiny::observe({
        rv_contr(contrRT())
      })

      reset_bucket <- function() {
        # update the rv values when the param changes
        rv$condition_start <- NULL
        rv$condition_group1 <- NULL
        rv$condition_group2 <- NULL
        shiny::updateTextInput(session, "newname", value = "")

        cond <- sel.conditions()
        if (length(cond) == 0 || is.null(cond)) {
          return(NULL)
        }
        ## items <- c("others"="<others>", sort(unique(cond)))
        items <- sort(unique(cond))
        rv$condition_start <- items
      }

      observeEvent(
        {
          list(input$param)
        },
        {
          reset_bucket()
        }
      )


      observeEvent(
        {
          list(phenoRT(), upload_wizard(), show_comparison_builder())
        },
        {
          req(
            upload_wizard() == "step_comparisons",
            show_comparison_builder() == TRUE
          )
          phenotypes <- c(sort(unique(colnames(phenoRT()))), "<samples>")
          phenotypes <- grep("_vs_", phenotypes, value = TRUE, invert = TRUE) ## no comparisons...
          psel <- c(grep("sample|patient|name|id|^[.]",
            phenotypes,
            value = TRUE, invert = TRUE
          ), phenotypes)[1]
          updateSelectInput(session, "param", choices = phenotypes, selected = psel)
        }
      )

      sel.conditions <- shiny::reactive({
        ## shiny::req( phenoRT(), countsRT())  ## shiny BUG if has NA!!!!
        shiny::req(dim(phenoRT()), dim(countsRT()))

        shiny::validate(shiny::need(
          length(input$param) > 0,
          "Please select at least one phenotype"
        ))

        ## create conditions from phenotype matrix
        df <- phenoRT()
        if ("<samples>" %in% input$param) {
          df <- cbind(df, "<samples>" = rownames(df))
        }
        df <- playbase::binarizeNumericalColumns(df)
        pp <- intersect(input$param, colnames(df))
        ss <- colnames(countsRT())
        df1 <- df[ss, pp, drop = FALSE]
        cond <- apply(df1, 1, paste, collapse = ".")
        names(cond) <- cond
        
        if (max(sapply(cond, function(x) nchar(x))) > 15) {
          shinyalert::shinyalert(
            title = "Warning",
            text = "Condition names are long (>15 characters). This may cause display issues.",
            type = "warning"
          )
        }

        ## get abbreviated phenotype, NOTE: commented out, some users complained.
        # minlen <- ifelse(length(pp) >= 2, 4, 8)
        # minlen <- ifelse(length(pp) >= 3, 3, minlen)
        # abv.df <- playbase::abbreviate_pheno(
        #   df,
        #   minlength = minlen, abbrev.colnames = FALSE
        # )
        # abv.df <- abv.df[ss, pp, drop = FALSE]
        # abv.cond <- apply(abv.df, 1, paste, collapse = ".")
        # names(cond) <- abv.cond

        cond
      })

      output$createcomparison <- shiny::renderUI({
        shiny::req(input$param)

        shiny::tagList(
          shiny::tags$head(shiny::tags$style(".default-sortable .rank-list-item {padding: 2px 15px;}")),
          shiny::tags$head(shiny::tags$style(".default-sortable.bucket-list-container {padding: 0px 0px;margin: -5px 0 0 -5px;}")),
          shiny::tags$head(shiny::tags$style(".default-sortable .rank-list-title {margin-bottom: 0px;}")),
          shiny::tags$head(shiny::tags$style(".default-sortable .rank-list-container {margin-top: 0px;}")),
          sortable::bucket_list(
            header = NULL,
            sortable::add_rank_list(
              input_id = ns("conditions"),
              text = "Conditions:",
              labels = rv$condition_start,
            ),
            sortable::add_rank_list(
              input_id = ns("group1"),
              text = "Main group:",
              labels = rv$condition_group1,
            ),
            sortable::add_rank_list(
              input_id = ns("group2"),
              text = "Control group:",
              labels = rv$condition_group2,
            ),
            group_name = "cmpbucket"
          )
        )
      })

      makebuttonInputs <- function(FUN, len, id, ...) {
        inputs <- character(len)
        for (i in seq_len(len)) {
          inputs[i] <- as.character(FUN(paste0(id, i), ...))
        }
        inputs
      }


      shiny::observeEvent(
        list(input$group1, input$group2, input$add_prefix),
        {
          shiny::updateTextInput(session, "newname", value = "")

          # make sure group1 and 2 are not NULL or ""
          req(input$group1, input$group2)

          # update rv
          rv$condition_group1 <- input$group1
          rv$condition_group2 <- input$group2

          # remove selected conditions from the start list
          rv$condition_start <- setdiff(rv$condition_start, c(input$group1, input$group2))

          # if some input$conditions are not in the start list, add them
          if (length(input$conditions) > 0 && any(!input$conditions %in% rv$condition_start)) {
            rv$condition_start <- sort(union(rv$condition_start, input$conditions))
          }

          ## create comparison name from group names
          cond <- sel.conditions()
          ## cond <- c("others" = "<others>", cond)
          g1 <- input$group1
          g2 <- input$group2
          ## if("<others>" %in% g2) g2 <- "<others>"
          g1 <- names(cond)[match(g1, cond)]
          g2 <- names(cond)[match(g2, cond)]
          g1 <- unique(unlist(strsplit(g1, split = "[.]")))
          g2 <- unique(unlist(strsplit(g2, split = "[.]")))
          g1 <- paste(g1, collapse = ".")
          g2 <- paste(g2, collapse = ".")

          if (is.null(g1) || length(g1) == 0) g1 <- ""
          if (is.null(g2) || length(g2) == 0) g2 <- ""
          if (is.na(g1)) g1 <- ""
          if (is.na(g2)) g2 <- ""
          if (g1 == "" && g2 == "") {
            tt <- ""
          } else {
            g1 <- substring(g1, 1, 20)
            g2 <- substring(g2, 1, 20)
            tt <- paste0(g1, "_vs_", g2)
            if (input$add_prefix) {
              abv.pheno <- abbreviate(colnames(phenoRT()), minlength = 10)
              sel <- match(input$param, colnames(phenoRT()))
              prm.name <- gsub("[-_.,<> ]", "", abv.pheno[sel])
              # prm.name <- stringr::str_to_title(trimws(prm.name))
              # prm.name[1] <- tolower(prm.name[1])
              prm.name <- paste(prm.name, collapse = ".")
              tt <- paste0(prm.name, ":", tt)
            }
          }

          shiny::updateTextInput(session, "newname", value = tt)
        }
      )

      shiny::observeEvent(c(input$group1, input$group2), {
        if (length(input$group1) && length(input$group2)) {
          shinyjs::removeClass(id = "addcontrast", class = "btn-outline-primary")
          shinyjs::addClass(id = "addcontrast", class = "btn-primary")
        } else {
          shinyjs::addClass(id = "addcontrast", class = "btn-outline-primary")
          shinyjs::removeClass(id = "addcontrast", class = "btn-primary")
        }
      })

      shiny::observeEvent(input$contrast_delete, {
        ## Observe if a contrast is to be deleted
        id <- as.numeric(gsub(".*_", "", input$contrast_delete))
        if (length(id) == 0) {
          return(NULL)
        }
        if (!is.null(rv_contr()) && NCOL(rv_contr()) <= 1) {
          rv_contr(rv_contr()[, 0, drop = FALSE])
        } else {
          rv_contr(rv_contr()[, -id, drop = FALSE])
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
        in.ref <- 1 * (cond %in% group2)
        # in.ref2 <- ("<others>" %in% group2) & (!cond %in% group1)
        # in.ref <- in.ref1 | in.ref2

        ct.name <- input$newname
        gr1 <- gsub(".*:|_vs_.*", "", ct.name) ## first is MAIN group!!!
        gr2 <- gsub(".*_vs_|@.*", "", ct.name)
        ctx <- c(NA, gr1, gr2)[1 + 1 * in.main + 2 * in.ref]

        if (sum(in.main) == 0 || sum(in.ref) == 0) {
          shinyalert::shinyalert("ERROR", "Both groups must have samples")
          return(NULL)
        }
        if (ct.name %in% c(NA, "", " ")) {
          shinyalert::shinyalert("ERROR", "You must give the comparison a name")
          return(NULL)
        }
        if (gr1 == gr2) {
          shinyalert::shinyalert("ERROR", "Invalid comparison name")
          shiny::updateTextInput(session, "newname", value = "")
          return(NULL)
        }
        if (!is.null(rv_contr()) && ct.name %in% colnames(rv_contr())) {
          shinyalert::shinyalert("ERROR", "Comparison name already exists.")
          shiny::updateTextInput(session, "newname", value = "")
          return(NULL)
        }
        if (!grepl("_vs_", ct.name)) {
          shinyalert::shinyalert("ERROR", "Comparison must include _vs_ in name")
          shiny::updateTextInput(session, "newname", value = "")
          return(NULL)
        }

        ## update reactive value
        samples <- colnames(countsRT())
        ctx1 <- matrix(ctx, ncol = 1, dimnames = list(samples, ct.name))
        if (is.null(rv_contr())) {
          rv_contr(ctx1)
        } else {
          rv_contr(cbind(rv_contr(), ctx1))
        }

        # reset text input
        reset_bucket()
      })

      shiny::observeEvent(autocontrast(), {
        shiny::req(phenoRT())
        df <- phenoRT()

        ctx <- NULL
        ct <- playbase::pgx.makeAutoContrasts(
          df,
          mingrp = 3, slen = 20, ref = NULL, fix.degenerate = FALSE
        )
        if (!is.null(ct)) {
          ctx <- ct$exp.matrix
          rownames(ctx) <- rownames(df)
        }
        if (is.null(ctx)) {
          return(NULL)
        }

        ## update reactive value
        ctx2 <- playbase::contrastAsLabels(ctx)
        if (!is.null(rv_contr())) {
          rv_contr(cbind(rv_contr(), ctx2))
        } else {
          rv_contr(ctx2)
        }
      })

      output$contrastTable <- DT::renderDataTable(
        {
          ct <- rv_contr()

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
            paste.max <- function(x, n = 5) {
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

            deleteButtons <- makebuttonInputs(
              FUN = actionButton,
              len = ncol(ct),
              id = paste0("contrast_delete_", sample(99999, 1), "_"), ## hack to allow double click
              label = "",
              width = "50px",
              inline = TRUE,
              icon = shiny::icon("trash-alt"),
              class = "btn-inline btn-outline-danger-hover",
              style = "padding:0px; margin:0px; font-size:95%; color: #B22222;",
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
            fillContainer = TRUE,
            rownames = FALSE,
            escape = c(-1),
            selection = "none",
            class = "compact",
            plugins = "scrollResize",
            options = list(
              dom = "t",
              scrollResize = TRUE,
              pageLength = 999,
              ## autoWidth = TRUE, ## scrollX=TRUE,

              scrollY = "25vh",
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

      output$contrastPlot <- renderPlot({
        ct <- rv_contr()
        shiny::req(NCOL(ct))
        plotPhenoDistribution(data.frame(ct))
      })

      ## pointer to reactive
      return(rv_contr)
    } ## end-of-server
  )
}
