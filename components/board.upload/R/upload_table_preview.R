# NOTES
#SAMPLES
#- all tabs should be the same size
#- convert column types from numeric to categorical
#- remove column
#- convert columns to categorical via binning and cuts
#
#CONTRASTS
#- make contrasts table editable
#
#- Save updated files
#
#- put the errors in the same tab
#- make the modal the whole screen and put the errors on the side
#

# no ui right now because the preview is in a modal dialog
upload_table_preview_ui <- function(id) {
  ns <- NS(id)
  tagList(

  )
}

upload_table_preview_server <- function(id, uploaded) {
  moduleServer(
    id,
    function(input, output, session) {

        rv_preview <- reactiveValues()
        # every time something is uploaded, it can be previewed
        observeEvent(uploaded$last_uploaded, {

            rv_preview$counts_approval <- 'Not Handled Yet'
            rv_preview$samples_approval <- 'Not Handled Yet'
            rv_preview$contrasts_approval <- 'Not Handled Yet'
            rv_preview$tab_order <- c()

            tabs <- list(
                id = session$ns('preview_panel')
            )
            has_counts <- 'counts.csv' %in% uploaded$last_uploaded
            has_samples <- 'samples.csv' %in% uploaded$last_uploaded
            has_contrasts <- ('contrasts.csv' %in% uploaded$last_uploaded) & (!is.null(uploaded$contrasts.csv))

            if (has_counts) {
                rv_preview$tab_order <- c(rv_preview$tab_order, 'Counts')
                tabs <- c(
                    tabs,
                    list(tabPanel(
                        'Counts',
                        fluidRow(
                            column(
                                width = 8,
                                div(
                                    shiny::actionButton(
                                        session$ns('convert_coltype'),
                                        label = 'Convert column type',
                                        style = 'display: inline-block',
                                        class = 'btn-secondary btn-sm'
                                    ),
                                    shiny::actionButton(
                                        session$ns('set_rownames'),
                                        label = 'Set column as rownames',
                                        style = 'display: inline-block',
                                        class = 'btn-secondary btn-sm'
                                    ),
                                    shiny::actionButton(
                                        session$ns('remove_col'),
                                        label = 'Remove column',
                                        style = 'display: inline-block',
                                        class = 'btn-secondary btn-sm'
                                    )
                                ),
                                br(),
                                DT::dataTableOutput(session$ns("counts_preview"))
                            ),
                            column(
                                width = 4,
                                'Summary:'
                            )
                        ),
                        br(),
                        div(
                            style = 'float: right',
                            shiny::actionButton(
                                session$ns('discard_counts'),
                                label = 'Discard counts',
                                style = 'display: inline-block; margin-left: 15px;',
                                class = 'btn-danger'
                            ),
                            shiny::actionButton(
                                session$ns('approve_counts'),
                                label = 'Approve counts',
                                style = 'display: inline-block',
                                class = 'btn-info'
                            )
                        )
                    ))
                )
            }
            if (has_samples) {
                rv_preview$tab_order <- c(rv_preview$tab_order, 'Samples')
                tabs <- c(
                    tabs,
                    list(tabPanel(
                        'Samples',
                        fluidRow(
                            column(
                                width = 8,
                                DT::dataTableOutput(session$ns("samples_preview"))
                            ),
                            column(
                                width = 4,
                                'Summary:'
                            )
                        ),
                        br(),
                        div(
                            style = 'float: right',
                            shiny::actionButton(
                                session$ns('discard_samples'),
                                label = 'Discard samples',
                                style = 'display: inline-block; margin-left: 15px;',
                                class = 'btn-danger'
                            ),
                            shiny::actionButton(
                                session$ns('approve_samples'),
                                label = 'Approve samples',
                                style = 'display: inline-block',
                                class = 'btn-info'
                            )
                        )
                    ))
                )
            }
            if (has_contrasts) {
                rv_preview$tab_order <- c(rv_preview$tab_order, 'Contrasts')
                tabs <- c(
                    tabs,
                    list(tabPanel(
                        'Contrasts',
                        fluidRow(
                            column(
                                width = 8,
                                DT::dataTableOutput(session$ns("contrasts_preview"))
                            ),
                            column(
                                width = 4,
                                'Summary:'
                            )
                        ),
                        fluidRow(
                           column(
                               width = 12,
                               br(),
                               div(
                                   style = 'float: right',
                                   shiny::actionButton(
                                       session$ns('discard_contrasts'),
                                       label = 'Discard contrasts',
                                       style = 'display: inline-block; margin-left: 15px;',
                                       class = 'btn-danger'
                                   ),
                                   shiny::actionButton(
                                       session$ns('approve_contrasts'),
                                       label = 'Approve contrasts',
                                       style = 'display: inline-block',
                                       class = 'btn-info'
                                   )
                               )
                           )
                        )

                    ))
                )
            }

            # add summary tab
            # this could include any diagnostics / errors for the uploaded data
            # and could actively check which datasets have been approved/discarded
            # so users would have to go here to close the modal
            summary_tab <- list()
            rv_preview$tab_order <- c(rv_preview$tab_order, 'Summary')

            if (has_counts) {
                summary_tab <- c(
                    summary_tab,
                    tagList(
                        'Counts:', shiny::uiOutput(session$ns('counts_approval')),
                        br()
                    )
                )
            }
            if (has_samples) {
                summary_tab <- c(
                    summary_tab,
                    tagList(
                        'Samples:', shiny::uiOutput(session$ns('samples_approval')),
                        br()
                    )
                )
            }
            if (has_contrasts) {
                summary_tab <- c(
                    summary_tab,
                    tagList(
                        'Contrasts:', shiny::uiOutput(session$ns('contrasts_approval')),
                        br()
                    )
                )
            }
            summary_tab <- do.call(tagList, summary_tab)
            tabs <- c(
                tabs,
                list(tabPanel(
                    'Summary',
                    'This is a summary of the uploaded data:',
                    br(), br(),
                    summary_tab,
                    shiny::actionButton(
                        session$ns('finish_preview'),
                        label = 'Finish Preview',
                        style = 'float: right;',
                        class = 'btn-success'
                    ) %>% shinyjs::disabled()
                ))
            )

            tab_panel <- do.call(tabsetPanel, tabs)
            shiny::showModal(
                shiny::modalDialog(
                    title = 'Data Upload Preview',
                    label = 'this is a label',
                    tab_panel,
                    footer = NULL,
                    easyClose = FALSE,
                    size = 'xl'
                ) %>%
                    tagAppendAttributes(
                        style = 'min-height: 90%; min-width: 90%',
                        .cssSelector = '.modal-dialog'
                    )
            )

        }, ignoreNULL = TRUE)

        # function to format datatable
        make_preview_table <- function(tbl) {

            # callback for highlighting column instead of row
            js <- c("table.on('mouseenter', 'td', function () {
                        // Remove highlight from all columns
                        table
                        .columns()
                        .nodes()
                        .flatten()  // Reduce to a 1D array
                        .to$()      // Convert to a jQuery object
                        .removeClass( 'highlight' );

                                // Add highlight to mouseover column
                        table
                        .column( this )
                        .nodes()
                        .to$()      // Convert to a jQuery object
                        .addClass( 'highlight' );
                    });",
                    "table.on('mouseleave', 'td', function () {
                        // Remove highlight from all columns
                        table
                        .columns()
                        .nodes()
                        .flatten()  // Reduce to a 1D array
                        .to$()      // Convert to a jQuery object
                        .removeClass( 'highlight' );
                    });")

            DT::datatable(tbl,
                          class = "compact",
                          rownames = TRUE,
                          options = list(
                              dom = "rtp",
                              pageLength = 20
                          ),
                          callback = DT::JS(js),
                          selection = list(target = 'column', mode = 'single')
            ) %>%
                DT::formatStyle(0, target = "row",
                                fontSize = "12px", lineHeight = "70%")
        }

        output$counts_preview <- DT::renderDataTable({
            make_preview_table(uploaded$counts.csv)
        })
        output$counts_approval <- shiny::renderUI({
            txt <- rv_preview$counts_approval
            if (txt == 'Discarded') {
                div(txt, style = 'color: red;')
            } else if (txt == 'Approved') {
                div(txt, style = 'color: green;')
            } else {
                div(txt, style = 'color: orange;')
            }
        })
        observeEvent(input$discard_counts, {
            rv_preview$counts_approval <- 'Discarded'
            shiny::hideTab(inputId = 'preview_panel', target = 'Counts')
            rv_preview$tab_order <- rv_preview$tab_order[rv_preview$tab_order != 'Counts']
            shiny::updateTabsetPanel(inputId = 'preview_panel', selected = rv_preview$tab_order[1])
        }, ignoreInit = TRUE)
        observeEvent(input$approve_counts, {
            rv_preview$counts_approval <- 'Approved'
            shiny::hideTab(inputId = 'preview_panel', target = 'Counts')
            rv_preview$tab_order <- rv_preview$tab_order[rv_preview$tab_order != 'Counts']
            shiny::updateTabsetPanel(inputId = 'preview_panel', selected = rv_preview$tab_order[1])
        }, ignoreInit = TRUE)

        output$samples_preview <- DT::renderDataTable({
            make_preview_table(uploaded$samples.csv)
        })
        output$samples_approval <- shiny::renderUI({
            txt <- rv_preview$samples_approval
            if (txt == 'Discarded') {
                div(txt, style = 'color: red;')
            } else if (txt == 'Approved') {
                div(txt, style = 'color: green;')
            } else {
                div(txt, style = 'color: orange;')
            }
        })
        observeEvent(input$discard_samples, {
            rv_preview$samples_approval <- 'Discarded'
            shiny::hideTab(inputId = 'preview_panel', target = 'Samples')
            rv_preview$tab_order <- rv_preview$tab_order[rv_preview$tab_order != 'Samples']
            shiny::updateTabsetPanel(inputId = 'preview_panel', selected = rv_preview$tab_order[1])
        }, ignoreInit = TRUE)
        observeEvent(input$approve_samples, {
            rv_preview$samples_approval <- 'Approved'
            shiny::hideTab(inputId = 'preview_panel', target = 'Samples')
            rv_preview$tab_order <- rv_preview$tab_order[rv_preview$tab_order != 'Samples']
            shiny::updateTabsetPanel(inputId = 'preview_panel', selected = rv_preview$tab_order[1])
        }, ignoreInit = TRUE)

        output$contrasts_preview <- DT::renderDataTable({
            make_preview_table(uploaded$contrasts.csv)
        })
        output$contrasts_approval <- shiny::renderUI({
            txt <- rv_preview$contrasts_approval
            if (txt == 'Discarded') {
                div(txt, style = 'color: red;')
            } else if (txt == 'Approved') {
                div(txt, style = 'color: green;')
            } else {
                div(txt, style = 'color: orange;')
            }
        })
        observeEvent(input$discard_contrasts, {
            rv_preview$contrasts_approval <- 'Discarded'
            shiny::hideTab(inputId = 'preview_panel', target = 'Contrasts')
            rv_preview$tab_order <- rv_preview$tab_order[rv_preview$tab_order != 'Contrasts']
            shiny::updateTabsetPanel(inputId = 'preview_panel', selected = rv_preview$tab_order[1])
        }, ignoreInit = TRUE)
        observeEvent(input$approve_contrasts, {
            rv_preview$contrasts_approval <- 'Approved'
            shiny::hideTab(inputId = 'preview_panel', target = 'Contrasts')
            rv_preview$tab_order <- rv_preview$tab_order[rv_preview$tab_order != 'Contrasts']
            shiny::updateTabsetPanel(inputId = 'preview_panel', selected = rv_preview$tab_order[1])
        }, ignoreInit = TRUE)


        # when the preview is finished, close the modal and NULLify any
        # datasets that were discarded
        observeEvent(input$finish_preview, {
            shiny::removeModal()

            if (rv_preview$counts_approval == 'Discarded') {
                uploaded$counts.csv <- NULL
            }
            if (rv_preview$samples_approval == 'Discarded') {
                uploaded$samples.csv <- NULL
            }
            if (rv_preview$contrasts_approval == 'Discarded') {
                uploaded$contrasts.csv <- NULL
            }
        })

        # only enable finish-preview button when all datasets have been handled
        observeEvent(c(
            rv_preview$counts_approval,
            rv_preview$samples_approval,
            rv_preview$contrasts_approval
        ), {
            has_counts <- 'counts.csv' %in% uploaded$last_uploaded
            has_samples <- 'samples.csv' %in% uploaded$last_uploaded
            has_contrasts <- ('contrasts.csv' %in% uploaded$last_uploaded) & (!is.null(uploaded$contrasts.csv))

            enable_button <- TRUE
            if (has_counts) {
                if (rv_preview$counts_approval == 'Not Handled Yet') {
                    enable_button <- FALSE
                }
            }
            if (has_samples) {
                if (rv_preview$samples_approval == 'Not Handled Yet') {
                    enable_button <- FALSE
                }
            }
            if (has_contrasts) {
                if (rv_preview$contrasts_approval == 'Not Handled Yet') {
                    enable_button <- FALSE
                }
            }

            if (enable_button) shinyjs::enable(id = 'finish_preview')
        }, ignoreInit = TRUE)


    }
  )
}