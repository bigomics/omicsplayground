
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

            rv_preview$counts_approval <- 'Not Approved'
            rv_preview$samples_approval <- 'Not Approved'
            rv_preview$contrasts_approval <- 'Not Approved'
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
                        DT::dataTableOutput(session$ns("counts_preview")),
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
                        DT::dataTableOutput(session$ns("samples_preview")),
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
                        DT::dataTableOutput(session$ns("contrasts_preview")),
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
                    )
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
                )
            )

        }, ignoreNULL = TRUE)

        output$counts_preview <- DT::renderDataTable({ uploaded$counts.csv })
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

        output$samples_preview <- DT::renderDataTable({ uploaded$samples.csv })
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

        output$contrasts_preview <- DT::renderDataTable({ uploaded$contrasts.csv })
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

        observeEvent(input$finish_preview, {
            shiny::removeModal()
        })


    }
  )
}