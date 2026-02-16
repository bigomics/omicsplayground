##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


admin_table_credentials_ui <- function(
  id,
  title = "User Credentials",
  height = c("100%", 800),
  info.text = "",
  caption = ""
) {
  ns <- shiny::NS(id)

  if (length(height) == 1) height <- c(height, 800)

  ifnotchar.int <- function(s) {
    suppressWarnings(
      ifelse(!is.na(as.integer(s)), paste0(as.integer(s), "px"), s)
    )
  }
  height.1 <- ifnotchar.int(height[1])

  header <- shiny::fillRow(
    flex = c(NA, 1, NA),
    class = "tablemodule-header",
    shiny::div(class = "tablemodule-title", title = title, title),
    "",
    DropdownMenu(
      shiny::div(
        class = "tablemodule-info",
        shiny::HTML(paste0("<b>", as.character(title), ".", "</b>", "&nbsp;", as.character(info.text)))
      ),
      width = "250px",
      size = "xs",
      icon = shiny::icon("info"),
      status = "default"
    )
  )

  bslib::card(
    class = "tablemodule",
    full_screen = FALSE,
    style = paste0("height:", height.1, ";overflow:visible;"),
    bslib::as.card_item(div(header)),
    bslib::card_body(
      DT::DTOutput(ns("credentials_tbl"), height = "100%") %>% bigLoaders::useSpinner(),
      shiny::div(
        style = "padding-top: 10px;",
        shiny::actionButton(ns("add_user"), "Add User", icon = shiny::icon("plus"), class = "btn-sm btn-primary"),
        shiny::actionButton(ns("delete_user"), "Delete Selected", icon = shiny::icon("trash"), class = "btn-sm btn-danger"),
        shiny::actionButton(ns("save_credentials"), "Save Changes", icon = shiny::icon("floppy-disk"), class = "btn-sm btn-success"),
        shiny::span(
          style = "margin-left: 15px;",
          shiny::textOutput(ns("status_msg"), inline = TRUE)
        )
      )
    ),
    bslib::card_body(
      class = "card-footer",
      div(class = "caption", shiny::HTML(paste0(
        "<b>", as.character(title), ".</b>",
        "&nbsp;", as.character(caption)
      )))
    )
  )
}

admin_table_credentials_server <- function(id, auth, credentials_file = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    cred_data <- shiny::reactiveVal(NULL)
    status <- shiny::reactiveVal("")
    render_version <- shiny::reactiveVal(0)

    input_style <- "height:20px;padding:0 2px;font-size:11px;border:1px solid #ccc;border-radius:2px;"

    ## Load credentials from file
    shiny::observe({
      shiny::req(isTRUE(auth$ADMIN))
      shiny::req(!is.null(credentials_file))
      shiny::req(file.exists(credentials_file))
      dbg("[admin_table_credentials] loading credentials from: ", credentials_file)
      df <- read.csv(credentials_file, colClasses = "character", stringsAsFactors = FALSE)
      cred_data(df)
      render_version(isolate(render_version()) + 1)
    })

    ## Build display data with only email, password, admin columns.
    ## Password and admin use custom HTML inputs; email is plain text.
    display_data <- shiny::reactive({
      render_version()
      df <- shiny::isolate(cred_data())
      shiny::req(df)
      current_email <- tolower(trimws(auth$email))

      display <- data.frame(
        email = df$email,
        password = sapply(seq_len(nrow(df)), function(i) {
          pw_escaped <- gsub('"', '&quot;', df$password[i])
          sprintf(
            paste0(
              '<input type="password" class="pw-input" data-row="%d" value="%s"',
              ' ondblclick="this.type=\'text\';this.select();"',
              ' onblur="this.type=\'password\';"',
              ' style="%swidth:100px;">'
            ),
            i, pw_escaped, input_style
          )
        }),
        admin = sapply(seq_len(nrow(df)), function(i) {
          val <- toupper(trimws(df$ADMIN[i]))
          is_self <- identical(tolower(trimws(df$email[i])), current_email)
          sprintf(
            paste0(
              '<select class="admin-select" data-row="%d" style="%s%s"%s>',
              '<option value="TRUE"%s>TRUE</option>',
              '<option value="FALSE"%s>FALSE</option>',
              '</select>'
            ),
            i,
            input_style,
            if (is_self) "opacity:0.5;" else "",
            if (is_self) " disabled" else "",
            if (val == "TRUE") " selected" else "",
            if (val != "TRUE") " selected" else ""
          )
        }),
        stringsAsFactors = FALSE
      )
      display
    })

    ## Render datatable (no DT inline editing; password and admin use custom HTML)
    output$credentials_tbl <- DT::renderDT(
      {
        shiny::req(display_data())
        DT::datatable(
          display_data(),
          class = "compact hover",
          rownames = FALSE,
          escape = 1,
          editable = FALSE,
          selection = list(mode = "single", target = "row"),
          callback = DT::JS(sprintf(
            "table.on('change', 'select.admin-select', function() {
               var row = parseInt($(this).data('row'));
               var val = $(this).val();
               Shiny.setInputValue('%s', {row: row, value: val}, {priority: 'event'});
             });
             table.on('change', 'input.pw-input', function() {
               var row = parseInt($(this).data('row'));
               var val = $(this).val();
               Shiny.setInputValue('%s', {row: row, value: val}, {priority: 'event'});
             });",
            ns("admin_select_change"),
            ns("password_change")
          )),
          options = list(
            dom = "lfrtip",
            scrollX = TRUE,
            scrollY = "calc(100vh - 400px)",
            deferRender = TRUE,
            pageLength = 50
          )
        ) %>%
          DT::formatStyle(0, target = "row", fontSize = "11px", lineHeight = "70%")
      },
      server = FALSE
    )

    ## Handle ADMIN dropdown changes
    shiny::observeEvent(input$admin_select_change, {
      shiny::req(isTRUE(auth$ADMIN))
      change <- input$admin_select_change
      df <- cred_data()
      df$ADMIN[change$row] <- change$value
      cred_data(df)
      status("Unsaved changes")
    })

    ## Handle password input changes
    shiny::observeEvent(input$password_change, {
      shiny::req(isTRUE(auth$ADMIN))
      change <- input$password_change
      df <- cred_data()
      df$password[change$row] <- change$value
      cred_data(df)
      status("Unsaved changes")
    })

    ## Add a new user row
    shiny::observeEvent(input$add_user, {
      shiny::req(isTRUE(auth$ADMIN))
      df <- cred_data()
      new_row <- data.frame(
        username = "new_user",
        email = "user@example.com",
        password = "changeme",
        expiry = format(Sys.Date() + 365, "%Y-%m-%d"),
        level = "free",
        limit = "3",
        ADMIN = "FALSE",
        stringsAsFactors = FALSE
      )
      cred_data(rbind(df, new_row))
      render_version(render_version() + 1)
      status("New user added (unsaved)")
    })

    ## Delete the selected user row
    shiny::observeEvent(input$delete_user, {
      shiny::req(isTRUE(auth$ADMIN))
      sel <- input$credentials_tbl_rows_selected
      if (is.null(sel) || length(sel) == 0) {
        status("No row selected for deletion")
        return()
      }
      df <- cred_data()
      df <- df[-sel, , drop = FALSE]
      cred_data(df)
      render_version(render_version() + 1)
      status("User deleted (unsaved)")
    })

    ## Save changes back to CSV
    shiny::observeEvent(input$save_credentials, {
      shiny::req(isTRUE(auth$ADMIN))
      shiny::req(!is.null(credentials_file))
      df <- cred_data()
      shiny::req(!is.null(df))
      tryCatch(
        {
          write.csv(df, file = credentials_file, row.names = FALSE, quote = TRUE)
          dbg("[admin_table_credentials] saved credentials to: ", credentials_file)
          status(paste("Saved successfully at", format(Sys.time(), "%H:%M:%S")))
        },
        error = function(e) {
          status(paste("Error saving:", e$message))
        }
      )
    })

    ## Render status message
    output$status_msg <- shiny::renderText({
      status()
    })
  }) ## end of moduleServer
}
