##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## --------------------------------------------------------------------------
## convert list of checks to html tags for display in the data preview modal
## --------------------------------------------------------------------------
check_to_html <- function(check, pass_msg = "", null_msg = "", false_msg = "",
                          details = TRUE) {
  error_list <- playbase::PGX_CHECKS
  tags <- NULL

  if (is.null(check)) {
    tags <- tagList(
      span(null_msg, style = "color: red"), br()
    )
  } else if (isFALSE(check)) {
    tags <- tagList(
      span(false_msg, style = "color: orange"), br()
    )
  } else {
    if (length(check) > 0) {
      hr1 <- shiny::hr(style = "border-top: 1px solid black;")
      if (!details) hr1 <- NULL

      tags <- tagList(
        lapply(1:length(check), function(idx) {
          error_id <- names(check)[idx]
          error_log <- check[[idx]]
          error_detail <- error_list[error_list$error == error_id, ]
          error_length <- length(error_log)
          ifelse(length(error_log) > 5, error_log <- error_log[1:5], error_log)
          if (error_detail$warning_type == "warning") {
            title_color <- "orange"
          } else if (error_detail$warning_type == "error") {
            title_color <- "red"
          }
          div(
            hr1,
            span(error_detail$title, style = paste("color:", title_color)),
            shiny::br(),
            ifelse(
              !details,
              "",
              paste(error_detail$message, "\n",
                paste(error_length, "case(s) identified, examples:"),
                paste(error_log, collapse = " "),
                sep = " "
              )
            ),
            shiny::br()
          )
        }),
        hr1
      )
    } else {
      if (pass_msg != "") {
        tags <- tagList(
          span(pass_msg, style = "color: green"), br()
        )
      }
    }
  }
  return(tags)
}

preview_module_legend <- shiny::div(
  class = "pt-4",
  style = "margin-top: 150px;",
  span(style = "color: green", "Green"),
  span("= data OK. "),
  br(),
  span(style = "color: orange", "Orange"),
  span("= warning but data will still be uploaded. "),
  br(),
  span(style = "color:red", "Red"),
  span("= error and data will not be uploaded.")
)

error_popup <- function(title, header, message, error, btn_id, onclick, show_consent = FALSE) {
  showModal(
    shiny::tagList(
      tags$div(
        id = "sendLogModal",
        class = "modal",
        style = "
              display: block;
              position: fixed;
              z-index: 1;
              left: 0;
              top: 0;
              width: 100%;
              height: 100%;
              overflow: auto;
              background-color: rgba(0,0,0,0.4);
            ",
        tags$div(
          class = "modal-content",
          style = "
              background-color: #fefefe;
              margin: 15% auto;
              padding: 20px;
              border: 1px solid #888;
              width: 45%;
              color: black;
              text-align: left;
              height: auto;
              overflow-y: auto;
            ",
          tags$button(
            class = "btn btn-danger", HTML("&times;"),
            onClick = "document.getElementById('sendLogModal').style.display = 'none';",
            style = "
                                  position: absolute;
                                  top: 5px;
                                  right: 5px;
                                  "
          ),
          shiny::tags$h1(
            title,
            style = "color:#BF616A;font-family:lato;"
          ),
          shiny::tags$h2(
            header,
            style = "color:#BF616A;font-family:lato;"
          ),
          shiny::br(),
          # add grey style to tag p, and corner edges
          tags$p(error, style = "font-size:12px; background-color: rgba(0,0,0,0.1); border-radius: 5px; padding: 10px; overflow: auto;"),
          shiny::p(message, style = "font-size:15px;"),
          if (show_consent) {
            tags$p(
              style = "font-size:13px; font-style: italic; color: #5E81AC; background-color: rgba(94, 129, 172, 0.1); border-left: 3px solid #5E81AC; padding: 10px; margin: 10px 0;",
              HTML("<strong>Note:</strong> By clicking \"Send data to customer support\", you authorize BigOmics to access your data solely for the purpose of fixing computation errors, if any.")
            )
          },
          div(
            tags$button(
              id = btn_id,
              class = "btn btn-danger", HTML("Send data to customer support"),
              onclick = onclick,
            )
          )
        )
      )
    )
  )
}

sendErrorMessageToCustomerSuport <- function(user_email, pgx_name, pgx_path, error, path_to_creds = "hubspot_creds") {
  if (!file.exists(path_to_creds)) {
    message("[sendErrorMessageToCustomerSuport] WARNING : ticket not opened. cannot get credential =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)

  message <- glue::glue(
    "The ds name is: {pgx_name}

          The ds path is: {pgx_path}

          The error is:

          {error}"
  )

  payload <- list(
    fields = list(
      list(
        objectTypeId = "0-1",
        name = "email",
        value = user_email
      ),
      list(
        objectTypeId = "0-5",
        name = "subject",
        value = "Technical crash"
      ),
      list(
        objectTypeId = "0-5",
        name = "content",
        value = message
      )
    )
  )

  json_payload <- jsonlite::toJSON(payload, auto_unbox = TRUE)

  # Send the POST request to HubSpot (send data to form api)
  response <- httr::POST(
    url = "https://api.hsforms.com/submissions/v3/integration/secure/submit/24974201/0597a423-f0e4-44b8-bbfc-48d9a6e6309a",
    httr::add_headers(
      `Content-Type` = "application/json",
      `Authorization` = paste("Bearer", readLines(path_to_creds))
    ),
    body = json_payload
  )
}

sendErrorMessageToUser <- function(user_email, pgx_name, error, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    message("[sendErrorMessageToUser] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)
  if (is.null(user_email) || user_email == "") {
    message("[sendSuccessMessageToUser] WARNING : mail not sent. invalid or empty email")
    return(NULL)
  }


  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello,

          We detected a problem in your dataset computation. We are sorry that this happened. For support, please contact us at support@bigomics.ch.

          Please find below the related logs.

          The dataset name is: {pgx_name}

          The error log is:

          {error}

          Yours,

          BigOmics Developers Team"
        )
      ),
      footer = blastula::md(
        glue::glue("Email sent on {blastula::add_readable_time()}.")
      )
    ),
    from = "bigomics.app@gmail.com",
    to = user_email,
    subject = paste("Omics Playground: Error when computing a dataset"),
    credentials = blastula::creds_file(path_to_creds)
  )
}

sendSuccessMessageToUser <- function(user_email, pgx_name, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    message("[sendSuccessMessageToUser] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)

  if (is.null(user_email) || user_email == "") {
    message("[sendSuccessMessageToUser] WARNING : mail not sent. invalid or empty email")
    return(NULL)
  }

  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello,

          Congratulations, the dataset {pgx_name} completed successfully!

          The omics revolution is one click away, simply <a href='https://auth.bigomics.ch/#!/login' target='_blank'>login here.</a>

          Yours,

          BigOmics Developers Team"
        )
      ),
      footer = blastula::md(
        glue::glue("Email sent on {blastula::add_readable_time()}.")
      )
    ),
    from = "bigomics.app@gmail.com",
    to = user_email,
    subject = paste("Omics Playground: Dataset computed successfully!"),
    credentials = blastula::creds_file(path_to_creds)
  )
}

isValidFileName <- function(name) {
  if (name == "") {
    return(FALSE)
  }
  pattern <- "/"
  if (grepl(pattern, name)) {
    return(FALSE)
  }
  return(TRUE)
}

write_check_output <- function(
    checks_list,
    file_type = c("SAMPLES", "COUNTS", "CONTRASTS", "SAMPLES_COUNTS", "SAMPLES_CONTRASTS"),
    raw_dir = raw_dir()) {
  file_type <- match.arg(file_type)
  # write date and hour and no error in cross_check samples counts
  date_hour <- paste0(Sys.time())
  # replece : by _
  date_hour <- gsub(":", "_", date_hour)
  if (length(checks_list) == 0 && !is.null(raw_dir)) {
    write(paste(date_hour, file_type, "PASS", sep = ":  "), file.path(raw_dir, "CHECKS_OUTPUT"), append = TRUE)
  } else if (length(checks_list) > 0 && !is.null(raw_dir)) {
    # Convert each element in the list to a line in the text file
    lines <- sapply(names(checks_list), function(name) {
      paste(date_hour, file_type, "FAILED", name, playbase::PGX_CHECKS[playbase::PGX_CHECKS$error == name, "checks"], paste(checks_list[[name]], collapse = "__"), sep = ": ")
    })
    write(unlist(lines), file.path(raw_dir, "CHECKS_OUTPUT"), append = TRUE)
  }
}
