##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

error_popup <- function(title, header, message, error, btn_id, onclick) {
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

sendErrorMessageToCustomerSuport <- function(user_email, pgx_name, pgx_path, error, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    info("[sendShareMessage] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)

  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello,

          The user <strong>{user_email}</strong> had a dataset that failed to compute and requested help.

          Please find below the log and data path.

          The ds name is: {pgx_name}

          The ds path is: {pgx_path}

          The error is:

          {error}

          Yours,

          Team Developers"
        )
      ),
      footer = blastula::md(
        glue::glue("Email sent on {blastula::add_readable_time()}.")
      )
    ),
    from = "bigomics.app@gmail.com",
    to = "support@bigomics.ch",
    subject = paste("Problem with dataset from", user_email),
    credentials = blastula::creds_file(path_to_creds)
  )
}
sendErrorMessageToUser <- function(user_email, pgx_name, error, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    info("[sendShareMessage] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)

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

          We are an open source company that value transparency. Please consider supporting us by purchasing a subscription at https://bigomics.ch/pricing/ and get access to advanced funcionalities!

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
    subject = paste("OmicsPlayground: Error when computing a dataset"),
    credentials = blastula::creds_file(path_to_creds)
  )
}

sendSuccessMessageToUser <- function(user_email, pgx_name, path_to_creds = "gmail_creds") {
  if (!file.exists(path_to_creds)) {
    info("[compute PGX process] WARNING : mail not sent. cannot get mail creds =", path_to_creds)
    return(NULL)
  }

  user_email <- trimws(user_email)

  blastula::smtp_send(
    blastula::compose_email(
      body = blastula::md(
        glue::glue(
          "Hello,

          Congratulation, the dataset {pgx_name} completed successfully!

          The omics revolution is one click away, go to https://bigomics.ch/ and find your dataset in Omics Playground Loading Module.

          We are an open source company that value transparency. Please consider supporting us by purchasing a subscription at https://bigomics.ch/pricing/ and get access to advanced funcionalities!

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
    subject = paste("OmicsPlayground: Dataset computed successfully!"),
    credentials = blastula::creds_file(path_to_creds)
  )
}
